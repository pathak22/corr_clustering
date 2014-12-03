%% Load MNIST data
clear all;
% addpath(genpath('/Users/wckuo/Documents/MATLAB/DeepLearning'));
addpath(genpath('../DeepLearning'));
% addpath(genpath('/Library/Frameworks/Python.framework/Versions/3.4/bin/python3'));
% addpath('/Applications/MATLAB_R2014a.app/toolbox/matlab/general');
addpath('../MNIST');
load('mnist_complete.mat');
% load('mnist_region_feat.mat');
train_x = double(reshape(train_x',28,28,60000))/255;
train_y = double(train_y');
test_x = double(reshape(test_x',28,28,10000))/255;
test_y = double(test_y');

digit_num = 4;%number of digits used in training
G0_digit_num = digit_num;%number of digits used in G0 
G0_samp_num = 20;%number of nodes per cluster in G0, has to be multiples of 20
sel_digit = randperm(10,digit_num);
fprintf('Graph size, num = %d, samp_num = %d\n',digit_num,G0_samp_num);

%% Generate the G0 map
train_sel_array = zeros(G0_samp_num,G0_digit_num);
train_label_array = ones(G0_samp_num,G0_digit_num);

%Randomly sample images with these labels
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(train_y(digit,:) == 1);
    perm_idx = randperm(numel(idx));
    samp_idx = perm_idx(1:G0_samp_num);
    train_sel_array(:,i) = idx(samp_idx);
    train_label_array(:,i) = (digit-1)*train_label_array(:,i);
end

% Generate all pair-wise similarity label by SVM
graph_sz = G0_samp_num*G0_digit_num;
train_sel_vect = train_sel_array(:);
train_label_vect = train_label_array(:);
N_d = graph_sz*(graph_sz-1)/2;

train_mat = zeros(28,56,N_d);
train_label = zeros(2,N_d);
cnt = 1;
for i = 1:graph_sz-1
    for j = i+1:graph_sz
        train_mat(:,:,cnt) = cat(2,train_x(:,:,train_sel_vect(i)),train_x(:,:,train_sel_vect(j)));
        label_train = train_label_vect(i) == train_label_vect(j);
        train_label(label_train+1,cnt) = 1;
        cnt = cnt+1;
    end
end


%% Setup CNN
rand('state',0)
cnn.layers = {
    struct('type', 'i') %input layer
    struct('type', 'c', 'outputmaps', 6, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %sub sampling layer
    struct('type', 'c', 'outputmaps', 12, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %subsampling layer
};

opts.alpha = 1;
opts.batchsize = digit_num*10;%assume nodes per cluster are multiples of 10
% opts.numepochs = graph_sz*(graph_sz-1)/(digit_num*10);
opts.numepochs = digit_num;
v_d = randperm(N_d);
train_mat = train_mat(:,:,v_d);
train_label = train_label(:,v_d);

cnn = cnnsetup(cnn, train_mat, train_label);
fprintf('start training');
cnn = cnntrain(cnn, train_mat, train_label, opts);

%% Test CNN
% G0_samp_num = G0_samp_num*10;
test_sel_array = zeros(G0_samp_num,G0_digit_num);
label_array = ones(G0_samp_num,G0_digit_num);

%Randomly sample images with these labels
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(test_y(digit,:) == 1);
    perm_idx = randperm(numel(idx));
    samp_idx = perm_idx(1:G0_samp_num);
    test_sel_array(:,i) = idx(samp_idx);
    label_array(:,i) = (digit-1)*label_array(:,i);
end

% Generate all pair-wise similarity label by CNN
graph_sz = G0_samp_num*G0_digit_num;
test_vect = test_sel_array(:);
label_vect = label_array(:);
G0 = zeros(graph_sz);
accum_er = 0;

fprintf('Start computing the similarity map...\n');
tic;
for i = 1:graph_sz-1
    fprintf('iteration i = %d\n',i);
    test_data = zeros(28,56,graph_sz-i);
    test_label = zeros(2,graph_sz-i);
    for j = i+1:graph_sz
        test_data(:,:,j-i) = cat(2,test_x(:,:,test_vect(i)),test_x(:,:,test_vect(j)));
        label_test = label_vect(i)==label_vect(j);
        test_label(label_test+1,j-i) = 1;
    end
    
    if i ~= graph_sz-1
        [er, bad] = cnntest(cnn, test_data, test_label);
        predict_label = test_label;
        predict_label(:,bad) = ~predict_label(:,bad);
        G0(i,i+1:end) = predict_label(2,:);
        accum_er = accum_er + er*(graph_sz-i);
    else
        [er, bad] = cnntest(cnn, repmat(test_data,1,1,2), repmat(test_label,1,2));
        if isempty(bad)
            predict_label = test_label;
        else
            predict_label = ~test_label;
        end
        G0(i,i+1:end) = predict_label(2,:);
        accum_er = accum_er + er*(graph_sz-i);

    end
end
toc;
sum_er = accum_er/(graph_sz*(graph_sz-1)/2-1);
display(sum_er);
fprintf('Finished computing the similarity map...\n');

G0 = G0+G0';

%% Write G0 to .txt file
fprintf('Writing files...\n');
fileID = fopen('../data/data_mnist.txt','w');
fprintf(fileID,'# This is MNIST graph. \n');
fprintf(fileID,'# Format is: Ground truth ID, Vertex ID, Neighbors. \n');
fprintf(fileID,'%d\n',graph_sz);
for i = 1:graph_sz
    fprintf(fileID,'%d %d',label_vect(i),i);
    for j = 1:graph_sz
        if (G0(i,j)==1)
            fprintf(fileID,' %d',j);
        end
    end
    fprintf(fileID,'\n');
end
fclose(fileID);


%% Compute metrics
for i = 1:graph_sz
    G0(i,i) = 1;
end
save(['G0_' num2str(digit_num) '_' num2str(G0_samp_num)],'G0','-v7.3');

fprintf('Now re-compute the solution_mnist in python ...\n');


