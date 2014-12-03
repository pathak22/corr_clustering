%% Prepare MNIST data
clear all;clc
addpath(genpath('../libsvm-master'));
addpath('../MNIST');
load('mnist_complete.mat');
load('mnist_region_feat.mat');
digit_num = 3;%number of digits used in training
G0_digit_num = digit_num;%number of digits used in G0 
G0_samp_num = 100;%number of nodes per cluster in G0
sel_digit = randperm(10,digit_num);
fprintf('Graph size, num = %d, samp_num = %d\n',digit_num,G0_samp_num);

%% Generate the G0 map
train_sel_array = zeros(G0_samp_num,G0_digit_num);
train_label_array = ones(G0_samp_num,G0_digit_num);

%Randomly sample images with these labels
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(train_y(:,digit) == 1);
    perm_idx = randperm(numel(idx));
    samp_idx = perm_idx(1:G0_samp_num);
    train_sel_array(:,i) = idx(samp_idx);
    train_label_array(:,i) = (digit-1)*train_label_array(:,i);
end

% Generate all pair-wise similarity label by SVM
graph_sz = G0_samp_num*G0_digit_num;
train_sel_vect = train_sel_array(:);
train_label_vect = train_label_array(:);
train_mat = zeros(graph_sz*(graph_sz-1)/2,28);
train_label = zeros(graph_sz*(graph_sz-1)/2,1);
cnt = 1;
for i = 1:graph_sz-1
    for j = i+1:graph_sz
        train_mat(cnt,:) = cat(2,train_feat(train_sel_vect(i),:),train_feat(train_sel_vect(j),:));
        train_label(cnt) = (train_label_vect(i) == train_label_vect(j));
        cnt = cnt+1;
    end
end

%% Libsvm training and testing
%Use default svm setting
fprintf('start svm training...\n');
tic;
model = svmtrain(train_label,train_mat,'-t 0');% Use Linear SVM only
toc;
fprintf('Training finished...\n');
% testout = test_pairwise(test_x,test_y,sel_digit,model,'full');


%% Generate the G0 map
test_sel_array = zeros(G0_samp_num,G0_digit_num);
label_array = ones(G0_samp_num,G0_digit_num);

%Randomly sample images with these labels
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(test_y(:,digit) == 1);
    perm_idx = randperm(numel(idx));
    samp_idx = perm_idx(1:G0_samp_num);
    test_sel_array(:,i) = idx(samp_idx);
    label_array(:,i) = (digit-1)*label_array(:,i);
end

% Generate all pair-wise similarity label by SVM
graph_sz = G0_samp_num*G0_digit_num;
test_vect = test_sel_array(:);
label_vect = label_array(:);
G0 = zeros(graph_sz);
accum_accuracy = 0;

fprintf('Start computing the similarity map...\n');
tic;
for i = 1:graph_sz-1
% for i = 1
    fprintf('iteration i = %d\n',i);
    test_data = zeros(graph_sz-i,28);
    test_label = zeros(graph_sz-i,1);
    for j = i+1:graph_sz
        test_data(j-i,:) = cat(2,test_feat(test_vect(i),:),test_feat(test_vect(j),:));
        test_label(j-i) = (label_vect(i)==label_vect(j));
    end
    [predict_label, accuracy, ~]= svmpredict(test_label, test_data, model);
    G0(i,i+1:end) = predict_label';
    accum_accuracy = accum_accuracy + accuracy*(graph_sz-i); 
%     clear test_data test_label
end
toc;
sum_accuracy = accum_accuracy/(graph_sz*(graph_sz-1)/2);
display(sum_accuracy(1));
fprintf('Finished computing the similarity map...\n');

G0 = G0+G0';

%% Write G0 to .txt file
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

% Let the diagonal of G0 be zeros
for i = 1:graph_sz
    G0(i,i) = 1;
end

save(['G0_' num2str(graph_sz)],'G0','-v7.3');
    
compute_metrics();

