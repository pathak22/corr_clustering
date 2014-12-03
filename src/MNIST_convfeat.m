clear all;
addpath(genpath('../DeepLearning'));
addpath('../MNIST');
load('mnist_complete.mat');
load('cnn_mnist_half.mat');

test_x = double(reshape(test_x',28,28,10000))/255;
test_y = double(test_y');

% Sample the digits 5-9 for testing
label_y = zeros(1,10000);
for i = 1:10000
    label_y(i) = find(test_y(:,i));
end
idx_half = label_y>5;
test_yh = test_y(:,idx_half);
test_xh = test_x(:,:,idx_half);
test_xh = test_xh(:,:,1:4800);
test_yh = test_yh(:,1:4800);

digit_num = 5;%number of digits used in testing
opt_thr = 0.2;
scale = 12;
G0_digit_num = digit_num;%number of digits used in G0 
G0_samp_num = 200;%number of nodes per cluster in G0, has to be multiples of 20
sel_digit = 6:10;

%% Prepare test data
test_sel_array = zeros(G0_samp_num,G0_digit_num);
test_label_array = ones(G0_samp_num,G0_digit_num);

%Randomly sample images with these labels
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(test_yh(digit,:) == 1);
    
    perm_idx = randperm(numel(idx));
    samp_idx = perm_idx(1:G0_samp_num);
    test_sel_array(:,i) = idx(samp_idx);
    test_label_array(:,i) = (digit-1)*test_label_array(:,i);
end

test_d = test_sel_array(:);
test_l = test_label_array(:);
n_test = numel(test_d);

%Prepare test data
test_data = zeros(28,28,n_test);
test_label = zeros(10,n_test);
for i = 1:n_test
    test_data(:,:,i) = test_xh(:,:,test_d(i));
    test_label(:,i) = test_yh(:,test_d(i));
end

%% Feature and pair-wise distance computation
% forward pass to obtain features
fv = cnntest_fv(cnn,test_data);

% Compute pairwise distance (Euclidean)
fv = fv';
fdist = real(pdist2(fv,fv,'euclidean'));

%visualize the distribution of pdist2
% nbins = 50;
% figure(1);hist(fdist(:),nbins);

% Threshold the distance
thr = prctile(fdist(:),(100/digit_num)*opt_thr*scale);

G0 = fdist<thr;
graph_sz = size(fdist,1);

%% Write G0 to .txt file 
% Make sure there is no self-edges when you run the file
for i = 1:graph_sz
    G0(i,i) = 0;
end

fprintf('Writing files...\n');
fileID = fopen('../data/data_mnist.txt','w');
fprintf(fileID,'# This is MNIST graph. \n');
fprintf(fileID,'# Format is: Ground truth ID, Vertex ID, Neighbors. \n');
fprintf(fileID,'%d\n',graph_sz);
for i = 1:graph_sz
    fprintf(fileID,'%d %d',test_l(i),i);
    for j = 1:graph_sz
        if (G0(i,j)==1)
            fprintf(fileID,' %d',j);
        end
    end
    fprintf(fileID,'\n');
end
fclose(fileID);

for i = 1:graph_sz
    G0(i,i) = 1;
end

save(['G0conv_' num2str(digit_num) '_' num2str(G0_samp_num)],'G0','-v7.3');





