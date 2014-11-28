%% Prepare MNIST data
clear all;clc
addpath(genpath('../libsvm-master'));
addpath('../MNIST');
load('mnist_complete.mat');
digit_num = 5;%number of digits used in training
G0_digit_num = digit_num;%number of digits used in G0 
G0_samp_num = 500;%number of nodes per cluster in G0

%% Random sampling from MNIST
%Randomly pick some digits
sel_digit = randperm(10,digit_num);

% Find max training set digit sample number
train_samp_num = Inf;
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(train_y(:,digit) == 1);
    train_samp_num = min(train_samp_num,numel(idx));    
end

train_array = zeros(train_samp_num,digit_num);
train_label = ones(train_samp_num,digit_num);
%Randomly sample images with these labels
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(train_y(:,digit) == 1);
    train_samp_idx = randperm(numel(idx),train_samp_num);
    train_array(:,i) = idx(train_samp_idx);
    train_label(:,i) = (digit-1)*train_label(:,i);%Change the sel_digit from indices to integer label
end

train_mat = double(train_x(train_array(:),:))/255;
train_label = double(train_label(:));


%% Libsvm training and testing
%Use default svm setting
fprintf('start svm training...\n');
model = svmtrain(train_label,train_mat,'-t 0');% Use Linear SVM only

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
test_data = zeros(graph_sz,784);
for i = 1:graph_sz
    test_data(i,:) = double(test_x(test_vect(i),:))/255;
end

p_label = svmpredict(label_vect, test_data, model);
for i = 1:graph_sz
    for j = i+1:graph_sz
        G0(i,j) = (p_label(i) == p_label(j));
    end
end

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
 
%% Backup code

% %% Generate data for pair-wise SVM training
% % Generate negative pairs
% neg_num = 4000;
% neg_pair = zeros(neg_num,2);
% for j = 1:neg_num
%    pair = randperm(digit_num,2);
%    id1 = sel_array(randperm(samp_num,1),pair(1));
%    id2 = sel_array(randperm(samp_num,1),pair(2));
%    neg_pair(j,:) = [id1,id2];
% end
% 
% % Generate positive pairs
% pos_num = 4000;
% pos_pair = zeros(pos_num,2);
% for j = 1:pos_num
%    cluster_id = randperm(digit_num,1);
%    ids = sel_array(randperm(samp_num,2),cluster_id);
%    pos_pair(j,:) = ids';
% end
% 
% % Construct the data matrix
% pair_mat1 = cat(1,neg_pair(:,1),pos_pair(:,1));
% pair_mat2 = cat(1,neg_pair(:,2),pos_pair(:,2));
% pair_num = size(pair_mat1,1);
% full_mat1 = zeros(pair_num,784);
% full_mat2 = zeros(pair_num,784);
% full_label1 = zeros(pair_num,1);
% full_label2 = zeros(pair_num,1);
% for i = 1:pair_num
%     full_mat1(i,:) = train_x(pair_mat1(i),:);
%     full_mat2(i,:) = train_x(pair_mat2(i),:);
%     full_label1(i) = find(train_y(pair_mat1(i),:),1,'first')-1;%find 1 index in each row
%     full_label2(i) = find(train_y(pair_mat2(i),:),1,'first')-1;%find 1 index in each row
% end
% %SVM pair-wise testing
% fprintf('start svm testing...\n');
% p_label1 = svmpredict(full_label1, full_mat1, model);
% p_label2 = svmpredict(full_label2, full_mat2, model);
% sim = (p_label1 == p_label2);
% fprintf('svm testing result: %f',sum(sim)/numel(p_label1));
    
