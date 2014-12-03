%% Prepare MNIST data
clear all;clc
addpath(genpath('../libsvm-master'));
addpath('../MNIST');
load('mnist_complete.mat');
digit_num = 3;%number of digits used in training
G0_digit_num = digit_num;%number of digits used in G0 
G0_samp_num = 30;%number of nodes per cluster in G0
sel_digit = randperm(10,digit_num);
fprintf('Graph size, num = %d, samp_num = %d\n',digit_num,G0_samp_num);

%% Generate the G0 map
train_sel_array = zeros(G0_samp_num,G0_digit_num);
train_label_array = ones(G0_samp_num,G0_digit_num);

%Randomly sample images with these labels
for i = 1:digit_num
    digit = sel_digit(i);
    idx = find(test_y(:,digit) == 1);
    perm_idx = randperm(numel(idx));
    samp_idx = perm_idx(1:G0_samp_num);
    train_sel_array(:,i) = idx(samp_idx);
    train_label_array(:,i) = (digit-1)*train_label_array(:,i);
end

% Generate all pair-wise similarity label by SVM
graph_sz = G0_samp_num*G0_digit_num;
train_sel_vect = train_sel_array(:);
train_label_vect = train_label_array(:);
train_mat = zeros(graph_sz*(graph_sz-1)/2,784);
train_label = zeros(graph_sz*(graph_sz-1)/2,1);
cnt = 1;
for i = 1:graph_sz-1
    for j = i+1:graph_sz
        train_mat(cnt,:) = abs(train_x(train_sel_vect(i),:)-train_x(train_sel_vect(j),:))/255;
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
% for i = 1:10
    fprintf('iteration i = %d\n',i);
    test_data = zeros(graph_sz-i,784);
    test_label = zeros(graph_sz-i,1);
    for j = i+1:graph_sz
        diff_vect = abs(test_x(test_vect(i),:)-test_x(test_vect(j),:));
        test_data(j-i,:) = double(diff_vect)/255;
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
 
%% Backup code
%% Random sampling from MNIST
% %Randomly pick some digits
% sel_digit = randperm(10,digit_num);
% train_samp_num = 60;
% train_array = zeros(train_samp_num,digit_num);
% train_label = ones(train_samp_num,digit_num);
% % Randomly sample images with these labels
% for i = 1:digit_num
%     digit = sel_digit(i);
%     idx = find(train_y(:,digit) == 1);
%     train_samp_idx = randperm(numel(idx),train_samp_num);
%     train_array(:,i) = idx(train_samp_idx);
%     train_label(:,i) = (digit-1)*train_label(:,i);%Change the sel_digit from indices to integer label
% end

%% Generate data for pair-wise SVM training
% Generate negative pairs
% neg_num = 4000;
% neg_pair = zeros(neg_num,2);
% for j = 1:neg_num
%    pair = randperm(digit_num,2);
%    id1 = train_array(ceil(rand(1)*train_samp_num),pair(1));
%    id2 = train_array(ceil(rand(1)*train_samp_num),pair(2));
%    neg_pair(j,:) = [id1,id2];
% end
% 
% % Generate positive pairs
% pos_num = 4000;
% pos_pair = zeros(pos_num,2);
% for j = 1:pos_num
%    cluster_id = randperm(digit_num,1);
%    ids = train_array(randperm(train_samp_num,2),cluster_id);
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
% 
% diff_fullmat = abs(full_mat1-full_mat2);
% train_mat = double(diff_fullmat)/255;
% train_label = double((full_label1 == full_label2));
    
