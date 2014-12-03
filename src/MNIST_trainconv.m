%% Training CNN on MNIST digits
% clear all
% addpath(genpath('../DeepLearning'));
% addpath('../MNIST');
% load('mnist_complete.mat');
% 
% train_x = double(reshape(train_x',28,28,60000))/255;
% train_y = double(train_y');
% 
% label_y = zeros(1,60000);
% for i = 1:60000
%     label_y(i) = find(train_y(:,i));
% end
% idx_half = label_y<6;
% train_yh = train_y(:,idx_half);
% train_xh = train_x(:,:,idx_half);
% train_xh = train_xh(:,:,1:30000);
% train_yh = train_yh(:,1:30000);
% 
% rand('state',0)
% 
% cnn.layers = {
%     struct('type', 'i') %input layer
%     struct('type', 'c', 'outputmaps', 6, 'kernelsize', 5) %convolution layer
%     struct('type', 's', 'scale', 2) %sub sampling layer
%     struct('type', 'c', 'outputmaps', 12, 'kernelsize', 5) %convolution layer
%     struct('type', 's', 'scale', 2) %subsampling layer
% };
% 
% opts.alpha = 1;
% opts.batchsize = 50;
% opts.numepochs = 20;
% 
% cnn = cnnsetup(cnn, train_xh, train_yh);
% fprintf('start training');
% cnn = cnntrain(cnn, train_xh, train_yh, opts);
% 
% save('cnn_mnist_half.mat','cnn','-v7.3');
load('cnn_mnist_half.mat');

%% Find best threshold for L2 distance
num = 6000;
digit_num = 5;
% forward pass to get features
fv = cnntest_fv(cnn,train_x(:,:,1:num));
fv = fv';
fdist = real(pdist2(fv,fv,'euclidean'));
c = 0.1:0.1:2;
Dist = zeros(1,numel(c));
label_v = zeros(1,num);
Gt = zeros(num);
dv = ones(1,num);
diagG = diag(dv);

%construct Gt
for i = 1:num
    label_v(i) = find(train_y(:,i));
end

for i = 1:num
    for j = 1:num
        Gt(i,j) = (label_v(i)==label_v(j));
    end
end

%construct G0 and store distance
for i = 1:numel(c)
    thr = prctile(fdist(:),100/digit_num*c(i));
    G0 = fdist<thr;
    G0(diagG>0) = 1;
    H0 = xor(G0,Gt);
    Dist(i) = sum(H0(:));
end

[~,idx]=min(Dist);
display(c(idx));

%10 digits num = 6000 => c=0.5;
%5 digits num = 6000 => c=0.2




