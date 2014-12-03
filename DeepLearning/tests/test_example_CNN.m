% function test_example_CNN
addpath(genpath('/Users/wckuo/Documents/MATLAB/DeepLearning'));
load mnist_uint8;

n = 0;% n = 3 err = 0.0515, n = 0 err = 0.1113
% train_x = mnist_augment(train_x,n);
% train_y = repmat(train_y,n+1,1);

train_x = double(reshape(train_x',28,28,60000*(n+1)))/255;
train_y = double(train_y');
test_x = double(reshape(test_x',28,28,10000))/255;
test_y = double(test_y');

%% ex1 Train a 6c-2s-12c-2s Convolutional neural network 
%will run 1 epoch in about 200 second and get around 11% error. 
%With 100 epochs you'll get around 1.2% error

rand('state',0)

cnn.layers = {
    struct('type', 'i') %input layer
    struct('type', 'c', 'outputmaps', 6, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %sub sampling layer
    struct('type', 'c', 'outputmaps', 12, 'kernelsize', 5) %convolution layer
    struct('type', 's', 'scale', 2) %subsampling layer
};


opts.alpha = 1;
opts.batchsize = 50;
opts.numepochs = 1;

cnn = cnnsetup(cnn, train_x, train_y);
fprintf('start training');
cnn = cnntrain(cnn, train_x, train_y, opts);

[er, bad] = cnntest(cnn, test_x, test_y);

%plot mean squared error
figure; plot(cnn.rL);
assert(er<0.12, 'Too big error');
