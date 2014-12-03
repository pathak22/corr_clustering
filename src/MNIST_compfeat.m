addpath('../MNIST');
load('mnist_complete.mat');
str = cell(13,1);
str{1} = 'Area';
str{2} = 'Centroid';
str{3} = 'MajorAxisLength';
str{4} = 'MinorAxisLength';
str{5} = 'Eccentricity';
str{6} = 'Orientation';
str{7} = 'ConvexArea';
str{8} = 'FilledArea';
str{9} = 'EulerNumber';
str{10} = 'EquivDiameter';
str{11} = 'Solidity';
str{12} = 'Extent';
str{13} = 'Perimeter';
str{14} = 'BoundingBox';

train_feat = [];
test_feat = [];
for i = 1:60000
    fprintf('iter:%d\n',i);
    ximg = reshape(train_x(i,:)',28,28)';
    xbin = ximg>mean(ximg(:));
    xprops = regionprops(xbin,'all');
    n = 1;
    while numel(xprops)>1
        se = strel('disk',n);
        xprops = regionprops(imclose(xbin,se),'all');
        n = n+1;
    end
          
    xfeat = [];
    for j = 1:13
        xstr = str{j};
        xfeat = cat(2,xfeat,xprops.(xstr));
    end
    train_feat = cat(1,train_feat,xfeat);
end

for i = 1:10000
    fprintf('iter:%d\n',i);
    ximg = reshape(test_x(i,:)',28,28)';
    xbin = ximg>mean(ximg(:));
    xprops = regionprops(xbin,'all');
    n = 1;
    while numel(xprops)>1
        se = strel('disk',n);
        xprops = regionprops(imclose(xbin,se),'all');
        n = n+1;
    end
          
    xfeat = [];
    for j = 1:13
        xstr = str{j};
        xfeat = cat(2,xfeat,xprops.(xstr));
    end
    test_feat = cat(1,test_feat,xfeat);
end

save('mnist_region_feat','train_feat','test_feat','-v7.3');