function data = mnist_augment(data,n)
    
%     load mnist_uint8;
%     n = 5;
%     id = 100;
%     data = double(train_x(id,:))/255;
%     disp(find(train_y(id,:)>0)-1);
    
    if (n == 0)
        return;
    end
    
    Ndata = size(data,1);
    for l = 1:n
        % Uniform random fields
        xrand = -1+2*rand(28);
        yrand = -1+2*rand(28);
%         fprintf('\n augment cycle %d',l);
        %Smooth out the random fields
        s = 3;
        h = fspecial('gaussian',6*s,s);
        xrand = imfilter(xrand,h,'conv','same');
        yrand = imfilter(yrand,h,'conv','same');
        
        for i = 1:Ndata
        
        if mod(i,5000) == 0
            fprintf('\n i=%d',i);
        end
          
        img = reshape(data(i,:),28,28);
        img = double(img)';

        % Apply image displacement field by bilinear interpolation
        xcoord = repmat(1:28,28,1)+xrand;
        ycoord = repmat((1:28)',1,28)+yrand;
        
        img_mask = img>0;
        se = strel('disk',2);
        img_mask = imdilate(img_mask,se);
        vimg = (img_mask>0);
        xd = zeros(28);
        yd = zeros(28);
        xd(vimg) = xcoord(vimg)-floor(xcoord(vimg));
        yd(vimg) = ycoord(vimg)-floor(ycoord(vimg));
        xcoord(xcoord<1) = 1;
        ycoord(ycoord<1) = 1;
        xcoord(xcoord>28) = 28;
        ycoord(ycoord>28) = 28;
        idx1 = sub2ind([28 28],floor(xcoord),floor(ycoord));
        p1 = img(idx1);
        idx2 = sub2ind([28 28],ceil(xcoord),floor(ycoord));
        p2 = img(idx2);
        idx3 = sub2ind([28 28],floor(xcoord),ceil(ycoord));
        p3 = img(idx3);
        idx4 = sub2ind([28 28],ceil(xcoord),ceil(ycoord));
        p4 = img(idx4);
        fX1 = xd.*p2+(1-xd).*p1;
        fX2 = xd.*p4+(1-xd).*p3;
        els_img = yd.*fX2+(1-yd).*fX1;       

        %Generate random field agumented data base
        data = cat(1,data,els_img(:)');

        end
    end
end
% figure(1);imagesc(reshape(data(1,:),28,28)');colormap(gray);
% figure(2);imagesc(reshape(data(3,:),28,28)');colormap(gray);
