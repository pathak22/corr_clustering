% This is the function to test SVM model
function testout=test_pairwise(test_x,test_y,sel_digit,model,opt)
digit_num = numel(sel_digit);
if strcmp(opt,'rand')
    %% Random sampling from MNIST test set
    % Find max training set digit sample number
    
    train_samp_num = Inf;
    for i = 1:digit_num
        digit = sel_digit(i);
        idx = find(test_y(:,digit) == 1);
        train_samp_num = min(train_samp_num,numel(idx));    
    end

    train_array = zeros(train_samp_num,digit_num);
    train_label = ones(train_samp_num,digit_num);
    % Randomly sample images with these labels
    for i = 1:digit_num
        digit = sel_digit(i);
        idx = find(test_y(:,digit) == 1);
        train_samp_idx = randperm(numel(idx),train_samp_num);
        train_array(:,i) = idx(train_samp_idx);
        train_label(:,i) = (digit-1)*train_label(:,i);%Change the sel_digit from indices to integer label
    end

    %% Generate data for pair-wise SVM training
    % Generate negative pairs
    neg_num = 1000;
    neg_pair = zeros(neg_num,2);
    for j = 1:neg_num
       pair = randperm(digit_num,2);
       id1 = train_array(ceil(rand(1)*train_samp_num),pair(1));
       id2 = train_array(ceil(rand(1)*train_samp_num),pair(2));
       neg_pair(j,:) = [id1,id2];
    end

    % Generate positive pairs
    pos_num = 1000;
    pos_pair = zeros(pos_num,2);
    for j = 1:pos_num
       cluster_id = randperm(digit_num,1);
       ids = train_array(randperm(train_samp_num,2),cluster_id);
       pos_pair(j,:) = ids';
    end

    % Construct the data matrix
    pair_mat1 = cat(1,neg_pair(:,1),pos_pair(:,1));
    pair_mat2 = cat(1,neg_pair(:,2),pos_pair(:,2));
    pair_num = size(pair_mat1,1);
    full_mat1 = zeros(pair_num,784);
    full_mat2 = zeros(pair_num,784);
    full_label1 = zeros(pair_num,1);
    full_label2 = zeros(pair_num,1);
    for i = 1:pair_num
        full_mat1(i,:) = test_x(pair_mat1(i),:);
        full_mat2(i,:) = test_x(pair_mat2(i),:);
        full_label1(i) = find(test_y(pair_mat1(i),:),1,'first')-1;%find 1 index in each row
        full_label2(i) = find(test_y(pair_mat2(i),:),1,'first')-1;%find 1 index in each row
    end

    diff_fullmat = abs(full_mat1-full_mat2);
    test_mat = double(diff_fullmat)/255;
    test_label = double((full_label1 == full_label2));
    fprintf('start svm testing...\n');
    testout = svmpredict(test_label, test_mat, model);
    
end
if strcmp(opt,'full')
    G0_digit_num = digit_num;
    G0_samp_num = 10;
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
        [predict_label, accuracy, ~] = svmpredict(test_label, test_data, model);
        G0(i,i+1:end) = predict_label';
        accum_accuracy = accum_accuracy + accuracy*(graph_sz-i); 
    %     clear test_data test_label
    end
    sum_accuracy = accum_accuracy/(graph_sz*(graph_sz-1)/2);
    display(sum_accuracy(1));
    testout = 0;
end
