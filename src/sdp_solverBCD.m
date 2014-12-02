%data.txt format:
%leader(bool) cluster_id(int) vertex_id(0-N-1) list of neighbors
%two lines of comments on top

% addpath('path/to/cvx/');
%addpath(genpath('/Users/wckuo/Documents/MATLAB/cvx'));
clear all;
str = '../data/data_leaderModel.txt';

%Open file and skip the first two lines of comments and the third lines
fileName = str;
inputfile = fopen(str);
fgetl(inputfile);
fgetl(inputfile);
data_info = fgetl(inputfile);
v_info = str2num(data_info);
num_nodes = v_info(1);
p = v_info(2);
epsilon = v_info(3);

% Read in the rest
tline = fgetl(inputfile);
i=1;
C=cell(num_nodes,1);
while ischar(tline)
    C{i} = str2num(tline);
    tline = fgetl(inputfile);
    i = i+1;
end

% Read the geometric random graph
% M = [p,e;e,p];%GRG
% M = zeros();%no knowledge about GRG

% Construct w_ij from G0
adj_m = zeros(num_nodes);
for i = 1:num_nodes
    conn_i = C{i}(4:end);
    adj_m(i,i)=0;%pij = 1 for node itself
    %     id = C{i}(1);
    for j = i+1:num_nodes
        %         scan_id = C{j}(1);
        if any(conn_i == j)
            adj_m(i,j) = 1;%1-p for existing edge in G0
        else
            adj_m(i,j) = -1;%p for non-existing edge in G0
        end
    end
end

W = adj_m+adj_m';
tic;

%% Solving SDP by Block coordinate Descent
% Variable : X
% Weight Matrix : W

% Initialization 
X = eye(num_nodes);
XHat = eye(num_nodes -1);

profile on
% Alternating optimization
for i = 1:2*num_nodes%num_nodes
    
    % Compute Permuted Matrix
    XHat(1: i-1, 1:i-1) = X(1: i-1, 1:i-1);
    XHat(:,i) = X(1:num_nodes-1, num_nodes);
    XHat(i,:) = X(num_nodes,1:num_nodes-1);
    XHat(i+1:num_nodes-1, i+1:num_nodes-1) = X(i+1:num_nodes-1, i+1:num_nodes-1);
    XHat(i,i) = 1;
    
    % Compute corresponding weight vector
    wBar = W(setdiff(1:num_nodes,i),i);
    
    % Optimizing CVX
    cvx_begin sdp
    
    variable u(num_nodes-1)
    minimize (u'*(1-wBar)+(1-u)'*wBar)
    
    subject to
    u >= 0                  % Elementwise positivity
    [XHat u; u' 1] >= 0     % Semidefinite Constraint
    
   for j = 1:num_nodes
        for k = 1:num_nodes
            if (j<i && k<i)
                u(j)-u(k) >= X(j,k)-1
                u(j)+u(k) <= X(j,k)+1
            elseif (j>i && k<i)
                u(j-1)-u(k) >= X(j,k)-1
                u(j-1)+u(k) <= X(j,k)+1
            elseif (j<i && k>i)
                u(j)-u(k-1) >= X(j,k)-1
                u(j)+u(k-1) <= X(j,k)+1
            elseif (j>i && k>i)
                u(j-1)-u(k-1) >= X(j,k)-1
                u(j-1)+u(k-1) <= X(j,k)+1
            end
        end
   end 
    cvx_end
    
    for j = 1:num_nodes
       if j < i
           X(i,j) = u(j);
           X(j,i) = u(j);
       elseif j > i
           X(i,j) = u(j-1);
           X(j,i) = u(j-1);
       end
    end
    fprintf('Overall Onjective Value: %d', sum(sum(X.*(1-W)+(1-X).*(1+W))));   
end

profile off


%% Clustering algorithm
max_num_cluster = num_nodes;
cluster = cell(max_num_cluster,1);
list = 1:num_nodes;
cnt = 0;
cl_id = 1;  %Cluster id

while any(list)
    v=list(list~=0);                                %   Find unchosen nodes
    perm_list = randperm(num_nodes-cnt);            %   Randomly permute the list
    sel_node = v(perm_list(1));                     %   Pick the first node from permuted list
    cluster{cl_id} = cat(1,cluster{cl_id},sel_node);%   Add that node to the cluster
    list(list==sel_node)=0;                         %   Remove node from the list
    cnt = cnt+1;                                    %   Update the count of chosen nodes
    v=list(list~=0);
    for i = 1:num_nodes-cnt
        vi = v(i);                                  %   Scan every unchosen node
        bintest = binornd(1,X(sel_node,vi));        %   Binormial sampling to include in the cluster
        if bintest
            cluster{cl_id} = cat(1,cluster{cl_id},vi); %    Add that node to the cluster
            list(list==vi)=0;                       %   Remove the node from list
            cnt = cnt+1;                            %   Update count
        end
    end
    cl_id = cl_id + 1;                              %   Update cluster id
end
toc;

%% Compute the Hamming distance between our clustering and ground truth
% Ground truth connection map
Gt = zeros(num_nodes);
for i = 1:num_nodes
    gt_id = C{i}(2);
    for j = 1:num_nodes
        node_id = C{j}(2);
        if (gt_id == node_id)
            Gt(i,j) = 1;
        end
    end
end

%Construct table for membership of nodes in my clustering
membership = zeros(2,num_nodes);
membership(1,:) = 1:num_nodes;
for i = 1:max_num_cluster
    list = cluster{i};
    membership(2,list) = i;
end


% My clustering connection map
myCl = zeros(num_nodes);
for i = 1:num_nodes
    myCl_id = membership(2,i);
    for j = 1:num_nodes
        node_id = membership(2,j);
        if (myCl_id == node_id)
            myCl(i,j) = 1;
        end
    end
end

% Compute Hanning distance
H_C = xor(myCl,Gt);
H_dist = sum(sum(H_C));

fprintf('SDP Hamming distance is %d\n',H_dist);
Gt_cost = sum(sum(Gt.*(1-W)+(1-Gt).*(1+W)));
fprintf('Ground truth cost is %d\n',Gt_cost);
% celldisp(cluster);





