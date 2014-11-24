%data.txt format: 
%leader(bool) cluster_id(int) vertex_id(0-N-1) list of neighbors
%two lines of comments on top

% addpath('path/to/cvx/');
addpath(genpath('/Users/wckuo/Documents/MATLAB/cvx'));
str = '../data/data_leaderModel.txt';

%Read the third line
[num_nodes, p, e] = textread(str,'%d %f %f', 3);

%Open file and skip the first two lines of comments and the third lines
fileName = str;
inputfile = fopen(str);
fgetl(inputfile);
fgetl(inputfile);
fgetl(inputfile);

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
    conn_i = C{i}(2:end);
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
        
w = adj_m+adj_m';
tic;
%% Solving SDP by CVX
% create and solve the problem
cvx_begin 
  % A is a PSD symmetric matrix (n-by-n)
  variable X(num_nodes,num_nodes) semidefinite;

  % constrained matrix entries.
  diag(X) == ones(num_nodes,1);
  X >= 0;
  for j = 1:num_nodes
      repmat(X(:,j),[1 num_nodes])+repmat(X(j,:),[num_nodes 1])<=ones(num_nodes,num_nodes)+X;
  end

  % find the solution to the problem
  minimize( sum(sum(X.*(1-w)+(1-X).*(1+w))))

cvx_end


%% Clustering algorithm
max_num_cluster = 20;
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
    gt_id = C{i}(1);
    for j = 1:num_nodes
        node_id = C{j}(1);
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

fprintf('SDP Hanning distance is %d\n',H_dist);

% celldisp(cluster);





