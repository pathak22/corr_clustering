% addpath('path/to/cvx/');
addpath(genpath('/Users/wckuo/Documents/MATLAB/cvx'));
str = '../data/data.txt';

%Read first line
[num_nodes, p, q, e] = textread(str, ...
'%d %f %f %f', 1);
%p: i=j, q: i!=j, E(i,j)=1, e: i!=j, E(i,j)=0

%Open file
fileName = str;
inputfile = fopen(str);

% Skip the first line
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
    id = C{i}(1);
    for j = i+1:num_nodes
        scan_id = C{j}(1);
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

  % find the solution to the problem
  minimize( sum(sum(X.*(1-w)+(1-X).*(1+w))))

cvx_end


%% Clustering algorithm
cluster = cell(20,1);
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

celldisp(cluster);





