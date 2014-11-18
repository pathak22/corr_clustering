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

% Construct p_ij table for G_0
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
v = diag(w);
Dw = diag(v);
w = w - Dw/2;

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

