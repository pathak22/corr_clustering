%Compute OPT and Hanning distance for greedy algorithm results against
%Ground truth
%Read first two lines
str = '../data/solution_leaderModel.txt';
inputfile = fopen(str);
fgetl(inputfile);
fgetl(inputfile);

%Read the rest into cell
C = textscan(inputfile,'%d %d %d');
fclose(inputfile);

%% Compute the Hamming distance between our clustering and ground truth
% Ground truth connection map
num_nodes = numel(C{1});
Gt = zeros(num_nodes);
for i = 1:num_nodes
    gt_id = C{2}(i);
    idx = (C{2}==gt_id);
    Gt(i,idx) = 1;
end   

% Obtained clustering connection map
oCl = zeros(num_nodes);
for i = 1:num_nodes
    oCl_id = C{3}(i);
    idx = (C{3}==oCl_id);
    oCl(i,idx) = 1;
end

% Compute Hanning distance
H_C = xor(oCl,Gt);
H_dist = sum(sum(H_C));

%% Compute the cost for our clustering and ground truth
str_data = '../data/data_leaderModel.txt';

%Open file and skip the first two lines of comments and the third lines
fileName = str_data;
inputfile = fopen(str_data);
fgetl(inputfile);
fgetl(inputfile);
fgetl(inputfile);

% Read in the rest
tline = fgetl(inputfile);
i=1;
C0=cell(num_nodes,1);
while ischar(tline)
    C0{i} = str2num(tline);
    tline = fgetl(inputfile);
    i = i+1;
end

% Construct w_ij from G0
adj_m = zeros(num_nodes);
for i = 1:num_nodes
    conn_i = C0{i}(2:end);
    adj_m(i,i)=0;%pij = 1 for node itself
    for j = i+1:num_nodes
        if any(conn_i == j)
            adj_m(i,j) = 1;%1-p for existing edge in G0
        else
            adj_m(i,j) = -1;%p for non-existing edge in G0
        end        
    end
end
        
w = adj_m+adj_m';

% Compute the cost
Gt_cost = sum(sum(Gt.*(1-w)+(1-Gt).*(1+w)));
oCl_cost = sum(sum(oCl.*(1-w)+(1-oCl).*(1+w)));
