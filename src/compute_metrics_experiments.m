function compute_metrics_experiments

% Download this package: http://masumhabib.com/publication-quality-graphs-matlab/
addpath('./PlotPub-1.3/lib');

niList = 10:10:100; 				% Number of nodes per cluster
pList = 0.0:0.05:0.45;   			% Intra-non-leader flipping probability
epsilonList = 0.0:0.05:0.4;         % Leader-neighbor flipping probability
        
% -----------------------------

k=20; p=0.4; epsilon=0.1; 
count=1;
dist_oCl=[]; dist_leaderCl=[]; dist_voteCl=[]; Gt_cost=[];
oCl_cost=[]; leaderCl_cost=[]; voteCl_cost=[];
for ni=niList
    str_data = sprintf('../data/experiments/data_leaderModel_k%d_ni%d_p%d_epsilon%d.txt',k,ni,single(p*100),single(epsilon*100));
    str_solution = sprintf('../data/experiments/solution_leaderModel_k%d_ni%d_p%d_epsilon%d.txt',k,ni,single(p*100),single(epsilon*100));
    [dist_oCl(count),dist_leaderCl(count),dist_voteCl(count),Gt_cost(count), ...
        oCl_cost(count),leaderCl_cost(count),voteCl_cost(count)] = ...
        compute_metrics(str_solution,str_data);
    count = count + 1;
end
figure;
plot(niList*k,voteCl_cost,niList*k,oCl_cost,niList*k,leaderCl_cost,niList*k,Gt_cost);
title('Objective v/s N at K:20, p:0.4, \epsilon:0.1','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Objective Cost';
opt.XLabel = 'Number of nodes, N';
opt.FileName = 'objVsN.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our','Ground-Truth'};
opt.LegendLoc = 'NorthWest';
setPlotProp(opt);

figure;
plot(niList*k,dist_voteCl,niList*k,dist_oCl,niList*k,dist_leaderCl);
title('Hamming Distance v/s N at K:20, p:0.4, \epsilon:0.1','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Hamming Distance';
opt.XLabel = 'Number of nodes, N';
opt.FileName = 'hammingVsN.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our'};
setPlotProp(opt);

% -----------------------------

k=20; ni=50; epsilon=0.1;
count=1;
dist_oCl=[]; dist_leaderCl=[]; dist_voteCl=[]; Gt_cost=[];
oCl_cost=[]; leaderCl_cost=[]; voteCl_cost=[];
for p=pList
    str_data = sprintf('../data/experiments/data_leaderModel_k%d_ni%d_p%d_epsilon%d.txt',k,ni,single(p*100),single(epsilon*100));
    str_solution = sprintf('../data/experiments/solution_leaderModel_k%d_ni%d_p%d_epsilon%d.txt',k,ni,single(p*100),single(epsilon*100));
    [dist_oCl(count),dist_leaderCl(count),dist_voteCl(count),Gt_cost(count), ...
        oCl_cost(count),leaderCl_cost(count),voteCl_cost(count)] = ...
        compute_metrics(str_solution,str_data);
    count = count + 1;
end
figure;
plot(pList,voteCl_cost,pList,oCl_cost,pList,leaderCl_cost,pList,Gt_cost);
title('Objective v/s p at K:20, N:1000, \epsilon:0.1','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Objective Cost';
opt.XLabel = 'Probability of Non-leader-edge flips, p';
opt.FileName = 'objVsP.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our','Ground-Truth'};
opt.LegendLoc = 'NorthWest';
setPlotProp(opt);

figure;
plot(pList,dist_voteCl,pList,dist_oCl,pList,dist_leaderCl);
title('Hamming Distance v/s p at K:20, N:1000, \epsilon:0.1','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Hamming Distance';
opt.XLabel = 'Probability of Non-leader-edge flips, p';
opt.FileName = 'hammingVsP.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our'};
setPlotProp(opt);

% -----------------------------

k=20; ni=50; p=0.4;
count=1;
dist_oCl=[]; dist_leaderCl=[]; dist_voteCl=[]; Gt_cost=[];
oCl_cost=[]; leaderCl_cost=[]; voteCl_cost=[];
for epsilon=epsilonList
    str_data = sprintf('../data/experiments/data_leaderModel_k%d_ni%d_p%d_epsilon%d.txt',k,ni,single(p*100),single(epsilon*100));
    str_solution = sprintf('../data/experiments/solution_leaderModel_k%d_ni%d_p%d_epsilon%d.txt',k,ni,single(p*100),single(epsilon*100));
    [dist_oCl(count),dist_leaderCl(count),dist_voteCl(count),Gt_cost(count), ...
        oCl_cost(count),leaderCl_cost(count),voteCl_cost(count)] = ...
        compute_metrics(str_solution,str_data);
    count = count + 1;
end
figure;
plot(epsilonList,voteCl_cost,epsilonList,oCl_cost,epsilonList,leaderCl_cost,epsilonList,Gt_cost);
title('Objective v/s \epsilon at K:20, N:1000, p:0.4','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Objective Cost';
opt.XLabel = 'Probability of Leader-edge flips, \epsilon';
opt.FileName = 'objVsEpsilon.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our','Ground-Truth'};
opt.LegendLoc = 'NorthWest';
setPlotProp(opt);

figure;
plot(epsilonList,dist_voteCl,epsilonList,dist_oCl,epsilonList,dist_leaderCl);
title('Hamming Distance v/s \epsilon at K:20, N:1000, p:0.4','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Hamming Distance';
opt.XLabel = 'Probability of Leader-edge flips, \epsilon';
opt.FileName = 'hammingVsEpsilon.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our'};
setPlotProp(opt);

% -----------------------------

end

function [dist_oCl,dist_leaderCl,dist_voteCl,Gt_cost,oCl_cost,leaderCl_cost,voteCl_cost] = compute_metrics(str_solution,str_data)

% Compute OPT and Hanning distance for greedy algorithm results against
% Ground truth
% Read first two lines
inputfile = fopen(str_solution);
fgetl(inputfile);
fgetl(inputfile);

% Read the rest into cell
C = textscan(inputfile,'%d %d %d %d %d');
fclose(inputfile);

% Compute the Hamming distance between our clustering and ground truth
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

% Obtained leader model clustering connection map
leaderCl = zeros(num_nodes);
for i = 1:num_nodes
    leaderCl_id = C{4}(i);
    idx = (C{4}==leaderCl_id);
    leaderCl(i,idx) = 1;
end

% Obtained vote clustering connection map
voteCl = zeros(num_nodes);
for i = 1:num_nodes
    voteCl_id = C{5}(i);
    idx = (C{5}==voteCl_id);
    voteCl(i,idx) = 1;
end

% Compute Hamming distance
H_oCl = xor(oCl,Gt);
dist_oCl = sum(sum(H_oCl))/(num_nodes^2)*100;
H_leader = xor(leaderCl,Gt);
dist_leaderCl = sum(sum(H_leader))/(num_nodes^2)*100;
H_vote = xor(voteCl,Gt);
dist_voteCl = sum(sum(H_vote))/(num_nodes^2)*100;


% Compute the cost for our clustering results and ground truth
% Open file and skip the first two lines of comments and the third line
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
leaderCl_cost = sum(sum(leaderCl.*(1-w)+(1-leaderCl).*(1+w)));
voteCl_cost = sum(sum(voteCl.*(1-w)+(1-voteCl).*(1+w)));

end
