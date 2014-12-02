% Plot comparison of plot 

addpath('./PlotPub-1.3/lib');

dist_sdp=[19.2 16.8 16.8 17.3 12.6 NaN];
sdp_cost=[231.86 1147.22 2909.17 5605.72 8574.68 NaN];

dist_lp=[25.6 22 23.36 24.94 25.34 NaN];
lp_cost=[191 834 1969 3664 5505 NaN];

load('greedy_sdp_vars.mat');

figure;
plot(niList*k,voteCl_cost,niList*k,oCl_cost,niList*k,leaderCl_cost,niList*k,Gt_cost,niList*k,sdp_cost,niList*k,lp_cost);
title('Objective v/s N at K:20, p:0.4, \epsilon:0.1','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Objective Cost';
opt.XLabel = 'Number of nodes, N';
opt.FileName = 'objVsN_sdp.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our','Ground-Truth','SDP','LP'};
opt.LegendLoc = 'NorthWest';
setPlotProp(opt);

figure;
plot(niList*k,dist_voteCl,niList*k,dist_oCl,niList*k,dist_leaderCl,niList*k,dist_sdp,niList*k,dist_lp);
title('Hamming Distance v/s N at K:20, p:0.4, \epsilon:0.1','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Hamming Distance';
opt.XLabel = 'Number of nodes, N';
opt.FileName = 'hammingVsN_sdp.png';
opt.Legend = {'Vote-Greedy','Pivot', 'Our','SDP','LP'};
opt.LegendLoc = 'NorthWest';
setPlotProp(opt);


% ---------------------------------------------------

load('mnist_exp_result.mat');
figure;
plot(thr,Vote_hamdist,thr,Pivot_hamdist,thr,Density_hamdist,thr,G0_hamdist);
title('Result for held-out 5 categories','FontSize',18,'FontWeight','normal');
opt = [];
opt.YLabel = 'Hamming Distance';
opt.XLabel = 'Threshold value';
opt.FileName = 'hamming_mnist.png';
opt.Legend = {'Vote-Greedy','Pivot','Our','5-Category Classifier'};
opt.LegendLoc = 'NorthWest';
setPlotProp(opt);