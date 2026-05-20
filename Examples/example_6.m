%In this example we analyse low resolution downsampled data set
%generated with the linear Ornstein-Uhlenbeck (OU) model
%Here we apply a spline method, but assume additive noise.
%see for details example_5.m
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

S = load(fullfile(rootDir, 'OUdata1D.mat'));
data = S.data; %load the data
%We downsample to illustrate the performance for low resolution data
data = data(1:100:end); % Only every 100 data points are considered
L = -2.5;
R = 2.5; %Since we have a  spline  model it is better to shrink the state space
data(data < L | data > R) = nan; %This is VERY important: In spline modeling if you consider a smaller range for your data then you must assign 'nan' to those few data points falling outside this range.
dt = 1;
nmu = 8; %number of knots for mu
nsigma = 1; %number of knots for sigma
result12 = euler_reconstruction(data, dt, 'nKnots', [nmu nsigma], 'spline', 'QQ', 'L', L, 'R', R, ...,
    'lb', [zeros(1, nmu) - 10, zeros(1, nsigma) + eps], 'ub', zeros(1, nmu + nsigma) + 10, 'solver', 'fmincon', 'search_agents', 5);
disp('Searching legitimate points')
legpoints = legitimate_points(data, dt, 'prev', result12, 'prev_range', 0.5, 'j', 3, 'k', 9);
result_her12 = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', 'fmincon', 'search_agents', 5);

mu = @(x,par)par(1).*x;
sigma = @(x,par)par(2)+0 .*x; %this is true model
par = zeros(2, 1);
par(1) = -1;
par(2) = 1; %true model parameters
xplot = linspace(L, R, 2000); % a dense mesh across the considered range
plot_results(result12, xplot, mu(xplot, par), sigma(xplot, par)); % Euler & true models
plot_results(result_her12, xplot, mu(xplot, par), sigma(xplot, par)); %Hermite & true models
