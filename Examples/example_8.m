thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

S = load(fullfile(rootDir, 'MayData1D.mat'));
data = S.data;
%We downsample to illustrate the performance for low resolution data
data = data(1:300:end);
L = 0;
R = max(data); %we consider the entire range of data
data(data < L | data > R) = nan;
dt = 3;
nmu = 8;
nsigma = 8;
result14 = euler_reconstruction(data, dt, 'nKnots', [nmu nsigma], 'spline', 'QQ', 'L', L, 'R', R, ...,
    'lb', [zeros(1, nmu) - 10, zeros(1, nsigma) + eps], 'ub', zeros(1, nmu + nsigma) + 10, 'solver', 'fmincon', 'search_agents', 5);
legpoints = legitimate_points(data, dt, 'prev', result14, 'prev_range', 0.5, 'j', 3, 'k', 9);
result_her14 = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', ...,
    'fmincon');
realmu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2 ./(x.^2+par(4).^2);
realsigma = @(x,par)par(5)+0 .*x; % true model
realpar = zeros(5, 1);
realpar(1) = 1.01;
realpar(2) = 10;
realpar(3) = 2.75;
realpar(4) = 1.6;
realpar(5) = 0.4; %true model parameters
xplot = linspace(L, R, 2000); % a dense mesh across the considered range
plot_results(result14, xplot, realmu(xplot, realpar), realsigma(xplot, realpar)); % Euler & true models
plot_results(result_her14, xplot, realmu(xplot, realpar), realsigma(xplot, realpar)); %Hermite & true models
