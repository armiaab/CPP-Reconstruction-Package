%In this example we analyse high resolution data set
%generated with a nonlinear model on high resolution data:
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

fprintf('The data set is full MayData1D.mat\n')
S = load(fullfile(rootDir, 'MayData1D.mat'));
data = S.data; %load the data
data = data(1:100000); %We only use the first third of the dataset
%%Data resolution

%%Augmented Dickey-Fuller test for data stationarity
[~, ~, ~, ~, reg] = adftest(data, 'model', 'ARD', 'lags', 0:20); % the input 0:20 is the number of lags we try in fitting an autoregressive model to data.
[~, lagndx] = min([reg(:).BIC]); % this tells us how many lags we need (lagndx is 1)
[h, Pvalue, ~] = adftest(data, 'model', 'ARD', 'lags', lagndx);
fprintf('Augmented Dickey-Fuller test for data stationarity\nOptimum lag=%d, h = %g, p = %g\n', lagndx, h, Pvalue);

%%Markovitiy Markov-Einstein timescale
order = 10; % order is the maximum AR order being considered (should be long enough)
%Markov-Einstein time scale
AR = ME_TimeScale(data, order);
AR1 = AR(AR > 0.1);
fprintf('AR model:%s\nThe Markov-Einstein scale is ca. %d\n', mat2str(AR, 5), numel(AR1));

%%data resolution
R = RelaxationTime(data);

fprintf('Data resolution\nRelaxation time = %g\n', R);
if R > 100
    disp('High resolution data')
elseif R > 50
    disp('Medium resolution data')
elseif R > 1
    disp('Low resolution data')
else
    disp('Extremely low resolution data');
end

%%parametric Euler reconstruction
L = min(data);
R = max(data);
dt = 0.01;
mu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2 ./(x.^2+par(4).^2);
sigma = @(x,par)par(5);
result4 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun', eulergrad(mu, sigma), ...
    'lb', zeros(1, 5) + eps, 'ub', 15 .* ones(1, 5), 'useparallel', true, 'search_agents', 5);
r = 1.01;
K = 10;
g = 2.75;
a = 1.6;
s = 0.4; % true parameter values
realpar = [r K g a s];
realmu = @(x,par)r.*x.*(1-x./K)-g.*x.^2 ./(x.^2+a.^2);
realsigma = @(x,par)s; % true model
xplot = linspace(L, R, 2000);

%%fitting a spline model
plot_results(result4, xplot, realmu(xplot, realpar), realsigma(xplot, realpar))
nmu = 8;
nsigma = 8; % A spline model with 8 knots for mu and 8 knots for sigma
result5 = euler_reconstruction(data, dt, 'nKnots', [nmu nsigma], 'spline', 'CC', 'L', ...
    L, 'R', R, 'lb', [zeros(1, nmu) - 10, zeros(1, nsigma) + eps], 'ub', zeros(1, nmu + nsigma) + 10, ...
    'solver', 'fmincon', 'search_agents', 1);
plot_results(result5, xplot, realmu(xplot, realpar), realsigma(xplot, realpar))
