%In this example we analyse high resolution data set
%generated with the linear Ornstein-Uhlenbeck (OU) model
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

fprintf('The data set is full OUdata1D.mat\n')
S = load(fullfile(rootDir, 'OUdata1D.mat'));
data = S.data;

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

%%Data resolution
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

%% fitting a parametric model on high resolution data
S = load(fullfile(rootDir, 'OUdata1D.mat'));
data = S.data; %load the data
data = data(1:20000); %This is a big timeseries with 10^6 data points. Here, we just use
%its first 20000 data points
dt = 0.01;
mu = @(x,par)par(1).*x;
sigma = @(x,par)par(2);
result1 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, ...
    'gradient_fun', eulergrad(mu, sigma), 'lb', [-200 eps], 'ub', [200 200]);
%model that generated this data set:
realpar = [-1 1];
realmu = @(x,par)par(1).*x;
realsigma = @(x,par)par(2);
x = [-2:0.01:2];
%plot comparing the fitted model with the true model
plot_results(result1, x, realmu(x, realpar), realsigma(x, realpar))


%%Fitting a slightly different parametric model, high resolution with some missing values:
S = load(fullfile(rootDir, 'OUdata1D.mat'));
data = S.data; %load the data
data = data(1:20000); %This is a big timeseries with 10^6 data points. Here, we just use
data(rand(size(data)) < 0.01) = nan; %this is to show you that the package works in the presence of NANs
dt = 0.01;
mu = @(x,par)par(1).*x+par(2).*x.^2;
sigma = @(x,par)par(3); %here, we have 3 %parameters
result2 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, ...
    'gradient_fun', eulergrad(mu, sigma), 'lb', [-200 -200 eps], 'ub', [200 200 200]);
realpar = [-1 1];
realmu = @(x,par)par(1).*x;
realsigma = @(x,par)par(2);
plot_results(result2, x, realmu(x, realpar), realsigma(x, realpar))


%% Fitting a spline model on high resolution data
S = load(fullfile(rootDir, 'OUdata1D.mat'));
data = S.data; %load the data
data = data(1:20000);
dt = 0.01;
L = -2;
R = 1.8; %since we have a 'spline' model it is better to shrink the state space
data(data < L | data > R) = nan; %This is VERY important: In spline modeling if you consider a smaller range for your data then you must assign 'nan' to those few data points falling outside this range.
nmu = 8;
nsigma = 8; %since nmu and nsigma are numbers this means that we want to consider %spline modeling with 8 knots for mu and 8 knots for sigma
result3 = euler_reconstruction(data, dt, 'nKnots', [nmu nsigma], 'spline', 'CC', 'L', ...
    L, 'R', R, 'lb', [zeros(1, nmu) - 10, zeros(1, nsigma) + eps], 'ub', zeros(1, nmu + nsigma) + 10, 'solver', 'fmincon'); %we have 8+8 parameters, so, 'lb' and 'ub' should have 16 %elements. The vector of lower bounds 'lb' has 8 lower bounds for mu (which are -10) and 8 %lower bounds for sigma (which are eps, i.e., infinitesimal). All 16 elements of 'ub' are 10
realpar = [-1 1];
realmu = @(x,par)par(1).*x;
realsigma = @(x,par)par(2);
plot_results(result3, x, realmu(x, realpar), realsigma(x, realpar))
