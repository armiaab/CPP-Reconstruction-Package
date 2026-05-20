%In this example we analyse low resolution downsampled data set
%generated with the nonlinear May model
%Here we apply  spline method, but assume additive noise.
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

S = load(fullfile(rootDir, 'MayData1D.mat'));
data = S.data; %load the data

%We downsample to illustrate the performance for low resolution data
data = data(1:300:end); %We consider every 300-th data point

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

dt = 3; % Since in the mother dataset dt=0.01 and here we considered every 300-th data points the actual time step is 3
%we standardize the data:
m = mean(data); s = std(data);
data = (data - m) ./ s;
mu = @(x,par)par(1).*x.^3+par(2).*x.^2+par(3).*x+par(4);
sigma = @(x,par)par(5);
result13 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun', eulergrad(mu, sigma), 'lb', [-5 .* ones(1, 4) 0], 'ub', [5 .* ones(1, 4) 5], 'useparallel', true, 'search_agents', 5); %Since data are standardized the vectors of lower and upper bounds for the drift part are better to be symmetrical about 0, which is [-5 5] here. For the diffusion parameters it should be [0 5]
%the original model
mu = @(x,par)s.*par(1).*((x-m)./s).^3+s.*par(2).*((x-m)./s).^2+s.*par(3).*(x-m)./s+s.*par(4);
sigma = @(x,par)s.*par(5);
%where m = mean(data) and s = std(data). Next, we get some legitimate points by the following command
legpoints = legitimate_points(data, dt, 'prev', result13, 'prev_range', 0.5, 'j', 3, 'k', 9);
%and finally we go for Hermite reconstruction as below
result_her13 = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', 'fmincon');

result = result13; %plot for Euler outcomes
S = load(fullfile(rootDir, 'MayData1D.mat'));
data = S.data;
m = mean(data);
s = std(data);
result.s = s;
result.m = m;
result.par_est = result.estimated_par;
result.mufun = @(x,par_est)s.*par_est(1).*((x-m)./s-par_est(2)).*((x-m)./s-par_est(3)).*((x-m)./s-par_est(4)); %back-transformed drift (subtract state by data mean and divide by data standard deviation. Finally multiply the whole by data standard deviation)
result.sigmafun = @(x,par_est)s.*par_est(5)+0 .*x; %back-transformed diffusion
mu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2 ./(par(4).^2+x.^2); %true drift
sigma = @(x,par)par(5)+0 .*x; %true diffusion
par = zeros(5, 1);
par(1) = 1.01;
par(2) = 10;
par(3) = 2.75;
par(4) = 1.6;
par(5) = 0.4; %true model parameters
xplot = linspace(0, 6, 2000);
plot_results(result, xplot, mu(xplot, par), sigma(xplot, par));

result = result_her13; %plot for Hermite outcomes
S = load(fullfile(rootDir, 'MayData1D.mat'));
data = S.data;
m = mean(data);
s = std(data);
result.s = s;
result.m = m;
result.par_est = result.estimated_par;
result.mufun = @(x,par_est)s.*par_est(1).*((x-m)./s-par_est(2)).*((x-m)./s-par_est(3)).*((x-m)./s-par_est(4));
result.sigmafun = @(x,par_est)s.*par_est(5)+0 .*x;
realmu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2 ./(par(4).^2+x.^2);
realsigma = @(x,par)par(5)+0 .*x;
realpar = zeros(5, 1);
realpar(1) = 1.01;
realpar(2) = 10;
realpar(3) = 2.75;
realpar(4) = 1.6;
realpar(5) = 0.4;
xplot = linspace(0, 6, 2000);
plot_results(result, xplot, realmu(xplot, realpar), realsigma(xplot, realpar));
