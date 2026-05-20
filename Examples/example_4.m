%In this example we analyse low resolution downsampled bdata set
%generated with the linear Ornstein-Uhlenbeck (OU) model
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

S = load(fullfile(rootDir, 'OUdata1D.mat'));
data = S.data;
%We downsample to illustrate the performance for low resolution data
data = data(1:100:end);

Rt = RelaxationTime(data);

fprintf('Data resolution\nRelaxation time = %g\n', Rt);
if Rt > 100
    disp('High resolution data')
elseif Rt > 50
    disp('Medium resolution data')
elseif Rt > 1
    disp('Low resolution data')
else
    disp('Extremely low resolution data');
end

dt = 1; % note that the mother dataset has the time step of dt=0.01 which is multiplied %by 100 to match the time scale of this sample
mu = @(x,par)par(1).*x; sigma = @(x,par)par(2);
result10 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun',  eulergrad(mu, sigma), 'lb', [-200 eps], 'ub', [200 200],'useparallel',true, 'solver', 'fmincon', 'search_agents', 5);

legpoints = legitimate_points(data, dt, 'prev', result10, 'prev_range', 0.5, 'j', 3, 'k', 9);
result_her10 = hermite_reconstruction(data, dt, 'prev', legpoints,'solver', 'fmincon');
mu = @(x,par)par(1).*x;
sigma = @(x,par)par(2)+0 .*x; %this is true model
par = zeros(2, 1);
par(1) = -1;
par(2) = 1; %true model parameters
L = min(data);
R = max(data);
xplot = linspace(L, R, 2000); % a dense mesh across the considered range
plot_results(result10, xplot, mu(xplot, par), sigma(xplot, par)); % Euler & true models
plot_results(result_her10, xplot, mu(xplot, par), sigma(xplot, par)); %Hermite & true models
