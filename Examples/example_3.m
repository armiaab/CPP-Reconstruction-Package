%Here, we apply a spline reconstruction to a univariate Cyanobacterial biome
%measured as phycocyanin concentrations in the Lake Mendota
%(Carpenter et al. 2020). This dataset has a high resolution, with measurements
%taken at minute intervals. We focus on a period during summer thermal
%stratifictiation in 2011, a period when Cyanobacterial blooms are common.
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

fprintf('The data set is downsampled Cyanobacterial data BGA_stdlevel_2011.csv\n')
data = readmatrix(fullfile(rootDir, 'BGA_stdlevel_2011.csv'));
data = data(:, 3);
data = data(1:3:end); % From Table1, we see that this dataset is not Markov. But, a rarified sample with every third data point is Markov
dt = 1; % This is completely arbitrary.

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
L = -6.5;
R = 6;
nmu = 8;
nsigma = 8; % A spline model with 8 knots for mu and 8 knots for sigma
result6 = euler_reconstruction(data, dt, 'nKnots', [nmu nsigma], 'spline', 'CC', 'L', ...
    L, 'R', R, 'lb', [zeros(1, nmu) - 10, zeros(1, nsigma) + eps], 'ub', zeros(1, nmu + nsigma) + 10, 'solver', 'fmincon', 'search_agents', 5);
xplot = linspace(L, R, 2000);
plot_results(result6, xplot)
