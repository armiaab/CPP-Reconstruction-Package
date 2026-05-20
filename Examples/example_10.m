%Handling large data sets
%When working with datasets containing millions of data points, the computational burden
%can be significant, leading us to consider using only a portion of the dataset.
%However, selecting the appropriate portion is crucial, as opting for the first 10%,
%last 10%, or middle portion can notably influence the final results,
%potentially introducing bias into the estimated parameters. Since diffusion
%models are Markovian, we can employ mini-batch optimization, where we sample a fraction of
%'data pairs' and solve the underlying optimization problem based on that fraction alone.
%Here, a 'data pair' refers to any consecutive pair (x_t,x_(t+1) ) across the data.
%By randomly selecting a sample comprising just 10% of all data pairs, we can conduct
%the analysis on this subset. This fraction is well-mixed across the entire dataset
%and provides a representative sample. To ensure an even more random selection
%compared to simple random sampling, we recommend and implement a 'stratified'
%random sampling of data pairs. This method offers an excellent representation
%of the entire dataset. After a random sample of data is obtained, we can follow
%either of Euler or Hermite reconstructions as explained in previous sections.
%implementing both Euler and Hermite reconstructions.
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

S = load(fullfile(rootDir, 'OUdata1D.mat'));
data = S.data;
dt = 0.01;
mu = @(x,par)par(1).*x;
sigma = @(x,par)par(2);
%'reconst_fraction' should have two values, the first is the number of
%strata and the second value is the probability for selecting a  data point
result16 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun', eulergrad(mu, sigma),  ...
    'reconst_fraction', [10 0.01], 'lb', [-200 eps], 'ub', [200 200], 'useparallel', true, 'solver', 'fmincon', 'search_agents', 5);
legpoints16 = legitimate_points(data, dt, 'prev', result16, 'prev_range', 0.5, 'j', 3, 'k', 4);
result_her16 = hermite_reconstruction(data, dt, 'prev', legpoints16, 'solver', 'fmincon');

realpar = [-1 1];
realmu = @(x,par)par(1).*x;
realsigma = @(x,par)par(2);
x = [-2:0.01:2];
%plots comparing the fitted model with the true model
plot_results(result16, x, realmu(x, realpar), realsigma(x, realpar))
plot_results(result_her16, x, realmu(x, realpar), realsigma(x, realpar))
