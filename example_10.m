here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

S    = load('OUdata1D.mat');
data = S.data;
dt   = 0.01;

L = min(data); R = max(data);
data(data < L | data > R) = nan;

nmu    = 8;
nsigma = 8;
lb = [zeros(1, nmu)-10,  zeros(1, nsigma)+eps];
ub =  zeros(1, nmu+nsigma)+10;

par_true = [-1 1];
realmu   = @(x, par) par(1) .* x;
realsig  = @(x, par) par(2) + zeros(size(x));
xplot    = linspace(L, R, 2000);

fprintf('Full OU dataset - 8+8 cubic spline (large data demonstration)\n');
result16 = euler_reconstruction(data, dt, ...
    'nKnots', [nmu nsigma], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb, 'ub', ub, ...
    'search_agents', 5);
fprintf('  NLL=%.4f  time=%.3f s  n_data=%d\n', ...
    result16.negloglik, result16.computation_time, sum(isfinite(data)));
plot_results(result16, xplot, realmu(xplot, par_true), realsig(xplot, par_true))
