here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

S    = load('OUdata1D.mat');
data = S.data(1:100:end);
dt   = 1;

L = -2.5; R = 2.5;
data(data < L | data > R) = nan;

par_true = [-1 1];
realmu   = @(x, par) par(1) .* x;
realsig  = @(x, par) par(2) + zeros(size(x));
xplot    = linspace(L, R, 2000);

fprintf('Low-resolution OU - 8+1 cubic + additive sigma\n');
lb = [zeros(1,8)-10,  eps];
ub = [zeros(1,8)+10, 10];
result10 = euler_reconstruction(data, dt, ...
    'nKnots', [8 1], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb, 'ub', ub, ...
    'search_agents', 5);
fprintf('  NLL=%.4f  time=%.3f s\n', result10.negloglik, result10.computation_time);
plot_results(result10, xplot, realmu(xplot, par_true), realsig(xplot, par_true))
