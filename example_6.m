here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

S    = load('OUdata1D.mat');
data = S.data(1:100:end);
L    = -2.5; R = 2.5;
data(data < L | data > R) = nan;
dt   = 1;

nmu    = 8;
nsigma = 1;
lb = [zeros(1, nmu)-10,  eps];
ub = [zeros(1, nmu)+10,  10];

par_true = [-1 1];
realmu   = @(x, par) par(1) .* x;
realsig  = @(x, par) par(2) + zeros(size(x));
xplot    = linspace(L, R, 2000);

fprintf('Low-resolution OU - 8+1 cubic + additive sigma\n');
result12 = euler_reconstruction(data, dt, ...
    'nKnots', [nmu nsigma], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb, 'ub', ub, ...
    'search_agents', 5);
fprintf('  NLL=%.4f  time=%.3f s\n', result12.negloglik, result12.computation_time);
plot_results(result12, xplot, realmu(xplot, par_true), realsig(xplot, par_true))
