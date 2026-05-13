here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

data = readmatrix(fullfile(here, 'NGRIP20.csv'));
data = data(2649:5081);
data = data(1:2:end);

L = -45.5; R = -38.2;
data(data < L | data > R) = nan;
dt   = 1;

nmu    = 7;
nsigma = 7;
lb = [zeros(1, nmu)-10,  zeros(1, nsigma)+eps];
ub =  zeros(1, nmu+nsigma)+10;

fprintf('NGRIP climate data - 7+7 cubic spline\n');
result15 = euler_reconstruction(data, dt, ...
    'nKnots', [nmu nsigma], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb, 'ub', ub, ...
    'search_agents', 5);
fprintf('  NLL=%.4f  time=%.3f s\n', result15.negloglik, result15.computation_time);

xplot = linspace(L, R, 2000);
plot_results(result15, xplot)
