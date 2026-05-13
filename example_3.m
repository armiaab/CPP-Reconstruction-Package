here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

data = readmatrix(fullfile(here, 'BGA_stdlevel_2011.csv'));
data = data(1:3:end, 3);
dt   = 1;

L   = -6.5;
R   = 6;
data(data < L | data > R) = nan;

nmu    = 8;
nsigma = 8;
lb = [zeros(1, nmu)-10,  zeros(1, nsigma)+eps];
ub =  zeros(1, nmu+nsigma)+10;

fprintf('Cyanobacterial data - 8+8 cubic spline\n');
result6 = euler_reconstruction(data, dt, ...
    'nKnots', [nmu nsigma], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb, 'ub', ub, ...
    'search_agents', 5);
fprintf('  NLL=%.4f  time=%.3f s\n', result6.negloglik, result6.computation_time);

xplot = linspace(L, R, 2000);
plot_results(result6, xplot)
