here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

S    = load('OUdata1D.mat');
data = S.data(1:20000);
dt   = 0.01;

realpar = [-1 1];
realmu  = @(x, par) par(1) .* x;
realsig = @(x, par) par(2) + zeros(size(x));
x_plot  = linspace(-2, 2, 500);

fprintf('Fit 1: 2-knot linear mu + additive sigma\n');
result1 = euler_reconstruction(data, dt, ...
    'nKnots', [2 1], 'spline', 'LC', ...
    'lb', [-200 -200 eps], 'ub', [200 200 200], ...
    'search_agents', 1);
fprintf('  NLL=%.4f  time=%.3f s\n', result1.negloglik, result1.computation_time);
plot_results(result1, x_plot, realmu(x_plot, realpar), realsig(x_plot, realpar))

fprintf('Fit 2: 4-knot cubic mu + additive sigma (missing data)\n');
data2 = data;
data2(rand(size(data2)) < 0.01) = nan;
lb2 = [repmat(-200, 1, 4), eps];
ub2 = repmat(200, 1, 5);
result2 = euler_reconstruction(data2, dt, ...
    'nKnots', [4 1], 'spline', 'CC', ...
    'lb', lb2, 'ub', ub2, ...
    'search_agents', 3);
fprintf('  NLL=%.4f  time=%.3f s\n', result2.negloglik, result2.computation_time);
plot_results(result2, x_plot, realmu(x_plot, realpar), realsig(x_plot, realpar))

fprintf('Fit 3: 8+8 cubic spline on bounded domain\n');
L = -2; R = 1.8;
data3 = S.data(1:20000);
data3(data3 < L | data3 > R) = nan;
lb3 = [zeros(1,8)-10,  zeros(1,8)+eps];
ub3 =  zeros(1,16)+10;
result3 = euler_reconstruction(data3, dt, ...
    'nKnots', [8 8], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb3, 'ub', ub3, ...
    'search_agents', 1);
fprintf('  NLL=%.4f  time=%.3f s\n', result3.negloglik, result3.computation_time);
plot_results(result3, x_plot, realmu(x_plot, realpar), realsig(x_plot, realpar))
