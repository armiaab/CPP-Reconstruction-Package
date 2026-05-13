here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

S    = load('MayData1D.mat');
data = S.data(1:100000);
dt   = 0.01;

r = 1.01; K = 10; g = 2.75; a = 1.6; s = 0.4;
realpar  = [r K g a s];
realmu   = @(x, par) par(1).*x.*(1-x./par(2)) - par(3).*x.^2./(x.^2+par(4).^2);
realsig  = @(x, par) par(5) + zeros(size(x));
L        = min(data);
R        = max(data);
xplot    = linspace(L, R, 2000);

fprintf('Spline fit on May model data (8+8 cubic)\n');
lb = [zeros(1,8)-10,  zeros(1,8)+eps];
ub =  zeros(1,16)+10;
result5 = euler_reconstruction(data, dt, ...
    'nKnots', [8 8], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb, 'ub', ub, ...
    'search_agents', 1);
fprintf('  NLL=%.4f  time=%.3f s\n', result5.negloglik, result5.computation_time);
plot_results(result5, xplot, realmu(xplot, realpar), realsig(xplot, realpar))
