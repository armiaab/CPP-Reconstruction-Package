here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(here);

S    = load('MayData1D.mat');
data = S.data(1:300:end);
dt   = 3;

data = (data - mean(data)) ./ std(data);
L    = min(data); R = max(data);
data(data < L | data > R) = nan;

lb = [zeros(1,8)-10, zeros(1,8)+eps];
ub =  zeros(1,16)+10;

r=1.01; K=10; g=2.75; a=1.6; s=0.4;
realpar = [r K g a s];
realmu  = @(x, par) par(1).*x.*(1-x./par(2)) - par(3).*x.^2./(x.^2+par(4).^2);
realsig = @(x, par) par(5) + zeros(size(x));
xplot   = linspace(0, 6, 2000);

fprintf('May model low-resolution - 8+8 cubic (standardised data)\n');
result13 = euler_reconstruction(data, dt, ...
    'nKnots', [8 8], 'spline', 'CC', ...
    'L', L, 'R', R, 'lb', lb, 'ub', ub, ...
    'search_agents', 5);
fprintf('  NLL=%.4f  time=%.3f s\n', result13.negloglik, result13.computation_time);
plot_results(result13, linspace(L, R, 2000))
