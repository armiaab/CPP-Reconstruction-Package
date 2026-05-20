S = load('OUdata1D.mat');data = S.data; %load the data
data = data(1:100:end); % Only every 100 data points are considered
L = -2.5;R = 2.5;%Since we have a spline model it is better to shrink the state space
data(data<L | data>R) = nan; 
dt = 1;
mu = 7; sigma = 7;
result7 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 'lb', [zeros(1, mu) - 10, zeros(, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);

