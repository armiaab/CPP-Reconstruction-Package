addpath('..\')

mu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2 ./(x.^2+par(4).^2);
sigma = @(x,par)par(5)+zeros(size(x));
realpar = [1.01 10 2.75 1.6 0.4];
S = load('MayData1.mat');
data = S.data;
data(rand(size(data)) < 0.01) = NaN;
dt = 0.01;
lb = [0 0 0 0 0] + eps;
ub = 100 .* ones(1, 5);
solvers = {'GWO', 'fmincon'};
Max_iter = 100;
%mask=true(1,5);
%mask(5)=false;
%sigma=@(x,par)fixpars(sigma,x,par,mask,0.4);
RelaxationTime(data)
for i = 1:length(solvers)
    res = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'lb', lb, 'ub', ub, 'solver', solvers{i}, 'maxiter', Max_iter, 'search_agents', 5);
    xdat = linspace(min(data), max(data), 500);
    plot_results(res, xdat, mu(xdat, realpar), sigma(xdat, realpar))
end
res = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'lb', lb, 'ub', ub, 'solver', 'fmincon', 'maxiter', Max_iter);
%realpar=[0.0345678      31.8749      66.6783      95.7615     0.400076]


%test ETA for Q
ETA = shared.create_ETA(1, 8, 2, 0, 'ETA');
d = num2cell(rand(1, 4));
assert(all(ETA_Q_1_8(d{:}) == ETA(d{:})), 'Error in ETA Q');
EZ = shared.create_EZ(3, 3, 1, 1, 'EZ');
d = num2cell(rand(1, 4));
assert(all(EZ_LL_3_3(1, d{:}) == EZ(1, d{:})), 'Error in EZ LL');
%test EZ for LL
EZ = shared.create_EZ(3, 3, 1, 1, 'EZ');
d = num2cell(rand(1, 4));
assert(all(EZ_LL_3_3(1, d{:}) == EZ(1, d{:})), 'Error in EZ LL');
%test EZ for QQ
EZ = shared.create_EZ(3, 3, 2, 2);
d = num2cell(rand(1, 6));
assert(all(EZ_QQ_3_3(1, d{:}) == EZ(1, d{:})), 'Error in EZ QQ');
%test EZ for CC
EZ = shared.create_EZ(3, 3, 3, 3);
d = num2cell(rand(1, 8));
assert(all(EZ_CC_3_3(1, d{:}) == EZ(1, d{:})), 'Error in EZ CC');
%test EZ for Q
EZ = shared.create_EZ(4, 4, 2, 0);
d = num2cell(rand(1, 4));
assert(all(EZ_Q_4_4(1, d{:}) == EZ(1, d{:})), 'Error in EZ Q');
%test EZ for L
EZ = shared.create_EZ(3, 4, 1, 0);
d = num2cell(rand(1, 3));
assert(all(EZ_L_3_4(1, d{:}) == EZ(1, d{:})), 'Error in EZ L');
%test EZ for C
EZ = shared.create_EZ(3, 2, 3, 0);
d = num2cell(rand(1, 5));
assert(all(EZ_C_3_2(1, d{:}) == EZ(1, d{:})), 'Error in EZ C');

%test EZ with zeros sigma
EZ = shared.get_hermite(1, 6, 3, 0);
d = num2cell(rand(1, 5));
EZ2 = shared.get_hermite(1, 6, 3, 3);
%check if we get the same answer with zeros filled in the more complex
%function
d2 = [d {0 0 0}];
assert(all(abs(EZ(0.4, d{:}) - EZ2(0.4, d2{:})) < 1E-15), 'Error in EZ 3,3');
delete('EZ.m');
delete('ETA.m');



%%Example data set Ornstad Ulenberg model
%we load the data OUdata1D and define

%load the data
S=load('OUdata1D');
data1=S.data;
data1(rand(size(data1)) < 0.01) = NaN;

%This data set is created the following mu and sigma functions 
mu = @(x,par)par(1)*x;
sigma = @(x,par)par(2)+zeros(size(x));
%the real parameters of the model
realpar = [-1 1];
%the mu and sigma on which the data is based (used for plotting later)
xtrue = linspace(min(data1), max(data1), 200);
mutrue = mu(xtrue, realpar);
sigmatrue = sigma(xtrue, realpar);



%Relaxation time, fit to autocorrelation function: acf=exp(-c*lag)
%Relaxation time is 1/c
fprintf('The relaxation time = %g\n', RelaxationTime(data1))


%%We assume additive noise and use cubic splines to describe the drift
%%function (mu)
%optionally you can use a number of different starts
result1 = euler_reconstruction(data1, dt, 'nKnots', [4 1], 'spline', 'CC',  ...
    'lb', [zeros(1, 4) - 10, zeros(1, 1)], 'ub', zeros(1, 4 + 1) + 10,  ...
    'solver', 'fmincon', 'search_agents', 1);


plot_results(result1, xtrue, mutrue, sigmatrue)

%result2 = hermite_reconstruction(data1,dt,'prev',result1,'prev_range',0.05,'solver','fmincon','search_agents',5)

%plot_results(result2, xtrue, mutrue, sigmatrue);


%% assume multiplicative noise and linear splines

result2 = euler_reconstruction(data1, dt, 'nknots', [4 4], 'spline', 'LL',  ...
    'lb', [zeros(1, 4) - 10, zeros(1, 4)], 'ub', zeros(1, 4 + 4) + 10);

plot_results(result2, xtrue, mutrue, sigmatrue)

%%parametric use with gradient, note that if we define mu and sigma it will
%be evaluated parametric

%In the example we use the fmincon solver. Here we illustrate how to
%adapt the solver options.
opts = optimoptions('fmincon', 'Display', 'iter');


result3 = euler_reconstruction(data1, dt, 'mu', mu, 'sigma', sigma,  ...
    'gradient_fun', eulergrad(mu, sigma), 'lb', [-200 eps], 'ub', [200 200],  ...
    'solver', 'fmincon', 'maxiter', 300, 'solveroptions', opts);

plot_results(result3, xtrue, mutrue, sigmatrue)

%%parametric without gradient

result4 = euler_reconstruction(data1, dt, 'mu', mu, 'sigma', sigma,  ...
    'lb', [-200 eps], 'ub', [200 200]);

plot_results(result4, xtrue, mutrue, sigmatrue);

%%parametric without gradient and with additive noise (empty assumes
%additive)

result5 = euler_reconstruction(data1, dt, 'mu', mu, 'sigma', [],  ...
    'lb', [-200 eps], 'ub', [200 200]);

plot_results(result5, xtrue, mutrue, sigmatrue);


legpoints2 = legitimate_points(data1, dt, 'prev', result5, 'prev_range', 0.5, 'j', 3, 'k', 5);
for i = 1:length(solvers)
    result4 = hermite_reconstruction(data1, dt, 'prev', legpoints2, 'solver', solvers{i}, 'search_agents', 10);

    plot_results(result4, xtrue, mutrue, sigmatrue);
end

example
example2


disp('all tests passed');

