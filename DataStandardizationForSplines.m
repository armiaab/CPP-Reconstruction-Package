
% X space (original space of data)
S = load('MayData1D.mat');data = S.data;   
data = data(1:100000);
L = 0;R = max(data); %we consider the entire range of data
dt = 0.01;
mu = 8; sigma = 8;  
result_x = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', L, 'R', R, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);

% Z space (standardized space)
m = mean(data);
s = std(data);
L_z = (L-m)/s;R_z = (R-m)/s;
data_z = (data-m)./s;
result_z = euler_reconstruction(data_z, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', L_z, 'R', R_z, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);

% Ploting the results in original space in two ways
%1) Direct
xplot = linspace(L,R,3000);
plot_results(result_x,xplot);

%2) Indirectly from Z-results (I need to do write lots of code lines as below. So, this should be done in the package) 
par_x = s.*result_z.estimated_par;
par_x_mu = par_x(1:mu);
par_x_sigma = par_x(mu+1:end);
knotsX_mu = result_x.knots{1};
knotsX_sigma = result_x.knots{2};
mufun = @(x)interp1(knotsX_mu,par_x_mu,x,'spline');
sigmafun = @(x)interp1(knotsX_sigma,par_x_sigma,x,'spline');
figure,plot(xplot,mufun(xplot),'-k');hold on;plot(xplot,0.*xplot,'-k');
figure,plot(xplot,sigmafun(xplot),'-k');ylim([0 0.45]);

% Conclusinon
% If the user wants to use spline modeling and if the data have very large
% scales (like climate data NGRIP20.csv) then it is better to first
% standardize the data. Otherwise, the package has to work with big numbers
% and I have experienced crazy results in such cases. However, if the scale
% of problem is not that big this is not necessary (like in this example)
% In spline modelin we can 
% 1) solve the problem in transformed space (knots should also be transformed when we solve optimization)
% 2) estimate the parameters in transformed space
% 3) parameters in original scale = s*parameters in transformed space. So, at the end
% we consider our ORIGINAL knots and only multiply the parameters by s.











