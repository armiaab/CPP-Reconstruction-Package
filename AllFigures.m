eu
% A short explanation for 3 data requirements

% There are three data requirements: 1) Data stationarity 2) Data
% Markovicity 3)Data resolution. If data are simulated then all we need to
% do is to check the third requirement since simulated data are stationary and
% Markov by their construction. However, if data are real then we have to check the
% three data requirements. 
% 
% If the first requirement fails, i.e, the dataset is not stationary then we should chunk data into smaller
% time periods (which can be overlapping) so that the smaller datasets over these periods are stationary and
% repeat reconstruction over each time period. This % is called 'moving-window' approach.
% 
% If the second requirement is not met, i.e., data are not Markov, then you should find a quntity
% called 'Markov-Einstein' (ME) time scale using the command 'ME_TimeScale'. If ME = 1 then your data is Markov and you
% an analyze the entire dataset. However, if ME>1 then your data is not Markov but a rarified sample with every ME-th
% data points is Markov and analysis must be performed on this rarified sample (or more rarified samples) not the entire dataset.

% The third requirement (data resolution):
% For all datasets (real or simulated) it is necessay to check data resolution and find to
% which category data resolution falls. As explained for the second requirement, you must check the
% resolution of a rarified sample which is Markov not the entire dataset. To investigate data resolution 
% you should find a quantity called the 'Relaxation time' (RT) using the command 'RelaxationTime'. Here, we
% have 3 categories as below:

% 1) 1< RT <= 100, which signifies that data falls in the category of 'low-resolution'
% 2) RT > 100, which signifies that data falls in the category of 'high-resolution'
% 3) RT < 1, which signifies that data falls in the category of 'bad-resolution' meaning that the data resolution
% is so low so that the true dynamics is likely not reflected in such a poor dataset

% NOTE 1
% If data resolution is high then you can follow 'Euler' reconstruction. Otherwise a 'Hermite' reconstruction must
% be followed. If your data is a real dataset it is recommended to try Hermite reconstruction irrespective of its
% resolution unless if its resolution is really high.

% NOTE 2
% To perform Hermite reconstruction we need to perform Euler reconstruction first

% We explain these data requirements as well as how to reconstruct a dataset in details across the 8 datasets
% corresponding with Figures 2-9 in the paper. Two datasets are real datasets from Ecology (Figure 5) and the
% climate(Figure 9). The details of these datasets are discussed in the paper.

%%
% Figure 2 (Analyzing a high-resolution dataset simulated from a linear model)

% This is a simulated dataset. So, we only check the third data requirement in below
S = load('OUdata1D.mat');data = S.data;  %load the data
data=data(1:20000); %This is a big timeseries with 10^6 data points. Here, we just use 
%its first 20000 data points
RT = RelaxationTime(data);

% Here, we find that RT = 94.4225 (50<RT<100) meaning that this dataset has a medium 
% (but close to high) resolution. So, we follow a Euler reconstruction in below

%******** A) Fitting a parametric model to data
% Reconstructing the data
dt=0.01;   % this is the time step we used to generate this dataset. Note that it can be arbitrary but we used the same dt for the reprodicibility of the true parameters
mu = @(x,par)par(1).*x;sigma = @(x,par)par(2);   % since mu and sigma are function handles this means that we want to consider parametric modeling
result1 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma,  ...
'gradient_fun', eulergrad(mu, sigma), 'lb', [-200 eps], 'ub', [200 200]);

%******* B) Fitting a spline model to data
L=-2;R=1.8; %since we have a ‘spline’ model it is better to shrink the state space
data(data<L | data>R) = nan;  %This is VERY important: In spline modeling if you consider a smaller range for your data then you must assign ‘nan’ to those few data points falling outside this range.
mu=8;sigma=8;  %since mu and sigma are numbers this means that we want to consider spline modeling with 8 knots for mu and 8 knots for sigma
result2 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', ...
L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon'); % 'CC' means that we want to consider cubic splines for both drift and diffusion functions

% Left panel
S = load('OUdata1D.mat');
data = S.data; 
data = data(1:20000);
t = linspace(0, 200, length(data));
plot(t, data, '-k');
ylim([-inf inf]);
xlabel('time [a.u]', 'FontSize', 16); 
ylabel('state [a.u]', 'FontSize', 16); 
title('Ornstein-Uhlenbeck model', 'FontSize', 16); 

% Right panel
mu = @(x, par) par(1).*x; 
sigma = @(x, par) par(2) + 0.*x; % this is true model  
par = zeros(2, 1); 
par(1) = -1; 
par(2) = 1; % true model parameters
xplot = linspace(L, R, 1000); % a dense mesh across the considered range
plot_results(result2, xplot, mu(xplot, par), sigma(xplot, par));
%%
% Figure 3 (Analyzing a small fraction of a high-resolution dataset)

% In this example we analyze a random stratified sample (with 10 strata) taken from a big dataset OUdata1D.mat where the sample size
% includes only 1% of the original dataset. This random sample is a Markov sample. Therefore, we do not need to check the first two data requirement 

% Checking the resolution of data
% First, note that the resolution of the original dataset and the smaller random sample are exactly the same. Here, the random
% sample is NOT a rarified sample with lower resolution. Rather, its resolution is the same as the original dataset but it just
% has a smaller size (smaller number of data points). For more technical details see section 11 of the tutorial. 

% The original dataset is the same dataset which was already analyzed in the previous section and we found that it had a
% high-resolution. So, any random sample of this dataset also has a high resolution. Therefore, we only apply Euler
% reconstruction to this synthetic dataset

%******** A) Fitting a parametric model to a small random fraction of data containing 1% of the entire data points
% Reconstructing the data
S = load('OUdata1D.mat');data = S.data;  
dt = 0.01;  % the time step remains unchanged (this should not be confused with data rarification)
mu = @(x,par)par(1).*x;sigma = @(x,par)par(2);  % since mu and sigma are function handles this means that we want to consider parametric modeling
result1 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun',  eulergrad(mu, sigma), ...
    'reconst_fraction', [10 0.01], 'lb', [-200 eps], 'ub', [200 200],'useparallel',true, 'solver', 'fmincon', 'search_agents', 5);  % 'reconst_fraction', [10 0.01] means that we want to analyze a stratified random sample with 10 strata containing 1% of the mother dataset

%******* B) Fitting a spline model to data
L=-2;R=1.8; % since we have a  spline’ model it is better to shrink the state space
data(data<L | data>R) = nan;  %This is VERY important: In spline modeling if you consider a smaller range for your data then you must assign ‘nan’ to those few data points falling outside this range.
mu=8;sigma=8;  %since mu and sigma are numbers this means that we want to consider spline modeling with 8 knots for mu and 8 knots for sigma

% If you use the following code lines then the package uses 1% of data for spline modeling. BUT this
% fraction will not be the same fraction we used in A

% result2 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', ...
% 'reconst_fraction', [10 0.01], L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon'); % 'CC' means that we want to consider cubic splines for both drift and diffusion functions

% However, to be able to compare parameteric and spline modelings it is better to use EXACTLY the same
% fraction of data points which was used in A. Bellow, are the needed code lines for this purpose 

sample_pairs = readmatrix('sample_pairs.csv');
n = length(sample_pairs);n = n+n/2-1;
data = nan(1, n);
data(1:3:end) = sample_pairs(1:2:end);
data(2:3:end) = sample_pairs(2:2:end);

result2 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', ...
L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon'); % 'CC' means that we want to consider cubic splines for both drift and diffusion functions

%********Plotting Figure3

% NOTE: You should not expect to get exactly the same result as in Figure 3 in the main text. This is because here we analyze a
% random sample of the original dataset which is different each time you run this section. However, the results are close to each other.

% Left panel
S = load('OUdata1D.mat');data = S.data;
t = linspace(0,10^4,length(data));
sample_pairs = readmatrix('sample_pairs.csv');
X0 = sample_pairs(1:2:end);
X = sample_pairs(2:2:end);
tolerance = 1e-10;  % to work well, the tollerance should not be so small or so big
idx_X0 = arrayfun(@(x) find(abs(data - x) < tolerance, 1), X0);
idx_X = arrayfun(@(x) find(abs(data - x) < tolerance, 1), X);
t_X0 = t(idx_X0);
t_X = t(idx_X);
plot(t,data,'-k');hold on;h1 = plot(t_X0,data(idx_X0), 'or');hold on;h2 = plot(t_X,data(idx_X), 'ob');
xlim([0 20]);  % the actual range of time is [0 10^4]. For the ease of illustration we just depict the initial part (otherwise, the plot gets messy). You can try any range you wish
ylim([min(data(1: 2000)) max(data(1: 2000))]);
xlabel('time [a.u]');
ylabel('state [a.u]');
title('Ornstein-Uhlenbeck model');
legend([h1 h2], {'First elements of data pairs','Second elements of data pairs'});

% Right panel
mu=@(x,par)par(1).*x;sigma=@(x,par)par(2)+0.*x; %this is true model  
par=zeros(2,1);par(1)=-1;par(2)=1; %true model parameters
xplot=linspace(L,R,1000); % a dense mesh across the considered range
plot_results(result2,xplot,mu(xplot,par),sigma(xplot,par));
%%
% Figure 4 (Analyzing a high-resolution dataset simulated from a nonlinear model)

% This is a simulated dataset. So, we only check the third data requirement (data resolution) in below
S = load('MayData1D.mat');data = S.data;  %load the data
data=data(1:100000); %This dataset has 3*10^5 data points. We only use the first third of the dataset 
RT = RelaxationTime(data);

% Here, we find that RT = 3.9785e+03 > 100 meaning that this dataset has a high resolution. So, we follow a Euler reconstruction only

%******** A) Fitting a parametric model to data
% NOTE: % Since this model is nonlinear you should run the following code lines a few times (Alternatively, you can consider increasing the
% 'search_agents' from 5 (default) to a bigger value) to make sure you do not miss the global solution. There is a possibility to get different
% outcomes each time. But, the solution with the lowest objective function (i.e., - sum of log-likelihoods) is the fitest solution. These information apear on the command window  

dt=0.01;  % This is the time step we used tosimulate this dataset (but it is arbitrary)
mu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2./(x.^2+par(4).^2);sigma=@(x,par)par(5);
result1 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun',  eulergrad(mu, sigma), ...
    'lb', zeros(1,5)+eps, 'ub', 15.*ones(1,5),'useparallel',true,'search_agents', 5);  

%******* B) Fitting a spline model to data
% NOTE: since spline models are linear in parameters (though nonlinear in state) we often need to execute the following code
% lines once. This is in contrast to parametric reconstruction in A where one run is likely not sufficient. This
% highlights the superiority and ease of using spline models
dt = 0.01;
L = min(data);R = max(data);  % here, the entire range of dataset is considered
mu = 8;sigma = 8; % A spline model with 8 knots for mu and 8 knots for sigma
result2 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', ...
L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);   % 'CC' means that we want to consider cubic splines for both drift and diffusion functions

%********Plotting Figure4
% Left panel
t = linspace(0,1000,length(data));
plot(t, data,'-k');
xlabel('time [a.u]');ylabel('biomass [a.u]');
title('Grazing model of May');

% Right panel
r=1.01;K=10;g=2.75;a=1.6;s=0.4;  % true parameter values
par = [r K g a s]; 
mu = @(x,par)r.*x.*(1-x./K)-g.*x.^2./(x.^2+a.^2);sigma=@(x,par)s;  % true model
xplot=linspace(L,R,2000);
plot_results(result2,xplot,mu(xplot,par),sigma(xplot,par));

%%
% Figure 5 (Analyzing a high-resolution ecological dataset)

% This is a real dataset (phycocyanin levels in Lake Mendota). So, we should check the three data requirement in below
data = readmatrix('BGA_stdlevel_2011.csv');data = data(:,3);  % Information aabout phycocyanin levels are contained in the third column of BGA_stdlevel_2011.csv dataset

%1) Test of stationarity (if the test result is 1, as it is, then data are stationary)
[~,~,~,~,reg] = adftest(data,'model','ARD','lags',0:15);  % the input 0:15 is the number of lags we try in fitting an autoregressive model to data (it should be big enough. See next command). 
[~,lagndx] = min([reg(:).BIC]);  % this tells us how many lags we need (lagndx is 8<15. So, considering 15 lags in the previous command was enough. Otherwise, increase it and repeat the process)
[h, Pvalue,~]=adftest(data,'model','ARD','lags',lagndx);   % h is the test result with a pvalue of Pvalue (here h = 1, so the dataset is stationary)

%2) Test of Markovicty
%NOTE: Here you need the 'Burg' folder. So, you should add its path to your
%working folder by typing the command 'addpath('The path of the Burg folder')' 
order = 10;  % Number of autoregressive lags to consider (should be big enough. Continue reading ... you will understand if this lag (10) is enough or not)
AR = ME_TimeScale(data,order);

% Here you get the following AR vector
%1.0000   -0.6261   -0.2177   -0.0840   -0.0339   -0.0120   -0.0068   -0.0009   -0.0080    0.0005   -0.0089

% Ignore the first element (which is always 1).
% This suggests that the dataset exhibits long-range correlations. However, the magnitudes of the AR coefficients beyond the fourth
% element (-0.0840) or the fifth element (-0.0339) are small and we can safely assume that the AR order is p = 4 or even p = 3. We
% consider the ME time scale of 4 (note that any time scale bigger than 4 is also Markov and we could consider it though this comes
% at the expense of reduced data resolution). This indicates that the dataset is not Markov but a data ratification by considering
% every forth data point data(1:4:end) is Markov and reconstruction should be performed on this smaller and coarser sample. 

% Checking the resolution of the Markov sample
% Note: We must investigate the resolution of a rarified sample which is Markov not the dataset. Here the Markov
% sample consist of every fourth data point data(1:4:end).
RT = RelaxationTime(data(1:3:end));   % Relaxation time of the Markov sample

% Here, we find that RT = 273.7843 > 100 meaning that the Markov sample has a high resolution. So, we apply a Euler reconstruction

%********************
% Performing reconstruction to the rarified Markov sample
% NOTE 1: This dataset has a high resolution. However, since it is a real dataset we recommend you to go for Hermite reconstruction
% NOTE 2: Since we want to applt Hermite reconstruction to this dataset we recommend you to follow 'quadratic splines' 
% (see the 'QQ' flag below). This significantly reduces the computational burden. However, at the end the package offers you a
% cubic spline approximation (cubic splines have nicer smoothing properties and we prefer them). In short, the package uses quadratic
% splines to reduce the computational time but at the end % delivers a cubic spline approximation

% Hermite reconstruction has two phases: in the first phase we still follow Euler reconstruction to get a rough result. In the second phase a Hermite
% reconstruction is performed to improve the Euler reconstruction outcomes

data = readmatrix('BGA_stdlevel_2011.csv');
data = data(1:4:end); % A rarified sample with every fourth data point is Markov and should be considered
dt = 1;  % This is completely arbitrary. 
L = -6.5;R = 6;
data(data<L | data>R) = nan;  % You must asign NAN to data points which fall outside the considered data range
mu = 8;sigma = 1; % An additive spline model with 8 knots for mu 
result_euler = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', ...
L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);
legpoints = legitimate_points(data, dt, 'prev', result_euler, 'prev_range', 0.5, 'j', 3, 'k', 12);
result_her = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', 'fmincon'); 

% Left panel
figure,
A = readmatrix('BGA_stdlevel_2011.csv');
T = 160:10:250;
T1 = zeros(1,length(T));
time = A(1:end-8,2);  % Information about time is contained in the second column of BGA_stdlevel_2011.csv dataset
for i=1:length(T)
    T1(i) = find(time>T(i),1);
end
my_color = [0 0.4470 0.7410];
plot(A(:, 3),'color',my_color);
xticks(T1);
xticklabels(T);
xlim([-inf inf]);
ylim([-inf inf]);
xlabel('Time (day number in 2011)');
ylabel('Phycocyanin level');
title('Lake Mendota');

% Right panel
figure,
par_euler = result_euler.estimated_par;mu_euler = result_euler.mufun;
par_her = result_her.estimated_par;mu_her = result_her.mufun;
knots = result_euler.knots{1};
xplot=linspace(L,R,2000);
h1 = plot(xplot,mu_euler(xplot,par_euler),'-b');hold on;plot(xplot,0.*xplot,'-k');hold on;
h2 = plot(xplot,mu_her(xplot,par_her),'-k');hold on;
plot(knots,mu_euler(knots,par_euler),'r*');hold on;
plot(knots,mu_her(knots,par_her),'r*');hold on;
legend([h1 h2],{'Euler','Hermite'});
xlabel('State (Phycocyanin level)');ylabel('Drift');
xlim([-inf inf]);
%%
% Figure 6 (Analyzing a low-resolution dataset simulated from a linear model)

% This is a simulated dataset. So, we only check the third data requirement in below
S = load('OUdata1D.mat');data = S.data;  %load the data
data = data(1:100:end); % Only every 100 data points are considered (this is to generate a dataset which is low-resolution and show that Euler reconstruction does not suffice and a Hermite reconstruction is needed)
RT = RelaxationTime(data);   % Here, we find that RT = 1.0074 ~ 1 meaning that this dataset has the lowest resolution we can, in theory, stil hope to perform a successful reconstruction. So, we must follow a Hermite reconstruction 

% Applying Hermite reconstruction using spline modeling 
% This has two phases: in the first phase we still follow Euler reconstruction to get a rough result. In the second phase a Hermite
% reconstruction is performed to improve the Euler reconstruction outcomes

%******** A) Fitting a parametric model to data
% We first go for Euler reconstruction
S = load('OUdata1D.mat');
data = S.data;  %load the data
data = data(1:100:end); % Only every 100 data points are considered
dt = 1; % note that the mother dataset has the time step of dt=0.01 which is multiplied %by 100 to match the time scale of this sample
mu = @(x,par)par(1).*x;sigma = @(x,par)par(2);
result1 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun',  eulergrad(mu, sigma), 'lb', [-200 eps], 'ub', [200 200],'useparallel',true, 'solver', 'fmincon', 'search_agents', 5);

% Then you get the estimates -0.63555 and 0.65456 for mu and sigma, respectively. Clerarly, these estimates deviate a lot from the true values and this necessitates to use Hermite reconstruction as the second

% Now, we go for Hermite reconstruction
legpoints = legitimate_points(data, dt, 'prev', result1, 'prev_range', 0.5, 'j', 3, 'k', 9);
result_her1 = hermite_reconstruction(data, dt, 'prev', legpoints,'solver', 'fmincon');

% Then you get the estimates -0.99683 and 0.99489 for mu and sigma which are wonderful results!

%******** B) Fitting a spline model to data
% We first go for Euler reconstruction

S = load('OUdata1D.mat');data = S.data;  %load the data
data = data(1:100:end); % Only every 100 data points are considered
L = -2.5;R = 2.5;   %Since we have a  spline  model it is better to shrink the state space
data(data<L | data>R) = nan;  %This is VERY important: In spline modeling if you consider a smaller range for your data then you must assign ‘nan’ to those few data points falling outside this range.
dt = 1;
mu = 7; sigma = 7; %In spline modeling mu and sigma are numbers 
result2 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);  % Here, 'QQ' means that we want to fit quadratic splines to both drift and diffusion functions. We recommend you not to use the flag 'CC' (which means cubic splines) for Hermite reconstruction

legpoints = legitimate_points(data, dt, 'prev', result2, 'prev_range', 0.5, 'j', 3, 'k', 9);
result_her2 = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', 'fmincon', 'search_agents', 5);   % The chosen search agent hear should match that used in the Euler reconstruction

% Left and right panels
mu = @(x,par)par(1).*x;sigma = @(x,par)par(2)+0.*x; %this is true model  
par = zeros(2,1);par(1) = -1;par(2) = 1; %true model parameters
xplot = linspace(L,R,2000); % a dense mesh across the considered range
plot_results(result2,xplot,mu(xplot,par),sigma(xplot,par)); % Euler & true models 
plot_results(result_her2,xplot,mu(xplot,par),sigma(xplot,par));%Hermite & true models 
%%
% Figure 7 (Analyzing a low-resolution dataset simulated from a nonlinear model)

% This is a simulated dataset. So, we only check the third data requirement in below
S = load('MayData1D.mat');data = S.data;  %load the data
data=data(1:300:end);  %We consider every 300-th data point (this is to generate a dataset which is low-resolution and show that Euler reconstruction does not suffice and a Hermite reconstruction is needed)
RT = RelaxationTime(data);

% Here, we find that 1 < RT = 12.7959< 100 meaning that this dataset has a low resolution. So, we must follow a Hermite reconstruction 

% Applying Hermite reconstruction using spline modeling 
% This has two phases: in the first phase we still follow Euler reconstruction to get a rough result. In the second phase a Hermite
% reconstruction is performed to improve the Euler reconstruction outcomes

% First and second phases (Euler and Hermite reconstructions)
S = load('MayData1D.mat');data = S.data;   
data=data(1:300:end); 
L = 0;R = max(data); %we consider the entire range of data
% data(data<L | data>R) = nan;  
dt = 3;  % Since in the mother dataset dt = 0.01 and here we considered every 300-th data points the actual time step is 300*0.01 = 3  
mu = 8; sigma = 8;  
% AN EXTREMELY IMPORTANT POINT: To improve execution time and also accuracy we strongly recommend to use 'quardatic' splines 
% using the flag 'QQ' when Hermite reconstruction is needed, as is the case here
result = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);
legpoints = legitimate_points(data, dt, 'prev', result, 'prev_range', 0.5, 'j', 3, 'k', 9);   % Note that for low-resolution data a K value between 6 and 12 is needed (here we used K = 9)
result_her = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver',...,
    'fmincon');

% Left panel
mu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2./(x.^2+par(4).^2);sigma=@(x,par)par(5)+0.*x;  % true model
par = zeros(5,1);par(1) = 1.01;par(2)= 10;par(3) = 2.75;par(4) = 1.6;par(5) = 0.4; %true model parameters
xplot=linspace(L,R,2000); % a dense mesh across the considered range
plot_results(result,xplot,mu(xplot,par),sigma(xplot,par)); % Euler & true models 

% Right panel
plot_results(result_her,xplot,mu(xplot,par),sigma(xplot,par));
%%

% Figure 7 (Analyzing a low-resolution and replicate dataset)

% This is a simulated dataset. So, we only check the third data requirement in below
S = load('MayData1D_Replicate.mat');
data = S.data; % replicate data should be supplied as a cell array  (this dataset has 3 replicates)
rarified_data = cellfun(@(x) x(end:-30:1), data, 'UniformOutput', false);   % the rarified sample is constructed by choosing every 30-th data points
rarified_data = prepare_replicateData(rarified_data); % to reconstruct replicate data we first need to use this function. The rest of calculations are similar to those for typical datasets
L = min(rarified_data); % L = min(data) is ok because the density of data near min(data) is high enough
R = max(rarified_data);R = 6.50;  % The 10 biggest data points are {6.9149,6.6910,6.6103,6.5917,6.5915,6.5668,6.5271,6.4973,6.4756,6.4597}. If you try R = max(data) you can get unreliable results as the density of data mear max(data) is too low. We considered a slightly smaller number. This requires a bit of experimentation
rarified_data(rarified_data<L | rarified_data>R) = nan;  %This is VERY important: In spline modeling if you consider a smaller range for your data then you must assign rnana to those few data points falling outside this range.

RT = RelaxationTime(rarified_data);
% Here, we find that RT = 9.5971 < 100 meaning that this dataset has a low resolution. So, we must follow a Hermite reconstruction 

dt = 3;  % this is the true resolution of replicates (the original data has dt = 0.1, we considered every 30-th data points. So, the rarified sample has dt = 30*0.1 = 3)
mu = 8; sigma = 1;  % this is an additive model because sigma = 1
result = euler_reconstruction(rarified_data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'useparallel', true, 'solver', 'fmincon', 'search_agents', 5);

legpoints = legitimate_points(rarified_data, dt, 'prev', result, 'prev_range', 0.5, 'j', 3, 'k', 9);   % Note that for low-resolution data a K value between 6 and 12 is needed (here we used K = 9)
result_her = hermite_reconstruction(rarified_data, dt, 'prev', legpoints, 'useparallel', true, 'solver',...,
    'fmincon');

% You might get slightly different outcomes each time for the Hermite reconstruction. Please run this part a few times. Note that the best
% solution is the one with the smalest objective value (as below):

% Estimated parameters : 
% 0.20682      0.1962   0.0030841    0.012677    0.093598    0.081107   -0.049113    -0.37251     0.38311
% - sum of log-likelihoods : 379.4046


% Top panel
S = load('MayData1D_Replicate.mat');
S = load('MayData1D_Replicate.mat');
data = S.data;
rep1_low = data{1}(end:-30:1);rep2_low= data{2}(end:-30:1);rep3_low = data{3}(end:-30:1);
t1_low = linspace(0, floor(length(rep1_low)*3), length(rep1_low));
t2_low = linspace(0, floor(length(rep2_low)*3), length(rep2_low));
t3_low = linspace(0, floor(length(rep3_low)*3), length(rep3_low));
h1 = plot(t1_low,fliplr(rep1_low'), '.-k', 'MarkerSize', 10);hold on;
h2 = plot(t2_low,fliplr(rep2_low'), '.-b', 'MarkerSize', 10);hold on;
h3 = plot(t3_low,fliplr(rep3_low'), '.-', 'MarkerSize', 10, 'color', [0.8500 0.3250 0.0980]);
xlim([0 750]);
xlabel('time [a.u]');ylabel('biomass [a.u]');
title('A small low-resolution replicate dataset');

% Bottom panels
r=1.01;K=10;g=2.75;a=1.6;s=0.4;  % true parameter values
par = [r K g a s]; 
mu = @(x,par)r.*x.*(1-x./K)-g.*x.^2./(x.^2+a.^2);sigma=@(x,par)s;  % true model
xplot=linspace(L,R,2000);
plot_results(result,xplot,mu(xplot,par),sigma(xplot,par));hold on;
plot_results(result_her,xplot,mu(xplot,par),sigma(xplot,par)); 
%%

% Figure 9

% This is a real dataset (a climate record covering the northern hemisphere climate from 70000 to 20000 years before the present time). So, we must check the three data requirement in below
data = readmatrix('NGRIP20.csv');   % Load the data
data=data(2649:5081);  % this corresponds with time period from 70000 to 20000 years before the present time (Note that the entire dataset does not meet stationarity (first data assumption). But, this portion of data is stationary. Note that reconstruction must be performed on stationary data. 

%1) Test of stationarity (if the test result is 1, as it is, then data are stationary)
[~,~,~,~,reg] = adftest(data,'model','ARD','lags',0:10);  % the input 0:10 is the number of lags we try in fitting an autoregressive model to data (it should be big enough. See next command). 
[~,lagndx] = min([reg(:).BIC]);  % this tells us how many lags we need (lagndx is 5<10. So, considering 10 lags in the previous command was enough. Otherwise, increase it and repeat the process)
[h, Pvalue,~]=adftest(data,'model','ARD','lags',lagndx);   % h is the test result with a pvalue of Pvalue (here h = 1, so the dataset is stationary)

%2) Test of Markovicty 
%NOTE: Here you need the 'Burg' folder. So, you should add its path to your
%working folder by typing the command 'addpath('The path of the Burg folder')' 
order = 10;  % Number of autoregressive lags to consider (should be big enough. Continue reading ... you will understand if this lag (10) is enough ro not)
AR = ME_TimeScale(data,order);

% Here you get the following AR vector
%1.0000   -0.4796   -0.2098   -0.1082   -0.0616   -0.0197   -0.0468

% Ignore the first element (which is always 1). The lag after which the
% rest of lements of AR are 0 is the 'Markov-Einstein' (ME) time scale. This real dataset shows long-range correlations. However, 
% after lag 3 (-0.1082) the correlations are rather samll and we can  roughly say that ME time scale is 3. Unfortunately, this dataset is quite
% small with only 2433 data points, so considering a rarified sample with every third data points severely shrinks the data. We, therefore, decided to
% consider every othewr data points data(1:2:end) for the analysis, so the ME = 2.  A dataset is Markov if it has a ME time scale of 1. Otherwise, the
% entire dataset is not Markov. However, in this dataset if we consider every other data point, data(1:2:end), then this rarified sample is nearly Markov. We call it Markov sample 

% Checking the resolution of Markov sample
% Note: We must investigate the resolution of a rarified sample which is Markov not the dataset. Here the Markov
% sample consist of every third data point data(1:3:end).
RT = RelaxationTime(data(1:2:end));   % Relaxation time of the Markov sample

% Here, we find that RT = 1 < 15.0951< 100 meaning that the Markov sample has a low-resolution. So, we apply a Hermite reconstruction

% Applying Hermite reconstruction using spline modeling 
% This has two phases: in the first phase we still follow Euler reconstruction to get a rough result. In the second phase a Hermite
% reconstruction is performed to improve the Euler reconstruction outcomes

% First and second phases (Euler and Hermite reconstructions)
data = readmatrix('NGRIP20.csv');
data=data(2649:5081);
data = data(1:2:end);  % Markov sample
L = -45.5;R = -38.2;  % We shrink the data range a bit since near data borders there is very few data points (which adversely impact the analysis)
dt = 1;  % This is completely arbitrary. Note that in spline modeling, (since splines are linear functions of parameters) if you, for instance
% would consider dt = 2 then drift parameters will be divided by 2 and diffusion parameters will be divided by sqert(2), so on.
data(data<L | data>R) = nan;
mu = 7; sigma = 7;  
% AN EXTREMELY IMPORTANT POINT: To improve execution time and also accuracy we strongly recommend to use 'quardatic' splines 
% using the flag 'QQ' when Hermite reconstruction is needed, as is the case here
result = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 20); % we used 20 search agents (default is 5)
legpoints = legitimate_points(data, dt, 'prev', result, 'prev_range', 0.5, 'j', 3, 'k', 9);   % Note that for low-resolution data a K value between 6 and 12 is needed (here we used K = 9)
result_her = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', 'fmincon');

% Top panel
data = readmatrix('NGRIP20.csv');
data=data(2649:5081);
data = data(1:2:end);
my_color = [0 0.4470 0.7410];
plot(2649:2:5081,data,'-','color',my_color);
age = linspace(2649,5081,6);
xticks(age);xticklabels([70 60 50 40 30 20]);
yticks(-47:2:-37);yticklabels(-47:2:-37);
ylabel('\delta^{18}O (Permille)');
xlabel('Age (ky before present)');
title('A low-resolution climate recrord');
xlim([-inf inf]);ylim([-47 -36]);

% bottom panels
xplot=linspace(L,R,2000); % a dense mesh across the considered range
plot_results(result,xplot);
plot_results(result_her,xplot,'eff_potential');  % Here, we can also ask the packge to plot the effective potential (it is similar to the usual concept of potential function).
