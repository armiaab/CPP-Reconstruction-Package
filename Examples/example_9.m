%Here we use a δ^18 Ο climate record, with a resolution of 20 years, extending from 70 to 20 thousand
%years before the present time from NGRIP. This is used as a proxy for the temperature of
%the northern hemisphere which shows that the northern hemisphere climate alternated
%between cold stadial and warmer interstadial alternative climate states.
%In this time period majority of Dansgaard-Oescher events occurred.
thisDir = fileparts(mfilename('fullpath'));
rootDir = fullfile(thisDir, '..');
addpath(rootDir);
addpath(fullfile(rootDir, 'Burg'));

data = readmatrix(fullfile(rootDir, 'NGRIP20.csv'));
%selection of the data range
data = data(2649:5081);
%%Markovitiy Markov-Einstein timescale
order = 10; % order is the maximum AR order being considered (should be long enough)
%Markov-Einstein time scale
AR = ME_TimeScale(data, order);
AR1 = AR(AR > 0.1);
fprintf('AR model:%s\nThe Markov-Einstein scale is ca. %d\n', mat2str(AR, 5), numel(AR1));

%Although this dataset is not Markov, a sample containing every other point exhibit
%Markov property approximately. Additionally, due to its low resolution,
%Hermite reconstruction is necessary. To reconstruct this dataset, use the
%following commands.
data = data(1:2:end);
%%Data resolution
R = RelaxationTime(data);

fprintf('Data resolution\nRelaxation time = %g\n', R);
if R > 100
    disp('High resolution data')
elseif R > 50
    disp('Medium resolution data')
elseif R > 1
    disp('Low resolution data')
else
    disp('Extremely low resolution data');
end
%we take the data range with the majority of data:
L = -45.5;
R = -38.2;
data(data < L | data > R) = nan;
dt = 1;
nmu = 7; %number of knots for mu
nsigma = 7;%number of knots for sigma
%Euler reconstruction using a spline model
result15 = euler_reconstruction(data, dt, 'nKnots', [nmu nsigma], 'spline', 'QQ', 'L', L, 'R', R, ...,
    'lb', [zeros(1, nmu) - 10, zeros(1, nsigma) + eps], 'ub', zeros(1, nmu + nsigma) + 10, 'solver', 'fmincon', 'search_agents', 20); % we used 20 search agents
%find legitimate points
disp('find legitimate points');
legpoints15 = legitimate_points(data, dt, 'prev', result15, 'prev_range', 0.5, 'j', 3, 'k', 9);
%Hermite reconstruction
result_her15 = hermite_reconstruction(data, dt, 'prev', legpoints15, 'solver', 'fmincon');
%plot the results
xplot = linspace(L,R,2000);
plot_results(result15,xplot);
plot_results(result_her15,xplot);
