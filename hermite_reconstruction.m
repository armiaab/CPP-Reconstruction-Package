%hermite_reconstruction - reconstruct Langevin equation using the refined Hermite
%method (Aït-Sahalia, 2002, Bakshi & Ju 2005). 
%
%
%usage:
%res=hermite_reconstruction(data,dt,'name',value,...)
%data - vector with the equally spaced time series
%dt   - the time interval between two data points
%name-value pairs:
%'prev' - struct with previous results of euler_reconstruction or
%         hermite_reconstruction
%'prev_range' - range around previous results as fraction of the 
%        previous parameters (in case of spline the range in parameters) - default 0.05 
%'j'   - the order of the Hermite expansion coeficients - default 3
%'k'   - the order of the temporal Taylor expansion - default 4
%'lb'  - vector with the lower bounds of all parameters/knots (by default Euler-prev_range))
%'ub'  - vector with the lower bounds of all parameters/knots (by default Euler+prev_range))
%'solver' solver for the log-likelihood (default fmincon)
%       fmincon, GWO and PS
%'solveroptions' object or struct with additional options dependent on the 
%      used solver (for instance made with optimoptions)
%'useparallel' use parallel computing (default false)
%'search_agents', number of search agents for GWO
%'maxiter' maximum number of iterations
%
% Spline (cannot be combined with parametric): 
%'nknots' the number of knots for mu and sigma use [n 1] for additive noise (default [8 8])
%'knots' (alternative for nknots) the values of the knots
%'spline' type of splines for mu and sigma (for instance LL for both
%         linear) or 'SCSL' ('SCS' for mu and 'L' for sigma)  default CC
%        'L' = linear interpolation
%        'C' = cubic (order 2) spline
%        'Q' = order 3 spline
%        'P' = pchip 
%        'SCS' = cubic smoothing spline
%        'Approximate' = SCS but uses L for fast fitting
%
%Parametric (cannot be combined with Spline):
%'mu'    parametric function handle for mu
%'sigma' parametric function handle for sigma (empty = additive)
%'npars' (optional) number of parameters default automatic
%'analytic_diffs' by default the sym toolbox is used to calculate
%derivatives, set this option to false to approximate them numerically
%
%References
%Aït‐Sahalia, Y. (2002) Maximum likelihood estimation of discretely sampled diffusions:
%     a closed‐form approximation approach. Econometrica 70, 223-262.
%Bakshi, G. & N. Ju (2005) A Refinement to Aït‐Sahalia's (2002)“Maximum Likelihood 
%     Estimation of Discretely Sampled Diffusions: A Closed‐Form Approximation
%     Approach”. The Journal of Business 78, 2037-2052.
%
%
%Author Babak M.S. Ariani (gmail: m.shojaeiarani@gmail.com)
function result = hermite_reconstruction(data, dt, varargin)
    defaults_param = struct('prev', [], 'j', 3, 'k', 4, 'l', [], 'r', [], 'refined', true, 'prev_range', 0.05, 'lb', [], 'ub', [], 'npars', [],  ...
        'mu', [], 'sigma', [], 'analytic_diffs', true, 'solver', 'fmincon', 'reconst_fraction', [], 'solveroptions', [], 'useparallel', false, 'search_agents', 5, 'maxiter', realmax);
    defaults_spline = struct('prev', [], 'j', 3, 'k', 4, 'l', [], 'r', [], 'refined', true, 'prev_range', 0.05, 'lb', [], 'ub', [], 'nknots', [8 8], 'knots', [], 'spline', 'CC',  ...
        'solver', 'fmincon', 'reconst_fraction', [], 'solveroptions', [], 'useparallel', false, 'search_agents', 5, 'maxiter', realmax);
    if ischar(data) && strcmp(data, '-d') %default options;
        result = struct('parametric', defaults_param, 'spline', defaults_spline);
        return;
    end

    %do initial checks and collect the used options
    [options, parametric, index_mu, knots] = shared.parse_options_reconst(data, varargin, defaults_param, defaults_spline);

    %replace points outside the domain borders (L and R) with NaN
    if ~isempty(options.l) && ~isnan(options.l)
        ndx = data < options.l;
    else
        ndx = false(size(data));
    end

    tim = tic;
    hfig = findobj(0, 'tag', 'optimdlg');
    if ~isempty(hfig)
        close(hfig);
    end
    if ~parametric
        [cost, mu, sigma] = shared.Create_Hermite_Cost_Spline(data, dt, options.j, options.k, knots, index_mu, options.spline{1}, options.spline{2}, options.refined, options.reconst_fraction);
    else
        mu = options.mu;
        sigma = options.sigma;
        cost = shared.Create_Hermite_Cost_Parametric(data, dt, options.j, options.k, mu, sigma, numel(options.lb), options.refined, options.analytic_diffs, options.reconst_fraction);
    end
    par0 = options.lb + (options.ub - options.lb) .* rand(1, numel(options.ub));
    dim = length(par0);
    switch lower(options.solver)
        case 'legpoints'
            dim = length(par0);
            index_sigma = ~(index_mu);
            % sp_type = options.spline;
            %             knots = options.prev.options.prev.knots
            % spline = isfield(options, 'nknots');
            options.Leg_Points = shared.LegitimateRegion(cost, dim, options.lb, options.ub, options.search_agents, index_sigma);
            estimated_par = [];
            f_best = [];
            disp('The lower bounds and upper bounds are adapted. No optimal solution found')
        case 'fmincon'
            %Optimization using fmincon
            try
                if isempty(options.solveroptions)
                    opts = optimoptions('fmincon', 'StepTolerance', 10^(-4), 'FiniteDifferenceStepSize', 10^(-10), 'BarrierParamUpdate', 'predictor-corrector', 'Display', 'iter', 'UseParallel', options.useparallel, 'MaxFunEvals', realmax, 'MaxIter', options.maxiter, 'OutputFcn', @shared.optimdlg);
                else
                    opts = optimoptions(options.solveroptions, 'StepTolerance', 10^(-4), 'FiniteDifferenceStepSize', 10^(-10), 'BarrierParamUpdate', 'predictor-corrector', 'UseParallel', options.useparallel, 'MaxIter', options.maxiter);
                end
            catch err
                if strcmp(err.identifier, 'optimlib:optimoptions:UnmatchedParameter') %BarrierParamUpdate not allowed
                    if isempty(options.solveroptions)
                        opts = optimoptions('fmincon', 'StepTolerance', 10^(-4), 'FiniteDifferenceStepSize', 10^(-10), 'Display', 'iter', 'UseParallel', options.useparallel, 'MaxFunEvals', realmax, 'MaxIter', options.maxiter, 'OutputFcn', @shared.optimdlg);
                    else
                        opts = optimoptions(options.solveroptions, 'StepTolerance', 10^(-4), 'FiniteDifferenceStepSize', 10^(-10), 'UseParallel', options.useparallel, 'MaxIter', options.maxiter);
                    end
                else
                    rethrow(err);
                end
            end
            
            fprintf('\n Solving maximum likelihood problem (Hermite)\n');
            opts = optimoptions(opts, 'display', 'none');
            lb = options.lb;
            ub = options.ub;

            LP = varargin{2}.options.Leg_Points;
            % We first apply surrogate optimization for a very short time
            IP = struct;
            IP.X = LP(:, 1:end - 1);
            IP.Fval = LP(:, end);
            MFE = [50 100 200];
            fval = zeros(3, 1);
            x0 = zeros(3, dim);
            for i = 1:3
                options1 = optimoptions('surrogateopt', 'Display', 'off', 'MaxFunctionEvaluations', MFE(i), 'InitialPoints', IP, 'PlotFcn', {});
                [x0(i, :), fval(i)] = surrogateopt(cost, lb, ub, [], options1);
            end
            A = sortrows(LP, dim + 1);
            f_LP = A(1, dim + 1);
            B = sortrows([x0 fval], dim + 1);
            f_surr = B(1, dim + 1);
           
            if f_surr < f_LP
                % Surrogate optimization could give us better solutions (we, then seek further improvements via fmincon)
                fprintf('\nSurrogate optimization could improve the legitimate solutions. The package attempts to improve this further using fmincon solver\n');
                x0 = B(1, 1:dim);
                [estimated_par, f_best] = fmincon(cost, x0, [], [], [], [], lb, ub, [], opts);
            else
                % Surrogate optimization could not give us better solutions (we directly go for fmincon)
                fprintf('\nSurrogate optimization failed to improve the legitimate solutions. The package attempts to improve the results using fmincon solver\n');
                LP = LP(:, 1:dim);
                Estimated_par = zeros(options.search_agents, dim);
                F_best = zeros(options.search_agents, 1);
                for n = 1:options.search_agents
                    [Estimated_par(n, :), F_best(n)] = fmincon(cost, LP(n, :), [], [], [], [], lb, ub, [], opts);
                end
                A = sortrows([Estimated_par, F_best], dim + 1);
                estimated_par = A(1, 1:dim);
                f_best = A(1, dim + 1);
            end
            disp('Estimated parameters : ');
            disp(num2str(estimated_par));
            disp(['- sum of log-likelihoods : ' num2str(f_best)]);
        case 'gwo'
            %Optimization using fmincon
            fprintf('\n Solving maximum likelihood problem (Hermite)\n');
            lb = options.lb;
            ub = options.ub;
            LP = varargin{2}.options.Leg_Points;
            % We first apply surrogate optimization for a very short time
            IP = struct;
            IP.X = LP(:, 1:end - 1);
            IP.Fval = LP(:, end);
            MFE = [50 100 200];
            fval = zeros(3, 1);
            x0 = zeros(3, dim);
            for i = 1:3
                options1 = optimoptions('surrogateopt', 'Display', 'off', 'MaxFunctionEvaluations', MFE(i), 'InitialPoints', IP, 'PlotFcn', {});
                [x0(i, :), fval(i)] = surrogateopt(cost, lb, ub, [], options1);
            end
            A = sortrows(LP, dim + 1);
            f_LP = A(1, dim + 1);
            B = sortrows([x0 fval], dim + 1);
            f_surr = B(1, dim + 1);

            %We now use Improved Grey-Wolf optimizer
            if parametric
                [~, mundx] = shared.getfundimension(options.mu, []);
                M = sum(mundx);
                [estimated_par, f_best] = shared.IGWO_Parametric(dim, options.search_agents, options.maxiter, options.lb, options.ub, cost, M, options.useparallel, @shared.optimdlg);
            else
                M = options.nknots(1);
                switch options.maxiter
                    case realmax
                        [estimated_par, f_best] = shared.IGWO_Spline(dim, M, options.search_agents, options.maxiter, options.lb, options.ub, cost, options.useparallel, LP, data, @shared.optimdlg); %Results are displayed only
                    otherwise %Results are neighter displayed nor ploted
                        [estimated_par, f_best] = shared.IGWO_Spline(dim, M, options.search_agents, options.maxiter, options.lb, options.ub, cost, options.useparallel, LP, data, @shared.optimdlg);
                        if isequal(f_best, realmax)
                            disp('The algorithm did not find a legitimate solution. Increase maximum iteration number ''maxiter'' and if this does not work consider simpler models with less parameters');
                        end
                end
            end
        case 'none'
            %return the cost function
            f_best = cost(par0);
            estimated_par = par0;
            disp('Solver=''none'': cost function returned');
        otherwise
            error('Solver ''%s'' is unknown, valid options are ''fmincon'',''fminunc'',''GWO'',''PS'', ''none''', options.solver);
    end

    hfig = findobj(0, 'tag', 'optimdlg');
    if ~isempty(hfig)
        close(hfig);
    end

    %   fine_mesh = linspace(min(data), max(data), 5000);
    tim2 = toc(tim);
    if ~parametric
        %  if any(sigma(fine_mesh, estimated_par) < 0)
        %      error('The fitted sigma function has negative points, please refit or use higher lower boundaries');
        %  end     
        if strcmp(options.spline{2}, 'Approximate')
            [costSCS, mu, sigma] = shared.Create_Hermite_Cost_Spline(data, dt, options.j, options.k, knots, index_mu, options.spline{1}, 'SCS', refined);
            f_best = costSCS(estimated_par);
            fprintf('f_best (real)= %g\n', f_best);
        end

        result = struct('estimated_par', estimated_par, 'negloglik', f_best, 'datarange', [min(data), max(data)], 'dt', dt,  ...
            'mufun', mu, 'sigmafun', sigma, 'costfun', cost, 'computation_time', tim2, 'knots', {{knots(index_mu), knots(~index_mu)}},  ...
            'options', options);
    else
        %paramtric has empty knots
        %        STD1 = options.sigma(fine_mesh, estimated_par);
        %        if any(STD1 < 0)
        %            error('The fitted sigma function has negative points, please refit or use higher lower boundaries');
        %        end
        %     if numel(STD1) == 1
        %if the sigma function returns a single value we adapt it to
        %return a vector of the same size as x
        %        options.sigma = @(x,par)zeros(size(x))+options.sigma(x,par);
        %    end
        result = struct('estimated_par', estimated_par, 'negloglik', f_best, 'datarange', [min(data), max(data)], 'dt', dt,  ...
            'mufun', options.mu, 'sigmafun', options.sigma, 'computation_time', tim2, 'knots', [],  ...
            'options', options);
    end

    if ~isempty(result.estimated_par)
        if ~isempty(result.knots)
            m = numel(result.knots{1});
            knots_mu = result.knots{1};
            if isequal(options.spline{1}, 'Q')
                if m == 1
                    mu = @(x,par)par(1: m)+zeros(size(x));
                else
                    mu = @(x,par)ppval(spline(knots_mu, par(1: m)),x);
                end
                result.mufun = mu;
            end
            if isequal(options.spline{2}, 'Q')
                knots_sigma = result.knots{2};
                if isscalar(result.knots{2})
                    sigma = @(x,par)par(m+1: end)+zeros(size(x));
                else
                    sigma = @(x,par)ppval(spline(knots_sigma, par(m+1: end)),x);
                end
                result.sigmafun = sigma;
            end
        end
    end
end

function [estimated_par, f_best] = SGD(cost_SGD, data, par0, knots, index_sigma, dt, N, dim, lb, ub, learning_rate, num_epochs, batch_size, opts)
estimated_par = par0;
    for epoch = 1:num_epochs % SGD loop
        shuffled_indices = randperm(N); % Shuffle dataset
        for i = 1:batch_size:N % Iterate over mini-batches
            batch_indices = shuffled_indices(i:min(i + batch_size - 1, N));
            X0 = data(batch_indices);
            X = data(batch_indices + 1);
            cost = @(par)cost_SGD(par,X0,X);

            % MATLAB way of estimating gradient
            opts = optimoptions(opts, 'display', 'none');
            [~, ~, ~, ~, ~, grad, ~] = fmincon(cost, estimated_par, [], [], [], [], lb, ub, [], opts);
            grad = grad(:).';

            % My approach to gradient (close to MATLAB)
            % Cost0 = cost(estimated_par);
            % G = zeros(size(estimated_par));
            % h = 10^(-10);
            % for k=1:length(estimated_par)
            %     A = estimated_par;
            %     A(k) = A(k)+h;
            %     G(k) = (cost(A)-Cost0)./h;
            % end
            % grad = G;

            grad = grad(:).';
            a = estimated_par - learning_rate .* grad;
            pp_sigma = interp1(knots, a(index_sigma), 'spline', 'pp');
            if fnmin(pp_sigma, [knots(1) knots(end)]) > 0
                estimated_par = a;
                estimated_par
            end
        end
    end
end

function res = normal_fun(aspline_derivs, x)
    s = aspline_derivs;
    res = spline_val(s{1}, x);
end

function nder = num_derivatives(aspline)
    switch aspline
        case {'L', 'Approximate'}
            nder = 1;
        case 'Q'
            nder = 2;
        case {'C', 'P', 'SCS'}
            nder = 3;
        case 'constant'
            nder = 0;
        otherwise
            nder = 3;
    end
end

function res = spline_derivs(aspline, nder)
    res = cell(nder + 1, 1);
    res{1} = aspline;
    if isstruct(aspline)
        for i = 1:nder
            aspline = fnder(aspline);
            if all(aspline.coefs == 0)
                res{i + 1} = 0;
            else
                res{i + 1} = aspline;
            end
        end
    else
        %constant
        for i = 1:nder
            res{i + 1} = 0;
        end
    end
end

function res = spline_val(aspline, x)
    if ~isstruct(aspline)
        res = zeros(size(x)) + aspline;
    elseif strcmp(aspline.form, 'pp')
        res = ppval(aspline, x);
    else
        res = fnval(aspline, x);
    end
end





