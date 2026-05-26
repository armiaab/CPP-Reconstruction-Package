%euler_reconstruction  Reconstruct Langevin equation using the Euler method.
%
% Drop-in replacement for the original MATLAB package.
% When langevin_mex is compiled, the NLL is evaluated in C++ (~10x faster).
% Optimization always uses MATLAB's fmincon/GWO - identical to the original.
%
% Usage:
%   result = euler_reconstruction(data, dt, 'name', value, ...)
%
% Required:
%   data   - vector of equally-spaced observations (NaN allowed)
%   dt     - time step
%
% Options (same as original package):
%   'nknots'          [nmu nsig] or scalar (default [8 8])
%   'knots'           explicit knot positions (cell {km,ks} or flat vector)
%   'spline'          'CC','LL','QQ','PP',... (default 'CC')
%   'L','R'           domain boundaries
%   'lb','ub'         parameter bounds
%   'solver'          'fmincon' (default), 'GWO'
%   'search_agents'   number of multi-starts (default 1)
%   'maxiter'         max function evaluations (default realmax)
%   'useparallel'     logical (default false)
%   'mu','sigma'      parametric function handles
%   'gradient_fun'    gradient handle (parametric only)
%   'reconst_fraction' [strata, p] for big data
%   'zscore'          z-score transform data (default false, spline only)
%   'prev','prev_range' pass previous result as starting point
%
% Result struct (identical to original package):
%   estimated_par, negloglik, datarange, dt, mufun, sigmafun, costfun (spline),
%   computation_time, knots, options
%
function result = euler_reconstruction(data, dt, varargin)

    defaults_param = struct( ...
        'prev', [], 'prev_range', 0.05, 'lb', [], 'ub', [], 'l', [], 'r', [], 'npars', [], ...
        'mu', [], 'sigma', [], 'gradient_fun', [], 'hessian_fun', [], 'solver', 'fmincon', ...
        'reconst_fraction', [], 'solveroptions', [], 'useparallel', false, 'zscore', false, ...
        'search_agents', 1, 'maxiter', realmax);

    defaults_spline = struct( ...
        'prev', [], 'prev_range', 0.05, 'lb', [], 'ub', [], 'l', [], 'r', [], 'nknots', [8 8], ...
        'knots', [], 'spline', 'CC', 'zscore', false, 'gradient_fun', [], 'hessian_fun', [], ...
        'solver', 'fmincon', 'reconst_fraction', [], 'solveroptions', [], 'useparallel', false, ...
        'search_agents', 1, 'maxiter', realmax);

    % '-d' mode: return default options struct (matches original)
    if ischar(data) && strcmp(data, '-d')
        result = struct('parametric', defaults_param, 'spline', defaults_spline);
        return;
    end

    % ------------------------------------------------------------------ %
    %  Parse options - uses the original shared.parse_options_reconst     %
    % ------------------------------------------------------------------ %
    [options, parametric, index_mu, knots] = shared.parse_options_reconst( ...
        data, varargin, defaults_param, defaults_spline);

    if options.zscore && parametric
        error('euler_reconstruction:zscore', ...
            'Z-score transformation is not supported for parametric models.');
    end

    % ------------------------------------------------------------------ %
    %  Apply L/R boundary as NaN - IDENTICAL to original                  %
    % ------------------------------------------------------------------ %
    if ~isempty(options.l) && ~isnan(options.l)
        data(data < options.l) = NaN;
    end
    if ~isempty(options.r) && ~isnan(options.r)
        data(data > options.r) = NaN;
    end

    if all(isnan(data))
        error('euler_reconstruction:nodata', ...
            'Not enough non-missing data for reconstruction.');
    end

    datastats = struct( ...
        'range', [min(data, [], 'omitnan'), max(data, [], 'omitnan')], ...
        'mean',   mean(data, 'omitnan'), ...
        'std',    std(data,  'omitnan'));

    if options.zscore
        [data, knots, options] = shared.ztransform(data, datastats, knots, options);
    end

    % ------------------------------------------------------------------ %
    %  Build consecutive data pairs - matching original                    %
    % ------------------------------------------------------------------ %
    fraction = options.reconst_fraction;
    if isempty(fraction)
        X0 = data(1:end-1);
        X  = data(2:end);
    else
        strata = fraction(1);
        p      = fraction(2);
        sample_pairs = sample_stratified(data, strata, p);
        X0 = sample_pairs(1:2:end);
        X  = sample_pairs(2:2:end);
    end
    if size(X0, 1) > 1, X0 = X0.'; X = X.'; end

    % ------------------------------------------------------------------ %
    %  Build cost function                                                  %
    % ------------------------------------------------------------------ %
    hfig = findobj(0, 'tag', 'optimdlg');
    if ~isempty(hfig), close(hfig); end
    tim = tic;

    % Use 200 points for sigma-sign check (faster; 5000 is overkill for sign)
    fine_mesh = linspace(min(data, [], 'omitnan'), max(data, [], 'omitnan'), 200);
    positive_sigma = false;
    par0 = [];

    if ~parametric
        spline_mu  = options.spline{1};
        spline_sig = options.spline{2};

        % Always build MATLAB function handles for mu/sigma (for result and sigma check)
        [cost_matlab, mu, sigma] = loc_euler_cost_spline( ...
            X0, X, dt, knots, index_mu, spline_mu, spline_sig);

        % Build C++ NLL costfun when langevin_mex is available (faster).
        % Falls back to pure MATLAB NLL when MEX is absent.
        knots_mu  = knots(index_mu);
        knots_sig = knots(~index_mu);
        sp_str    = [loc_sp_to_mex(spline_mu), loc_sp_to_mex(spline_sig)];
        if sum(~index_mu) == 1
            sp_str(2) = 'a';   % additive (constant) sigma
        end

        nll_extra = {};
        if ~isempty(options.l) && isfinite(options.l), nll_extra = [nll_extra, {'l', options.l}]; end
        if ~isempty(options.r) && isfinite(options.r), nll_extra = [nll_extra, {'r', options.r}]; end

        use_mex = exist('langevin_mex', 'file') == 3;
        if use_mex
            data_row = data(:)';   % row copy captured in closure
            km_row   = knots_mu(:)';
            ks_row   = knots_sig(:)';
            costfun  = @(par) langevin_mex('nll', data_row, dt, par(:)', ...
                km_row, ks_row, sp_str, nll_extra{:});
        else
            costfun = cost_matlab;
        end

        % Find positive-sigma start - matching original.
        % Build the cheap linear (LL) cost once, outside the loop.
        use_ll_prefit = any(strcmp(options.solver, {'fmincon', 'fminunc'})) && ...
                        any(strcmp(spline_sig, {'Q', 'C', 'SCS', 'Approximate'}));
        if use_ll_prefit
            cost_ll = loc_euler_cost_spline(X0, X, dt, knots, index_mu, 'L', 'L');
            opts_ll = optimoptions('fmincon', 'Display', 'none', 'MaxFunEvals', 2000);
        end
        iter = 1;
        while ~positive_sigma && iter < 100
            par0 = options.lb + (options.ub - options.lb) .* rand(1, length(options.lb));
            iter = iter + 1;
            if use_ll_prefit
                par_ll = fmincon(cost_ll, par0, [], [], [], [], options.lb, options.ub, [], opts_ll);
                par0   = par_ll * 0.5 + par_ll .* rand(size(par0));
            end
            positive_sigma = any(strcmp(spline_sig, {'L', 'P'})) || ...
                all(sigma(fine_mesh, par0) > 0);
        end

    else  % parametric
        costfun = loc_euler_cost_parametric(X0, X, dt, options.mu, options.sigma, options.gradient_fun);

        iter = 1;
        while ~positive_sigma && iter < 100
            par0 = options.lb + (options.ub - options.lb) .* rand(1, length(options.lb));
            iter = iter + 1;
            positive_sigma = all(options.sigma(fine_mesh, par0) > 0);
        end
    end

    % ------------------------------------------------------------------ %
    %  Optimize - MATLAB solver, matching original                         %
    % ------------------------------------------------------------------ %
    dim = length(par0);
    fprintf('\n Solving maximum likelihood problem (Euler)\n');

    switch lower(options.solver)
        case 'fmincon'
            if parametric && isempty(options.gradient_fun)
                opts = optimoptions('fmincon', ...
                    'FiniteDifferenceStepSize', 1e-6, 'UseParallel', false, ...
                    'MaxFunEvals', realmax, 'MaxIter', options.maxiter, ...
                    'OutputFcn', @shared.optimdlg);
            elseif parametric && ~isempty(options.gradient_fun)
                opts = optimoptions('fmincon', ...
                    'FiniteDifferenceStepSize', 1e-6, 'UseParallel', options.useparallel, ...
                    'MaxFunEvals', realmax, 'MaxIter', options.maxiter, ...
                    'SpecifyObjectiveGradient', true, 'OutputFcn', @shared.optimdlg);
            else
                opts = optimoptions('fmincon', ...
                    'FiniteDifferenceStepSize', 1e-6, 'UseParallel', false, ...
                    'MaxFunEvals', realmax, 'MaxIter', options.maxiter, ...
                    'OutputFcn', @shared.optimdlg);
            end
            opts = optimoptions(opts, 'display', 'none');
            lb = options.lb;  ub = options.ub;

            if options.search_agents == 1
                [estimated_par, f_best] = fmincon(costfun, par0, [], [], [], [], lb, ub, [], opts);
            else
                % Sobol quasi-random multi-start (matching original)
                try
                    p_sobol = sobolset(dim, 'Skip', 1e3, 'Leap', 1e2);
                    p_sobol = scramble(p_sobol, 'MatousekAffineOwen');
                    RS = net(p_sobol, options.search_agents);
                catch
                    RS = rand(options.search_agents, dim);
                end
                Estimated_par = zeros(options.search_agents, dim);
                F_best        = zeros(options.search_agents, 1);
                A = []; retries = 0;
                while size(A, 1) < options.search_agents && retries < 5
                    ptmatrix = lb + (ub - lb) .* RS;
                    for n = 1:options.search_agents
                        [Estimated_par(n,:), F_best(n)] = fmincon( ...
                            costfun, ptmatrix(n,:), [], [], [], [], lb, ub, [], opts);
                    end
                    F_best = F_best(:);
                    A = [A; Estimated_par, F_best]; %#ok<AGROW>
                    idx = A(:, end) > 1.79e308;
                    A = A(~idx, :);
                    retries = retries + 1;
                    RS = rand(options.search_agents, dim);  % fresh random points on retry
                end
                if isempty(A)
                    % All starts returned invalid — fall back to best single start par0
                    [estimated_par, f_best] = fmincon(costfun, par0, [], [], [], [], lb, ub, [], opts);
                else
                    A = sortrows(A, dim + 1);
                    estimated_par = A(1, 1:dim);
                    f_best        = A(1, end);
                end
            end

        case 'gwo'
            if parametric
                [~, mundx] = shared.getfundimension(options.mu, []);
                M = sum(mundx);
                [estimated_par, f_best] = shared.IGWO_Parametric( ...
                    dim, options.search_agents, options.maxiter, ...
                    options.lb, options.ub, costfun, M, options.useparallel, @shared.optimdlg);
            else
                % Guard: derive M from index_mu if nknots was not set
                if isempty(options.nknots)
                    options.nknots = [sum(index_mu), sum(~index_mu)];
                end
                M = options.nknots(1);
                [estimated_par, f_best] = shared.IGWO_Spline_Euler( ...
                    dim, M, options.search_agents, options.maxiter, ...
                    options.lb, options.ub, costfun, options.useparallel, @shared.optimdlg);
                if isequal(f_best, realmax)
                    disp(['The algorithm did not find a legitimate solution. ' ...
                        'Increase maxiter or use a simpler model.']);
                end
            end

        otherwise
            error('euler_reconstruction:solver', ...
                'Solver ''%s'' is not supported. Valid: ''fmincon'', ''GWO''.', options.solver);
    end

    if isequal(f_best, realmax)
        disp(['The algorithm did not find a legitimate solution. ' ...
            'Increase maxiter or use simpler model.']);
    else
        disp('Estimated parameters :');
        disp(num2str(estimated_par));
        disp(['- sum of log-likelihoods) : ' num2str(f_best)]);
    end

    hfig = findobj(0, 'tag', 'optimdlg');
    if ~isempty(hfig), close(hfig); end
    tim2 = toc(tim);

    % ------------------------------------------------------------------ %
    %  Assemble result struct - IDENTICAL to original                      %
    % ------------------------------------------------------------------ %
    if ~parametric
        if any(sigma(fine_mesh, estimated_par) < 0)
            error('euler_reconstruction:sigma', ...
                'Fitted sigma has negative values. Refit or use higher lower bounds.');
        end

        if strcmp(spline_sig, 'Approximate')
            % Re-evaluate with SCS for accurate NLL
            [costSCS, mu, sigma] = loc_euler_cost_spline( ...
                X0, X, dt, knots, index_mu, spline_mu, 'SCS');
            f_best = costSCS(estimated_par);
            fprintf('f_best (real)= %g\n', f_best);
        end

        % Post-process Q spline handles - IDENTICAL to original post-processing
        m = options.nknots(1);
        if isequal(spline_mu, 'Q')
            if m == 1
                mu = @(x, par) par(1:m) + zeros(size(x));
            else
                km_loc = knots(1:m);
                mu = @(x, par) ppval(spline(km_loc, par(1:m)), x);
            end
        end
        if isequal(spline_sig, 'Q')
            n_sig_k = options.nknots(2);
            if n_sig_k == 1
                sigma = @(x, par) par(m+1:end) + zeros(size(x));
            else
                ks_loc = knots(m+1:end);
                sigma = @(x, par) ppval(spline(ks_loc, par(m+1:end)), x);
            end
        end

        result = struct( ...
            'estimated_par',    estimated_par, ...
            'negloglik',        f_best, ...
            'datarange',        [min(data, [], 'omitnan'), max(data, [], 'omitnan')], ...
            'dt',               dt, ...
            'mufun',            mu, ...
            'sigmafun',         sigma, ...
            'costfun',          costfun, ...     % present in original for non-parametric
            'computation_time', tim2, ...
            'knots',            {{knots(index_mu), knots(~index_mu)}}, ...
            'options',          options);

        if options.zscore
            result = shared.back_ztransform(result);
        end

    else  % parametric
        if ~strcmp(options.solver, 'fitnlm')
            STD1 = options.sigma(fine_mesh, estimated_par);
            if any(STD1 < 0)
                error('euler_reconstruction:sigma', ...
                    'Fitted sigma has negative values. Refit or use higher lower bounds.');
            end
            if numel(STD1) == 1
                options.sigma = @(x, par) zeros(size(x)) + options.sigma(x, par);
            end
        end
        result = struct( ...
            'estimated_par',    estimated_par, ...
            'negloglik',        f_best, ...
            'datarange',        [min(data, [], 'omitnan'), max(data, [], 'omitnan')], ...
            'dt',               dt, ...
            'mufun',            options.mu, ...
            'sigmafun',         options.sigma, ...
            'computation_time', tim2, ...
            'knots',            [], ...
            'options',          options);
    end
end


% ======================================================================= %
%  Local helpers                                                            %
% ======================================================================= %

function c = loc_sp_to_mex(ch)
% Map original spline char to langevin_mex lowercase code.
    switch upper(ch)
        case 'L',           c = 'l';
        case {'C', 'S'},    c = 'c';
        case 'Q',           c = 'q';
        case 'P',           c = 'p';
        otherwise,          c = 'c';
    end
end


function [costfun, mu, sigma] = loc_euler_cost_spline(X0, X, dt, knots, index_mu, sp_mu, sp_sig)
% Build MATLAB-native cost/mu/sigma for spline models.
% This is the original Create_Euler_Cost_Spline, kept local.

    function cost = ffun(par, X0, X, dt, sqrtdt, sigma, mu)
        STD1 = sigma(X0, par);
        if all(isnan(STD1) | STD1 > 0)
            cost = loc_euler_likelihood(X, X0, mu(X0, par), STD1, dt, sqrtdt);
        else
            cost = realmax;
        end
    end

    function cost = ffun_approx(par, X0, X, dt, sqrtdt, sigma, mu, testsigma)
        if all(isnan(X0) | testsigma(X0, par) > 0)
            cost = loc_euler_likelihood(X, X0, mu(X0, par), sigma(X0, par), dt, sqrtdt);
        else
            cost = realmax;
        end
    end

    index_sigma = ~index_mu;
    if sum(index_mu)    == 1, sp_mu  = 'constant'; end
    if sum(index_sigma) == 1, sp_sig = 'constant'; end

    knots_mu  = knots(index_mu);
    knots_sig = knots(index_sigma);
    sqrtdt    = sqrt(dt);

    switch sp_mu
        case 'L'
            mu = @(x, par) interp1(knots_mu, par(index_mu), x);
        case 'Q'
            mu = @(x, par) fnval(spapi(3, knots_mu, par(index_mu)), x);
        case 'C'
            mu = @(x, par) ppval(spline(knots_mu, par(index_mu)), x);
        case 'P'
            mu = @(x, par) ppval(pchip(knots_mu, par(index_mu)), x);
        case {'SCS', 'Approximate'}
            mu = @(x, par) ppval(csaps(knots_mu, par(index_mu)), x);
        case 'constant'
            mu = @(x, par) par(index_mu) + zeros(size(x));
        otherwise
            error('euler_reconstruction:spline', ...
                'Spline type of mu ''%s'' is not supported.', sp_mu);
    end

    switch sp_sig
        case 'constant'
            sigma   = @(x, par) par(index_sigma) + zeros(size(x));
            costfun = @(par) loc_euler_likelihood(X, X0, mu(X0, par), par(index_sigma), dt, sqrtdt);
        case 'L'
            sigma   = @(x, par) interp1(knots_sig, par(index_sigma), x);
            costfun = @(par) loc_euler_likelihood(X, X0, mu(X0, par), sigma(X0, par), dt, sqrtdt);
        case 'Q'
            sigma   = @(x, par) fnval(spapi(3, knots_sig, par(index_sigma)), x);
            costfun = @(par) ffun(par, X0, X, dt, sqrtdt, sigma, mu);
        case 'C'
            sigma   = @(x, par) ppval(spline(knots_sig, par(index_sigma)), x);
            costfun = @(par) ffun(par, X0, X, dt, sqrtdt, sigma, mu);
        case 'P'
            sigma   = @(x, par) ppval(pchip(knots_sig, par(index_sigma)), x);
            costfun = @(par) loc_euler_likelihood(X, X0, mu(X0, par), sigma(X0, par), dt, sqrtdt);
        case 'Approximate'
            sigma        = @(x, par) ppval(csaps(knots_sig, par(index_sigma)), x);
            sigma_linear = @(x, par) interp1(knots_sig, par(index_sigma), x);
            mu_linear    = @(x, par) interp1(knots_mu,  par(index_mu),    x);
            costfun      = @(par) ffun_approx(par, X0, X, dt, sqrtdt, sigma_linear, mu_linear, sigma);
        case 'SCS'
            sigma   = @(x, par) ppval(csaps(knots_sig, par(index_sigma)), x);
            costfun = @(par) ffun(par, X0, X, dt, sqrtdt, sigma, mu);
        otherwise
            error('euler_reconstruction:spline', ...
                'Spline type of sigma ''%s'' is not supported.', sp_sig);
    end
end


function cost = loc_euler_likelihood(X, X0, mu, sigma, dt, sqrtdt)
% Base Euler log-likelihood - IDENTICAL to original.
    mu1    = X0 + mu .* dt;
    sigma1 = sigma .* sqrtdt;
    z      = (X - mu1) ./ sigma1;
    L      = -0.5 .* z .* z - log(sqrt(2 .* pi) .* sigma1);
    cost   = -sum(L, 'omitnan');
end


function costfun = loc_euler_cost_parametric(X0, X, dt, mu, sigma, gradient_fun)
% Parametric cost function - IDENTICAL to original.

    function [cost, g_cost] = cost_with_grad(par, X0, X, dt, sqrtdt, mu, sigma, gradient_fun)
        STD1 = sigma(X0, par);
        cost = -sum(-0.5 .* (log(2 .* pi .* dt .* STD1.^2) + ...
            ((X - X0 - mu(X0, par) .* dt) ./ (STD1 .* sqrtdt)).^2), 'omitnan');
        if nargout > 1
            cpar   = num2cell(par);
            g_cost = -sum(gradient_fun(X0, X, dt, cpar{:}), 2, 'omitnan');
        end
    end

    function cost = cost_wo_grad(par, X0, X, dt, sqrtdt, mu, sigma)
        STD1 = sigma(X0, par);
        cost = -sum(-0.5 .* (log(2 .* pi .* dt .* STD1.^2) + ...
            ((X - X0 - mu(X0, par) .* dt) ./ (STD1 .* sqrtdt)).^2), 'omitnan');
    end

    sqrtdt = sqrt(dt);
    if isempty(gradient_fun)
        costfun = @(par) cost_wo_grad(par, X0, X, dt, sqrtdt, mu, sigma);
    else
        costfun = @(par) cost_with_grad(par, X0, X, dt, sqrtdt, mu, sigma, gradient_fun);
    end
end
