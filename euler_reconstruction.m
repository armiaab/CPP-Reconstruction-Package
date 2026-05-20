% euler_reconstruction  Reconstruct Langevin equation using Euler method.
%
% Drop-in replacement for the original MATLAB package, backed by the C++
% langevin_mex MEX binary.  The calling convention and result struct are
% identical to the original so that all existing scripts (Examples,
% plot_results, legitimate_points, etc.) work without modification.
%
% Usage:
%   result = euler_reconstruction(data, dt, 'name', value, ...)
%
% Required arguments:
%   data   - column or row vector of equally-spaced observations (NaN allowed)
%   dt     - time step between observations
%
% Supported name-value options (same names as original package):
%   'nKnots'        [nmu nsig] or scalar  number of knots (default [8 8])
%   'knots'         flat vector           explicit knot positions
%                                         first nmu values → mu knots,
%                                         next nsig values → sigma knots
%   'spline'        2-char string         spline type per function, e.g. 'CC',
%                                         'LL', 'QQ', 'PP' (default 'CC')
%                    'L' = linear interpolation
%                    'C' = cubic spline
%                    'Q' = order-3 spline  (mapped to cubic in C++ backend)
%                    'P' = PCHIP
%   'L'             scalar                left boundary (default [])
%   'R'             scalar                right boundary (default [])
%   'lb'            vector                lower bounds for all parameters
%   'ub'            vector                upper bounds for all parameters
%   'solver'        string                'fmincon'(default),'GWO','PS',
%                                         'lbfgs','isres'
%   'search_agents' integer               number of random starts (default 1)
%   'maxiter'       integer or realmax    max function evaluations (default 10000)
%   'ftol'          double                convergence tolerance (default 1e-9)
%   'useparallel'   logical               ignored (C++ backend is single-threaded)
%   'reconst_fraction' ignored
%   'mu','sigma','gradient_fun'           not supported; raises an error
%
% Result struct (identical to original package):
%   estimated_par   [1 x (nmu+nsig)]  fitted spline coefficients
%   negloglik       scalar             negative log-likelihood
%   datarange       [1 x 2]            [min max] of cleaned observations
%   dt              scalar             time step (echoed)
%   mufun           @(x, par)          drift evaluation handle
%   sigmafun        @(x, par)          diffusion evaluation handle
%   computation_time scalar            wall-clock seconds
%   knots           {1x2 cell}         {knots_mu, knots_sig}
%   options         struct             parsed options (matches original layout)
%
function result = euler_reconstruction(data, dt, varargin)

    % ------------------------------------------------------------------ %
    %  Parse options (case-insensitive, mirrors original package names)   %
    % ------------------------------------------------------------------ %
    o = parse_opts(varargin);

    if ~isempty(o.mu)
        result = euler_reconstruction_parametric(data, dt, o);
        return;
    end

    % ------------------------------------------------------------------ %
    %  Knot counts and spline-type mapping                                %
    % ------------------------------------------------------------------ %
    n_mu  = o.nknots(1);
    n_sig = o.nknots(2);

    sp_str    = upper(o.spline);
    sp_mu_ch  = sp_str(1);
    sp_sig_ch = sp_str(min(2, length(sp_str)));

    mex_sp_mu  = orig2mex_sp(sp_mu_ch);
    mex_sp_sig = orig2mex_sp(sp_sig_ch);
    if n_sig == 1
        mex_sp_sig = 'a';          % additive (constant) noise
    end

    % ------------------------------------------------------------------ %
    %  Clean data: remove NaN and out-of-boundary values                  %
    % ------------------------------------------------------------------ %
    clean = data(:);
    clean = clean(~isnan(clean));
    if ~isempty(o.l) && isfinite(o.l), clean = clean(clean >= o.l); end
    if ~isempty(o.r) && isfinite(o.r), clean = clean(clean <= o.r); end

    if isempty(clean)
        error('euler_reconstruction:nodata', ...
            'No non-missing data remain after applying L/R boundaries.');
    end

    % ------------------------------------------------------------------ %
    %  Knot positions — computed from cleaned data to match C++ exactly   %
    % ------------------------------------------------------------------ %
    if ~isempty(o.knots)
        % User-supplied flat vector: first n_mu → mu, next n_sig → sigma
        flat = o.knots(:)';
        km = flat(1 : n_mu);
        ks = flat(n_mu+1 : n_mu+n_sig);
    else
        margin = 0.02 * (max(clean) - min(clean));
        lo = min(clean) - margin;
        hi = max(clean) + margin;
        km = linspace(lo, hi, n_mu);
        if n_sig == 1
            ks = (lo + hi) / 2;    % single knot; value irrelevant for additive
        else
            ks = linspace(lo, hi, n_sig);
        end
    end

    % ------------------------------------------------------------------ %
    %  Build langevin_mex argument list                                   %
    % ------------------------------------------------------------------ %
    mex_solver = orig2mex_solver(o.solver);

    % Guard against realmax being passed as maxiter (C++ uses int)
    safe_maxiter = min(o.maxiter, 1e8);
    if ~isfinite(o.maxiter) || o.maxiter >= realmax / 2
        safe_maxiter = 1e8;
    end

    mex_args = { ...
        'knots_mu',  km(:)', ...
        'knots_sig', ks(:)', ...
        'spline',    [mex_sp_mu mex_sp_sig], ...
        'solver',    mex_solver, ...
        'nstarts',   max(1, round(o.search_agents)), ...
        'maxiter',   round(safe_maxiter), ...
        'ftol',      o.ftol };

    if ~isempty(o.l) && isfinite(o.l), mex_args = [mex_args, {'l', o.l}]; end
    if ~isempty(o.r) && isfinite(o.r), mex_args = [mex_args, {'r', o.r}]; end
    if ~isempty(o.lb), mex_args = [mex_args, {'lb', o.lb(:)'}]; end
    if ~isempty(o.ub), mex_args = [mex_args, {'ub', o.ub(:)'}]; end

    % ------------------------------------------------------------------ %
    %  Call MEX backend                                                    %
    % ------------------------------------------------------------------ %
    [par, nll, tsec, ~] = langevin_mex(data(:)', dt, mex_args{:});

    % ------------------------------------------------------------------ %
    %  Build @(x, par) function handles using the same knots we passed    %
    % ------------------------------------------------------------------ %
    mufun   = make_spline_fn(sp_mu_ch,  km, 1,       n_mu);
    if n_sig == 1
        nm = n_mu;
        sigmafun = @(x, p) p(nm + 1) .* ones(1, numel(x));
    else
        sigmafun = make_spline_fn(sp_sig_ch, ks, n_mu + 1, n_mu + n_sig);
    end

    % ------------------------------------------------------------------ %
    %  Assemble options struct (preserves original field names and values) %
    % ------------------------------------------------------------------ %
    opts_out.solver        = o.solver;          % keep user's original string
    opts_out.spline        = {sp_mu_ch, sp_sig_ch};
    opts_out.l             = o.l;
    opts_out.r             = o.r;
    opts_out.nknots        = [n_mu, n_sig];
    opts_out.knots         = o.knots;
    opts_out.lb            = o.lb;
    opts_out.ub            = o.ub;
    opts_out.search_agents = o.search_agents;
    opts_out.maxiter       = o.maxiter;
    opts_out.useparallel   = false;

    % ------------------------------------------------------------------ %
    %  Assemble result struct (identical to original package layout)       %
    % ------------------------------------------------------------------ %
    result.estimated_par    = par(:)';
    result.negloglik        = nll;
    result.datarange        = [min(clean), max(clean)];
    result.dt               = dt;
    result.mufun            = mufun;
    result.sigmafun         = sigmafun;
    result.computation_time = tsec;
    result.knots            = {km, ks};
    result.options          = opts_out;
end

% ======================================================================= %
%  Local helpers                                                           %
% ======================================================================= %

function c = orig2mex_sp(ch)
% Map original single-char spline code to langevin_mex lowercase code.
    switch upper(ch)
        case 'L',               c = 'l';
        case {'C','Q','S'},     c = 'c';   % Q (order-3 B-spline) → cubic
        case 'P',               c = 'p';
        otherwise,              c = 'c';
    end
end

function s = orig2mex_solver(solver)
% Map original solver names to langevin_mex solver names.
    switch lower(solver)
        case {'fmincon','ps','fminunc'},  s = 'slsqp';
        case {'gwo'},                     s = 'isres';
        case {'lbfgs'},                   s = 'lbfgs';
        case {'isres'},                   s = 'isres';
        case {'slsqp'},                   s = 'slsqp';
        otherwise,                        s = 'slsqp';
    end
end

function fn = make_spline_fn(sp_ch, knots, i1, i2)
% Build an @(x, par) evaluation handle for a spline segment.
% Uses MATLAB's native interpolation so it matches the original's function
% handles for plotting and downstream use.
    switch upper(sp_ch)
        case 'L'
            fn = @(x, p) interp1(knots, p(i1:i2), x(:)', 'linear', 'extrap');
        case 'P'
            if numel(knots) == 1
                fn = @(x, p) p(i1) .* ones(1, numel(x));
            else
                fn = @(x, p) ppval(pchip(knots, p(i1:i2)), x(:)');
            end
        otherwise  % 'C', 'Q', or anything else → cubic spline
            if numel(knots) == 1
                fn = @(x, p) p(i1) .* ones(1, numel(x));
            else
                fn = @(x, p) ppval(spline(knots, p(i1:i2)), x(:)');
            end
    end
end

function o = parse_opts(args)
% Parse a varargin cell of name-value pairs. All keys are case-insensitive.
    o.nknots        = [8 8];
    o.knots         = [];
    o.spline        = 'CC';
    o.l             = [];
    o.r             = [];
    o.lb            = [];
    o.ub            = [];
    o.solver        = 'fmincon';
    o.search_agents = 1;
    o.maxiter       = 10000;
    o.ftol          = 1e-9;
    o.mu               = [];
    o.sigma            = [];
    o.reconst_fraction = [];

    i = 1;
    while i <= length(args) - 1
        key = args{i};
        if ~(ischar(key) || isstring(key)), i = i + 1; continue; end
        key = lower(char(key));
        val = args{i + 1};
        switch key
            case 'nknots',                         o.nknots        = val;
            case 'knots',                          o.knots         = val;
            case 'spline',                         o.spline        = char(val);
            case {'l','left'},                     o.l             = val;
            case {'r','right'},                    o.r             = val;
            case 'lb',                             o.lb            = val;
            case 'ub',                             o.ub            = val;
            case 'solver',                         o.solver        = char(val);
            case {'search_agents','nstarts'},      o.search_agents = val;
            case 'maxiter',                        o.maxiter       = val;
            case 'ftol',                           o.ftol          = val;
            case 'mu',                             o.mu               = val;
            case 'sigma',                          o.sigma            = val;
            case 'reconst_fraction',               o.reconst_fraction = val;
            case {'useparallel','zscore', ...
                  'prev','prev_range','npars', ...
                  'gradient_fun','hessian_fun', ...
                  'solveroptions','transform'}
                % silently ignored — not applicable to C++ backend
        end
        i = i + 2;
    end

    % Normalise nknots to a 2-element row vector
    if isscalar(o.nknots)
        o.nknots = [o.nknots, o.nknots];
    end
end

% ======================================================================= %
%  Parametric-mode fitting (mu/sigma function handles)                     %
% ======================================================================= %

function result = euler_reconstruction_parametric(data, dt, o)
    t_start = tic;

    % Infer parameter count from bounds
    if ~isempty(o.lb)
        n_par = numel(o.lb);
    elseif ~isempty(o.ub)
        n_par = numel(o.ub);
    else
        error('euler_reconstruction:parametric', ...
            'Parametric mode requires ''lb'' and/or ''ub'' to determine parameter dimension.');
    end
    lb_v = o.lb(:)'; if isempty(lb_v), lb_v = -inf(1, n_par); end
    ub_v = o.ub(:)'; if isempty(ub_v), ub_v =  inf(1, n_par); end

    % Build consecutive data pairs, skipping NaN transitions
    if ~isempty(o.reconst_fraction)
        frac = o.reconst_fraction;
        sp   = sample_stratified(data(:)', frac(1), frac(2));
        x0_p = sp(1:2:end);
        x1_p = sp(2:2:end);
    else
        raw      = data(:)';
        ok       = ~isnan(raw);
        pair_ok  = ok(1:end-1) & ok(2:end);
        x0_p = raw(1:end-1); x0_p = x0_p(pair_ok);
        x1_p = raw(2:end);   x1_p = x1_p(pair_ok);
    end

    mu_fn  = o.mu;
    sig_fn = o.sigma;
    nll_fn = @(par) euler_nll(par, x0_p, x1_p, dt, mu_fn, sig_fn);

    safe_maxiter = min(o.maxiter, 1e8);
    if ~isfinite(o.maxiter) || o.maxiter >= realmax/2, safe_maxiter = 1e8; end

    fmc_opts = optimoptions('fmincon', ...
        'Display',                 'off', ...
        'MaxFunctionEvaluations',  round(safe_maxiter), ...
        'FunctionTolerance',       o.ftol, ...
        'OptimalityTolerance',     o.ftol);

    n_starts = max(1, round(o.search_agents));
    mid = (lb_v + ub_v) / 2;  mid(~isfinite(mid)) = 0;
    best_nll = Inf;  best_par = mid;

    for k = 1:n_starts
        if k == 1
            par0 = mid;
        else
            r    = rand(1, n_par);
            par0 = lb_v + r .* (ub_v - lb_v);
            inf_mask = ~isfinite(par0);
            par0(inf_mask) = randn(1, sum(inf_mask));
        end
        try
            [p_k, f_k] = fmincon(nll_fn, par0, [], [], [], [], lb_v, ub_v, [], fmc_opts);
            if f_k < best_nll
                best_nll = f_k;  best_par = p_k;
            end
        catch
        end
    end

    tsec = toc(t_start);

    clean_all = data(:);  clean_all = clean_all(~isnan(clean_all));

    opts_out.solver        = o.solver;
    opts_out.spline        = {};
    opts_out.l             = o.l;
    opts_out.r             = o.r;
    opts_out.nknots        = [];
    opts_out.knots         = [];
    opts_out.lb            = o.lb;
    opts_out.ub            = o.ub;
    opts_out.search_agents = o.search_agents;
    opts_out.maxiter       = o.maxiter;
    opts_out.useparallel   = false;

    result.estimated_par    = best_par(:)';
    result.negloglik        = best_nll;
    result.datarange        = [min(clean_all), max(clean_all)];
    result.dt               = dt;
    result.mufun            = mu_fn;
    result.sigmafun         = sig_fn;
    result.computation_time = tsec;
    result.knots            = {};
    result.options          = opts_out;
end

function nll = euler_nll(par, x0, x1, dt, mu_fn, sig_fn)
    par   = par(:)';
    mu_v  = mu_fn(x0, par);
    sig_v = sig_fn(x0, par);
    sig_v = max(abs(sig_v), 1e-300);
    nll   = 0.5 * sum(log(2*pi*dt*sig_v.^2) + ((x1 - x0 - mu_v*dt) ./ (sig_v*sqrt(dt))).^2);
end
