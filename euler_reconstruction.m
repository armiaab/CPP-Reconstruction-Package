function result = euler_reconstruction(data, dt, varargin)
%EULER_RECONSTRUCTION  C++ MEX wrapper matching the original API.
%
%  result = euler_reconstruction(data, dt, 'name', value, ...)
%
%  Supported options (spline model only):
%    'nKnots'        [n_mu n_sig]          default [8 8]
%    'knots'         explicit knot vector  (alternative to nKnots)
%    'spline'        two-char string       default 'CC'
%                    'L'=linear  'C'=cubic  'P'=pchip
%                    'Q'/'SCS'/'Approximate' mapped to cubic
%    'L'             left domain bound
%    'R'             right domain bound
%    'lb'            lower bounds vector
%    'ub'            upper bounds vector
%    'search_agents' number of multi-starts   default 1
%    'maxiter'       max optimizer iterations default 10000
%    'solver'        'slsqp'(default) | 'lbfgs' | 'isres'

    opts = parse_opts(varargin);

    if ~isempty(opts.mu)
        error('euler_reconstruction (C++ wrapper): parametric models are not supported. Use nKnots + spline instead.');
    end

    n_mu  = opts.nknots(1);
    n_sig = opts.nknots(end);

    sp_str    = upper(opts.spline);
    sp_mu_ch  = sp_str(1);
    sp_sig_ch = sp_str(min(2, length(sp_str)));
    sp_mu_mex = matlab2mex_sp(sp_mu_ch);
    if n_sig == 1
        sp_mex = sp_mu_mex;
    else
        sp_mex = [sp_mu_mex matlab2mex_sp(sp_sig_ch)];
    end

    mex_args = {'nknots', [n_mu n_sig], 'spline', sp_mex, ...
                'nstarts', max(1, opts.search_agents), ...
                'maxiter', min(round(opts.maxiter), 1000000), ...
                'solver',  opts.solver};

    if ~isempty(opts.l) && isfinite(opts.l)
        mex_args = [mex_args, {'l', opts.l}];
    end
    if ~isempty(opts.r) && isfinite(opts.r)
        mex_args = [mex_args, {'r', opts.r}];
    end
    if ~isempty(opts.lb)
        mex_args = [mex_args, {'lb', opts.lb(:)'}];
    end
    if ~isempty(opts.ub)
        mex_args = [mex_args, {'ub', opts.ub(:)'}];
    end

    if ~isempty(opts.knots)
        km = opts.knots(1:n_mu);
        ks = opts.knots(n_mu+1:end);
        mex_args = [mex_args, {'knots_mu', km, 'knots_sig', ks}];
    end

    [par, nll, tsec] = langevin_mex(data, dt, mex_args{:});

    if isempty(opts.knots)
        km = auto_knots(data, opts.l, opts.r, n_mu);
        if n_sig == 1
            ks = km(1);
        else
            ks = auto_knots(data, opts.l, opts.r, n_sig);
        end
    end

    mu_fn  = make_spline_fn(sp_mu_ch,  km,  1,        n_mu);
    if n_sig == 1
        sig_fn = @(x, p) p(n_mu+1) .* ones(size(x));
    else
        sig_fn = make_spline_fn(sp_sig_ch, ks, n_mu+1, n_mu+n_sig);
    end

    sp_label_mu  = spline_label(sp_mu_ch);
    sp_label_sig = spline_label(sp_sig_ch);
    if n_sig == 1, sp_label_sig = 'constant'; end

    opts_out.solver        = ['C++ ' upper(opts.solver)];
    opts_out.spline        = {sp_label_mu, sp_label_sig};
    opts_out.l             = opts.l;
    opts_out.r             = opts.r;
    opts_out.nknots        = [n_mu n_sig];
    opts_out.knots         = [];
    opts_out.lb            = opts.lb;
    opts_out.ub            = opts.ub;
    opts_out.search_agents = opts.search_agents;
    opts_out.maxiter       = opts.maxiter;
    opts_out.useparallel   = false;

    d_clean = data(isfinite(data));
    result.estimated_par    = par(:)';
    result.negloglik        = nll;
    result.datarange        = [min(d_clean), max(d_clean)];
    result.dt               = dt;
    result.mufun            = mu_fn;
    result.sigmafun         = sig_fn;
    result.computation_time = tsec;
    result.knots            = {km, ks};
    result.options          = opts_out;
end

function c = matlab2mex_sp(ch)
    switch upper(ch)
        case 'L',                       c = 'l';
        case {'C','Q','SCS','A'},        c = 'c';
        case 'P',                       c = 'p';
        otherwise,                      c = 'c';
    end
end

function lbl = spline_label(ch)
    switch upper(ch)
        case 'L',  lbl = 'L';
        case 'C',  lbl = 'C';
        case 'P',  lbl = 'P';
        case 'Q',  lbl = 'C';
        otherwise, lbl = 'C';
    end
end

function fn = make_spline_fn(sp_ch, knots, i1, i2)
    switch upper(sp_ch)
        case 'L'
            fn = @(x, p) interp1(knots, p(i1:i2), x, 'linear', 'extrap');
        case 'P'
            fn = @(x, p) ppval(pchip(knots, p(i1:i2)), x);
        otherwise
            fn = @(x, p) ppval(spline(knots, p(i1:i2)), x);
    end
end

function km = auto_knots(data, L, R, nk)
    if ~isempty(L) && isfinite(L) && ~isempty(R) && isfinite(R)
        lo = L; hi = R;
    else
        d  = data(isfinite(data));
        lo = min(d); hi = max(d);
    end
    m  = (hi - lo) * 0.02;
    lo = lo - m; hi = hi + m;
    km = linspace(lo, hi, nk);
end

function opts = parse_opts(args)
    opts.nknots        = [8 8];
    opts.knots         = [];
    opts.spline        = 'CC';
    opts.l             = [];
    opts.r             = [];
    opts.lb            = [];
    opts.ub            = [];
    opts.search_agents = 1;
    opts.maxiter       = 10000;
    opts.solver        = 'slsqp';
    opts.mu            = [];
    opts.sigma         = [];

    i = 1;
    while i <= length(args)
        key = args{i};
        if ~ischar(key) && ~isstring(key), i = i+1; continue; end
        key = lower(char(key));
        val = args{i+1};
        switch key
            case 'nknots',         opts.nknots        = val;
            case 'knots',          opts.knots         = val;
            case 'spline',         opts.spline        = char(val);
            case {'l','left'},     opts.l             = val;
            case {'r','right'},    opts.r             = val;
            case 'lb',             opts.lb            = val;
            case 'ub',             opts.ub            = val;
            case 'search_agents',  opts.search_agents = val;
            case 'maxiter',        opts.maxiter       = val;
            case 'solver',         opts.solver        = lower(char(val));
            case 'mu',             opts.mu            = val;
            case 'sigma',          opts.sigma         = val;
        end
        i = i + 2;
    end

    if length(opts.nknots) == 1
        opts.nknots = [opts.nknots opts.nknots];
    end
    if isempty(opts.l),  opts.l = NaN; end
    if isempty(opts.r),  opts.r = NaN; end
end
