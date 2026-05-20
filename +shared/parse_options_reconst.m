function [options, parametric, index_mu, knots] = parse_options_reconst(data, args, defaults_param, defaults_spline)
% Parse name-value options for hermite_reconstruction / legitimate_points.
% Handles the three-step chain: euler -> legitimate_points -> hermite (legpoints) -> hermite (fmincon).

    % --- parse raw key-value pairs from args ---
    kv = struct();
    i = 1;
    while i <= numel(args) - 1
        key = args{i};
        if ischar(key) || isstring(key)
            kv.(lower(char(key))) = args{i+1};
        end
        i = i + 2;
    end

    prev       = [];
    prev_range = defaults_param.prev_range;
    if isfield(kv, 'prev'),       prev       = kv.prev;       end
    if isfield(kv, 'prev_range'), prev_range = kv.prev_range; end

    % --- build base options ---
    if ~isempty(prev) && isstruct(prev) && isfield(prev, 'options')
        % Inherit from a previous result's options struct (steps 2 and 3 of the chain)
        options = prev.options;
        % Apply any explicit overrides from the current call
        override = rmfield_safe(kv, {'prev', 'prev_range'});
        fnames = fieldnames(override);
        for fi = 1:numel(fnames)
            options.(fnames{fi}) = override.(fnames{fi});
        end
        % If the previous result has fitted parameters, tighten lb/ub around them
        if isfield(prev, 'estimated_par') && ~isempty(prev.estimated_par)
            ep = prev.estimated_par(:)';
            new_lb = ep - prev_range * abs(ep);
            new_ub = ep + prev_range * abs(ep);
            if isfield(options, 'lb') && ~isempty(options.lb)
                new_lb = max(new_lb, options.lb(:)');
            end
            if isfield(options, 'ub') && ~isempty(options.ub)
                new_ub = min(new_ub, options.ub(:)');
            end
            options.lb = new_lb;
            options.ub = new_ub;
        end
    else
        % First call: pick appropriate defaults then apply all kv
        if isfield(kv, 'mu') && ~isempty(kv.mu)
            options = defaults_param;
        elseif isfield(kv, 'nknots') || isfield(kv, 'knots')
            options = defaults_spline;
        elseif ~isempty(prev) && isstruct(prev) && isfield(prev, 'knots') ...
                && iscell(prev.knots) && ~isempty(prev.knots)
            options = defaults_spline;
        elseif ~isempty(prev) && isstruct(prev) && isfield(prev, 'options') ...
                && isfield(prev.options, 'mu') && ~isempty(prev.options.mu)
            options = defaults_param;
        else
            options = defaults_spline;
        end
        fnames = fieldnames(kv);
        for fi = 1:numel(fnames)
            f = fnames{fi};
            if ~ismember(f, {'prev', 'prev_range'})
                options.(f) = kv.(f);
            end
        end
        % Inherit knot positions and function handles from prev result
        if ~isempty(prev) && isstruct(prev)
            if isfield(prev, 'options')
                po = prev.options;
                if isfield(po, 'mu')    && ~isempty(po.mu),    options.mu    = po.mu;    end
                if isfield(po, 'sigma') && ~isempty(po.sigma), options.sigma = po.sigma; end
                if isfield(po, 'spline') && ~isempty(po.spline) && ~isfield(kv,'spline')
                    options.spline = po.spline;
                end
                if isfield(po, 'nknots') && ~isfield(kv,'nknots'), options.nknots = po.nknots; end
                if isfield(po, 'l') && ~isempty(po.l) && ~isfield(kv,'l'), options.l = po.l; end
                if isfield(po, 'r') && ~isempty(po.r) && ~isfield(kv,'r'), options.r = po.r; end
            end
            if isfield(prev, 'dt') && ~isfield(kv,'dt'), options.dt = prev.dt; end
            % lb/ub from prev estimated_par (first call, prev is euler result)
            if isfield(prev, 'estimated_par') && ~isempty(prev.estimated_par)
                ep = prev.estimated_par(:)';
                base_lb = [];  base_ub = [];
                if isfield(prev, 'options') && isfield(prev.options, 'lb'), base_lb = prev.options.lb(:)'; end
                if isfield(prev, 'options') && isfield(prev.options, 'ub'), base_ub = prev.options.ub(:)'; end
                new_lb = ep - prev_range * abs(ep);
                new_ub = ep + prev_range * abs(ep);
                if ~isempty(base_lb), new_lb = max(new_lb, base_lb); end
                if ~isempty(base_ub), new_ub = min(new_ub, base_ub); end
                options.lb = new_lb;
                options.ub = new_ub;
            end
            % Store knot positions from euler result
            if isfield(prev, 'knots') && iscell(prev.knots) && numel(prev.knots) == 2
                options.knots_mu  = prev.knots{1}(:)';
                options.knots_sig = prev.knots{2}(:)';
            end
        end
    end

    % --- normalise spline field to a 2-element cell ---
    if isfield(options, 'spline')
        sp = options.spline;
        if ischar(sp) && ~isempty(sp)
            options.spline = {sp(1), sp(min(2,end))};
        end
    else
        options.spline = {'C', 'C'};
    end

    % --- ensure required scalar fields exist ---
    if ~isfield(options,'j'),               options.j               = 3;      end
    if ~isfield(options,'k'),               options.k               = 4;      end
    if ~isfield(options,'refined'),         options.refined         = true;   end
    if ~isfield(options,'analytic_diffs'),  options.analytic_diffs  = true;   end
    if ~isfield(options,'reconst_fraction'),options.reconst_fraction= [];     end
    if ~isfield(options,'solver'),          options.solver          = 'fmincon'; end
    if ~isfield(options,'search_agents'),   options.search_agents   = 5;      end
    if ~isfield(options,'maxiter'),         options.maxiter         = realmax; end
    if ~isfield(options,'useparallel'),     options.useparallel     = false;   end
    if ~isfield(options,'solveroptions'),   options.solveroptions   = [];     end
    if ~isfield(options,'l'),               options.l               = [];     end
    if ~isfield(options,'r'),               options.r               = [];     end
    if ~isfield(options,'mu'),              options.mu              = [];     end
    if ~isfield(options,'sigma'),           options.sigma           = [];     end

    % --- determine parametric vs spline ---
    parametric = ~isempty(options.mu);

    % --- compute knots vector and index_mu ---
    if parametric
        n_total = numel(options.lb);
        n_mu_dim = shared.getfundimension(options.mu, []);
        if n_mu_dim >= n_total, n_mu_dim = n_total; end
        index_mu = [true(1, n_mu_dim), false(1, n_total - n_mu_dim)];
        knots = [];
    else
        % Spline mode — recover knot positions
        if isfield(options,'knots_mu') && ~isempty(options.knots_mu)
            km = options.knots_mu(:)';
            ks = options.knots_sig(:)';
        elseif ~isempty(prev) && isstruct(prev) && isfield(prev,'knots') ...
                && iscell(prev.knots) && numel(prev.knots)==2
            km = prev.knots{1}(:)';
            ks = prev.knots{2}(:)';
            options.knots_mu  = km;
            options.knots_sig = ks;
        else
            % Generate from data
            n_mu  = options.nknots(1);
            n_sig = options.nknots(2);
            cdata = data(:);
            cdata = cdata(~isnan(cdata));
            if ~isempty(options.l) && isfinite(options.l), cdata = cdata(cdata >= options.l); end
            if ~isempty(options.r) && isfinite(options.r), cdata = cdata(cdata <= options.r); end
            mg = 0.02 * (max(cdata) - min(cdata));
            lo = min(cdata)-mg;  hi = max(cdata)+mg;
            km = linspace(lo, hi, n_mu);
            ks = linspace(lo, hi, n_sig);
            options.knots_mu  = km;
            options.knots_sig = ks;
        end
        knots = [km, ks];
        n_mu  = numel(km);
        n_sig = numel(ks);
        index_mu = [true(1, n_mu), false(1, n_sig)];
        if ~isfield(options,'nknots') || isempty(options.nknots)
            options.nknots = [n_mu, n_sig];
        end
    end
end

function s = rmfield_safe(s, fnames)
    for i = 1:numel(fnames)
        if isfield(s, fnames{i}), s = rmfield(s, fnames{i}); end
    end
end
