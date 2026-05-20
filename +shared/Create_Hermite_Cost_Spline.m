function [cost, mu_fn, sigma_fn] = Create_Hermite_Cost_Spline(data, dt, J, K, knots, index_mu, sp_mu, sp_sig, refined, reconst_fraction)
% Build the Hermite NLL cost function for spline models.

    km = knots(index_mu);
    ks = knots(~index_mu);
    n_mu  = numel(km);
    n_sig = numel(ks);
    additive = (n_sig == 1);

    mu_fn    = make_sp_fn(sp_mu, km, 1, n_mu);
    if additive
        sigma_fn = @(x, par) par(n_mu+1) * ones(1, numel(x));
    else
        sigma_fn = make_sp_fn(sp_sig, ks, n_mu+1, n_mu+n_sig);
    end

    ndiff_mu  = sp2ndiff(sp_mu);
    ndiff_sig = sp2ndiff(sp_sig);
    if additive, ndiff_sig = 0; end

    dmu_fn  = make_sp_deriv_fn(sp_mu, km, 1, n_mu, ndiff_mu);
    if additive
        dsig_fn = @(x, par) [par(n_mu+1)*ones(1,numel(x))];
    else
        dsig_fn = make_sp_deriv_fn(sp_sig, ks, n_mu+1, n_mu+n_sig, ndiff_sig);
    end

    % Build data pairs
    raw     = data(:)';
    ok      = ~isnan(raw);
    pair_ok = ok(1:end-1) & ok(2:end);
    x0_all  = raw(1:end-1);  x0_all = x0_all(pair_ok);
    x1_all  = raw(2:end);    x1_all = x1_all(pair_ok);

    if ~isempty(reconst_fraction)
        frac = reconst_fraction;
        sp   = sample_stratified([x0_all; x1_all], frac(1), frac(2));
        x0_all = sp(1, :);
        x1_all = sp(2, :);
    end

    EZ_fn = load_EZ_fn(ndiff_mu, ndiff_sig, J, K);

    cost = @(par) hermite_nll_spline(par, x0_all, x1_all, dt, ...
                    dmu_fn, dsig_fn, additive, EZ_fn, ndiff_mu, ndiff_sig, J, refined);
end

% -----------------------------------------------------------------------

function nll = hermite_nll_spline(par, x0, x1, dt, dmu_fn, dsig_fn, additive, EZ_fn, ndiff_mu, ndiff_sig, J, refined)
    par  = par(:)';
    n    = numel(x0);
    Dmu  = dmu_fn(x0, par);   % (ndiff_mu+1) x n
    if additive
        sig0 = dsig_fn(x0, par);  % 1 x n
        sig0 = sig0(1,:);
        Dsig = [sig0; zeros(ndiff_sig, n)];
    else
        Dsig = dsig_fn(x0, par);  % (ndiff_sig+1) x n
        sig0 = Dsig(1,:);
    end
    sig0 = max(abs(sig0), 1e-300);
    z    = (x1 - x0) ./ (sig0 * sqrt(dt));

    EZ   = compute_EZ_vec(EZ_fn, dt, Dmu, Dsig, ndiff_mu, ndiff_sig, J, n);
    corr = hermite_correction(z, EZ, J, refined);
    corr = max(corr, 1e-300);

    nll  = sum(0.5*z.^2 + log(sig0*sqrt(dt)) + 0.5*log(2*pi) - log(corr));
end

function EZ = compute_EZ_vec(EZ_fn, dt, Dmu, Dsig, ndiff_mu, ndiff_sig, J, n)
    EZ = zeros(J, n);
    if ~isempty(EZ_fn)
        dm_args = cell(1, ndiff_mu+1);
        for di = 0:ndiff_mu
            dm_args{di+1} = Dmu(di+1, :);
        end
        if ndiff_sig == 0
            ds_args = {Dsig(1,:)};
        else
            ds_args = cell(1, ndiff_sig+1);
            for di = 0:ndiff_sig
                ds_args{di+1} = Dsig(di+1, :);
            end
        end
        try
            ez_mat = EZ_fn(dt, dm_args{:}, ds_args{:});
            EZ = ez_mat(1:J, :);
            return;
        catch
        end
    end
    % Fallback: leading-order approximation
    sqdt = sqrt(dt);
    m0   = Dmu(1,:);
    s0   = max(abs(Dsig(1,:)), 1e-300);
    EZ(1,:) = m0 .* sqdt ./ s0;
    if J >= 2, EZ(2,:) = 1 + (m0./s0).^2 * dt; end
    if J >= 3, EZ(3,:) = 3 * m0 .* dt ./ (s0 * sqdt); end
end

function c = hermite_correction(z, EZ, J, refined)
    c   = ones(1, numel(z));
    ez1 = EZ(1,:);
    ez2 = EZ(2,:);
    ez3 = EZ(3,:);
    c1  = ez1;
    c2  = ez2 - 1;
    c3  = ez3 - 3*ez1;
    He1 = z;
    He2 = z.^2 - 1;
    He3 = z.^3 - 3*z;
    if refined
        c = 1 + c1.*He1 + (c2/2).*He2 + (c3/6).*He3;
    else
        c = 1 + (c3/6).*He3;
    end
end

function nd = sp2ndiff(sp)
    switch upper(sp)
        case {'L','APPROXIMATE'}, nd = 1;
        case 'Q',                 nd = 2;
        otherwise,                nd = 3;  % C, P, SCS, CC, QQ, etc.
    end
end

function EZ_fn = load_EZ_fn(ndiff_mu, ndiff_sig, J, K)
    if ndiff_sig == 0
        fname = sprintf('EZ_m%d_s0_%d_%d', ndiff_mu, J, K);
    else
        fname = sprintf('EZ_m%d_s%d_%d_%d', ndiff_mu, ndiff_sig, J, K);
    end
    full_name = ['coefficients.EZ.' fname];
    if ~isempty(which(full_name))
        EZ_fn = str2func(full_name);
    else
        EZ_fn = [];
    end
end

function fn = make_sp_fn(sp, knots_k, i1, i2)
    switch upper(sp)
        case 'L'
            fn = @(x,p) interp1(knots_k, p(i1:i2), x(:)', 'linear', 'extrap');
        case 'P'
            if numel(knots_k)==1
                fn = @(x,p) p(i1)*ones(1,numel(x));
            else
                fn = @(x,p) ppval(pchip(knots_k, p(i1:i2)), x(:)');
            end
        otherwise
            if numel(knots_k)==1
                fn = @(x,p) p(i1)*ones(1,numel(x));
            else
                fn = @(x,p) ppval(spline(knots_k, p(i1:i2)), x(:)');
            end
    end
end

function fn = make_sp_deriv_fn(sp, knots_k, i1, i2, ndiff)
    switch upper(sp)
        case 'L'
            fn = @(x,p) linear_deriv_matrix(knots_k, p(i1:i2), x(:)', ndiff);
        case 'P'
            fn = @(x,p) pchip_deriv_matrix(knots_k, p(i1:i2), x(:)', ndiff);
        otherwise
            fn = @(x,p) spline_deriv_matrix(knots_k, p(i1:i2), x(:)', ndiff);
    end
end

function D = spline_deriv_matrix(knots_k, vals, x, ndiff)
    n  = numel(x);
    D  = zeros(ndiff+1, n);
    if numel(knots_k) == 1
        D(1,:) = vals(1);
        return;
    end
    pp = spline(knots_k, vals);
    D(1,:) = ppval(pp, x);
    for d = 1:ndiff
        pp = fnder_pp(pp);
        if isstruct(pp)
            D(d+1,:) = ppval(pp, x);
        else
            D(d+1,:) = pp * ones(1, n);
        end
    end
end

function D = linear_deriv_matrix(knots_k, vals, x, ndiff)
    n  = numel(x);
    D  = zeros(ndiff+1, n);
    D(1,:) = interp1(knots_k, vals, x, 'linear', 'extrap');
    if ndiff >= 1 && numel(knots_k) > 1
        slopes = diff(vals(:)') ./ diff(knots_k(:)');
        for xi = 1:n
            seg = find(knots_k <= x(xi), 1, 'last');
            if isempty(seg), seg = 1; end
            seg = min(seg, numel(slopes));
            D(2, xi) = slopes(seg);
        end
    end
end

function D = pchip_deriv_matrix(knots_k, vals, x, ndiff)
    n  = numel(x);
    D  = zeros(ndiff+1, n);
    if numel(knots_k) == 1
        D(1,:) = vals(1);
        return;
    end
    pp = pchip(knots_k, vals);
    D(1,:) = ppval(pp, x);
    for d = 1:ndiff
        pp = fnder_pp(pp);
        if isstruct(pp)
            D(d+1,:) = ppval(pp, x);
        else
            D(d+1,:) = pp * ones(1, n);
        end
    end
end

function pp2 = fnder_pp(pp)
    if ~isstruct(pp), pp2 = 0; return; end
    [breaks, coefs, l, k, d] = unmkpp(pp);
    if k <= 1, pp2 = 0; return; end
    new_coefs = coefs(:, 1:k-1) .* repmat(k-1:-1:1, l, 1);
    pp2 = mkpp(breaks, new_coefs, d);
end
