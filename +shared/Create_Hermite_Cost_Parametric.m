function cost = Create_Hermite_Cost_Parametric(data, dt, J, K, mu_fn, sigma_fn, n_par, refined, ~, reconst_fraction)
% Build Hermite NLL cost function for parametric (function-handle) models.

    if ~isempty(reconst_fraction)
        frac = reconst_fraction;
        sp   = sample_stratified(data(:)', frac(1), frac(2));
        x0_all = sp(1:2:end);
        x1_all = sp(2:2:end);
    else
        raw     = data(:)';
        ok      = ~isnan(raw);
        pair_ok = ok(1:end-1) & ok(2:end);
        x0_all  = raw(1:end-1);  x0_all = x0_all(pair_ok);
        x1_all  = raw(2:end);    x1_all = x1_all(pair_ok);
    end

    h_num = 1e-5;
    cost  = @(par) hermite_nll_param(par, x0_all, x1_all, dt, J, K, mu_fn, sigma_fn, refined, h_num);
end

function nll = hermite_nll_param(par, x0, x1, dt, J, K, mu_fn, sigma_fn, refined, h)
    par  = par(:)';
    n    = numel(x0);

    mu0  = mu_fn(x0, par);
    sig0 = sigma_fn(x0, par);
    if isscalar(mu0),  mu0  = mu0  * ones(1, n); end
    if isscalar(sig0), sig0 = sig0 * ones(1, n); end
    sig0 = max(abs(sig0), 1e-300);

    z = (x1 - x0) ./ (sig0 * sqrt(dt));

    % Numerical derivatives of mu and sigma at x0 (up to 3rd order)
    dmu_fn  = @(x) mu_fn(x, par);
    dsig_fn = @(x) sigma_fn(x, par);
    Dmu  = num_derivs(dmu_fn,  x0, h, 3);   % 4 x n
    Dsig = num_derivs(dsig_fn, x0, h, 3);   % 4 x n

    % Compute EZ moments — try precomputed files first, fall back to leading-order
    additive = all(abs(Dsig(2,:)) < 1e-10 * max(1, max(abs(Dsig(2,:)))));
    EZ = compute_EZ(dt, Dmu, Dsig, J, K, additive);

    corr = hermite_correction(z, EZ, J, refined);
    corr = max(corr, 1e-300);

    nll  = sum(0.5*z.^2 + log(sig0*sqrt(dt)) + 0.5*log(2*pi) - log(corr));
end

function D = num_derivs(fn, x, h, nder)
    n  = numel(x);
    D  = zeros(nder+1, n);
    f0 = fn(x);
    if isscalar(f0), f0 = f0 * ones(1, n); end
    D(1,:) = f0;
    if nder < 1, return; end
    fp = fn(x+h);  fm = fn(x-h);
    if isscalar(fp), fp = fp*ones(1,n); end
    if isscalar(fm), fm = fm*ones(1,n); end
    D(2,:) = (fp - fm) / (2*h);
    if nder < 2, return; end
    D(3,:) = (fp - 2*f0 + fm) / (h^2);
    if nder < 3, return; end
    fp2 = fn(x+2*h);  fm2 = fn(x-2*h);
    if isscalar(fp2), fp2 = fp2*ones(1,n); end
    if isscalar(fm2), fm2 = fm2*ones(1,n); end
    D(4,:) = (fp2 - 2*fp + 2*fm - fm2) / (2*h^3);
end

function EZ = compute_EZ(dt, Dmu, Dsig, J, K, additive)
    n    = size(Dmu, 2);
    EZ   = zeros(J, n);
    sqdt = sqrt(dt);
    m0   = Dmu(1,:);
    s0   = max(abs(Dsig(1,:)), 1e-300);

    if additive
        for ndiff = [3 2 1]
            fname = sprintf('coefficients.EZ.EZ_m%d_s0_%d_%d', ndiff, J, K);
            try
                EZ_fn  = str2func(fname);
                ndiff2 = min(ndiff, size(Dmu,1)-1);
                dm_args = arrayfun(@(k) Dmu(k,:), 1:ndiff2+1, 'UniformOutput', false);
                ez_mat = EZ_fn(dt, dm_args{:}, s0);
                EZ     = ez_mat(1:J, :);
                return;
            catch
            end
        end
    else
        for ndiff = [3 2]
            fname = sprintf('coefficients.EZ.EZ_m%d_s%d_%d_%d', ndiff, ndiff, J, K);
            try
                EZ_fn  = str2func(fname);
                ndiff2 = min(ndiff, size(Dmu,1)-1);
                nds    = min(ndiff, size(Dsig,1)-1);
                dm_args = arrayfun(@(k) Dmu(k,:),  1:ndiff2+1, 'UniformOutput', false);
                ds_args = arrayfun(@(k) Dsig(k,:), 1:nds+1,   'UniformOutput', false);
                ez_mat  = EZ_fn(dt, dm_args{:}, ds_args{:});
                EZ      = ez_mat(1:J, :);
                return;
            catch
            end
        end
    end

    % Leading-order fallback
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
