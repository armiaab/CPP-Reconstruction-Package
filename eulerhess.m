function h_LL = eulerhess(mu, sigma, dim, filename)

    if nargin < 4
        filename = 'h_LL';
    end
    if nargin == 3 && ischar(dim)
        filename = dim;
        dim = [];
    end
    
    if nargin < 3 || isempty(dim)
       dim=shared.getfundimension(mu, sigma);
    end

    syms x0 x dt
    par = sym('par', [1, dim]);
    LL = -1 / 2 * (log(2 * pi * dt * sigma(x0, par)^2) + ((x - x0 - mu(x0, par) * dt) / (sigma(x0, par) * sqrt(dt)))^2);
    h_LL = hessian(LL, par);
    h_LL = reshape(h_LL, [], 1);
    
    comment = {'Hessian matrix of log-likelihood for Euler reconstruction',  'based on the following equations:'...
        sprintf( '   mu=%s', func2str(mu)), sprintf('   sigma=%s', func2str(sigma))};

    %VERY IMPORTANT: the reason we SAVE LL and g_LL in below is that this way matlab optimizes them which is very important to reduce computational burden
    h_LL = matlabFunction(h_LL, 'File', filename, 'vars', [x0 x dt par], 'Outputs', {'h_ll'}, 'comment', comment); % Hessian matrix of log-likelihood
end








