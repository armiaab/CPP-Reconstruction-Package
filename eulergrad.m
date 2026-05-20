function g_LL = eulergrad(mu, sigma, dim, filename)
    %eulergrad - create gradient vector of the log-likelihood for a Langevin equation
    %This function can derive a gradient of the cost function for Euler Reconstruction of
    %Langevin equations:
    % dx= mu(x)*dt + sigma(x)*dW
    %
    %Usage:
    %h=eulergrad(mu,sigma,dim,filename)
    %mu - function handle of the function of the deterministic part:
    %     mu(x,par) should take the x and a vector of all parameters
    %sigma - function handle of the function of the stochastic part
    %    sigma(x,par)  should take the x and a vector of all parameters
    %dim - number of parameters
    %filename - name of the file for the gradient function (default 'g_LL'
    %
    %the resulting gradient function has the following arguments:
    %g_LL(x,x0,dt,par)
    if nargin < 4
        filename = 'g_LL';
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
    g_LL = gradient(LL, par);


    comment = {'Gradient vector of log-likelihood for Euler reconstruction',  'based on the following drift and diffusion functions:'...
        sprintf( '   mu=%s', func2str(mu)), sprintf('   sigma=%s', func2str(sigma))};


    %VERY IMPORTANT: the reason we SAVE LL and g_LL in below is that this way matlab optimizes them which is very important to reduce computational burden
    %LL = matlabFunction(LL, 'File', 'LL', 'vars', [x x0 par], 'Outputs', {'ll'}); %Log-Likelihood
    g_LL = matlabFunction(g_LL, 'File', filename, 'vars', [x0 x dt par], 'Outputs', {'g_ll'}, 'comment', comment); % Gradient of log-likelihood

end

