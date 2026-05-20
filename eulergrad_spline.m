function G = eulergrad_spline(mu, sigma, L, R, spline)

knots_mu = linspace(L,R,mu); 
knots_sigma = linspace(L,R,sigma);
par_mu = rand(1,mu);
par_sigma = rand(1,sigma);
mesh = linspace(L, R, 5000);

spline_mu = spline(1);
spline_sigma = spline(2);

% Gradient vector
h=10^(-6);h=1;  % for splne types other than 'pchip' h can be any value but for 'pchip' it should be small
G_mu=zeros(length(mesh),mu);
for i=1:mu
    par=par_mu;
    par(i)=par(i)+h;   % perturb the i-th parameter (i-th knot value)
    g=@(x)(interp1(knots_mu,par,x,spline_mu)-interp1(knots_mu,par_mu,x,spline_mu))./h;
    G_mu(:,i)=g(mesh);
end

G_sigma=zeros(length(mesh),sigma);
for i=1:sigma
    par=par_sigma;
    par(i)=par(i)+h;   % perturb the i-th parameter (i-th knot value)
    g=@(x)(interp1(knots_sigma,par,x,spline_sigma)-interp1(knots_sigma,par_sigma,x,spline_sigma))./h;
    G_sigma(:,i)=g(mesh);
end

G=[G_mu G_sigma];
G = @(x) interp1(mesh,G,x);
end
