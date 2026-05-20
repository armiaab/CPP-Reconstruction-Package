function H = eulerhess_spline(mu, sigma, L, R, spline)

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

% Hessian matrix
H_mumu=zeros(mu,mu,length(mesh));
H_musigma=zeros(mu,sigma,length(mesh));
H_sigmamu=zeros(sigma,mu,length(mesh));
H_sigmasigma=zeros(sigma,sigma,length(mesh));

for i=1:mu
    for j=1:mu
        H_mumu(i,j,:)=G_mu(:,i).*G_mu(:,j);
    end
end
for i=1:mu
    for j=1:sigma
        H_musigma(i,j,:)=G_mu(:,i).*G_sigma(:,j);
    end
end
for i=1:sigma
    for j=1:mu
        H_sigmamu(i,j,:)=G_sigma(:,i).*G_mu(:,j);
    end
end
for i=1:sigma
    for j=1:sigma
        H_sigmasigma(i,j,:)=G_sigma(:,i).*G_sigma(:,j);
    end
end

H=[H_mumu H_musigma;H_sigmamu H_sigmasigma];%Q=H;
H=reshape(H,[(mu+sigma)^2 length(mesh)]);H=H.';
H = @(x) interp1(mesh,H,x);

% A=H(mesh);
% size(A)
% plot(mesh,A(:,7*16+7),'-k');hold on;plot(mesh,squeeze(Q(7,8,:)),'-r');

end
