function x=simulate(varargin)

% this code simulates from a multi-variate and multi-noise source
% stochastic differential equation dX=mu(X)dt+sigmas(X)dW where
% where mu and sigma are the deterministic and stochastic components of
% this system.


ModelType=lower(varargin{1});
switch ModelType
    case 'spline'
        SplineType=varargin{2};estimated_par=varargin{3};L=varargin{4};R=varargin{5};knots=varargin{6};dt=varargin{7};x0=varargin{8};T=varargin{9};
        M=length(knots);
        L1=knots(1);R1=knots(end);
        mesh_fine=linspace(L1,R1,ceil(R1-L1)*300);   %We consider 300 points for an interval of unit size. This should be enough for the entire interval [L R]

        if length(estimated_par)==2*M && isequal(SplineType,'SCS'), SplineType='SCS1';end
        if length(estimated_par)==2*M && isequal(SplineType,'Approximate'), SplineType='Approximate1';end

        switch SplineType
            case {'L','Approximate'}
                mu=@(x)interp1(knots,estimated_par(1:M),x,'linear','extrap');sigma=@(x)estimated_par(M+1)+0.*x;
            case 'Q'
                sp_mu=spapi(3,knots,estimated_par(1:M));y=fnval(sp_mu,mesh_fine);mu=@(x)interp1(mesh_fine,y,x,'spline','extrap');sigma=@(x)estimated_par(M+1)+0.*x;
            case 'C'
                mu=@(x)interp1(knots,estimated_par(1:M),x,'spline','extrap');sigma=@(x)estimated_par(M+1)+0.*x;
            case 'P'
                mu=@(x)interp1(knots,estimated_par(1:M),x,'pchip','extrap');sigma=@(x)estimated_par(M+1)+0.*x;
            case 'SCS'
                pp_mu=csaps(knots,estimated_par(1:M));y=ppval(pp_mu,mesh_fine);mu=@(x)interp1(mesh_fine,y,x,'spline','extrap');sigma=@(x)estimated_par(M+1)+0.*x;
            case {'LL','Approximate1'}
                mu=@(x)interp1(knots,estimated_par(1:M),x,'linear','extrap');sigma=@(x)interp1(knots,estimated_par(M+1:end),x,'linear',0);   %We do not extrapolate sigma outside its domain (we assign 0 for sigma outside its domain)
            case 'QQ'
                sp_mu=spapi(3,knots,estimated_par(1:M));y1=fnval(sp_mu,mesh_fine);mu=@(x)interp1(mesh_fine,y1,x,'spline','extrap');
                sp_sigma=spapi(3,knots,estimated_par(M+1:end));y2=fnval(sp_sigma,mesh_fine);sigma=@(x)interp1(mesh_fine,y2,x,'spline',0);    %We do not extrapolate sigma outside its domain (we assign 0 for sigma outside its domain)
            case 'CC'
                mu=@(x)interp1(knots,estimated_par(1:M),x,'spline','extrap');
                sigma=@(x)interp1(knots,estimated_par(M+1:end),x,'spline',0);   %We do not extrapolate sigma outside its domain (we assign 0 for sigma outside its domain)
            case 'SCS1'
                pp_mu=csaps(knots,estimated_par(1:M));y1=ppval(pp_mu,mesh_fine);mu=@(x)interp1(mesh_fine,y1,x,'spline','extrap');
                pp_sigma=csaps(knots,estimated_par(M+1:end));y2=ppval(pp_sigma,mesh_fine);sigma=@(x)interp1(mesh_fine,y2,x,'spline',0);      %We do not extrapolate sigma outside its domain (we assign 0 for sigma outside its domain)
            case 'CP'
                mu=@(x)interp1(knots,estimated_par(1:M),x,'spline','extrap');sigma=@(x)interp1(knots,estimated_par(M+1:end),x,'pchip',0);    %We do not extrapolate sigma outside its domain (we assign 0 for sigma outside its domain)
            case 'PP'
                mu=@(x)interp1(knots,estimated_par(1:M),x,'pchip','extrap');sigma=@(x)interp1(knots,estimated_par(M+1:end),x,'pchip',0);     %We do not extrapolate sigma outside its domain (we assign 0 for sigma outside its domain)
        end
    case 'parametric'
        L=varargin{2};R=varargin{3};mu=varargin{4};sigma=varargin{5};dt=varargin{6};x0=varargin{7};T=varargin{8};
end
if isempty(L)
    L=-realmax;
end

if isempty(R)
    R=realmax;
end

x0=x0(:);
dim=length(x0);
x=zeros(dim,T);x(:,1)=x0;

m=size(sigma(x0),2);   %m is the number of noise sources
Noise=randn(m,T);

for i=2:T
    a=x(:,i-1)+mu(x(:,i-1)).*dt+sqrt(dt).*sigma(x(:,i-1))*Noise(:,i);
    idx=a<L;a(idx)=L(idx);
    idx=a>R;a(idx)=R(idx);
    x(:,i)=a;
end
x=x.';
end