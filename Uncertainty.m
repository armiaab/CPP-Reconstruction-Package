function [hess,err_hess,err_par]=Uncertainty(varargin)

ModelType=varargin{1};
switch ModelType(1)
    case 'Spline'
        estimated_par=varargin{2};h=varargin{3};data=varargin{4};J=varargin{5};K=varargin{6};dt=varargin{7};knots=varargin{8};mesh_fine=varargin{9};SplineType=varargin{10};
    case 'Parametric'
        estimated_par=varargin{2};h=varargin{3};data=varargin{4};J=varargin{5};K=varargin{6};dt=varargin{7};par=varargin{8};mu=varargin{9};sigma=varargin{10};
        assume(par,'real');
end

% NOTE1: The typical data (meaning a single dataset) should be supplied as an N*d matrix where d is the spatial (or state) dimension and N is the temporal length. Clearly, N>d almost always. So, we interpret max(size(data)) as the temporal length and main(size(data)) as the state length (numbreer of state variables). If data is a square matrix we consider size(data,1) as the temporal length  
% NOTE2: The replicate data (meaning multiple datasets together) should be suplied as a cell array. The format of each cell (i.e., each replicate) should follow the points in NOTE1 

if iscell(data)
    A=data{1};
    D=min(size(A));  %state (or spatial) dimension
else
    D=min(size(data));  
end

if D==1 && isequal(ModelType(1),'Parametric')
    x=argnames(mu);assume(x,'real');
end

if D>1 && isequal(ModelType(1),'Parametric')
    V=argnames(mu);
    assume(V,'real');
    d=length(V);   %state (or spatial) dimension
end

if (isequal(ModelType(1),'Spline') || isequal(ModelType(1),'Parametric')) && D==1   %the user wants to use the whole data (Univarate)
    if iscell(data)   %Replicate datasets
        data1=[];
        for i=1:length(data)
            data1=[data1 data{i}(:).' nan];
        end
        data=data1;
    else   %Typical datasets
        data=data(:).';
    end
    N=length(data);
    X0=data(1:end-1);X=data(2:end);
    diff_data=diff(data);
end

if D>1 && nargin==10  %the user wants to use the whole data (Multivariate)
    if iscell(data)   %Replicate datasets
        data1=[];
        for i=1:length(data)
            if size(data{i},1)<size(data{i},2)
                data{i}=data{i}.';
            end
            data1=[data1;data{i};nan(1,d)];
        end
        data=data1;
        N=size(data,1);   %N is temporal length
        XX0=data(1:N-1,:);XX=data(2:N,:);
    else   %Typical datasets
        if size(data,1)<size(data,2)
            data=data.';
        end
        N=size(data,1);
        XX0=data(1:N-1,:);XX=data(2:N,:);
    end
end

if (isequal(ModelType(1),'Spline') && nargin==11 && D==1)  %the user wants to use a (random) fraction of data for a spline model (Univariate)
    sample_pairs=varargin{11};
    X0=sample_pairs(1:2:end);
    X=sample_pairs(2:2:end);
    N=length(X0)+1;
    diff_data=X-X0;
end
if (isequal(ModelType(1),'Parametric') && nargin==12 && D==1)  %the user wants to use a (random) fraction of data for a parametric model (Univariate)
    sample_pairs=varargin{12};
    X0=sample_pairs(1:2:end);
    X=sample_pairs(2:2:end);
    N=length(X0)+1;
    diff_data=X-X0;
end

if D>1 && nargin==12   %the user wants to use a (random) fraction of data for a parametric model (Multivariate)
    sample_pairs=varargin{12};
    XX0=sample_pairs(1:2:end,:);
    XX=sample_pairs(2:2:end,:);
    N=size(XX0,1)+1;
end


switch ModelType(1)
    case 'Spline'
        M=length(knots);
        switch ModelType(2)
            case 'Hermite'
                switch ModelType(3)
                    case 'Additive noise1'
                        [~,ETA,H,~,~]=Coefficients_Spline(J,K,ModelType(3),SplineType);
                        
                        %Calculating the objective function
                        cost=@(par)Cost_Spline1(par,dt,diff_data,N,X0,knots,M,H,ETA,SplineType);
                    case 'Additive noise2'
                        if J<3,disp('J should be equal or bigger than 3');return;end
                        [EZ,~,H,H_coeff,Coeff]=Coefficients_Spline(J,K,ModelType(3),SplineType);
                        
                        %Calculating the objective function
                        cost=@(par)Cost_Spline2(par,dt,J,N,X0,X,EZ,H,H_coeff,Coeff,knots,M,SplineType);
                    case 'Multiplicative noise'
                        if J<3,disp('J should be equal or bigger than 3');return;end
                        [EZ,~,H,H_coeff,Coeff]=Coefficients_Spline(J,K,ModelType(3),SplineType);                        
                        
                        %Calculating the objective function
                        switch SplineType
                            case {'CP','PP','LL'}
                                cost=@(par)Cost_Spline3(par,dt,J,N,X0,X,EZ,H,H_coeff,Coeff,knots,mesh_fine,M,SplineType);
                            case {'CC','SCS','QQ'}
                                cost=@(par)Cost_Spline4(par,dt,J,N,X0,X,EZ,H,H_coeff,Coeff,knots,mesh_fine,M,SplineType,[]);
                            case 'Approximate'
                                dx=knots(2)-knots(1);
                                E=linspace(knots(1)-dx/2,knots(end)+dx/2,M+1);
                                [w,~]=histcounts(data,E);
                                cost=@(par)Cost_Spline4(par,dt,J,N,X0,X,EZ,H,H_coeff,Coeff,knots,mesh_fine,M,SplineType,w);
                        end
                end
            case 'Euler'
                %Calculating the objective function
                cost=@(par)Cost_Spline_Euler(par,dt,X0,X,knots,M,ModelType(3),SplineType);
        end
    case 'Parametric'
        switch ModelType(2)
            case 'Hermite'
                switch ModelType(3)
                    case 'Simple noise'
                        [d_muY,Gam,H,eta,sigma,~]=Coefficients_Parametric(J,K,dt,mu,sigma,par,ModelType(3));
                        
                        %Calculating the objective function
                        cost=@(par)Cost_Simple(par,data,X,N,dt,sigma,d_muY,Gam,H,eta);
                    case 'Complex noise'
                        [d_mu,sigma,d_sigma,Z,EZ,H,H_coeff,Coeff,~]=Coefficients_Parametric(J,K,dt,mu,sigma,par,ModelType(3));
                        
                        %Calculating the objective function
                        cost=@(par)Cost_Complex(par,J,Coeff,X0,X,N,dt,sigma,d_mu,d_sigma,EZ,H,H_coeff,Z);

                    case 'Additive noise'
                        [d_mu,EZ,H,H_coeff,Coeff,~]=Coefficients_Parametric(J,K,dt,mu,sigma,par,ModelType(3));
                        
                        %Calculating the objective function
                        cost=@(par)Cost_Additive(par,J,Coeff,X0,X,N,dt,d_mu,EZ,H,H_coeff);
                end
            case 'Euler'
                syms x0
                LL=-1/2*(log(2*pi*dt*sigma(x0)^2)+((x-x0-mu(x0)*dt)/(sigma(x0)*sqrt(dt)))^2);  %Log-Likelihood
                %VERY IMPORTANT: the reason we SAVE LL in bellow is that this way matlab optimizes it which is very important to reduce computational burden
                LL=matlabFunction(LL,'File','LL','vars',[x x0 par],'Outputs',{'ll'});   
                
                %Calculating the objective function
                cost=@(par)Cost_Parametric_Euler(par,X,X0,LL);   
            case 'MultivariateEuler'
                assume(par,'real');
                syms X X0 [1 length(V)]
                p=sym2cell(X0);
                
                mu_X0=subs(mu,V,X0).';mu_X0=mu_X0(p{:});
                sigma_X0=subs(sigma,V,X0);sigma_X0=sigma_X0(p{:});
                D=sigma_X0*sigma_X0.';
                det_D=simplify(det(D));  
                inv_D=simplify(inv(D));
                
                %Log-likelihood function 
                S=1/dt*(X-X0-mu_X0*dt)*inv_D*(X-X0-mu_X0*dt).';
                LL=-1/2*(d*log(2*pi*dt)+log(det_D)+S);
                
                %VERY IMPORTANT: the reason we SAVE the following quntities in bellow is that this way matlab optimizes them which is very important to reduce computational burden
                LL=matlabFunction(LL,'File','LL','vars',{X,X0,par},'Outputs',{'ll'});   %Log-Likelihood
                
                %Calculating the objective function
                cost=@(par)Cost_Parametric_Euler_MultiVariate(par,XX,XX0,LL);     
        end
end

%Calculating the Hessian matrix at estimated values of parameters
fprintf('\nOnly in case of ''Hermite reconstruction'': If this lasts rather long with no results then you should re-run this code using a smaller step size h \n');
fprintf('Calculating the Hessian matrix (observed Fisher information matrix) at estimated parameters...\n\n');
[hess,err_hess]=hessian(cost,estimated_par,h);

%the uncertainty of the parameters in terms of standard deviation
err_par=sqrt(diag(inv(hess)));

if (norm(imag(err_par))<0.001)   %Note: for systems with small noise intensities we might get complex uncertainties in which the imaginary part is super small (this often happens in deterministic systems)
    err_par=real(err_par);
end

if isreal(err_par)
    fprintf('Observed Fisher information\n');
    disp(num2str(hess));
    fprintf('\nError of observed Fisher information\n');
    disp(num2str(err_hess))
    fprintf('\nUncertainty of the parameters (in terms of standard deviation)\n');
    disp(num2str(err_par(:).'))
else
    fprintf('\nTrying to estimate the uncertainty failed. Possible reasons are:\n');
    fprintf('\n 1) You were not careful to use exactly the same settings you used in the process of reconstruction');
    fprintf('\n 2) The parameters are badly estimated leading to the failure in estimating the uncertainties\n');
    fprintf('\nIf for any reason it is not possible to fix the problem then the following steps could be helpful:');
    fprintf('\n 1) At least you can get the uncertainties corresponding with the Euler reconstruction (this always works)');
    fprintf('\n 2) You have a range of simpler alternatives to try. You can try the following options:\n');
    fprintf('\n 1) You can consider an additive modeling strategy if your model is multiplicative');
    fprintf('\n 2) You can use a range of spline models if your current model is parametric: try splines with ''C'' flags, if this fails then try splines with ''Q'' flags, if this one also fails then try ''Approximative'' flags\n');
end   
end



%*************************************************************************************
%************************************SUBROUTINES**************************************
%*************************************************************************************
function cost=Cost_Spline1(par,dt,diff_data,N,X0,knots,M,H,ETA,SplineType)
p_end=par(end);
switch SplineType
    case {'L','Approximate'}
        pp_mu=interp1(knots,par(1:M),'linear','pp');
        pp_dmu=fnder(pp_mu);
        D0=ppval(pp_mu,X0)./p_end;D1=ppval(pp_dmu,X0);
        ETA=ETA(dt,D0,D1);
    case 'Q'
        sp_mu=spapi(3,knots,par(1:M));
        sp_dmu=fnder(sp_mu);sp_d2mu=fnder(sp_dmu);
        D0=fnval(sp_mu,X0)./p_end;D1=fnval(sp_dmu,X0);D2=fnval(sp_d2mu,X0).*p_end;
        ETA=ETA(dt,D0,D1,D2);
    case 'C'
        pp_mu=spline(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0=ppval(pp_mu,X0)./p_end;D1=ppval(pp_dmu,X0);D2=ppval(pp_d2mu,X0).*p_end;D3=ppval(pp_d3mu,X0).*p_end.^2;
        ETA=ETA(dt,D0,D1,D2,D3);
    case 'P'
        pp_mu=pchip(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0=ppval(pp_mu,X0)./p_end;D1=ppval(pp_dmu,X0);D2=ppval(pp_d2mu,X0).*p_end;D3=ppval(pp_d3mu,X0).*p_end.^2;
        ETA=ETA(dt,D0,D1,D2,D3);
    case 'SCS'
        pp_mu=csaps(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0=ppval(pp_mu,X0)./p_end;D1=ppval(pp_dmu,X0);D2=ppval(pp_d2mu,X0).*p_end;D3=ppval(pp_d3mu,X0).*p_end.^2;
        ETA=ETA(dt,D0,D1,D2,D3);
end

Z=diff_data./(p_end*sqrt(dt));
A=1./(sqrt(dt).*p_end).*normpdf(Z);
B=[ones(1,N-1);ETA];
C=[ones(1,N-1);H(Z)];
D=sum(B.*C);
L=A.*D;   % L is the vector of likelihoods
L=L(~isnan(L));   %this is to account for replicate and missing data

%The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
if ~isreal(L)
    cost=realmax;
elseif min(L)>0
    cost=-sum(log(L));   % objective function is the negative of sum of log-likelihoods
else
    cost=realmax;
end
end

function cost=Cost_Spline2(par,dt,J,N,X0,X,EZ,H,H_coeff,Coeff,knots,M,SplineType)
switch SplineType
    case {'L','Approximate'}
        pp_mu=interp1(knots,par(1:M),'linear','pp');
        pp_dmu=fnder(pp_mu);
        D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);
        EZ=EZ(dt,D0_mu,D1_mu,par(end));
    case 'Q'
        sp_mu=spapi(3,knots,par(1:M));
        sp_dmu=fnder(sp_mu);sp_d2mu=fnder(sp_dmu);
        D0_mu=fnval(sp_mu,X0);D1_mu=fnval(sp_dmu,X0);D2_mu=fnval(sp_d2mu,X0);
        EZ=EZ(dt,D0_mu,D1_mu,D2_mu,par(end));
    case 'C'
        pp_mu=spline(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);D2_mu=ppval(pp_d2mu,X0);D3_mu=ppval(pp_d3mu,X0);
        EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D3_mu,par(end));
    case 'P'
        pp_mu=pchip(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);D2_mu=ppval(pp_d2mu,X0);D3_mu=ppval(pp_d3mu,X0);
        EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D3_mu,par(end));
    case 'SCS'
        pp_mu=csaps(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);D2_mu=ppval(pp_d2mu,X0);D3_mu=ppval(pp_d3mu,X0);
        EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D3_mu,par(end)); 
end

VAR=EZ(2,:)-EZ(1,:).^2;
if min(VAR)>0
    rho=VAR.^(-1./2);
    Z=(X-X0)./(par(end).*sqrt(dt));
    Z1=rho.*(Z-EZ(1,:));
    EZh=zeros(J,N-1);EZh(2,:)=1;
    for j=3:J
        EZh(j,:)=rho.^j.*sum([EZ(j:-1:1,:).' ones(N-1,1)].*Coeff{j}(EZ(1,:).'),2).';
    end
    
    %Hermite expansion coeficients of transition density for arbitrary number of coefficients (J) and arbitrary
    %temporal tylor expansion order (K) for each coefficient
    ETA=zeros(J,N-1);
    for j=3:2:J
        ETA(j,:)=1/factorial(j).*sum(H_coeff{j}.*EZh(j:-2:1,:));
    end
    for j=4:2:J
        ETA(j,:)=1/factorial(j).*sum(H_coeff{j}.*[EZh(j:-2:1,:);ones(1,N-1)]);
    end
    ETA=ETA(3:end,:);   % We do not need the first 2 rows of ETA as they are 0
    
    A=rho./(sqrt(dt).*par(end)).*normpdf(Z1);
    B=[ones(1,N-1);ETA];
    C=[ones(1,N-1);H(Z1)];
    D=sum(B.*C);
    L=A.*D;   % L is the vector of likelihoods
    L=L(~isnan(L));   %this is to account for replicate and missing data
    
    %The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
    if ~isreal(L)
        cost=realmax;
    elseif min(L)>0
        cost=-sum(log(L));   % objective function is the negative of sum of log-likelihoods
    else
        cost=realmax;
    end
else
    cost=realmax;
end
end


function cost=Cost_Spline3(par,dt,J,N,X0,X,EZ,H,H_coeff,Coeff,knots,mesh_fine,M,SplineType)
% pp_sigma=pchip(knots,par(M+1:end));pp_sigmaInv=interp1(mesh_fine,1./ppval(pp_sigma,mesh_fine),'pchip','pp');    %we could also use this command and later use the commands pp_Z=fnint(pp_sigmaInv);Z=(ppval(pp_Z,X)-ppval(pp_Z,X0)); but this is computationally a bit more expensive
switch SplineType
    case 'LL'
        pp_mu=interp1(knots,par(1:M),'linear','pp');
        pp_dmu=fnder(pp_mu);
        D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);
        pp_sigma=interp1(knots,par(M+1:end),'linear','pp');
        sigma_meshfine=ppval(pp_sigma,mesh_fine);
        pp_dsigma=fnder(pp_sigma);
        D0_sigma=ppval(pp_sigma,X0);D1_sigma=ppval(pp_dsigma,X0);
        EZ=EZ(dt,D0_mu,D1_mu,D0_sigma,D1_sigma);
    case 'CP'
        pp_mu=spline(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);D2_mu=ppval(pp_d2mu,X0);D3_mu=ppval(pp_d3mu,X0);
        pp_sigma=pchip(knots,par(M+1:end));
        sigma_meshfine=ppval(pp_sigma,mesh_fine);
        pp_dsigma=fnder(pp_sigma);pp_d2sigma=fnder(pp_dsigma);pp_d3sigma=fnder(pp_d2sigma);
        D0_sigma=ppval(pp_sigma,X0);D1_sigma=ppval(pp_dsigma,X0);D2_sigma=ppval(pp_d2sigma,X0);D3_sigma=ppval(pp_d3sigma,X0);
        EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D3_mu,D0_sigma,D1_sigma,D2_sigma,D3_sigma);
    case 'PP'
        pp_mu=pchip(knots,par(1:M));
        pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
        D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);D2_mu=ppval(pp_d2mu,X0);D3_mu=ppval(pp_d3mu,X0);
        pp_sigma=pchip(knots,par(M+1:end));
        sigma_meshfine=ppval(pp_sigma,mesh_fine);
        pp_dsigma=fnder(pp_sigma);pp_d2sigma=fnder(pp_dsigma);pp_d3sigma=fnder(pp_d2sigma);
        D0_sigma=ppval(pp_sigma,X0);D1_sigma=ppval(pp_dsigma,X0);D2_sigma=ppval(pp_d2sigma,X0);D3_sigma=ppval(pp_d3sigma,X0);
        EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D3_mu,D0_sigma,D1_sigma,D2_sigma,D3_sigma);
end

VAR=EZ(2,:)-EZ(1,:).^2;
if min(VAR)>0
    rho=VAR.^(-1./2);
    z=cumtrapz(mesh_fine,1./sigma_meshfine);
    z=@(x)interp1(mesh_fine,z,x);
    Z=(z(X)-z(X0))./sqrt(dt);
    Z1=rho.*(Z-EZ(1,:));
    EZh=zeros(J,N-1);EZh(2,:)=1;
    for j=3:J
        EZh(j,:)=rho.^j.*sum([EZ(j:-1:1,:).' ones(N-1,1)].*Coeff{j}(EZ(1,:).'),2).';
    end
    
    %Hermite expansion coeficients of transition density for arbitrary number of coefficients (J) and arbitrary
    %temporal tylor expansion order (K) for each coefficient
    eta=zeros(J,N-1);
    for j=3:2:J
        eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*EZh(j:-2:1,:));
    end
    for j=4:2:J
        eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*[EZh(j:-2:1,:);ones(1,N-1)]);
    end
    eta=eta(3:end,:);   % We do not need the first 2 rows of ETA as they are 0
    
    A=rho./(sqrt(dt).*ppval(pp_sigma,X)).*normpdf(Z1);
    B=[ones(1,N-1);eta];
    C=[ones(1,N-1);H(Z1)];
    D=sum(B.*C);
    L=A.*D;   % L is the vector of likelihoods
    L=L(~isnan(L));   %this is to account for replicate and missing data
    
    %The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
    if ~isreal(L)
        cost=realmax;
    elseif min(L)>0
        cost=-sum(log(L));   % objective function is the negative of sum of log-likelihoods
    else
        cost=realmax;
    end
else
    cost=realmax;
end
end

function cost=Cost_Spline4(par,dt,J,N,X0,X,EZ,H,H_coeff,Coeff,knots,mesh_fine,M,SplineType,w)
% sp_sigma=spline(knots,par(M+1:end));pp_sigmaInv=interp1(mesh_fine,1./ppval(sp_sigma,mesh_fine),'pchip','pp');   %we could also use this command and later use the commands pp_Z=fnint(sp_sigmaInv);Z=(ppval(pp_Z,X)-ppval(pp_Z,X0)); but this is computationally a bit more expensive
switch SplineType
    case 'QQ'
        sp_sigma=spapi(3,knots,par(M+1:end));
    case 'CC'
        sp_sigma=spline(knots,par(M+1:end));   %this is 'pp' but to make the code coherent I chose the name 'sp'
    case 'SCS'
        sp_sigma=csaps(knots,par(M+1:end));   %this is 'pp' but to make the code coherent I chose the name 'sp'
    case 'Approximate'
        sp_sigma=csaps(knots,par(M+1:end),[],[],w);
end
if fnmin(sp_sigma)>0
    switch SplineType
        case 'QQ'
            sp_mu=spapi(3,knots,par(1:M));
            sp_dmu=fnder(sp_mu);sp_d2mu=fnder(sp_dmu);
            D0_mu=fnval(sp_mu,X0);D1_mu=fnval(sp_dmu,X0);D2_mu=fnval(sp_d2mu,X0);
            sigma_meshfine=fnval(sp_sigma,mesh_fine);
            sigma_X=fnval(sp_sigma,X);
            sp_dsigma=fnder(sp_sigma);sp_d2sigma=fnder(sp_dsigma);
            D0_sigma=fnval(sp_sigma,X0);D1_sigma=fnval(sp_dsigma,X0);D2_sigma=fnval(sp_d2sigma,X0);
            EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D0_sigma,D1_sigma,D2_sigma);
        case 'CC'
            pp_mu=spline(knots,par(1:M));
            pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
            D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);D2_mu=ppval(pp_d2mu,X0);D3_mu=ppval(pp_d3mu,X0);
            sigma_meshfine=ppval(sp_sigma,mesh_fine);
            sigma_X=ppval(sp_sigma,X);
            pp_dsigma=fnder(sp_sigma);pp_d2sigma=fnder(pp_dsigma);pp_d3sigma=fnder(pp_d2sigma);
            D0_sigma=ppval(sp_sigma,X0);D1_sigma=ppval(pp_dsigma,X0);D2_sigma=ppval(pp_d2sigma,X0);D3_sigma=ppval(pp_d3sigma,X0);
            EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D3_mu,D0_sigma,D1_sigma,D2_sigma,D3_sigma);
        case 'SCS'
            pp_mu=csaps(knots,par(1:M));
            pp_dmu=fnder(pp_mu);pp_d2mu=fnder(pp_dmu);pp_d3mu=fnder(pp_d2mu);
            D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);D2_mu=ppval(pp_d2mu,X0);D3_mu=ppval(pp_d3mu,X0);
            sigma_meshfine=ppval(sp_sigma,mesh_fine);
            sigma_X=ppval(sp_sigma,X);
            pp_dsigma=fnder(sp_sigma);pp_d2sigma=fnder(pp_dsigma);pp_d3sigma=fnder(pp_d2sigma);
            D0_sigma=ppval(sp_sigma,X0);D1_sigma=ppval(pp_dsigma,X0);D2_sigma=ppval(pp_d2sigma,X0);D3_sigma=ppval(pp_d3sigma,X0);
            EZ=EZ(dt,D0_mu,D1_mu,D2_mu,D3_mu,D0_sigma,D1_sigma,D2_sigma,D3_sigma);
        case 'Approximate'   %here, we just apply the procedure we used for the SplineType with flag 'LL' in the previous cost function
            pp_mu=interp1(knots,par(1:M),'linear','pp');
            pp_dmu=fnder(pp_mu);
            D0_mu=ppval(pp_mu,X0);D1_mu=ppval(pp_dmu,X0);
            pp_dsigma=fnder(sp_sigma);
            D0_sigma=ppval(sp_sigma,X0);D1_sigma=ppval(pp_dsigma,X0);
            sigma_meshfine=ppval(sp_sigma,mesh_fine);
            sigma_X=ppval(sp_sigma,X);
            EZ=EZ(dt,D0_mu,D1_mu,D0_sigma,D1_sigma);    
    end
    
    VAR=EZ(2,:)-EZ(1,:).^2;
    if min(VAR)>0
        rho=VAR.^(-1./2);
        z=cumtrapz(mesh_fine,1./sigma_meshfine);
        z=@(x)interp1(mesh_fine,z,x);
        Z=(z(X)-z(X0))./sqrt(dt);
        Z1=rho.*(Z-EZ(1,:));
        EZh=zeros(J,N-1);EZh(2,:)=1;
        for j=3:J
            EZh(j,:)=rho.^j.*sum([EZ(j:-1:1,:).' ones(N-1,1)].*Coeff{j}(EZ(1,:).'),2).';
        end
        
        %Hermite expansion coeficients of transition density for arbitrary number of coefficients (J) and arbitrary
        %temporal tylor expansion order (K) for each coefficient
        eta=zeros(J,N-1);
        for j=3:2:J
            eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*EZh(j:-2:1,:));
        end
        for j=4:2:J
            eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*[EZh(j:-2:1,:);ones(1,N-1)]);
        end
        eta=eta(3:end,:);   % We do not need the first 2 rows of ETA as they are 0
        
        A=rho./(sqrt(dt).*sigma_X).*normpdf(Z1);
        B=[ones(1,N-1);eta];
        C=[ones(1,N-1);H(Z1)];
        D=sum(B.*C);
        L=A.*D;   % L is the vector of likelihoods
        L=L(~isnan(L));   %this is to account for replicate and missing data
        
        %The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
        if ~isreal(L)
            cost=realmax;
        elseif min(L)>0
            cost=-sum(log(L));   % objective function is the negative of sum of log-likelihoods
        else
            cost=realmax;
        end
    else
        cost=realmax;
    end
else
    cost=realmax;
end
end

function cost=Cost_Spline_Euler(par,dt,X0,X,knots,M,ModelType,SplineType)
switch ModelType
    case 'Additive noise'
        STD=par(end);
        switch SplineType
            case 'L'
                mu=@(x)interp1(knots,par(1:M),x);
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);
            case 'Q'
                sp_mu=spapi(3,knots,par(1:M));
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-fnval(sp_mu,X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);
            case 'C'
                mu=@(x)interp1(knots,par(1:M),x,'cubic');
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);
            case 'P'
                mu=@(x)interp1(knots,par(1:M),x,'pchip');
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);
            case 'SCS'
                pp_mu=csaps(knots,par(1:M));
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);       
        end
    case 'Multiplicative noise'
        switch SplineType
            case 'LL'
                mu=@(x)interp1(knots,par(1:M),x);sigma=@(x)interp1(knots,par(M+1:end),x);
                STD=sigma(X0);
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);
            case 'QQ'
                sp_sigma=spapi(3,knots,par(M+1:end));
                if fnmin(sp_sigma,[knots(1) knots(end)])>0
                    sp_mu=spapi(3,knots,par(1:M));
                    STD=fnval(sp_sigma,X0);
                    LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-fnval(sp_mu,X0).*dt)./(STD.*sqrt(dt))).^2);
                    LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                    cost=-sum(LL);
                else
                    cost=realmax;
                end
            case 'CC'
                pp_sigma=spline(knots,par(M+1:end));
                if fnmin(pp_sigma,[knots(1) knots(end)])>0
                    pp_mu=spline(knots,par(1:M));
                    STD=ppval(pp_sigma,X0);
                    LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-ppval(pp_mu,X0).*dt)./STD.*sqrt(dt)).^2);
                    LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                    cost=-sum(LL);
                else
                    cost=realmax;
                end
            case 'SCS'
                pp_sigma=csaps(knots,par(M+1:end));
                if fnmin(pp_sigma,[knots(1) knots(end)])>0
                    pp_mu=csaps(knots,par(1:M));
                    STD=ppval(pp_sigma,X0);
                    LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-ppval(pp_mu,X0).*dt)./STD.*sqrt(dt)).^2);
                    LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                    cost=-sum(LL);
                else
                    cost=realmax;
                end
            case 'Approximate'
                pp_sigma1=csaps(knots,par(M+1:end));
                if fnmin(pp_sigma1,[knots(1) knots(end)])>0
                     mu=@(x)interp1(knots,par(1:M),x);sigma=@(x)interp1(knots,par(M+1:end),x);                   
                     STD=sigma(X0);
                     LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                     LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                     cost=-sum(LL);
                else
                    cost=realmax;
                end
            case 'CP'
                pp_mu=spline(knots,par(1:M));pp_sigma=pchip(knots,par(M+1:end));
                STD=ppval(pp_sigma,X0);
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);
            case 'PP'
                pp_mu=pchip(knots,par(1:M));pp_sigma=pchip(knots,par(M+1:end));
                STD=ppval(pp_sigma,X0);
                LL=-1./2.*(log(2.*pi.*dt.*STD.^2)+((X-X0-mu(X0).*dt)./(STD.*sqrt(dt))).^2);
                LL=LL(~isnan(LL));   %this is to account for replicate and missing data
                cost=-sum(LL);
        end          
end
end

function cost=Cost_Simple(par,data,X,N,dt,sigma,d_muY,Gam,H,eta)
p=num2cell(par);
Y=Gam(data,p{:});
Z=diff(Y)./sqrt(dt);

D_muY=d_muY(Y(1:end-1),p{:});
G1=D_muY.';
G2=cellfun(@transpose,num2cell(G1,1),'uniform',0);

A=1./(sqrt(dt).*sigma(X,p{:})).*normpdf(Z);
B=[ones(1,N-1);eta(G2{:})];
C=[ones(1,N-1);H(Z)];
D=sum(B.*C); 
L=A.*D;   % L is the vector of likelihoods
L=L(~isnan(L));   %this is to account for replicate and missing data

%The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
if ~isreal(L) 
    cost=realmax;   % objective function is the negative of sum of log-likelihoods
elseif min(L)>0
    cost=-sum(log(L));
else
    cost=realmax;
end
end

function cost=Cost_Additive(par,J,Coeff,X0,X,N,dt,d_mu,EZ,H,H_coeff)

p=num2cell(par);
Dmu=d_mu(X0,p{:});
G1=[Dmu.' p{end}.*ones(N-1,1)];
G2=cellfun(@transpose, num2cell(G1,1),'uniform',0);
EZ=EZ(G2{:});

VAR=EZ(2,:)-EZ(1,:).^2;
if min(VAR)>=0
    rho=VAR.^(-1./2);
    Z=(X-X0)./(p{end}.*sqrt(dt));
    Z1=rho.*(Z-EZ(1,:));
    
    EZh=zeros(J,N-1);EZh(2,:)=1;
    for j=3:J
        EZh(j,:)=rho.^j.*sum([EZ(j:-1:1,:).' ones(N-1,1)].*Coeff{j}(EZ(1,:).'),2).';
    end

    %Hermite expansion coeficients of transition density for arbitrary number of coefficients (J) and arbitrary
    %temporal tylor expansion order (K) for each coefficient
    eta=zeros(J,N-1);
    for j=3:2:J
        eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*EZh(j:-2:1,:));
    end
    for j=4:2:J
        eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*[EZh(j:-2:1,:);ones(1,N-1)]);
    end
    eta=eta(3:end,:);   % We do not need the first 2 rows of eta as they are 0

    A=rho./(sqrt(dt).*p{end}).*normpdf(Z1);
    B=[ones(1,N-1);eta];
    C=[ones(1,N-1);H(Z1)];
    D=sum(B.*C);
    L=A.*D;   % L is the vector of likelihoods
    L=L(~isnan(L));   %this is to account for replicate and missing data

    % The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening due to the rather 'far' 
    % distance between the current parameter values from the optimal parameters 
    if ~isreal(L) 
        cost=realmax;   
    elseif min(L)>0
        cost=-sum(log(L));   % objective function is the negative of sum of log-likelihoods
    else
        cost=realmax;
    end
else
    cost=realmax;
end
end

function cost=Cost_Complex(par,J,Coeff,X0,X,N,dt,sigma,d_mu,d_sigma,EZ,H,H_coeff,Z)
p=num2cell(par);
Dmu=d_mu(X0,p{:});
Dsigma=d_sigma(X0,p{:});
G1=[Dmu.' Dsigma.'];
G2=cellfun(@transpose, num2cell(G1,1),'uniform',0);

EZ=EZ(G2{:});
Z1=Z(X0,X,p{:});

VAR=EZ(2,:)-EZ(1,:).^2;
if min(VAR)>=0
    rho=VAR.^(-1./2);
    Z2=rho.*(Z1- EZ(1,:));
    EZh=zeros(J,N-1);EZh(2,:)=1;
    for j=3:J
        EZh(j,:)=rho.^j.*sum([EZ(j:-1:1,:).' ones(N-1,1)].*Coeff{j}(EZ(1,:).'),2).';
    end

    %Hermite expansion coeficients of transition density for arbitrary number of coefficients (J) and arbitrary
    %temporal tylor expansion order (K) for each coefficient
    eta=zeros(J,N-1);
    for j=3:2:J
        eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*EZh(j:-2:1,:));
    end
    for j=4:2:J
        eta(j,:)=1/factorial(j).*sum(H_coeff{j}.*[EZh(j:-2:1,:);ones(1,N-1)]);
    end
    eta=eta(3:end,:);   % We do not need the first 2 rows of eta as they are 0

    A=rho./(sqrt(dt).*sigma(X,p{:})).*normpdf(Z2);
    B=[ones(1,N-1);eta];
    C=[ones(1,N-1);H(Z2)];
    D=sum(B.*C);
    L=A.*D;   % L is the vector of likelihoods
    L=L(~isnan(L));   %this is to account for replicate and missing data

    %The following is a 'death penalty' implementation of the constraint of not geting illegitimate likelihoods happening when parameters are rather 'far' from the optimal parameters
    if ~isreal(L) 
        cost=realmax;   
    elseif min(L)>0
        cost=-sum(log(L));   % objective function is the negative of sum of log-likelihoods
    else
        cost=realmax;
    end
else
    cost=realmax;
end
end

function cost=Cost_Parametric_Euler(par,X,X0,LL)
p=num2cell(par);
LL=LL(X,X0,p{:});
cost=-sum(LL(~isnan(LL)));   % 'isnan' is used to account for replicate and missing data
end

function cost=Cost_Parametric_Euler_MultiVariate(par,XX,XX0,LL)
LL=LL(XX,XX0,par);
cost=-sum(LL(~isnan(LL)));   % 'isnan' is used to account for replicate and missing data
end