function [ar,ASAsellog,ASAcontrol] = sig2ar_s(sig,cand_order,last,estKeff_model,estKeff_mean)
%SIG2AR AR model identification from multiple segments
%   [AR,SELLOG] = SIG2AR(SIG) estimates autoregressive models from the 
%   data matrix/vector SIG and selects a model with optimal predictive 
%   qualities. The selected model is returned in the parameter vector AR. 
%   The structure SELLOG provides additional information on the selection 
%   process.
%
%   If SIG is a matrix, its columns are regarded as different segments of
%   the same process. These segments do not have to be independent, and
%   there is an effective number of segments both for estimation of
%   the model parameters and for estimation of the overall mean of the
%   data (These effective numbers of segments need not be the same). With
%   estKeff_model/estKeff_mean we can indicate whether we want to estimate
%   these effective numbers of segments. If so, Keff(model) is estimated
%   from the variance of the insignificant reflection coefficients in high-
%   order AR models for the individual segments, and obtained from all
%   segments at once. Keff(model) is used to adapt the order selection
%   criteria.
%   Keff(mean) is estimated from the average correlation function of the
%   segments and influences the estimation of the prediction error.
%   With estKeff_model we can indicate whether we want to estimate
%   Keff(model) here, or use a value chosen by the user (provided by
%   armasel.m, for example).
%   If estKeff_model<=0, we estimate Keff(model) from the data.
%   If estKeff_model>0, we use its value as Keff(model) (but force it to
%   be >= 1 and smaller than the number of segments).
%   If ASAglob_subtr_mean=0, Keff(mean) is only estimated when
%   ASAglob_mean_adj=1. If the latter is also 0, Keff(mean) is taken as the
%   number of segments. If ASAglob_subtr_mean=1 (the default), then the
%   first element of estKeff_mean is used as before:
%   If estKeff_mean(1)<=0, we estimate Keff(mean) from the data.
%   If estKeff_mean(1)>0, we use its value as Keff(mean) (but force it to
%   be >= 1).
%   The second element, estKeff_mean(2), indicate the subsampling of the
%   correlation function we are going to use (this stems from our
%   application where we look at images with circularly symmetric
%   correlation functions. Each segment is then a column or row of the
%   image, but these do not have to be adjacent. We could be scanning every
%   second line for example.
%   Note that if we have more than 1 segment, the overall mean is
%   subtracted from each segment (that is, the mean of the segment mean
%   values).
%   
%   SIG2AR(SIG,CAND_ORDER) selects only from candidate models whose 
%   orders are entered in CAND_ORDER. CAND_ORDER must either be a row of 
%   ascending orders, or a single order (in which case no true order 
%   selection is performed).
%   
%   Without user intervention, the mean of SIG is subtracted from the 
%   data. To control the subtraction of the mean, see the help topics on 
%   ASAGLOB_SUBTR_MEAN and ASAGLOB_MEAN_ADJ.
%     
%   SIG2AR is an ARMASA main function.
%   
%   See also: SIG2MA, SIG2ARMA, ARMASEL.

%   References: P. M. T. Broersen, Facts and Fiction in Spectral
%               Analysis, IEEE Transactions on Instrumentation and
%               Measurement, Vol. 49, No. 4, August 2000, pp. 766-772.

%Header
%=========================================================================

%Declaration of variables
%------------------------

%Declare and assign values to local variables
%according to the input argument pattern
switch nargin
case 1 
   if isa(sig,'struct'), ASAcontrol=sig; sig=[];
   else, ASAcontrol=[];
   end
   cand_order=[];
case 2 
   if isa(cand_order,'struct'), ASAcontrol=cand_order; cand_order=[]; 
   else, ASAcontrol=[]; 
   end
case 3 
   if isa(last,'struct'), ASAcontrol=last;
   else, error(ASAerr(39))
   end
case 4
   %pass arguments without modification   
   %ASAcontrol=[];
case 5
   %pass arguments without modification 
   %ASAcontrol=[];
otherwise
   error(ASAerr(1,mfilename))
end

if isequal(nargin,1) & ~isempty(ASAcontrol)
      %ASAcontrol is the only input argument
   ASAcontrol.error_chk = 0;
   ASAcontrol.run = 0;
end

%Declare ASAglob variables 
ASAglob = {'ASAglob_subtr_mean';'ASAglob_mean_adj'; ...
      'ASAglob_rc';'ASAglob_ar'};

%Assign values to ASAglob variables by screening the
%caller workspace
for ASAcounter = 1:length(ASAglob)
   ASAvar = ASAglob{ASAcounter};
   eval(['global ' ASAvar]);
   if evalin('caller',['exist(''' ASAvar ''',''var'')'])
      eval([ASAvar '=evalin(''caller'',ASAvar);']);
   else
      eval([ASAvar '=[];']);
   end
end

%ARMASA-function version information
%-----------------------------------

%This ARMASA-function is characterized by
%its current version,
ASAcontrol.is_version = [2000 12 30 20 0 0];
%and its compatability with versions down to,
ASAcontrol.comp_version = [2000 12 30 20 0 0];

%This function calls other functions of the ARMASA
%toolbox. The versions of these other functions must
%be greater than or equal to:
ASAcontrol.req_version.burg = [2000 12 30 20 0 0];
ASAcontrol.req_version.cic = [2000 12 30 20 0 0];
ASAcontrol.req_version.rc2arset = [2000 12 30 20 0 0];

%Checks
%------

if ~any(strcmp(fieldnames(ASAcontrol),'error_chk')) | ASAcontrol.error_chk
      %Perform standard error checks
   %Input argument format checks
   ASAcontrol.error_chk = 1;
   if ~isnum(sig)
      error(ASAerr(11,'sig'))
%   elseif ~isavector(sig)
%      error([ASAerr(14) ASAerr(15,'sig')])
%   elseif size(sig,2)>1
%      sig = sig(:);
%      warning(ASAwarn(25,{'row';'sig';'column'},ASAcontrol))
   end
   if ~isempty(cand_order)
      if ~isnum(cand_order) | ~isintvector(cand_order) |...
            cand_order(1)<0 | ~isascending(cand_order)
         error(ASAerr(12,{'candidate';'cand_order'}))
      elseif size(cand_order,1)>1
         cand_order = cand_order';
         warning(ASAwarn(25,{'column';'cand_order';'row'},ASAcontrol))
      end
   end
   
   %Input argument value checks
   if ~isreal(sig)
      error(ASAerr(13))
   end
   if max(cand_order) > length(sig)-1
      error(ASAerr(21))
   end
end

if ~any(strcmp(fieldnames(ASAcontrol),'version_chk')) | ASAcontrol.version_chk
      %Perform version check
   ASAcontrol.version_chk = 1;
      
   %Make sure the requested version of this function
   %complies with its actual version
   ASAversionchk(ASAcontrol);
   
   %Make sure the requested versions of the called
   %functions comply with their actual versions
   burg(ASAcontrol);
   cic(ASAcontrol);
   rc2arset(ASAcontrol);
end

if ~any(strcmp(fieldnames(ASAcontrol),'run')) | ASAcontrol.run
      %Run the computational kernel
   ASAcontrol.run = 1;
   ASAcontrol.version_chk = 0;
   ASAcontrol.error_chk = 0;
   ASAtime = clock;
   ASAdate = now;

%Main   
%=====================================================
  
%Initialization of variables
%---------------------------

n_obs = size(sig,1);n_segs=size(sig,2);
% In case the signal is 1 row vector, turn it into 1 column vector
if (n_obs==1 || n_segs==1) && n_segs>n_obs
    sig=sig(:);
    n_obs=n_segs;n_segs=1;
end
if isempty(ASAglob_subtr_mean) | ASAglob_subtr_mean
    ASAglob_subtr_mean=1;
   %Subtract mean per segment
   %sig  = sig - ones(n_obs,1)*mean(sig);
   %Subtract global mean
   sig =sig - mean(mean(sig));
   if isempty(ASAglob_mean_adj)
      ASAglob_mean_adj = 1;
   end
elseif isempty(ASAglob_mean_adj)
   ASAglob_mean_adj = 0;
end

%Estimation of the effective number of segments Keff for model estimation.
%-------------------------------------------------------------------------

%Minimum number of reflection coefficients required to estimate Keff
numKmin=1;
%Maximum AR order we consider
Pmax = min([n_obs-5 fix(200*log10(n_obs)) 1000]);
% finite-sample variance coefficients for Burg's method
t=1:n_obs;
vi=1./(n_obs+1-t);clear t
if estKeff_model(1)>0
    Keff_model=max(estKeff_model,1);
elseif n_segs>1
   %find maximum+1 of selected orders of all segments
   arorders=zeros(1,n_segs);
   RCs=zeros(n_segs,Pmax);
   %To avoid having to run sig2ar twice for every segment, we will
   %compute the average correlation function already here (only needed
   %for Keff_mean) when necessary
%    if estKeff_mean(1)<=0 & (ASAglob_subtr_mean==1 || ASAglob_mean_adj==1)
%        DoCORest=1;
%        COR=zeros(n_segs,n_segs);
%        %sub sample factor
%        S=estKeff_mean(2);
%    else
       DoCORest=0;
%    end
   for k=1:n_segs
       [par,sellog]=sig2ar(sig(:,k),0:Pmax);
       arorders(k)=length(par);
       RCs(k,:)=sellog.rcarlong(2:Pmax+1);
       if DoCORest
          cor=arma2cor(par,1,(n_segs-1)*S);
          % Average with gain
          COR(k,:)=std(sig(:,k))^2*cor(1:S:end);
       end
   end
   ASAglob_rc=[];
   Pmin=max(max(arorders),fix(n_obs/2)+1);
   if Pmax-Pmin<numKmin-1
       Pmin=max(arorders);
   end
   disp([Pmin Pmax])
   if Pmax-Pmin<numKmin-1
       Keff_model=1
   else
       %Estimate AR model of order Pmax from all segments at once
       rc=burg_s_modified(sig,Pmax);%NOTE: rc(1)=1
       rc=rc(2:end);
       %Keff_model= mean(mean(RCs(:,Pmin:Pmax).^2))/mean(rc(Pmin:Pmax).^2);
       Keff_model=1/mean(rc(Pmin:Pmax).^2./vi(Pmin:Pmax));
       Keff_model=min(max(Keff_model,1),n_segs);
   end
else
    Keff_model=n_segs;
end

%Estimation of the effective number of segments Keff for mean estimation.
%-------------------------------------------------------------------------
% Since it is not clear how exactly we should correct for the subtraction
% of the mean (does it depend on the correlation in the data?), we use
% Keff_mean=Keff_model.

DoCORest=0;
if estKeff_mean(1)>0
    Keff_mean=max(estKeff_mean(1),1);
elseif n_segs>1 && DoCORest
    COR=mean(COR);COR=COR/COR(1);
    t=2:n_segs;
    Keff_mean=n_segs/(1+2*sum((1-(t-1)/n_segs).*COR(t)))
else
    Keff_mean=n_segs;
end
%Keff_mean=Keff_model;
%disp('Keff_model Keff_mean')
%disp([Keff_model Keff_mean])

%Determination of the maximum candidate AR order
%-----------------------------------------------

if ~isempty(cand_order)
   max_order = cand_order(end);
else
   max_order = ...
      min(fix(n_obs/2),fix(200*log10(n_obs)));
   if max_order > 1000; 
      max_order = 1000;
   end
end


%Estimation procedure
%--------------------

[rc,varsig] = burg_s_modified(sig,max_order,ASAcontrol);

%AR model order selection
%------------------------

rc(1) = 0;
res = varsig*cumprod(1-rc.^2);
rc(1) = 1;
%[cicar,pe_est] = cic(res,n_obs,cand_order,ASAcontrol);
[cicar,pe_est] = cic_s(res,n_obs,Keff_model,Keff_mean,cand_order,ASAcontrol);
[min_value,sel_location] = min(cicar);
if isempty(cand_order)
   cand_order = 0:max_order;
end
sel_order = cand_order(sel_location);

%Arranging output arguments
%--------------------------

ar = rc2arset(rc(1:sel_order+1),ASAcontrol);

ASAglob_rc = rc;
ASAglob_ar = ar;

if nargout>1
   ASAsellog.funct_name = mfilename;
   ASAsellog.funct_version = ASAcontrol.is_version;
   ASAsellog.date_time = ...
      [datestr(ASAdate,8) 32 datestr(ASAdate,0)];
   ASAsellog.comp_time = etime(clock,ASAtime);
   ASAsellog.ar = ar;
   ASAsellog.rcarlong = rc;
   ASAsellog.mean_adj = ASAglob_mean_adj;
   ASAsellog.cand_order = cand_order;
   ASAsellog.cic = cicar;
   ASAsellog.pe_est = pe_est;
end

%Footer
%=====================================================

else %Skip the computational kernel
   %Return ASAcontrol as the first output argument
   if nargout>1
      warning(ASAwarn(9,mfilename,ASAcontrol))
   end
   ar = ASAcontrol;
   ASAcontrol = [];
end

function Neff=NewtonOptNeff(n_obs,n_segs,Pmin,Pmax,Sk,stepsize,numiter)
% N: initial condition
Neff=n_segs;
Pmin
Pmax
disp((Pmax-Pmin+1)/sum(Sk)/n_obs)
p=Pmin:Pmax;
% To detect sign change in dfdN
s=1;
for iter=1:numiter
    %vi=1./(n_obs+1-p)/Neff;
    % vi's voor 1 segment, Burg
    vi=1./(n_obs+1-p);
    % criterion sum{sum(vi)-sum(rc_i ^2)}^2
    dfdN=2/Neff^2*(sum(Sk)-sum(vi)/Neff)*sum(vi);
    d2fdN2=6/Neff^4*sum(vi)^2-4/Neff^3*sum(Sk)*sum(vi);
    % criterion sum{sqrt(vi)-abs(rc_i)}^2
    %dfdN=1/Neff^1.5*sum(sqrt(vi).*sqrt(Sk))-1/Neff^2*sum(vi);
    %d2fdN2=2/Neff^3*sum(vi)-1.5/Neff^2.5*sum(sqrt(vi).*sqrt(Sk));
    if sign(dfdN)~=s,stepsize=stepsize/2;s=sign(dfdN);
    else
        %Neff=Neff-stepsize*dfdN/d2fdN2;
        %steepest descend, dfdN is order 1/N^3
        if dfdN<0, Neff=min(Neff*(1+stepsize),5*n_segs);else Neff=max(Neff/(1+stepsize),1);end
    end
    %Neff=Neff-stepsize*dfdN*min(N,Neff)^3/length(p);
    if ~rem(iter,pi)
        disp(' '),disp(Neff),disp([dfdN d2fdN2]),disp(s),disp([sum(Sk) sum(vi)])
        figure(1),t=1:length(vi);plot(t,sqrt(Sk),t,sqrt(vi/Neff))
        pause
    end
end


%Program history
%======================================================================
%
% Version                Programmer(s)          E-mail address
% -------                -------------          --------------
% former versions        P.M.T. Broersen        p.m.t.broersen@tudelft.nl
% [2000 12 30 20 0 0]    W. Wunderink           wwunderink01@freeler.nl

% Try to use 1-p/N for determination of Neff:
%         %If P<0, P=0.1*N is used, if P==0, Neff=N is assumed.
%         P=acfmode(2);
%         if P<0, P=floor(0.1*n_obs);end
%         if P
%             % First estimate residual variance from product (1-rc(i)^2) for
%             % each segment separately and average these.
%             resvar=0;
%             for k=1:n_segs
%                 [rc,varsig]=burg_modified(sig(:,k),P);
%                 %resvar=resvar+varsig*cumprod(1-rc(2:end).^2)/n_segs;
%                 par=rc2poly(rc(2:end));
%                 res=filter(par,1,sig(:,k));
%                 resvar=resvar+mean(res(P+1:end).^2)/n_segs;
%             end
%             % Now do the same for all segments at once
%             [rc,varsig]=burg_s_modified(sig,P);
%             %resvars=varsig*cumprod(1-rc(2:end).^2);
%             par=rc2poly(rc(2:end));
%             res=filter(par,1,sig);
%             resvars=mean(mean(res(P+1:end,:).^2));
%             %figure(1),clf,plot(resvar),hold on,plot(resvars,'g--'),keyboard
%             resvar=resvar(end);resvars=resvars(end);
%             disp([resvar resvars])
%             % Estimate Neff
%             Neff=P/(1-(resvars/resvar)*(1-P/n_obs))
%             %Neff=0.5*(P+1)/(1-(resvars/resvar)*(1-(P+1)/2/n_obs));
%          else%P==0
%             Neff=n_obs
%          end
