function [ar,ma,ASAsellog,ASAcontrol]=sig2arma_s(sig,cand_ar_order,arma_order_diff,last,estKeff_model,estKeff_mean)
%SIG2ARMA ARMA model identification from multiple segments
%   [AR,MA,SELLOG] = SIG2ARMA(SIG) estimates autoregressive moving 
%   average models from the data vector SIG and selects a model with 
%   optimal predictive qualities. Only ARMA(P,P-1) models are considered 
%   with AR order P being greater than the MA order by one. The selected 
%   model is returned in the parameter vectors AR and MA. The structure 
%   SELLOG provides additional information on the selection process.
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
%   SIG2ARMA(SIG,CAND_AR_ORDER,ARMA_ORDER_DIFF) selects only from 
%   candidate models ARMA(CAND_AR_ORDER,CAND_AR_ORDER - ARMA_ORDER_DIFF). 
%   CAND_AR_ORDER must either be a row of ascending orders, or a single 
%   order (in which case no true order selection is performed). 
%   ARMA_ORDER_DIFF, a scalar, is the difference between AR and MA orders 
%   being constant during selection. Only a positive difference greater 
%   than 0 is allowed, that forms valid pairs of orders in combination 
%   with CAND_AR_ORDER. As an exception to the rule, the first element of 
%   CAND_AR_ORDER may always be chosen 0, to include analysis of the 
%   zero-order ARMA(0,0) model.
%   
%   CAND_AR_ORDER or ARMA_ORDER_DIFF may be passed as empty arguments. In 
%   case of empty ARMA_ORDER_DIFF, the difference between orders will be 
%   chosen 1. In case of empty CAND_AR_ORDER, an appropriate set of 
%   candidate orders will be chosen depending on the value of 
%   ARMA_ORDER_DIFF and the number of observations.
%   
%   Without user intervention, the mean of SIG is subtracted from the 
%   data. To control the subtraction of the mean, see the help topics on 
%   ASAGLOB_SUBTR_MEAN and ASAGLOB_MEAN_ADJ.
%     
%   SIG2ARMA is an ARMASA main function.
%   
%   See also: SIG2AR, SIG2MA, ARMASEL.

%   References: P. M. T. Broersen, Autoregressive Model Orders for
%               Durbin's MA and ARMA estimators, IEEE Transactions on
%               Signal Processing, Vol. 48, No. 8, August 2000,
%               pp. 2454-2457.

%Header
%===================================================================================

%Declaration of variables
%------------------------

%Declare and assign values to local variables
%according to the input argument pattern
switch nargin
case 1 
   if isa(sig,'struct'), ASAcontrol=sig; sig=[];
   else, ASAcontrol=[];
   end
   cand_ar_order=[]; arma_order_diff=[];
case 2 
   if isa(cand_ar_order,'struct'), ASAcontrol=cand_ar_order; cand_ar_order=[]; 
   else, ASAcontrol=[]; 
   end
   arma_order_diff=[];
case 3
   if isa(arma_order_diff,'struct'), ASAcontrol=arma_order_diff; arma_order_diff=[]; 
   else, ASAcontrol=[]; 
   end
case 4 
   if isa(last,'struct'), ASAcontrol=last;
   else, error(ASAerr(39))
   end
case 5
   %pass arguments without modification
   %ASAcontrol=[];
case 6
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
      'ASAglob_rc';'ASAglob_ar';'ASAglob_final_f'; ...
      'ASAglob_final_b';'ASAglob_ar_cond'};

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
ASAcontrol.is_version = [2003 4 2 15 0 0];
%and its compatability with versions down to,
ASAcontrol.comp_version = [2000 12 30 20 0 0];

%This function calls other functions of the ARMASA
%toolbox. The versions of these other functions must
%be greater than or equal to:
ASAcontrol.req_version.burg = [2000 12 30 20 0 0];
ASAcontrol.req_version.cic = [2000 12 30 20 0 0];
ASAcontrol.req_version.rc2arset = [2000 12 30 20 0 0];
ASAcontrol.req_version.ar2arset = [2000 12 30 20 0 0];
ASAcontrol.req_version.cov2arset = [2000 12 30 20 0 0];
ASAcontrol.req_version.armafilter = [2000 12 12 14 0 0];
ASAcontrol.req_version.convol = [2000 12 6 12 17 20];
ASAcontrol.req_version.convolrev = [2000 12 6 12 17 20];
ASAcontrol.req_version.deconvol = [2000 12 12 12 0 0];

%Checks
%------

if ~any(strcmp(fieldnames(ASAcontrol),'error_chk')) | ASAcontrol.error_chk
      %Perform standard error checks
   %Input argument format checks
   ASAcontrol.error_chk = 1;
   if ~isnum(sig)
      error(ASAerr(11,'sig'))
   end
%    if ~isavector(sig)
%       error([ASAerr(14) ASAerr(15,'sig')])
%    elseif size(sig,2)>1
%       sig = sig(:);
%       warning(ASAwarn(25,{'row';'sig';'column'},ASAcontrol))
%    end
   if ~isempty(cand_ar_order)
      if ~isnum(cand_ar_order) | ~isintvector(cand_ar_order) |...
            cand_ar_order(1)<0 | ~isascending(cand_ar_order)
         error(ASAerr(12,{'candidate';'cand_ar_order'}))
      elseif size(cand_ar_order,1)>1
         cand_ar_order = cand_ar_order';
         warning(ASAwarn(25,{'column';'cand_ar_order';'row'},ASAcontrol))
      end
   end
   if ~isempty(arma_order_diff) & ...
         (~isnum(arma_order_diff) | ...
         ~isintscalar(arma_order_diff) |...
         arma_order_diff<0)
      error(ASAerr(17,'arma_order_diff'))
   end
   
   %Input argument value checks
   if ~isreal(sig)
      error(ASAerr(13))
   end
   if ~isempty(cand_ar_order) & ...
         ~isempty(arma_order_diff)
      if cand_ar_order(1)~=0 & ...
            (arma_order_diff < 1 | ...
            arma_order_diff > cand_ar_order(1))
         error(ASAerr(18,{'arma_order_diff';'1';...
               num2str(cand_ar_order(1))}))
      elseif length(cand_ar_order)>1 & ...
            (arma_order_diff < 1 | ...
            arma_order_diff > cand_ar_order(2))
         error(ASAerr(18,{'arma_order_diff';'1';...
               num2str(cand_ar_order(2))}))
      end
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
   ar2arset(ASAcontrol);
   cov2arset(ASAcontrol);
   armafilter(ASAcontrol);
   convol(ASAcontrol);
   convolrev(ASAcontrol);
   deconvol(ASAcontrol);
end

if ~any(strcmp(fieldnames(ASAcontrol),'run')) | ASAcontrol.run
      %Run the computational kernel
   ASAcontrol.run = 1;
   ASAcontrol.version_chk = 0;
   ASAcontrol.error_chk = 0;
   ASAtime = clock;
   ASAdate = now;
   
%Main   
%==================================================================================================
  
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
%n_obs = size(sig,1);
ar_orig = cell(5,1);
ar_entry = ones(1,5);
rc = [];
ar_sel = 1;
ma_sel = 1;

warn_state = warning;
if ischar(warn_state)
    warn_recov_comm = ['warning ' warn_state];
else
    warn_recov_comm = ['warning(warn_state);'];
end

%Combined determination of maximum candidate AR
%and MA orders
%----------------------------------------------

def_max_ar_order = min(fix(n_obs/10),fix(40*log10(n_obs)));
if def_max_ar_order > 200; 
   def_max_ar_order = 200;
end

if isempty(cand_ar_order)
   if isempty(arma_order_diff)
      cand_ar_order = 0:def_max_ar_order;   
   elseif arma_order_diff < 1 | arma_order_diff > def_max_ar_order
      error(ASAerr(18,{'arma_order_diff';'1';...
         [num2str(def_max_ar_order) ...
         ' (== max. candidate AR order, selected by default)']}))
   else
      cand_ar_order = [0 arma_order_diff:def_max_ar_order];
   end
end
if isempty(arma_order_diff)
   arma_order_diff = 1;
end
max_ar_order = cand_ar_order(end);
max_ma_order = max_ar_order-arma_order_diff;

if max_ar_order <= def_max_ar_order
   max_slid_ar_order = fix(5*def_max_ar_order);
else
   max_slid_ar_order = fix(5*max_ar_order);
   if max_slid_ar_order > n_obs-1
      max_slid_ar_order = n_obs-1;
      if max_slid_ar_order-max_ar_order < max_ma_order
         max_ar_order = fix((arma_order_diff+max_slid_ar_order)/2);
         max_ma_order = max_ar_order - arma_order_diff;
         warning(ASAwarn(16,{num2str(max_ar_order);num2str(max_ma_order)},ASAcontrol));
      end
      if ~any(max_ar_order==cand_ar_order)
         error(ASAerr(38))
      end
   end
end

min_ar_order = arma_order_diff;
min_ma_order = 0;

%Preparations for the estimation procedure
%-----------------------------------------
penalty_factor=3;
l_cand_ar_order = length(cand_ar_order);
ma = zeros(1,max_ma_order);
ar = zeros(1,max_ar_order);
gic3 = zeros(1,l_cand_ar_order);
pe_est = zeros(1,l_cand_ar_order);
test = cand_ar_order(1:min(end,2));
zero_incl = sum(0==test);
min_incl = sum(min_ar_order==test);
%var = sig'*sig/n_obs;
var=mean(mean(sig.^2));

req_counter = 0;
if zero_incl
   req_counter = 1;
   %gic3(1) = log(var)+3/n_obs/Keff_model;
   %gic3(1) = log(var)+3/Neff;
   if ASAglob_mean_adj
      pe_est(1) = var*(n_obs*Keff_mean+1)/(n_obs*Keff_mean-1);
      gic3(1) = log(var)+penalty_factor/n_obs/Keff_mean;
   else
      pe_est(1) = var;
      gic3(1) = log(var);
   end
end

if max_ar_order > 0
   %Conditioning AR orders to the previously selected AR model
   if isequal(ASAglob_ar_cond,1) & ~isempty(ASAglob_ar)
      ar = ASAglob_ar;
      sel_ar_order = length(ar)-1;
      if  3*sel_ar_order+max_ar_order+max_ma_order < max_slid_ar_order
         max_slid_ar_order = 3*sel_ar_order+max_ar_order+max_ma_order;
      end
   end 

   %AR model estimation
   l_rc = length(ASAglob_rc);
   if l_rc>1
      rc = ASAglob_rc;
      if (l_rc < max_slid_ar_order+1)
         if isempty(ASAglob_final_f)
            ar_det = rc2arset(ASAglob_rc(1:end-1),ASAcontrol);
            ASAglob_final_f = convol(sig,ar_det,l_rc-1,n_obs,ASAcontrol);
            ASAglob_final_b = convolrev(ar_det,sig,l_rc-1,n_obs,ASAcontrol);
         end
         rc = [ASAglob_rc burg(ASAglob_final_f,...
               ASAglob_final_b,max_slid_ar_order+1-l_rc,ASAcontrol)];
      end
   else
      rc = burg_s_modified(sig,max_slid_ar_order,ASAcontrol);
      %par=arburg(sig,max_slid_ar_order);
      %rc=[1 poly2rc(par)'];clear par
   end
   
   %AR model order selection
   if ~isequal(ASAglob_ar_cond,1) | isempty(ASAglob_ar)
      rc(1) = 0;
      res = var*cumprod(1-rc(1:max_slid_ar_order+1).^2);
      rc(1) = 1;
      %[min_value,sel_location] = min(cic(res,n_obs,ASAcontrol));
      [min_value,sel_location] = min(cic_s(res,n_obs,Keff_model,Keff_mean,ASAcontrol));
      sel_ar_order = sel_location-1;
   end
   
   pred_ar_order = min(3*sel_ar_order+min(9,1+fix(n_obs/10)),max_slid_ar_order);
   
   if l_cand_ar_order>zero_incl+min_incl
      try_ar_order = cand_ar_order(zero_incl+min_incl+1);
      try_slid_ar_order = min(3*sel_ar_order+2*try_ar_order-arma_order_diff,max_slid_ar_order);
      min_slid_ar_order = min(3*sel_ar_order+2*(min_ar_order+1)-arma_order_diff,max_slid_ar_order);
   else
      try_ar_order = 0;
      try_slid_ar_order = 0;
      min_slid_ar_order = 0;
   end
   
   %Determine a minimum set of AR parameter vectors, as needed for the preparations
   [cand_ar_orig_order,ar_entry] = ...
      sort([sel_ar_order pred_ar_order min_ar_order min_slid_ar_order try_slid_ar_order]);
   equal_entry = zeros(1,5);
   [dummy,redirect] = sort(ar_entry);
   equal_counter = 0;
   for i = 2:5
      if isequal(cand_ar_orig_order(i),cand_ar_orig_order(i-1));
         equal_counter = equal_counter+1;
         equal_entry(i) = i;
         index = find(max(0,redirect-i+equal_counter));
         redirect(index) = redirect(index)-1;
      end
   end
   cand_ar_orig_order(find(equal_entry)) = [];
   ar_entry = redirect;
   ar_orig = rc2arset(rc,cand_ar_orig_order,ASAcontrol);
   
   %Compute the first ARMA model, which equals an
   %AR model of the minimum candidate AR order
   sel_index = 1;
   ar = ar_orig{ar_entry(3)};
   ma = 1;
   res = var*prod(1-rc(2:min_ar_order+1).^2);
   if min_incl
      req_counter = req_counter+1;
      if ASAglob_mean_adj
        %gic3_temp = log(res)+3*(min_ma_order+min_ar_order+1)/n_obs;
        gic3_temp = log(res)+penalty_factor*(min_ma_order+min_ar_order+Keff_model/Keff_mean)/n_obs/Keff_model;
      else
          gic3_temp = log(res)+penalty_factor*(min_ma_order+min_ar_order)/n_obs/Keff_model;
      end
      if req_counter==2 & gic3_temp<gic3(sel_index)
         sel_index = req_counter;
         ar_sel = ar;
         ma_sel = ma;
      elseif req_counter==1
         ar_sel = ar;
         ma_sel = ma;
      end
      gic3(req_counter) = gic3_temp;
      if ASAglob_mean_adj
         %pe_est(req_counter) = res*...
         %   (n_obs+min_ar_order+min_ma_order+1)/(n_obs-min_ar_order-min_ma_order-1);
         pe_est(req_counter) = res*...
            (n_obs*Keff_model+min_ar_order+min_ma_order+Keff_model/Keff_mean)/(n_obs*Keff_model-min_ar_order-min_ma_order-Keff_model/Keff_mean);
      else
         %pe_est(req_counter) = res*...
         %   (n_obs+min_ar_order+min_ma_order)/(n_obs-min_ar_order-min_ma_order);
         pe_est(req_counter) = res*...
            (n_obs*Keff_model+min_ar_order+min_ma_order)/(n_obs*Keff_model-min_ar_order-min_ma_order);
      end
   end
   
   ar_pred = ar_orig{ar_entry(2)};
   l_ar_pred = length(ar_pred);
   % We could put all the equations in the initial estimation
   % below that originate from the different segments into 1 matrix. It is 
   % simpler to calculate the correlation matrices for each
   % segment separately and average them.
   Pred_Sig=zeros(fix(n_obs/2),n_segs);
   E=zeros(size(sig));
   RXX=zeros(max_ar_order+1,n_segs);
   REE=zeros(max_ar_order+1,n_segs);
   RXE=zeros(2*max_ar_order+1,n_segs);
   for k=1:n_segs
       pred_sig_rev = armafilter(zeros(fix(n_obs/2),1),ar_pred,1,...
      sig(end:-1:1,k),convolrev(sig(:,k),ar_pred,1,l_ar_pred,ASAcontrol),ASAcontrol);
       pred_sig = pred_sig_rev(end:-1:1);
       Pred_Sig(:,k)=pred_sig;
       pred_e = armafilter(pred_sig,1,ar_orig{ar_entry(2)},ASAcontrol);
       e = armafilter(sig(:,k),1,ar_orig{ar_entry(2)},pred_e,pred_sig,ASAcontrol);
       E(:,k)=e;
       RXX(:,k)=convolrev(sig(:,k),max_ar_order,ASAcontrol);
       REE(:,k) = convolrev(e,max_ar_order,ASAcontrol);
       RXE(:,k) = convolrev(sig(:,k),e,n_obs-max_ar_order,n_obs+max_ar_order,ASAcontrol);
   end
   s_R = 2*max_ar_order-arma_order_diff;
   
   %Try Durbin's first method to find an initial AR estimate at the start order
   if try_ar_order
      req_counter = req_counter+1;
      ar_order_start = try_ar_order;
      ph = ar_order_start+1;
      Lrtr=2*ar_order_start-arma_order_diff;
      RTR=zeros(Lrtr,Lrtr,2);
      RTb=zeros(Lrtr,1,2);
      for k=1:n_segs
        [rtr rtb]= durbinprep(sig(:,k),E(:,k),ar_order_start,ar_order_start-arma_order_diff,RXX(:,k),REE(:,k),RXE(:,k));
        RTR(:,:,1)=RTR(:,:,1)+rtr;
        RTb(:,:,1)=RTb(:,:,1)+rtb;
      end
      lastwarn('')
      warning off
      arma_ini = (RTR(:,:,1)\RTb(:,:,1))';
      eval(warn_recov_comm);
      if lastwarn
         type_ar_ini = 1;
         ar_order_start = min_ar_order+1;
         ar_slid = ar_orig{ar_entry(4)};
         clear('RTR','RTb')
      else
         type_ar_ini = 0;
         ar_slid = ar_orig{ar_entry(5)};
         ar_ini = [1 arma_ini(1:ar_order_start)];
         [ar_ini rc_ini] = ar2arset(ar_ini,ASAcontrol);
         if any(abs(rc_ini)>1)
            type_ar_ini = 1;
            ar_order_start = min_ar_order+1;
            ar_slid = ar_orig{ar_entry(4)};
            clear('RTR','RTb')
         end
      end
      n_iter = 1;
   else
      ar_order_start = max_ar_order+1;
   end
   
   %Estimation loop initializations
   reset_type=1;
   c = 0;
   ma_order=ar_order_start-arma_order_diff;
   slid_ar_order = 3*sel_ar_order+ar_order_start+ma_order;
   if slid_ar_order > max_slid_ar_order
      slid_ar_order = max_slid_ar_order;
   end
   
   %Estimation procedure
   %--------------------
   
   for ar_order = ar_order_start:max_ar_order  
      arma_order = ar_order+ma_order;
      
      %The first method of computing an initial
      %estimate of the AR parameters of the ARMA model
      if type_ar_ini == 1 
         if c, co = 2; cn = 1;
         else co = 1; cn = 2;
         end        
         if ar_order==ar_order_start
             Lrtr=ar_order+ma_order;
             RTR(:,:,cn)=zeros(Lrtr,Lrtr);
             RTb(:,:,cn)=zeros(Lrtr,1);
             for k=1:n_segs
                [rtr,rtb] = durbinprep(sig(:,k),E(:,k),ar_order,ma_order,RXX(:,k),REE(:,k),RXE(:,k));
                RTR(:,:,cn)=RTR(:,:,cn)+rtr;
                RTb(:,:,cn)=RTb(:,:,cn)+rtb;
             end
         else
            ph = ar_order+1;
            pn = ar_order;
            po = ar_order-1;
            qn = ma_order;
            qo = ma_order-1;
            pqo = po+qo;
            pqn = pn+qo;
            RTR(1:s_R,1:s_R,cn) = zeros(s_R,s_R);
            RTb(1:s_R,1,cn) = zeros(s_R,1);
            RTR(1:po,1:po,cn) = RTR(1:po,1:po,co);
            RTR(ph:pqn,ph:pqn,cn) = RTR(pn:pqo,pn:pqo,co);
            RTR(1:po,ph:pqn,cn) = RTR(1:po,pn:pqo,co); 
            RTR(ph:pqn,1:po,cn) = RTR(pn:pqo,1:po,co);
            for k=1:n_segs
                [RTR_add,RTb_add] = durbinprep(sig(:,k),E(:,k),ar_order,ma_order,RXX(:,k),REE(:,k),RXE(:,k),po,qo);
                RTR(1:arma_order,1:arma_order,cn) = RTR(1:arma_order,1:arma_order,cn)+RTR_add;
                RTb(1:arma_order,1,cn) = RTb(1:arma_order,1,cn)+RTb_add;
            end
         end
         c=~c;
         lastwarn('')
         warning off
         arma_ini = (RTR(1:arma_order,1:arma_order,cn)\RTb(1:arma_order,1,cn))';
         eval(warn_recov_comm);
         if lastwarn
            type_ar_ini = 2;
            reset_type = 0;
         else
            ar_ini = [1 arma_ini(1:ar_order)];
            [ar_ini rc_ini] = ar2arset(ar_ini,ASAcontrol);
            if any(abs(rc_ini)>1)
               type_ar_ini = 2;
               reset_type = 1;
            end
         end
         n_iter = 1;
      end
      
      %The second method to compute an initial
      %estimate of the AR parameters of the ARMA model
      if type_ar_ini == 2
          m=zeros(2,2);y=zeros(2,1);
          for k=1:n_segs
            pred_x = armafilter(Pred_Sig(:,k),ma,ar,ASAcontrol);
            x = armafilter(sig(:,k),ma,ar,pred_x,Pred_Sig(:,k),ASAcontrol);
            [mk,yk] = durbinprep(x,E(:,k),1,1,convolrev(x,1,ASAcontrol),...
            REE(1:2,k),convolrev(x,E(:,k),n_obs-1,n_obs+1,ASAcontrol));
            m=m+mk;y=y+yk;
          end
         arma_1_1 = (m\y)';
         ar_1 = [1 arma_1_1(1)];
         if any(abs(ar_1)>1)
            ar_ini = ar;
         else
            ar_ini = convol(ar,ar_1,ASAcontrol);
            rc_ini=poly2rc(ar_ini);
            if any(~isfinite(rc_ini)) || any(abs(rc_ini)>1)
                ar_ini=ConvolImpRes2Par(ar,ar_1,length(ar_ini)-1,1024);
            end
         end
         n_iter = 3+fix(log10(n_obs));
      end
      
      %Definitive ARMA model estimation
      k=1;
      while k <= n_iter
         ar_i = ar_ini;
         ar_interm = deconvol(ar_slid,ar_i,ASAcontrol);
         rc_interm=poly2rc(ar_interm);
         if any(~isfinite(rc_interm)) || any(abs(rc_interm)>1)
            ar_interm=DeconvolImpRes2Par(ar_slid,ar_i,length(ar_interm)-1,1024);
         end
         % fall back onto second method?
%          if any(~isfinite(ar_interm))
%              m=zeros(2,2);y=zeros(2,1);
%              for k=1:n_segs
%                 pred_x = armafilter(Pred_Sig(:,k),ma,ar,ASAcontrol);
%                 x = armafilter(sig(:,k),ma,ar,pred_x,Pred_Sig(:,k),ASAcontrol);
%                 [mk,yk] = durbinprep(x,E(:,k),1,1,convolrev(x,1,ASAcontrol),...
%                 REE(1:2,k),convolrev(x,E(:,k),n_obs-1,n_obs+1,ASAcontrol));
%                 m=m+mk;y=y+yk;
%              end
%              arma_1_1 = (m\y)';
%              ar_1 = [1 arma_1_1(1)];
%              if any(abs(ar_1)>1)
%                 ar_ini = ar;
%              else
%                 ar_ini = convol(ar,ar_1,ASAcontrol);
%              end
%              n_iter = 3+fix(log10(n_obs))+1;
%              ar_i = ar_ini;
%              ar_interm = deconvol(ar_slid,ar_i,ASAcontrol);
%              ma_i= arburg([zeros(1,ma_order) ar_interm/max(ar_interm)],ma_order);
%              k=k+1;
%              continue
%          end
         %ar_corr = convolrev(ar_interm,ma_order,ASAcontrol);
         %ma_i = cov2arset(ar_corr,ASAcontrol)
         ma_i= arburg([zeros(1,ma_order) ar_interm/max(ar_interm)],ma_order);
         %rc=poly2rc(ma_i);
         %if any(abs(rc)>1)
         %    ma_i=arburg([zeros(1,ma_order) ar_interm],ma_order);
         %end
         ar_i = convol(ar_slid,ma_i,ASAcontrol);
         ar_i = ar2arset(ar_i,ar_order,ASAcontrol);
         ar_i = ar_i{1};
         rc_i=poly2rc(ar_i);
         if any(~isfinite(rc_i)) || any(abs(rc_i)>1)
            ar_i=ConvolImpRes2Par(ar_slid,ma_i,ar_order,1024);
         end
         ar_ini = ar_i;
         k=k+1;
      end
      
      ar = ar_i;
      ma = ma_i;
      
      %Order selection
      if cand_ar_order(req_counter)==ar_order
          res=0;
          for k=1:n_segs
             e_interm = armafilter(sig(:,k),ma,ar,filter(ar,ma,Pred_Sig(:,k)),Pred_Sig(:,k),ASAcontrol);
             res = res+ e_interm'*e_interm/n_obs/n_segs;
          end
         if ASAglob_mean_adj
            gic3_temp = log(res)+penalty_factor*(ar_order+ma_order+Keff_model/Keff_mean)/n_obs/Keff_model;
         else
            gic3_temp = log(res)+penalty_factor*(ar_order+ma_order)/n_obs/Keff_model;
         end
         gic3(req_counter) = gic3_temp;
         if gic3_temp <= gic3(sel_index)
            sel_index = req_counter;
            ar_sel = ar;
            ma_sel = ma;
         end
         if ASAglob_mean_adj
            pe_est(req_counter) = res*(n_obs*Keff_model+ar_order+ma_order+Keff_model/Keff_mean)/(n_obs*Keff_model-ar_order-ma_order-Keff_model/Keff_mean);
         else
            pe_est(req_counter) = res*(n_obs*Keff_model+ar_order+ma_order)/(n_obs*Keff_model-ar_order-ma_order);
         end
         req_counter = req_counter+1;            
      end
      
      %Determine the sliding AR parameter vector
      for k=[1 2]
         if slid_ar_order <= max_slid_ar_order-1
            slid_ar_order = slid_ar_order+1;
            rc_temp = rc(slid_ar_order+1);
            ar_slid(2:slid_ar_order) = ...
               ar_slid(2:slid_ar_order)+rc_temp*ar_slid(slid_ar_order:-1:2);
            ar_slid(slid_ar_order+1) = rc_temp;
         end
      end
      
      ma_order = ma_order+1;
      
      if reset_type
         type_ar_ini = 1;
      end
   end
end
  
%Arranging output arguments
%--------------------------

%Retrieve the parameters of the proper model
ar = ar_sel;
ma = ma_sel;

%Assign reflectioncoefficients to ASAglob_rc, in order
%to make them available for other ARMASA functions
if ~isempty(rc)
   ASAglob_rc = rc;
end

%Generate a structure variable ASAsellog to report
%the selection process
if nargout>2
   ASAsellog.funct_name = mfilename;
   ASAsellog.funct_version = ASAcontrol.is_version;
   ASAsellog.date_time = [datestr(ASAdate,8) 32 datestr(ASAdate,0)];
   ASAsellog.comp_time = etime(clock,ASAtime);
   ASAsellog.ar = ar;
   ASAsellog.ma = ma;
   ASAsellog.ar_sel = ar_orig{ar_entry(1)};
   ASAsellog.mean_adj = ASAglob_mean_adj;
   ASAsellog.cand_ar_order = cand_ar_order;
   ASAsellog.arma_order_diff = arma_order_diff;
   ASAsellog.gic3 = gic3;
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

%Helper functions
%=================================================================================================

function par=ConvolImpRes2Par(ar1,ar2,P,L)
%p1=length(ar1);p2=length(ar2);
% Impulsresponses of ar1 and ar2
h1=filter(1,ar1,[1 zeros(1,L)]);
h2=filter(1,ar2,[1 zeros(1,L)]);
h12=conv(h1,h2);
par=arburg([zeros(1,P) h12],P);

function par=DeconvolImpRes2Par(ar1,ar2,P,L)
%p1=length(ar1);p2=length(ar2);
% Impulsresponses of ar1 and ar2
h1=filter(1,ar1,[1 zeros(1,L)]);
%h2=filter(1,ar2,[1 zeros(1,L)]);
h12=conv(h1,ar2);
par=arburg([zeros(1,P) h12],P);
      

function [RTR,RTb,covx,covepshat,kcov] = durbinprep(x,epshat,k,l,covx,covepshat,kcov,ko,lo)

%function [RTR,RTb,covx,covepshat,kcov] = durbinprep(x,epshat,k,l,covx,covepshat,kcov,ko,lo)
%
% Prepares the estimation of the ARMA-parameters with Durbin's method.
%
% Durbinprep has 2 basic modes:
%
% 1] Calculation of RTR from scratch
% Call: [RTR,RTb,covx,covepshat,kcov] = durbinprep(x,epshat,k,l,covx,covepshat,kcov)
%       The RTR and RTb are calculated.
%
% 2] Increment to results from Durbinprep from an ARMA(ko,lo)-model.
% Call: [RTR,RTb,covx,covepshat,kcov] = durbinprep(x,epshat,k,l,covx,covepshat,kcov,ko,lo)
%       Increments to the old RTR and RTb, RTRo and RTbo are calculated.
%
% Used in: sig2arma.

% Background:
%
% Durbin's ARMA method is the least-sqares solution for 
% x[n] + a(1)x[n-1] ...+a(k)x[n-k] = e[n] + b(1)e[n-1] ... + a(p)e[n-l]
%
% or
%
% a(1)x[n-1] ...+a(k)x[n-k] - b(1)e[n-1] ... - a(p)e[n-l] = x[n] - e[n]
%
% This can be re-written as:
% AO                        - BE                          = b  (least-squares)
%
% Combing the AR and the MA paramters in the vector pareps = {A | B) and
% the Observations O and Epsilons E in R = [0 E], this yields
%
% R*pareps = b (least-squares)
%
% So the parameters are found with pareps=RTR\RTb (see SIG2ARMA).
%

%programmer: S. de Waele

%Initialization
nobs = length(x);
m = max(k,l); afm = k+l;
xa = [zeros(m,1); x; zeros(m,1)];
epshata = [zeros(m,1); epshat; zeros(m,1)];
mmm = length(kcov); mmm = floor(mmm/2);

if nargin==9, adding = 1; mo = max(ko,lo); else adding = 0; end 

i1 = 0:m-1; i2 = 0:-1:1-m;
ce = toeplitz(epshata(i1+m),epshata(i2+m));				co = toeplitz(xa(i1+m),xa(i2+m));
de = toeplitz(epshata(nobs+i1+m),epshata(nobs+i2+m));	do = toeplitz(xa(nobs+i1+m),xa(nobs+i2+m));

if ~adding,
   %--------------------------------------------------------------------------
   %Calculation of RTR from scratch
   %--------------------------------------------------------------------------
	%epshat's with x
   ctot = ce'*co+de'*do; ctot = ctot(1:l,1:k);
   ETO = -(toeplitz(kcov((0:-1:-(l-1))+mmm+1),kcov((1:k)+mmm))-ctot);
   
	%x with x
   co = co(:,1:k); do = do(:,1:k); 
   ctot = co'*co+do'*do;
	OTO = toeplitz(covx(1:k))-ctot;
   
	%epshats with epshat
	ce = ce(:,1:l); de = de(:,1:l); 
	ctot = ce'*ce+de'*de;
	ETE = toeplitz(covepshat(1:l))-ctot;   
   
else
   %--------------------------------------------------------------------------   
   %Calculation of increments on RTR from previous model
   %--------------------------------------------------------------------------
   % epshat's with x
   ETO = zeros(l,k);
	i1 = mo:(m-1); i2 = mo:-1:mo-(ko-1);
	toevo = co(mo+1:m,1:mo);
	toeve = ce(mo+1:m,1:mo);   
   toev = toeve'*toevo;
   ETO(1:lo,1:ko) = +toev(1:lo,1:ko);
   
   cet = ce(:,lo+1:l); det = de(:,lo+1:l);
   ctot = cet'*co+det'*do; ctot = ctot(1:l-lo,1:k);
   if l-lo == 1,
      ETO(lo+1:l,:) = -(kcov(-toeplitz(lo:l-1,lo:-1:lo-(k-1))+mmm+1)' -ctot);
   else
      ETO(lo+1:l,:) = -(kcov(-toeplitz(lo:l-1,lo:-1:lo-(k-1))+mmm+1) -ctot);
   end
   cot = co(:,ko+1:k); dot = do(:,ko+1:k);
   ctot = ce'*cot+de'*dot;	ctot = ctot(1:lo,:);
   ETO(1:lo,ko+1:k) = -(toeplitz(kcov((ko:-1:(ko-lo+1))+mmm+1),kcov((ko:k-1)+mmm+1))-ctot);
   
   % x with x
   c = co(:,1:k); d = do(:,1:k); 
	OTO = zeros(k);
	toev = c(mo+1:m,1:ko);
   OTO(1:ko,1:ko) = -toev'*toev;
   ct = c(:,ko+1:k); dt = d(:,ko+1:k);
   ctot = ct'*c+dt'*d;
   if k-ko == 1,
   	OTO(ko+1:k,:) = covx(abs(toeplitz(ko:k-1,ko:-1:ko-(k-1)))+1)'-ctot;
   else
  		OTO(ko+1:k,:) = covx(abs(toeplitz(ko:k-1,ko:-1:ko-(k-1)))+1)-ctot;
   end
   OTO(:,ko+1:k) = OTO(ko+1:k,:)';
   
	% epshats with epshat
	c = ce(:,1:l); d = de(:,1:l); 
   ETE = zeros(l);
	toev = ce(mo+1:m,1:lo);
   ETE(1:lo,1:lo) = -toev'*toev;
   ct = c(:,lo+1:l); dt = d(:,lo+1:l);
	ctot = ct'*c+dt'*d;
   if l-lo == 1,
   	ETE(lo+1:l,:) = covepshat(abs(toeplitz(lo:l-1,lo:-1:lo-(l-1)))+1)'-ctot;
   else
   	ETE(lo+1:l,:) = covepshat(abs(toeplitz(lo:l-1,lo:-1:lo-(l-1)))+1)-ctot;
   end
   ETE(:,lo+1:l) = ETE(lo+1:l,:)';
end

%RTR 
RTR = [OTO ETO'; ...
       ETO ETE];
    
%RTb     
c = zeros(k+l,1);
for s = 1:k,
   t = 1:(m-s);
   c(s) = -(sum(x(t).*x(t+s))-sum(x(t).*epshat(t+s)));
end
for s = 1:l,
   t = 1:(m-s);
   c(k+s) = sum(epshat(t).*x(t+s))-sum(epshat(t).*epshat(t+s));   
end 
RTb = [-covx(2:k+1)+kcov((2:k+1)+mmm); ...
		  kcov((-1:-1:-l)+mmm+1)-covepshat(2:l+1)];
RTb = RTb-c;

%Program history
%======================================================================
%
% Version                Programmer(s)          E-mail address
% -------                -------------          --------------
% former versions        P.M.T. Broersen        p.m.t.broersen@tn.tudelft.nl
%                        S. de Waele            waele@tn.tudelft.nl
% [2000 12 30 20 0 0]    W. Wunderink           wwunderink01@freeler.nl
% [2001  1  7 12 0 0]    W. Wunderink           wwunderink01@freeler.nl
%                        S. de Waele            waele@tn.tudelft.nl
% [2003  4  2 15 0 0]    W. Wunderink           wwunderink01@freeler.nl