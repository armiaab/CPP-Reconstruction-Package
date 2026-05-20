function [ma,ASAsellog,ASAcontrol] = sig2ma_s(sig,cand_order,last,estKeff_model,estKeff_mean)
%SIG2MA MA model identification from multiple segments
%   [MA,SELLOG] = SIG2MA(SIG) estimates moving average models from the 
%   data vector SIG and selects a model with optimal predictive 
%   qualities. The selected model is returned in the parameter vector MA. 
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
%   SIG2MA(SIG,CAND_ORDER) selects only from candidate models whose 
%   orders are entered in CAND_ORDER. CAND_ORDER must either be a row of 
%   ascending orders, or a single order (in which case no true order 
%   selection is performed).
%   
%   Without user intervention, the mean of SIG is subtracted from the 
%   data. To control the subtraction of the mean, see the help topics on 
%   ASAGLOB_SUBTR_MEAN and ASAGLOB_MEAN_ADJ.
%     
%   SIG2MA is an ARMASA main function.
%   
%   See also: SIG2AR, SIG2ARMA, ARMASEL.

%   References: P. M. T. Broersen, Autoregressive Model Orders for
%               Durbin's MA and ARMA estimators, IEEE Transactions on
%               Signal Processing, Vol. 48, No. 8, August 2000,
%               pp. 2454-2457.

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
ASAcontrol.is_version = [2000 12 30 20 0 0];
%and its compatability with versions down to,
ASAcontrol.comp_version = [2000 12 30 20 0 0];

%This function calls other functions of the ARMASA
%toolbox. The versions of these other functions must
%be greater than or equal to:
ASAcontrol.req_version.burg = [2000 12 30 20 0 0];
ASAcontrol.req_version.cic = [2000 12 30 20 0 0];
ASAcontrol.req_version.rc2arset = [2000 12 30 20 0 0];
ASAcontrol.req_version.cov2arset = [2000 12 30 20 0 0];
ASAcontrol.req_version.armafilter = [2000 12 12 14 0 0];
ASAcontrol.req_version.convol = [2000 12 6 12 17 20];
ASAcontrol.req_version.convolrev = [2000 12 6 12 17 20];

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
   cov2arset(ASAcontrol);
   armafilter(ASAcontrol);
   convol(ASAcontrol);
   convolrev(ASAcontrol);
end

if ~any(strcmp(fieldnames(ASAcontrol),'run')) | ASAcontrol.run
      %Run the computational kernel
   ASAcontrol.run = 1;
   ASAcontrol.version_chk = 0;
   ASAcontrol.error_chk = 0;
   ASAtime = clock;
   ASAdate = now;

%Main   
%================================================================================================
  
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

%n_obs = length(sig);
ar_stack = cell(4,1);
ar_entry = ones(1,4);
rc = [];
ma_sel = 1;

%Combined determination of the maximum candidate MA
%order and the max. candidate sliding AR order 
%--------------------------------------------------

def_max_ma_order = min(fix(n_obs/5),fix(80*log10(n_obs)));
if def_max_ma_order > 400; 
   def_max_ma_order = 400;
end

if isempty(cand_order)
   cand_order = 0:def_max_ma_order;
end
max_ma_order = cand_order(end);

if max_ma_order <= def_max_ma_order
   max_slid_ar_order = fix(2.5*def_max_ma_order);
else
   max_slid_ar_order = fix(2.5*max_ma_order);
   if max_slid_ar_order > n_obs-1
      max_slid_ar_order = n_obs-1;
   end
end

%Preparations for the estimation procedure
%-----------------------------------------
penalty_factor=3;
l_cand_order = length(cand_order);
ma = zeros(1,max_ma_order);
%var = sig'*sig/n_obs;
var=mean(mean(sig.^2));
if cand_order(1)==0
   zero_incl = 1;
   gic3 = zeros(1,l_cand_order);
   %gic3(1) = log(var)+3/n_obs/Keff_model;
   %gic3(1) = log(var)+3/n_obs/Keff;
   pe_est = zeros(1,l_cand_order);
   if ASAglob_mean_adj
      pe_est(1) = var*(n_obs*Keff_mean+1)/(n_obs*Keff_mean-1);
      gic3(1) = log(var)+penalty_factor/n_obs/Keff_mean;
   else
      pe_est(1) = var;
      gic3(1) = log(var);
   end
else
   zero_incl = 0;
   cand_order = [0 cand_order];
   gic3 = zeros(1,l_cand_order+1);
   gic3(1) = inf;
   pe_est = zeros(1,l_cand_order+1);
end

if max_ma_order > 0
   %Conditioning AR orders to the previously selected AR model
   if isequal(ASAglob_ar_cond,1) & ~isempty(ASAglob_ar)
      ar = ASAglob_ar;
      sel_ar_order = length(ar)-1;
      if  2*sel_ar_order+max_ma_order < max_slid_ar_order
         max_slid_ar_order = 2*sel_ar_order+max_ma_order;
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
         rc = [ASAglob_rc burg(ASAglob_final_f, ...
               ASAglob_final_b,max_slid_ar_order+1-l_rc,ASAcontrol)];
      end
   else
      rc = burg_s_modified(sig,max_slid_ar_order,ASAcontrol);
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
   
   min_ma_order = max(1,cand_order(1));
   
   slid_ar_order = 2*sel_ar_order+min_ma_order;
   if slid_ar_order > max_slid_ar_order
      slid_ar_order = max_slid_ar_order;
   elseif slid_ar_order < 3
      slid_ar_order = min(3,max_slid_ar_order);
   end
   
   pred_ar_order = min(3*sel_ar_order+min(9,1+fix(n_obs/10)),max_slid_ar_order);
   
   %Determine a minimum set of AR parameter vectors, as needed for the preparations
   [cand_ar_order,ar_entry] = sort([sel_ar_order pred_ar_order slid_ar_order max_slid_ar_order]);
   equal_entry = zeros(1,4);
   [~,redirect] = sort(ar_entry);
   equal_counter = 0;
   for i = 2:4
      if isequal(cand_ar_order(i),cand_ar_order(i-1));
         equal_counter = equal_counter+1;
         equal_entry(i) = i;
         index = find(max(0,redirect-i+equal_counter));
         redirect(index) = redirect(index)-1;
      end
   end
   cand_ar_order(find(equal_entry)) = [];
   ar_entry = redirect;
   ar_stack = rc2arset(rc,cand_ar_order,ASAcontrol);
   
   ar_pred = ar_stack{ar_entry(2)};
   l_ar_pred = length(ar_pred);
   l_pred_sig = fix(n_obs/2);
   pred_sig_rev=zeros(l_pred_sig,n_segs);
   for k=1:n_segs
       pred_sig_rev(:,k) = armafilter(zeros(l_pred_sig,1),ar_pred,1,...
      sig(end:-1:1,k),convolrev(sig(:,k),ar_pred,1,l_ar_pred,ASAcontrol),ASAcontrol);
   end
   pred_sig = pred_sig_rev(end:-1:1,:);
   
   counter = 2;
   req_counter = 2;
   sel_index = 1;
   ar_slid = zeros(1,max_slid_ar_order+1); 
   ar_slid(1:slid_ar_order+1) = ar_stack{ar_entry(3)};
   
   %Estimation procedure and model order selection
   %----------------------------------------------
   
   for order = min_ma_order:max_ma_order
      if cand_order(req_counter)==order
         %ar_corr = convolrev(ar_slid(1:slid_ar_order+1),order,ASAcontrol);
         %ma = cov2arset(ar_corr,ASAcontrol)
         ma = arburg([zeros(1,order) ar_slid(1:slid_ar_order+1)]',order);
         res=0;
         for k=1:n_segs
            e = armafilter(sig(:,k),ma,1,filter(1,ma,pred_sig(:,k)),pred_sig(:,k),ASAcontrol);
            res = res+e'*e/n_obs/n_segs;
         end
         if ASAglob_mean_adj
            %gic3_temp = log(res)+3*(order+1)/n_obs;
            gic3_temp = log(res)+penalty_factor*(order+Keff_model/Keff_mean)/n_obs/Keff_model;
         else
            gic3_temp = log(res)+penalty_factor*order/n_obs/Keff_model;
         end
         gic3(req_counter) = gic3_temp;
         if gic3_temp < gic3(sel_index)
           sel_index = req_counter;
           ma_sel = ma;
         end
         if ASAglob_mean_adj
            pe_est(req_counter) = res*(n_obs*Keff_model+order+Keff_model/Keff_mean)/(n_obs*Keff_model-order-Keff_model/Keff_mean);
         else
            pe_est(req_counter) = res*(n_obs*Keff_model+order)/(n_obs*Keff_model-order);
         end
         req_counter = req_counter+1;            
      end
      
      if slid_ar_order < max_slid_ar_order
         slid_ar_order = slid_ar_order+1;
         rc_temp = rc(slid_ar_order+1);
         ar_slid(2:slid_ar_order) = ar_slid(2:slid_ar_order)+rc_temp*ar_slid(slid_ar_order:-1:2);
         ar_slid(slid_ar_order+1) = rc_temp;
      end
      
      counter = counter+1;
   end
end

%Arranging output arguments
%--------------------------

ma = ma_sel;

if ~zero_incl
   gic3(1) = [];
   pe_est(1) = [];
   cand_order(1) =[];
end

if ~isempty(rc)
   ASAglob_rc = rc;
end

if nargout>1
   ASAsellog.funct_name = mfilename;
   ASAsellog.funct_version = ASAcontrol.is_version;
   ASAsellog.date_time = [datestr(ASAdate,8) 32 datestr(ASAdate,0)];
   ASAsellog.comp_time = etime(clock,ASAtime);
   ASAsellog.ma = ma_sel;
   ASAsellog.ar_sel = ar_stack{ar_entry(1)};
   ASAsellog.mean_adj = ASAglob_mean_adj;
   ASAsellog.cand_order = cand_order;
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
   ma = ASAcontrol;
   ASAcontrol = [];
end

%Program history
%======================================================================
%
% Version                Programmer(s)          E-mail address
% -------                -------------          --------------
% former versions        P.M.T. Broersen        p.m.t.broersen@tudelft.nl
% [2000 12 30 20 0 0]    W. Wunderink           wwunderink01@freeler.nl
