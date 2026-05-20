% Perform simulations comparing sig2ar, sig2ma, sig2arma, armasel and their
% versions for segments on pure ar, ma, and arma processes. We filter the
% data in 2 directions, to test our estimator for the efffective number of
% segments as well.
% Define model_type in the workspace to choose between ar, ma, or arma.
%model_type=1;%1=AR, 2=MA, 3=ARMA
% Penalty_factor used in versions for segments
penalty_factor=3;
% Reflection coefficients are generated with k(m)=b^m.
% Vector of b-values used in the experiments
bvec=-0.8:0.4:0.8;
ARMAorder=5;%Simulate ar(p), ma(p), arma(p,p-3)
t1=1:ARMAorder;
t2=1:ARMAorder-3;
% Segment lenghts
%N=[30 75 200];% 250 500];
N=[75 200];
% To let initial conditions decay
Ninit=100;
% Number of segments
%S=[1 2 5 10 20];% 20 50];
S=[1 5 20];
Ln=length(N);
Ls=length(S);
Lb=length(bvec);
% number of simulation runs to determine averages
if model_type==1
    iter=1000;%ik had er al 900 gedaan voor AR
else
    iter=1000;
end
% Tables with average model errors
% 1-segment versions sig2ar, sig2ma, and sig2arma
ARtable_1=zeros(Lb,Ln,Ls);
MAtable_1=ARtable_1;
ARMAtable_1=ARtable_1;
% Selecting between AR/MA/ARMA
SELtable_1=ARtable_1;
% For the 1-segment versions, the autocorrelation functions of segment
% models are averaged. Model Error is calculated by transforming to a long
% AR model of order long_order.
long_order=100;
% multiple segment version (uses different penalty factors for order
% selection)
ARtable_s=zeros(Lb,Ln,Ls);
MAtable_s=zeros(Lb,Ln,Ls);
ARMAtable_s=zeros(Lb,Ln,Ls);
SELtable_s=zeros(Lb,Ln,Ls);
for b=1:Lb
    rcar=bvec(b).^t1;rcma=(-bvec(b)).^t2;
    parar=rc2poly(rcar);parma=rc2poly(rcma);
    if model_type==1
        truear=parar;truema=1;
    elseif model_type==2
        truear=1;truema=parar;
    else
        truear=parar;truema=parma;
    end
    for n=1:Ln
        %tic
        Nn=N(n);
        data=zeros(Nn,max(S));
        for k=1:iter
            % Generate segments of ar, ma, or arma
            e=randn(max(S)+Ninit,Nn+Ninit);
            for k2=1:max(S)
                x1=filter(truema,truear,e);
                data=filter(truema,truear,x1');
                data=data(Ninit+1:end,Ninit+1:end);%Nn x max(S)
            end%k2=1:max(S)
            % versions on segments
            for s=1:Ls
                % Initialize covariance functions
                cov_ar=zeros(1,long_order+1);cov_ma=cov_ar;
                cov_arma=cov_ar;cov_armasel=cov_ar;
                % 1 segment versions: average segment covariance functions
                for k2=1:S(s)
                    [armasel_ar_1,armasel_ma_1,sellog]=armasel_s(data(:,k2),[],[],[],[],[],1,1);
                    ar_1=sellog.ar.ar;
                    ma_1=sellog.ma.ma;
                    arma_ar_1=sellog.arma.ar;arma_ma_1=sellog.arma.ma;
                    vars=var(data(:,k2));
                    cov_ar=cov_ar+vars*arma2cor(ar_1,1,long_order);
                    cov_ma=cov_ma+vars*arma2cor(1,ma_1,long_order);
                    cov_arma=cov_arma+vars*arma2cor(arma_ar_1,arma_ma_1,long_order);
                    cov_armasel=cov_armasel+vars*arma2cor(armasel_ar_1,armasel_ma_1,long_order);
                end
                % Compute model errors
                ARtable_1(b,n,s)=ARtable_1(b,n,s)+moderr_cov2ar(cov_ar,truear,truema,Nn,long_order);
                MAtable_1(b,n,s)=MAtable_1(b,n,s)+moderr_cov2ar(cov_ma,truear,truema,Nn,long_order);
                ARMAtable_1(b,n,s)=ARMAtable_1(b,n,s)+moderr_cov2ar(cov_arma,truear,truema,Nn,long_order);
                SELtable_1(b,n,s)=SELtable_1(b,n,s)+moderr_cov2ar(cov_armasel,truear,truema,Nn,long_order);
                % All segments at once
                [armasel_ar_s,armasel_ma_s,sellog]=armasel_s(data(:,1:S(s)),[],[],[],[],[],[0 1],0);    
                ar_s=sellog.ar.ar;
                ma_s=sellog.ma.ma;
                arma_ar_s=sellog.arma.ar;arma_ma_s=sellog.arma.ma;
                % Model errors
                ARtable_s(b,n,s)=ARtable_s(b,n,s)+moderr(ar_s,1,truear,truema,Nn);
                MAtable_s(b,n,s)=MAtable_s(b,n,s)+moderr(1,ma_s,truear,truema,Nn);
                ARMAtable_s(b,n,s)=ARMAtable_s(b,n,s)+moderr(arma_ar_s,arma_ma_s,truear,truema,Nn);
                SELtable_s(b,n,s)=SELtable_s(b,n,s)+moderr(armasel_ar_s,armasel_ma_s,truear,truema,Nn);
            end%s=1:Ls
        end%k=1:iter
        %toc,disp([b n])
    end%n=1:Ln
    disp(b)
end%b=1:Lb
ARtable_1=ARtable_1/iter;
MAtable_1=MAtable_1/iter;
ARMAtable_1=ARMAtable_1/iter;
SELtable_1=SELtable_1/iter;
ARtable_s=ARtable_s/iter;
MAtable_s=MAtable_s/iter;
ARMAtable_s=ARMAtable_s/iter;
SELtable_s=SELtable_s/iter;
clear data t b n s k k2 e x1
% cmt='Deze versie gebruikt armasel ipv sig2ar for Seff schatten.';
% if model_type==1
%     eval(['save Compare_s3_' num2str(model_type) 'armasel_rest'])
% else
%     eval(['save Compare_s3_' num2str(model_type) 'armasel'])
% end
