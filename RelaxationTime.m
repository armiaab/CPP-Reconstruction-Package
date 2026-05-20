function R=RelaxationTime(data)

% This code estimates the relaxation time of the data. To this end we need to know the autocorrelation function of data for the first few lags (we use minimum 20 lags). Note that here we do not need a high accuracy. In
% the case of replicate data we use the 'Burg algorithm for segments' which is one of the most accurate approaches base on the following reference (you should add the MATLAB package ARMASA into your working directory)

% Broersen, P. M. (2006). Automatic autocorrelation and spectral analysis. Springer Science & Business Media.

% Implemented in Matlab by Babak M.S. Ariani (gmail: m.shojaeiarani@gmail.com)


% NOTE1: The typical data (meaning a single dataset) should be supplied as an N*d matrix where d is the spatial (or state) dimension and N is the temporal length. Clearly, N>d almost always. So, we interpret max(size(data)) as the temporal length and main(size(data)) as the state length (numbreer of state variables). If data is a square matrix we consider size(data,1) as the temporal length
% NOTE2: The replicate data (meaning multiple datasets together) should be suplied as a cell array. The format of each cell (i.e., each replicate) should follow the points in NOTE1

if iscell(data)
    A=data{1};
    D=min(size(A));  %state (or spatial) dimension
else
    D=min(size(data));
end

switch D
    case 1
        if iscell(data)   %Replicate datasets (Univariate)
            if size(data,2)>size(data,1)
                data = data';
            end
            N=length(data);
            for i=1:N
                data{i}=data{i}(:);
            end
            [~,idx]=sort(cellfun(@length,data));
            data=data(idx);
            L=cellfun(@length,data);
            U=L.*(N:-1:1)';
            m=find(U==max(U),1,'last');
            data1=cell(1,N-m+1);
            for i=m:N
                data1{i-m+1}=data{i}(1:L(m));
            end
            data=cell2mat(data1);
            [ar_s,ma_s]=armasel_s(data,[],[],[],[],1,size(data,2),size(data,2));
            acf=arma2cor(ar_s,ma_s);
        else
            data=data(:);
            % In the following autocorrelation is calculated based on the Burg algorithm which is very accurate but slow. If you like you can use the following (but, we do not care about high accuracy here)
            %     [ar,ma,~,~]=armasel(data,[],[],[],[]);
            %     acf=arma2cor(ar,ma);
            acf=autocorr(data);
        end
        lag=0:1:(length(acf)-1);
        acf=acf(:)';
        [fit1,~]=createFit(lag,acf);
        R=1/fit1.c;
    otherwise
        if iscell(data)   %Replicate datasets (Multuvariate)
            N=length(data);
            for i=1:N
                if size(data{i},1)<size(data{i},2)
                    data{i}=data{i}.';
                end
            end
            
        else
            R=zeros(1,D);
            if size(data,1)<size(data,2)
                data=data.';
            end
            
            M=xcorr(data,'normalized');   % Matrix whose columns are auto- and cross-correlations
            idx=(size(M,1)-1)/2+1;   %this corresponds with the first lag we should start with
            M=M(idx:end,:);

            idx=1:D+1:D^2;   % index of corresponding auto-correlation columns in matrix M
            AC=M(:,idx);      % Mtrix whose columns are data autocorrelations
            for i=1:size(AC,2)
                idx=find(AC(:,i)<=0,1);
                if ~isempty(idx)
                    ac=AC(1:idx,i);
                    lag=(0:idx-1).';
                else
                    ac=AC(1:end,i);
                    lag=(0:length(ac)-1).';
                end
                [fit,~]=createFit(lag,ac);
                R(i)=1/fit.c;
            end
        end
end
end

%*************************************************************************************
%************************************SUBROUTINES**************************************
%*************************************************************************************
function [fitresult, gof] = createFit(lag, acf)
[xData, yData] = prepareCurveData( lag, acf );

% Set up fittype and options.
ft = fittype( 'exp(-c*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Lower', 0, 'Upper', 1);
opts.Display = 'Off';
opts.Lower = 0;
opts.StartPoint = 0.5;  % Any number in the interval (0,1) should work here

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
end