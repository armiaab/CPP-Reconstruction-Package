function EZ=Coefficients_Spline(varargin)

J=varargin{1};K=varargin{2};ModelType=varargin{3};SplineType=varargin{4};

%Classical Hermite polynomials
syms x dt
H=sym('H',[J 1]);
for j=1:J
    H(j)=horner(simplify(exp(x^2/2)*diff(exp(-x^2/2),j)),x);
end

%Coefficients of H
H_coeff=cell(1,J);
for j=1:J
    H_coeff{j}=nonzeros(sym2poly(H(j)));
end
if ~strcmp(ModelType,'Additive noise1')
    H=H(3:end);   % We do not need the first two elements of H
    H=matlabFunction(H);
end

%The binomial-like coefficient matrix Coeff
syms w
Coeff=cell(1,J);
Coeff{1}=[];Coeff{2}=[];
for j=3:J
    A=zeros(1,j+1);
    for k=1:j+1
        A(k)=(-1)^(k-1)*nchoosek(j,k-1);
    end
    C=A.* w.^(0:j);
    Coeff{j}=MyMatlabFunction(C,w,w);
end
switch ModelType
    case 'Additive noise1'
        if J<=5 && K<=30
            switch SplineType
                case 'L'
                    filename=sprintf('ETA_L_%d_%d',J,K);ETA=str2func(filename);EZ=[];
                case 'Q'
                    filename=sprintf('ETA_Q_%d_%d',J,K);ETA=str2func(filename);EZ=[];
                case {'C','P','SCS'}
                    filename=sprintf('ETA_C_%d_%d',J,K);ETA=str2func(filename);EZ=[];
            end
            H=matlabFunction(H);
        else
            syms y y0 positive
            syms muY(y)
            Hz=subs(H,x,(y-y0)/dt^(1/2));
            ETA=sym('ETA',[J 1]);
            for j=1:J
                A0=Hz(j);S=A0;
                for k=1:K
                    A0=diff(A0,y);
                    A0=simplifyFraction(muY(y)*A0+1/2*diff(A0,y));
                    switch SplineType
                        case {'L','Approximate'}
                            A0=subs(A0,{diff(muY(y),2),diff(muY(y),3)},{0,0});   %Since we consider linear spline all derivatives of muY of order>=2 are 0
                        case 'Q'
                            A0=subs(A0,{diff(muY(y),3),diff(muY(y),4)},{0,0});   %Since we consider quadratic spline all derivatives of muY of order>=3 are 0
                        case {'C','P','SCS'}
                            A0=subs(A0,{diff(muY(y),4),diff(muY(y),5)},{0,0});   %Since we consider cubic spline all derivatives of muY of order>=4 are 0
                    end
                    S=S+A0*dt^k/factorial(k);
                end
                ETA(j)=1/factorial(j)*limit(S,y,y0);
            end
            ETA=simplifyFraction(ETA);
            H=matlabFunction(H);
            syms D0
            switch SplineType
                case {'L','Approximate'}
                    D=sym('D',[1 1]);ETA=subs(ETA,{diff(muY(y0)),muY(y0)},{D(1),D0});filename=strcat(['ETA_L_' num2str(J) '_' num2str(K)]);
                case 'Q'
                    D=sym('D',[1 2]);ETA=subs(ETA,{diff(muY(y0),2),diff(muY(y0),1),muY(y0)},{D(2),D(1),D0});filename=strcat(['ETA_Q_' num2str(J) '_' num2str(K)]);
                case {'C','P','SCS'}
                    D=sym('D',[1 3]);ETA=subs(ETA,{diff(muY(y0),3),diff(muY(y0),2),diff(muY(y0),1),muY(y0)},{D(3),D(2),D(1),D0});filename=strcat(['ETA_C_' num2str(J) '_' num2str(K)]);
            end
            [n,d]=numden(ETA);ETA=horner(n,D0)./horner(d,D0);
            I=hasSymType(ETA,'variable');ETA(~I)=ETA(~I)+10^(-100)*D0;   %this is just to turn constant component of ETA into a vector by adding an infinitesimal vector
            D=[dt D0 D];
            ETA=matlabFunction(ETA,'File',filename,'Vars',D,'Outputs',{'ETA'});  %The reason we SAVE ETA as an m-file is that this way matlab optimizes ETA which is very important to reduce computational burden
            EZ=[];
        end
    case 'Additive noise2'
        if J<3,disp('J should be equal or bigger than 3');return;end
        if J<=5 && K<=30
            switch SplineType
                case {'L','Approximate'}
                    filename=sprintf('EZ_L_%d_%d',J,K);EZ=str2func(filename);ETA=[];
                case 'Q'
                    filename=sprintf('EZ_Q_%d_%d',J,K);EZ=str2func(filename);ETA=[];
                case {'C','P','SCS'}
                    filename=sprintf('EZ_C_%d_%d',J,K);EZ=str2func(filename);ETA=[];
            end
        else
            %Calculating expectations of Z using the infinitesimal generator of the model
            EZ=sym('EZ',[J 1]);
            syms u m(u) s0 x x0
            for j=1:J
                B0=((x-x0)/s0)^j;
                S=B0;
                for k=1:K
                    B0=diff(B0,x);
                    B0=simplifyFraction(m(x)*B0+1/2*s0^2*diff(B0,x));
                    switch SplineType
                        case {'L','Approximate'}
                            B0=subs(B0,{diff(m(x),2),diff(m(x),3)},{0,0});   %Since we consider linear splines all derivatives of m(x) of order>=2 are 0
                        case 'Q'
                            B0=subs(B0,{diff(m(x),3),diff(m(x),4)},{0,0});   %Since we consider quadratic splines all derivatives of m(x) of order>=3 are 0
                        case {'C','P','SCS'}
                            B0=subs(B0,{diff(m(x),4),diff(m(x),5)},{0,0});   %Since we consider cubic splines all derivatives of m(x) of order>=4 are 0
                    end
                    S=S+B0*dt^k/factorial(k);
                end
                EZ(j)=simplifyFraction(dt^(-j/2)*limit(S,x,x0));
            end
            syms Dm0
            switch SplineType
                case {'L','Approximate'}
                    Dm=sym('Dm',[1 1]);EZ=subs(EZ,{diff(m(x0)),m(x0)},{Dm(1),Dm0});
                    filename=strcat(['EZ_L_' num2str(J) '_' num2str(K)]);
                case 'Q'
                    Dm=sym('Dm',[1 2]);EZ=subs(EZ,{diff(m(x0),2),diff(m(x0),1),m(x0)},{Dm(2),Dm(1),Dm0});
                    filename=strcat(['EZ_Q_' num2str(J) '_' num2str(K)]);
                case {'C','P','SCS'}
                    Dm=sym('Dm',[1 3]);EZ=subs(EZ,{diff(m(x0),3),diff(m(x0),2),diff(m(x0),1),m(x0)},{Dm(3),Dm(2),Dm(1),Dm0});
                    filename=strcat(['EZ_C_' num2str(J) '_' num2str(K)]);
            end
            %Unfortunately, here we canot benifit from the horner scheme when J and K are big (MATLAB does not throw any error message but later the generated m.file might throw
            %the error message: "Error: Nesting of {, [, and ( cannot exceed a depth of 32"
            [n,d]=numden(EZ);EZ=horner(n,s0)./horner(d,s0);   %%based on previous comment you should avoid this line in case of a possible error message
            I=hasSymType(EZ, 'variable');EZ(~I)=EZ(~I)+10^(-100)*s0;
            Dm=[Dm0 Dm];
            EZ=matlabFunction(EZ,'File',filename,'Vars',[dt Dm s0],'Outputs',{'EZ'});
            ETA=[];
        end
    case 'Multiplicative noise'
        if J<3,disp('J should be equal or bigger than 3');return;end
        I=0;
        switch SplineType
            case {'LL','Approximate'}
                if J<=4 && K<=20,I=1;end
            case 'QQ'
                if (J==3 && K<=15) || (J==4 && K<=14),I=1;end
            case {'CC','PP','CP','SCS'}
                if J<=4 && K<=11,I=1;end
        end
        I=0;
        if I
            switch SplineType
                case {'LL','Approximate'}
                    filename=sprintf('EZ_LL_%d_%d',J,K);EZ=str2func(filename);ETA=[];
                case 'QQ'
                    filename=sprintf('EZ_QQ_%d_%d',J,K);EZ=str2func(filename);ETA=[];
                case {'CC','PP','CP','SCS'}
                    filename=sprintf('EZ_CC_%d_%d',J,K);EZ=str2func(filename);ETA=[];
            end
        else
            %Calculating expectations of Z using the infinitesimal generator of the model
            EZ=sym('EZ',[J 1]);
            syms u m(u) s(u) x x0
            B=int(1/s,u,x0,x);B^3
            for j=1:J
                B0=B^j;
                S=B0;
                for k=1:K
                    B0=diff(B0,x);
                    B0=simplifyFraction(m(x)*B0+1/2*s(x)^2*diff(B0,x));
                    switch SplineType
                        case {'LL','Approximate'}
                            B0=subs(B0,{diff(m(x),2),diff(m(x),3),diff(s(x),2),diff(s(x),3)},{0,0,0,0});   %Since we consider linear splines all derivatives of m(x) and s(x) of order>=2 are 0
                        case 'QQ'
                            B0=subs(B0,{diff(m(x),3),diff(m(x),4),diff(s(x),3),diff(s(x),4)},{0,0,0,0});   %Since we consider quadratic splines all derivatives of m(x) and s(x) of order>=3 are 0
                        case {'CC','PP','CP','SCS'}
                            B0=subs(B0,{diff(m(x),4),diff(m(x),5),diff(s(x),4),diff(s(x),5)},{0,0,0,0});   %Since we consider cubic splines all derivatives of m(x) and s(x) of order>=4 are 0
                    end
                    S=S+B0*dt^k/factorial(k);
                end
                simplifyFraction(dt^(-j/2)*S)
                EZ(j)=simplifyFraction(dt^(-j/2)*limit(S,x,x0));
            end
            syms Dm0 Ds0
            switch SplineType
                case {'LL','Approximate'}
                    Dm=sym('Dm',[1 1]);Ds=sym('Ds', [1 1]);EZ=subs(EZ,{diff(m(x0)),m(x0),diff(s(x0)),s(x0)},{Dm(1),Dm0,Ds(1),Ds0});filename=strcat(['EZ_LL_' num2str(J) '_' num2str(K)]);
                case 'QQ'
                    Dm=sym('Dm',[1 2]);Ds=sym('Ds', [1 2]);EZ=subs(EZ,{diff(m(x0),2),diff(m(x0),1),m(x0),diff(s(x0),2),diff(s(x0),1),s(x0)},{Dm(2),Dm(1),Dm0,Ds(2),Ds(1),Ds0});
                        filename=strcat(['EZ_QQ_' num2str(J) '_' num2str(K)]);
                case {'CC','PP'}
                    Dm=sym('Dm',[1 3]);Ds=sym('Ds', [1 3]);EZ=subs(EZ,{diff(m(x0),3),diff(m(x0),2),diff(m(x0),1),m(x0),diff(s(x0),3),diff(s(x0),2),diff(s(x0),1),s(x0)},{Dm(3),Dm(2),Dm(1),Dm0,Ds(3),Ds(2),Ds(1),Ds0});
                    filename=strcat(['EZ_CC_' num2str(J) '_' num2str(K)]);
            end
            EZ
            %Unfortunately, here we canot benifit from the horner scheme when J and K are big (MATLAB does not throw any error message but later the generated m.file might throw
            %the error message: "Error: Nesting of {, [, and ( cannot exceed a depth of 32"

            % [n,d]=numden(EZ);EZ=horner(n,Ds0)./horner(d,Ds0);   %based on previous comment you should avoid this line in case of a possible error message
            % I=hasSymType(EZ, 'variable');EZ(~I)=EZ(~I)+10^(-100)*Ds0;   % this is just to turn constant component of EZ into a vector by adding an infinitesimal vector
            % Dm=[Dm0 Dm];Ds=[Ds0 Ds];
            % %The reason we SAVE EZ as an m-file in bellow is that this way matlab optimizes EZ which is very important to reduce computational burden
            % EZ=matlabFunction(EZ,'File',filename,'Vars',[dt Dm Ds],'Outputs',{'EZ'});
            % ETA=[];
        end
end
end
