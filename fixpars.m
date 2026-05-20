function res=fixpars(func,x,fitpars,mask,fixedvals)
pars=zeros(size(mask));
pars(mask)=fitpars;
pars(~mask)=fixedvals;
res=func(x,pars);