function F=MyMatlabFunction(f,Vars,var)
Vars=char(Vars);
if Vars(1)=='['
    Vars=Vars(2:end-1);
end
syms HELLO1
Idx=has(f,var);
if length(Idx)>1
    f(Idx==0)=f(Idx==0)+HELLO1;
elseif Idx==0
        f=f+HELLO1;
end
f=char(f);
f=vectorize(f);
f=replace(f,'HELLO1',strcat(['0.*' char(var)]));
f=append(strcat(['@(' Vars ')']),f);
F=str2func(f);
end