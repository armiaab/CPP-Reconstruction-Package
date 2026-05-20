function sample_pairs=sample_stratified(data,N_strata,p)

if iscell(data)
    data1=[];
    for i=1:length(data)
        data1=[data1 data{i}(:).' nan];
    end
    data=data1;  
else
    data=data(:).';
end

N=length(data);
idx_not=find(isnan(data));
idx_not=union(idx_not,idx_not-1);
idx_ok=setdiff(1:N,idx_not);

% Stratified random sampling (we stratify time not data range)
data_ok=data(idx_ok);
N=length(data_ok);
N1=floor(N/N_strata);
m=floor(p*N1);

idx_X0=[];
for i=1:N_strata
    A=randperm(N1,m);
    idx_X0=[idx_X0 A+(i-1)*N1];
end

sample_pairs=zeros(1,2*length(idx_X0));
sample_pairs(1:2:end)=data(idx_X0);
sample_pairs(2:2:end)=data(idx_X0+1);
end