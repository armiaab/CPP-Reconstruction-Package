function out=ME_TimeScale(data,p)

if iscell(data)
    N=length(data);
    for i=1:N
        data{i}=data{i}(:);
    end
    [~,idx]=sort(cellfun(@length,data));
    data=data(idx);
    L=cellfun(@length,data);
    V=L.*(N:-1:1);
    m=find(V==max(V),1,'last');

    data1=cell(1,N-m+1);
    for i=m:N
        data1{i-m+1}=data{i}(1:L(m));
    end
    data=cell2mat(data1);
    [out,~]=armasel_s(data,1:p,1,1,[],1,size(data,2),size(data,2));
else
    data=data(:);
    [out,~,~,~]=armasel(data,1:p,1,1,[]);
end

end
