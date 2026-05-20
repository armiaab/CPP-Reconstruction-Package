function data = prepare_replicateData(data)

% Handling replicate data
if iscell(data)  % replicate data
    data1=[];
    for i=1:length(data)
        data1=[data1 data{i}(:).' nan];
    end
    data=data1;
    X0=data(1:end-1);X=data(2:end);
else % typical data
    data=data(:).';
    X0=data(1:end-1);X=data(2:end);
end
end