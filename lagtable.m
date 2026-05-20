function tbl = lagtable(data, maxlags)
    res1 = cell(maxlags, 1);
    for lag = 1:maxlags
        res1{lag} = nan(size(data));
        res1{lag}(1 + lag:end) = data(1:end - lag);
    end
    tbl = table(res1{:}, data); %do not change the order of arguments, the last variable should be the response variable!
end