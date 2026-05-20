function [Coeff, bestorder, models] = armodel(data, maxorder, includeintercept)
    %find the optimal AR model with multiple linear regression
    %fit all orders from 1 to maxorder.
    %we use the AIC criterion to find the optimal parsimonious model (one
    %could also use the BIC criterion)
    %data can have NaN values or be a cell with different samples
    %
    if nargin < 2
        maxorder = 10;
    end
    if nargin < 3
        includeintercept = false;
    end
    %create table with lagged data (1:maxorder)
    %last variable should be the original data
    if iscell(data) && ~isempty(data)
        maxlength = 0;
        tbl = lagtable(data{1}(:), maxorder);
        for i = 2:numel(data)
            if length(data{i}) > maxlength
                maxlength = length(data{i});
            end
            tbl = [tbl; lagtable(data{i}(:), maxorder)]; %#ok<AGROW>
        end
    else
        maxlength = length(data);
        tbl = lagtable(data(:), maxorder);
    end
    if maxorder > maxlength
        error('armodel:toofewdata', 'Too few data for determining the maximum order');
    end
    %rowwise delete if any of the columns has nan values
    %we can also opt for pairwise delete of nan, but maybe it is better 
    %to base all models on the same number of data
    data
    tbl
    ndx1 = ~any(isnan(table2array(tbl)), 2);
    tbl = tbl(ndx1, :)
 
    ndx = true(size(tbl, 2), 1)
    aic = zeros(size(tbl, 2) - 1, 1);
    models = cell(maxorder, 1);
    for order = maxorder: -1:1
        %fit linear model store result in an object
        order
        tbl(:, ndx)
        models{order} = fitlm(tbl(:, ndx), 'intercept', includeintercept);
        aic(order) = models{order}.ModelCriterion.AIC;
        %deselect a column
        ndx(order) = false;
    end
    bestorder = find(min(aic) == aic); %minimum aic
    bestorder = bestorder(1); %take the first if there are draws
    %I add 1 to have similar output as armasel, maybe better to remove this
    Coeff = [1, models{bestorder}.Coefficients.Estimate'];
end

function tbl = lagtable(data, maxlags)
    res1 = cell(maxlags, 1);
    for lag = 1:maxlags
        res1{lag} = nan(size(data));
        res1{lag}(1 + lag:end) = data(1:end - lag);
    end
    tbl = table(res1{:}, data); %do not change the order of arguments, the last variable should be the response variable!
end
