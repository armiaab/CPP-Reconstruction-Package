function LP = LegitimateRegion(cost, dim, lb, ub, n_agents, index_sigma)
% Sample random points in [lb,ub], evaluate cost, return rows with finite cost.
% LP(:, 1:dim) = parameter vectors, LP(:, dim+1) = cost values.
% index_sigma: logical mask — sigma parameters must stay positive.
    lb = lb(:)';  ub = ub(:)';
    n_try   = max(500, n_agents * 50);
    LP_rows = zeros(0, dim + 1);

    for attempt = 1:n_try
        par = lb + rand(1, dim) .* (ub - lb);
        % Enforce sigma positivity at sampled point
        if any(index_sigma) && any(par(index_sigma) <= 0)
            continue;
        end
        try
            f = cost(par);
            if isfinite(f) && isreal(f)
                LP_rows(end+1, :) = [par, f]; %#ok<AGROW>
                if size(LP_rows,1) >= n_agents * 10, break; end
            end
        catch
        end
    end

    if isempty(LP_rows)
        % Fallback: use midpoints with small perturbations
        mid = (lb + ub) / 2;
        mid(~isfinite(mid)) = 0;
        LP_rows = zeros(n_agents, dim + 1);
        for k = 1:n_agents
            par = mid + 0.01*(ub-lb).*randn(1,dim);
            par = max(lb, min(ub, par));
            try
                f = cost(par);
                if ~isfinite(f) || ~isreal(f), f = realmax; end
            catch
                f = realmax;
            end
            LP_rows(k,:) = [par, f];
        end
    end

    % Return at most n_agents rows, sorted by cost
    LP_rows = sortrows(LP_rows, dim + 1);
    n_keep  = min(n_agents, size(LP_rows,1));
    LP = LP_rows(1:n_keep, :);
end
