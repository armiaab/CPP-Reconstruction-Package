function [dim, mundx] = getfundimension(mu, sigma)
    for dim = 1:200
        par = zeros(1, dim);
        try
            if ~isempty(mu),    mu(0, par);    end
            if ~isempty(sigma), sigma(0, par); end
            mundx = true(1, dim);
            return;
        catch
        end
    end
    error('shared:getfundimension', ...
        'Cannot determine parameter dimension from mu/sigma handles (tried up to 200).');
end
