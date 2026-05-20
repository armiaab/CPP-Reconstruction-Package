function result = legitimate_points(data, dt, varargin)
%adapt the upper bounds and lower bounds of the search to avoid failure of
%fmincon
%
%usage: 
%   result=legitimate_points(data,dt, 'parameter', value ...)
%   for parameter-value pairs see hermite_reconstruction
%
    defaults=hermite_reconstruction('-d');   
    %do initial checks and collect the used options
    [options] = shared.parse_options_reconst(data, varargin, defaults.parametric, defaults.spline);
    prev_results=struct('options',options,'estimated_par',[]);
    %run hermite_reconstruction with solver 'legpoints' to adapt the ub and
    %lb to avoid failure of the fmincon solver
    prev_results.options.solver='legpoints';
    result=hermite_reconstruction(data,dt,'prev',prev_results);
    result.options.solver=options.solver;
end
