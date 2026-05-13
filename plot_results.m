function U_eff = plot_results(res, x, y, y2)

U_eff = [];
%   struct with fields:
%
%        estimated_par: [1.6676 -1.8720 0.6711 0.6360]
%            negloglik: 9.9498e+03
%            datarange: [-2.6289 2.9383]
%                   dt: 1
%                mufun: @(x,par)interp1(knots_mu,par(ndx_mu),x)
%             sigmafun: @(x,par)interp1(knots_sigma,par(ndx_sigma),x)
%     computation_time: 56.5869
%                knots: {[-2.6289 2.9383]  [-2.6289 2.9383]}
%              options: [1×1 struct]
%
% result2.options
%
% ans =
%
%   struct with fields:
%
%                lb: [-10 -10 1.0000e-10 1.0000e-10]
%                ub: [10 10 10 10]
%            nknots: [2 2]
%             knots: []
%            spline: {'L'  'L'}
%            solver: 'fmincon'
%       useparallel: 0
%     search_agents: 5
%           maxiter: 1.7977e+308

if ~(nargin ==3 && (ischar(y) || issring(y)))

    if nargin>2
        if length(x)~=length(y)
            y=y+0.*x;
        end
        if length(x)~=length(y2)
            y2=y2+0.*x;
        end
    end

    if nargin < 2
        ran=res.datarange;
        if ~isempty(res.options.l)&&~isnan(res.options.l)
            ran(1)=res.options.l;
        end
        if ~isempty(res.options.r)&&~isnan(res.options.r)
            ran(2)=res.options.r;
        end
        x = linspace(ran(1), ran(2), 1000);
    end
    if nargin < 3
        y = [];
    end
    if nargin < 4
        y2 = [];
    end

    figure
    subplot(2, 1, 1)
    if ~isempty(y)
        plot(x, y, '-k', 'DisplayName', 'True model');
    end
    hold on;
    h = yline(0, '-k');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
    if ~isempty(res.knots)
        if numel(res.knots{1}) == 1
            plot(x, res.mufun(x, res.estimated_par), '-.b', 'DisplayName', 'constant');
        else
            plot(x, res.mufun(x, res.estimated_par), '-.b', 'DisplayName', sprintf('''%s'' spline', res.options.spline{1}));
            plot(res.knots{1}, res.mufun(res.knots{1}, res.estimated_par), '*r', 'DisplayName', 'knots'); %Estimated drift at knots
        end
    else
        plot(x, res.mufun(x, res.estimated_par), '-.b', 'DisplayName', 'mu');
    end
    if isfield(res.options, 'j')
        title(['Hermite Reconstruction - ' res.options.solver])
    else
        title(['Euler Reconstruction - ' res.options.solver])
    end
    xlim([-inf inf]);
    ylim([-inf inf]);
    %  xlabel('State x');
    ylabel('Drift \mu (x)');


    legend;
    subplot(2, 1, 2)
    if ~isempty(y2)
        plot(x, y2, '-k', 'DisplayName', 'True model');
    end
    hold on;
    h = yline(0, '-k');
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
    if ~isempty(res.knots)
        if numel(res.knots{2}) == 1
            plot(x, res.sigmafun(x, res.estimated_par), '-.b', 'DisplayName', 'constant');
        else
            plot(x, res.sigmafun(x, res.estimated_par), '-.b', 'DisplayName', sprintf('''%s'' spline', res.options.spline{2}));
            plot(res.knots{2}, res.sigmafun(res.knots{2}, res.estimated_par), '*r', 'DisplayName', 'knots'); %Estimated drift at knots
        end
    else
        plot(x, res.sigmafun(x, res.estimated_par), '-.b', 'DisplayName', 'sigma');
    end
    xlim([-inf inf]);
    ylim([-inf inf]);
    xlabel('State x');
    ylabel('Diffusion \sigma (x)');
    legend

else
    figure
    subplot(3, 1, 1)
    h = yline(0, '-k');hold on;
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
    if ~isempty(res.knots)
        if numel(res.knots{1}) == 1
            plot(x, res.mufun(x, res.estimated_par), '-.b', 'DisplayName', 'constant');
        else
            plot(x, res.mufun(x, res.estimated_par), '-.b', 'DisplayName', sprintf('''%s'' spline', res.options.spline{1}));hold on;
            plot(res.knots{1}, res.mufun(res.knots{1}, res.estimated_par), '*r', 'DisplayName', 'knots'); %Estimated drift at knots
        end
    else
        plot(x, res.mufun(x, res.estimated_par), '-.b', 'DisplayName', 'mu');
    end
    if isfield(res.options, 'j')
        title(['Hermite Reconstruction - ' res.options.solver])
    else
        title(['Euler Reconstruction - ' res.options.solver])
    end
    xlim([-inf inf]);
    ylim([-inf inf]);
    %  xlabel('State x');
    ylabel('Drift \mu (x)');


    legend;
    subplot(3, 1, 2)
    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
    if ~isempty(res.knots)
        if numel(res.knots{2}) == 1
            plot(x, res.sigmafun(x, res.estimated_par), '-.b', 'DisplayName', 'constant');
        else
            plot(x, res.sigmafun(x, res.estimated_par), '-.b', 'DisplayName', sprintf('''%s'' spline', res.options.spline{2}));hold on;
            plot(res.knots{2}, res.sigmafun(res.knots{2}, res.estimated_par), '*r', 'DisplayName', 'knots'); %Estimated drift at knots
        end
    else
        plot(x, res.sigmafun(x, res.estimated_par), '-.b', 'DisplayName', 'sigma');
    end
    xlim([-inf inf]);
    ylim([0 inf]);
    %     xlabel('State x');
    ylabel('Diffusion \sigma (x)');
    legend

    subplot(3,1,3)
    A1 = res.mufun(x, res.estimated_par);
    A2 = res.sigmafun(x, res.estimated_par);
    U = -2 .* cumtrapz(x, A1 ./ A2.^2) + 2 .* log(A2);
    U_eff = U - min(U) + 0.5;

    set(get(get(h, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off')
    if ~isempty(res.knots)
        h1=plot(x, U_eff, '-.b', 'DisplayName', sprintf('''%s'' spline', res.options.spline{2}));hold on;
        if isequal(res.knots{1},res.knots{2})
            U1 = @(z)interp1(x,U_eff,z);
            U_knots = U1(res.knots{1});hold on;
            h2=plot(res.knots{2}, U_knots, '*r', 'DisplayName', 'knots');hold on;
            xlim([-inf inf]);
            xlabel('State x');
            ylabel({'Effective potential'; 'U_{eff} (x)'});
        end
    end
end
end
