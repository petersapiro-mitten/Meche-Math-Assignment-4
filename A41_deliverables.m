function A41_deliverables()

close all; clc;

% Orbit parameters
orbit_params = struct();
orbit_params.m_sun = 1;
orbit_params.m_planet = 1;
orbit_params.G = 40;

V0 = [8; 0; 0; 1.5];
gravity_rate = @(t,V) gravity_rate_func(t,V,orbit_params);

% Runge-Kutta Butcher tableaux
ForwardEuler.A = 0;
ForwardEuler.B = 1;
ForwardEuler.C = 0;

Midpoint.A = [0 0; 0.5 0];
Midpoint.B = [0 1];
Midpoint.C = [0 0.5];

RK4.A = [0 0 0 0;
         0.5 0 0 0;
         0 0.5 0 0;
         0 0 1 0];
RK4.B = [1/6 1/3 1/3 1/6];
RK4.C = [0 0.5 0.5 1];

%% 1. Example integration: numerical solutions vs exact solution

tspan = [0 30];
h = 0.01;

[t_FE,X_FE,h_FE,evals_FE] = ...
    explicit_RK_fixed_step_integration( ...
    gravity_rate,tspan,V0,h,ForwardEuler);

[t_MP,X_MP,h_MP,evals_MP] = ...
    explicit_RK_fixed_step_integration( ...
    gravity_rate,tspan,V0,h,Midpoint);

[t_RK4,X_RK4,h_RK4,evals_RK4] = ...
    explicit_RK_fixed_step_integration( ...
    gravity_rate,tspan,V0,h,RK4);

X_true = compute_planetary_motion(t_FE,V0,orbit_params);

% Plot every 50th numerical marker to reduce clutter.
% All numerical integration points are still calculated.
marker_idx = 1:50:length(t_FE);

figure('Color','w','Position',[100 100 1050 750]);

% x(t)
subplot(2,1,1);

plot(t_FE,X_true(:,1),'k-','LineWidth',0.8);
hold on;

plot(t_FE(marker_idx),X_FE(marker_idx,1),'o', ...
    'LineStyle','none','MarkerSize',4);

plot(t_MP(marker_idx),X_MP(marker_idx,1),'s', ...
    'LineStyle','none','MarkerSize',4);

plot(t_RK4(marker_idx),X_RK4(marker_idx,1),'^', ...
    'LineStyle','none','MarkerSize',4);

xlabel('$t$','Interpreter','latex','FontSize',15);
ylabel('$x(t)$','Interpreter','latex','FontSize',15);
title('$x(t)$: Exact and Numerical Solutions', ...
    'Interpreter','latex','FontSize',16);

legend({'Exact solution','Forward Euler', ...
    'Explicit Midpoint','Classical RK4'}, ...
    'Interpreter','latex','FontSize',10, ...
    'Location','eastoutside');

set(gca,'FontSize',12,'TickLabelInterpreter','latex');
grid on;
box on;
xlim(tspan);

% y(t)
subplot(2,1,2);

plot(t_FE,X_true(:,2),'k-','LineWidth',0.8);
hold on;

plot(t_FE(marker_idx),X_FE(marker_idx,2),'o', ...
    'LineStyle','none','MarkerSize',4);

plot(t_MP(marker_idx),X_MP(marker_idx,2),'s', ...
    'LineStyle','none','MarkerSize',4);

plot(t_RK4(marker_idx),X_RK4(marker_idx,2),'^', ...
    'LineStyle','none','MarkerSize',4);

xlabel('$t$','Interpreter','latex','FontSize',15);
ylabel('$y(t)$','Interpreter','latex','FontSize',15);
title('$y(t)$: Exact and Numerical Solutions', ...
    'Interpreter','latex','FontSize',16);

legend({'Exact solution','Forward Euler', ...
    'Explicit Midpoint','Classical RK4'}, ...
    'Interpreter','latex','FontSize',10, ...
    'Location','eastoutside');

set(gca,'FontSize',12,'TickLabelInterpreter','latex');
grid on;
box on;
xlim(tspan);

sgtitle('Runge-Kutta Example Integration', ...
    'Interpreter','latex','FontSize',18);


%% 2. Global truncation error vs step size

h_values = [0.5 0.25 0.125 0.0625 0.03125 0.015625];

error_FE = zeros(size(h_values));
error_MP = zeros(size(h_values));
error_RK4 = zeros(size(h_values));

for i = 1:length(h_values)

    h_test = h_values(i);

    [t_FE_test,X_FE_test,~,~] = ...
        explicit_RK_fixed_step_integration( ...
        gravity_rate,tspan,V0,h_test,ForwardEuler);

    [t_MP_test,X_MP_test,~,~] = ...
        explicit_RK_fixed_step_integration( ...
        gravity_rate,tspan,V0,h_test,Midpoint);

    [t_RK4_test,X_RK4_test,~,~] = ...
        explicit_RK_fixed_step_integration( ...
        gravity_rate,tspan,V0,h_test,RK4);

    % Exact solution evaluated at each method's time points
    true_FE = compute_planetary_motion( ...
        t_FE_test,V0,orbit_params);

    true_MP = compute_planetary_motion( ...
        t_MP_test,V0,orbit_params);

    true_RK4 = compute_planetary_motion( ...
        t_RK4_test,V0,orbit_params);

    % Euclidean state error at each time; use maximum over time
    err_FE_each = vecnorm(X_FE_test-true_FE,2,2);
    err_MP_each = vecnorm(X_MP_test-true_MP,2,2);
    err_RK4_each = vecnorm(X_RK4_test-true_RK4,2,2);

    error_FE(i) = max(err_FE_each);
    error_MP(i) = max(err_MP_each);
    error_RK4(i) = max(err_RK4_each);

end

% Linear fits in log-log space
fit_FE = polyfit(log(h_values),log(error_FE),1);
fit_MP = polyfit(log(h_values),log(error_MP),1);
fit_RK4 = polyfit(log(h_values),log(error_RK4),1);

% Measured convergence orders
p_FE = fit_FE(1);
p_MP = fit_MP(1);
p_RK4 = fit_RK4(1);

% Fitted error curves
fit_error_FE = exp(polyval(fit_FE,log(h_values)));
fit_error_MP = exp(polyval(fit_MP,log(h_values)));
fit_error_RK4 = exp(polyval(fit_RK4,log(h_values)));

% Plot global truncation error
figure('Color','w','Position',[150 150 900 650]);

% Data points: larger markers
loglog(h_values,error_FE,'o', ...
    'LineStyle','none','MarkerSize',9);
hold on;

loglog(h_values,error_MP,'s', ...
    'LineStyle','none','MarkerSize',9);

loglog(h_values,error_RK4,'^', ...
    'LineStyle','none','MarkerSize',9);

% Fitted convergence lines: thicker dashed lines
loglog(h_values,fit_error_FE,'--','LineWidth',1.8);
loglog(h_values,fit_error_MP,'--','LineWidth',1.8);
loglog(h_values,fit_error_RK4,'--','LineWidth',1.8);

xlabel('$h$ (step size)', ...
    'Interpreter','latex','FontSize',15);

ylabel('$E(h)$ (global truncation error)', ...
    'Interpreter','latex','FontSize',15);

title('Global Truncation Error vs. Step Size', ...
    'Interpreter','latex','FontSize',17);

legend({'Forward Euler data', ...
    'Explicit Midpoint data', ...
    'Classical RK4 data', ...
    sprintf('Forward Euler fit: p = %.2f',p_FE), ...
    sprintf('Midpoint fit: p = %.2f',p_MP), ...
    sprintf('RK4 fit: p = %.2f',p_RK4)}, ...
    'Interpreter','none','FontSize',10, ...
    'Location','southwest');

set(gca,'FontSize',12,'TickLabelInterpreter','latex');
grid on;
box on;


%% 3. Convergence order table

Method = {'Forward Euler'; ...
          'Explicit Midpoint'; ...
          'Classical RK4'};

Theoretical_Order = [1; 2; 4];
Measured_Order = [p_FE; p_MP; p_RK4];

Convergence_Table = table( ...
    Method,Theoretical_Order,Measured_Order);

disp('Global truncation error convergence results:');
disp(Convergence_Table);

fprintf('Forward Euler:     p = %.4f\n',p_FE);
fprintf('Explicit Midpoint: p = %.4f\n',p_MP);
fprintf('Classical RK4:     p = %.4f\n',p_RK4);

fprintf(['Rate function evaluations: FE=%d, ', ...
    'Midpoint=%d, RK4=%d\n'], ...
    evals_FE,evals_MP,evals_RK4);

end