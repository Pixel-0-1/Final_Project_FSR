close all
clear
clc

% Apply default settings for all plots
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultFigureColor', 'w')
lw = 2; % linewidth

q_i = [0 0 0]'; % q_i=[x_i y_i theta_i]^T
random_vector = rand(1,3);
%q_f = 200*(random_vector/ norm(random_vector))'; % q_f=[x_f y_f theta_f]^T
q_f = [163.8064 113.4084 17.4931]'; % q_f used for the final project

% Display configurations
disp('Initial configuration:');
disp(q_i);
disp('Final configuration:');
disp(q_f);

x_i = q_i(1);
y_i = q_i(2);
theta_i = q_i(3);
x_f = q_f(1);
y_f = q_f(2);
theta_f = q_f(3);

% Novelty : Generate intermediate waypoints
num_waypoints = 3; % N° of waypoints
waypoints_x = zeros(num_waypoints+2, 1);
waypoints_y = zeros(num_waypoints+2, 1);
waypoints_theta = zeros(num_waypoints+2, 1);

waypoints_x(1) = x_i;
waypoints_y(1) = y_i;
waypoints_theta(1) = theta_i;
waypoints_x(end) = x_f;
waypoints_y(end) = y_f;
waypoints_theta(end) = theta_f;

% Generate a spline trajectory
for i = 2:num_waypoints+1

    progress = (i-1)/(num_waypoints+1);
    
    base_x = x_i + progress * (x_f - x_i);
    base_y = y_i + progress * (y_f - y_i);
    
    deviation_magnitude = 30; % Amplitude of deviation
    if mod(i, 2) == 0
        lateral_deviation = deviation_magnitude;
    else
        lateral_deviation = -deviation_magnitude;
    end
    
    % Normal vector wrt the main direction
    main_direction = [x_f - x_i, y_f - y_i];
    main_direction = main_direction / norm(main_direction);
    perpendicular = [-main_direction(2), main_direction(1)];
    
    waypoints_x(i) = base_x + lateral_deviation * perpendicular(1);
    waypoints_y(i) = base_y + lateral_deviation * perpendicular(2);
    
    % Compute the next orientation
    if i < num_waypoints+2
        next_x = waypoints_x(i+1);
        next_y = waypoints_y(i+1);
    else
        next_x = x_f;
        next_y = y_f;
    end
    waypoints_theta(i) = atan2(next_y - waypoints_y(i), next_x - waypoints_x(i));
end

% Spline
s_waypoints = linspace(0, 1, length(waypoints_x));
s_fine = 0:0.001:1; % Risoluzione più fine per la spline

x_spline = spline(s_waypoints, waypoints_x, s_fine);
y_spline = spline(s_waypoints, waypoints_y, s_fine);
theta_spline = spline(s_waypoints, waypoints_theta, s_fine);

% To avoid discontinuity
theta_spline = unwrap(theta_spline);

k = 9; %k_i=k_f=k>0

alpha_x = k*cos(theta_spline(end)) - 3*x_spline(end);
alpha_y = k*sin(theta_spline(end)) - 3*y_spline(end);
beta_x = k*cos(theta_spline(1)) + 3*x_spline(1);
beta_y = k*sin(theta_spline(1)) + 3*y_spline(1);

T_initial = 30;
step = 0.01;

fprintf('Initial Time = %.2f s\n\n',T_initial);

a_0 = 0; 
a_1 = 0; 
a_2 = 3/(T_initial^2); 
a_3 = -2/(T_initial^3);

t_initial = 0:step:T_initial;
s_initial = a_0 + a_1*t_initial + a_2*t_initial.^2 + a_3*t_initial.^3;

x_s_initial = interp1(s_fine, x_spline, s_initial, 'spline');
y_s_initial = interp1(s_fine, y_spline, s_initial, 'spline');

dt = step;
x_s_dot_initial = gradient(x_s_initial, dt) ./ gradient(s_initial, dt);
y_s_dot_initial = gradient(y_s_initial, dt) ./ gradient(s_initial, dt);

x_s_ddot_initial = gradient(x_s_dot_initial, dt) ./ gradient(s_initial, dt);
y_s_ddot_initial = gradient(y_s_dot_initial, dt) ./ gradient(s_initial, dt);

v_tilde_initial = sqrt(x_s_dot_initial.^2 + y_s_dot_initial.^2);
w_tilde_initial = (y_s_ddot_initial.*x_s_dot_initial - x_s_ddot_initial.*y_s_dot_initial)./(x_s_dot_initial.^2 + y_s_dot_initial.^2);

s_dot_initial = a_1 + 2*a_2*t_initial + 3*a_3*t_initial.^2;
initial_v = v_tilde_initial .* s_dot_initial;
initial_w = w_tilde_initial .* s_dot_initial;

v_max_initial = max(abs(initial_v));
w_max_initial = max(abs(initial_w));

fprintf('initial v = %.2f m/s\n',v_max_initial);
fprintf('initial w = %.2f rad/s\n',w_max_initial);

T = T_initial;
bounds_unsatisfied = true;
final_scaling_factor = 1;

while bounds_unsatisfied
    
    %Cubic polynomial s(t) = a_0 + a_1*T + a_2*T^2 + a_3*T^3
    a_0 = 0;
    a_1 = 0;
    a_2 = 3/(T^2);
    a_3 = -2/(T^3);

    t = 0:step:T;
    s = a_0 + a_1*t + a_2*t.^2 + a_3*t.^3;
    s_dot = a_1 + 2*a_2*t + 3*a_3*t.^2;

    x_s = interp1(s_fine, x_spline, s, 'spline');
    y_s = interp1(s_fine, y_spline, s, 'spline');
    
    % Derivatives
    x_s_dot = gradient(x_s, step) ./ gradient(s, step);
    y_s_dot = gradient(y_s, step) ./ gradient(s, step);
    
    x_s_ddot = gradient(x_s_dot, step) ./ gradient(s, step);
    y_s_ddot = gradient(y_s_dot, step) ./ gradient(s, step);
     
    theta = atan2(y_s_dot, x_s_dot);
    v_tilde = sqrt(x_s_dot.^2 + y_s_dot.^2);
    w_tilde = (y_s_ddot.*x_s_dot - x_s_ddot.*y_s_dot)./(x_s_dot.^2 + y_s_dot.^2);
    v = v_tilde .* s_dot;
    w = w_tilde .* s_dot;

    % Velocity Bounds
    v_max = max(abs(v));
    w_max = max(abs(w));

    if (v_max <= 10 && w_max <= 2)
        fprintf('\nVelocity bounds satisfied:\n|v(t)| <= 10 m/s\n|w(t)| <= 2 rad/s\n\n');
        if final_scaling_factor~=1.00
            fprintf('Scaling factor = %.2f\n',final_scaling_factor);
            fprintf('Time scaled = %.2f s\n',T);
            fprintf('final v = %.2f m/s\n',v_max);
            fprintf('final w = %.2f rad/s\n',w_max);
        end
        bounds_unsatisfied = false;
    else
        scaling_factor_v = v_max / 10;
        scaling_factor_w = w_max / 2;
        scaling_factor = max(scaling_factor_v, scaling_factor_w);

        T = T * scaling_factor;
        final_scaling_factor = final_scaling_factor * scaling_factor;
    end
end

% Time derivatives
x_dot = gradient(x_s, step);  % dx/dt
y_dot = gradient(y_s, step);  % dy/dt

% Computing desired velocities
v_d_l_nl = sqrt(x_dot.^2 + y_dot.^2);  % Desired Linear Velocity
w_d_l_nl = gradient(theta, step);      % Desired Angular Velocity (derivative of theta)

fprintf('\n=== Desired Velocities ===\n');
fprintf('Max v_d = %.2f m/s\n', max(abs(v_d_l_nl)));
fprintf('Max w_d = %.2f rad/s\n', max(abs(w_d_l_nl)));

L = 3.5; % Wheel Base

x_d = timeseries(x_s,t);
y_d = timeseries(y_s,t);
theta_d = timeseries(theta,t);
v_d= timeseries(v_d_l_nl, t);     
w_d = timeseries(w_d_l_nl, t);    

% Position 
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
subplot(2,2,1);
plot(t, x_s, 'r', t, y_s, 'b', 'LineWidth', lw);
legend('x_s(t)', 'y_s(t)', 'Location', 'best', 'Orientation', 'vertical');
title('Position');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 14);
xlim([t(1) t(end)]);

% Orientation
subplot(2,2,2);
plot(t, theta, 'k', 'LineWidth', lw);
legend('theta(t)', 'Location', 'best', 'Orientation', 'horizontal');
title('Orientation');
xlabel('t [s]');
ylabel('Orientation [rad]');
grid on;
box on;
set(gca, 'FontSize', 14);
xlim([t(1) t(end)]);

% Heading Velocity
subplot(2,2,3);
if final_scaling_factor~=1.00
    initial_v_interp = interp1(t_initial, initial_v, t, 'linear', 'extrap');
    plot(t, v, 'r', t, initial_v_interp, 'b', 'LineWidth', lw);
else
    plot(t, v, 'k', 'LineWidth', lw);
end

yline(10, '--k', 'LineWidth', 1.5);
yline(-10, '--k', 'LineWidth', 1.5);
title('Heading Velocity v(t)');
xlabel('t [s]');
ylabel('v [m/s]');

if final_scaling_factor~=1.00
    legend('scaled v(t)','initial v(t)','+v_{max}','-v_{max}','Location','best', 'Orientation', 'vertical');
else
    legend('v(t)','+v_{max}','-v_{max}','Location','best', 'Orientation', 'vertical');
end

grid on;
box on;
set(gca, 'FontSize', 14);
xlim([t(1) t(end)]);

% Angular Velocity
subplot(2,2,4);
if final_scaling_factor~=1.00
    initial_w_interp = interp1(t_initial, initial_w, t, 'linear', 'extrap');
    plot(t, w, 'r', t, initial_w_interp, 'b', 'LineWidth', lw);
else
    plot(t, w, 'k', 'LineWidth', lw);
end

title('Angular Velocity w(t)');
xlabel('t [s]');
ylabel('w [rad/s]');
yline(2.0, '--k', 'LineWidth', 1.5);
yline(-2.0, '--k', 'LineWidth', 1.5);

if final_scaling_factor~=1.00
    legend('scaled w(t)','initial w(t)','+w_{max}','-w_{max}','Location','best', 'Orientation', 'vertical');
else
    legend('w(t)','+w_{max}','-w_{max}','Location','best', 'Orientation', 'horizontal');
end

grid on;
box on;
set(gca, 'FontSize', 14);
xlim([t(1) t(end)]);

% Time Laws 
figure('Renderer', 'painters', 'Position', [10 10 900 500]);
subplot(2,1,1);
plot(t, s, 'k', 'LineWidth', lw);
legend('s(t)', 'Location', 'best', 'Orientation', 'horizontal');
xlabel('t [s]');
ylabel('Value [m]');
title('Curvilinear Abscissa');
grid on;
box on;
set(gca, 'FontSize', 20);
xlim([t(1) t(end)]);

subplot(2,1,2);
plot(t, s_dot, 'k', 'LineWidth', lw);
legend('ds/dt(t)', 'Location', 'best', 'Orientation', 'horizontal');
xlabel('t [s]');
ylabel('Value [m/s]');
title('Time Derivative');
grid on;
box on;
set(gca, 'FontSize', 20);
xlim([t(1) t(end)]);

% Trajectory + Animation
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
hold on;
axis equal;
grid on;
box on;
xlabel('x [m]');
ylabel('y [m]');
title('Trajectory with Waypoints');
set(gca, 'FontSize', 20);

% Limits
padding = max(abs([x_s, y_s])) * 0.1;
xlim([min(x_s)-padding, max(x_s)+padding]);
ylim([min(y_s)-padding, max(y_s)+padding]);

% Plot waypoints
plot(waypoints_x, waypoints_y, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');  % Waypoints
plot(x_s, y_s, 'k--', 'LineWidth', lw);  % Path
plot(x_i, y_i, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Initial point
plot(x_f, y_f, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');  % Final point

% Bicycle
robot_plot = plot(x_s(1), y_s(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
direction_arrow = quiver(x_s(1), y_s(1), cos(theta(1))*10, sin(theta(1))*10, ...
                         'b', 'LineWidth', lw, 'MaxHeadSize', 2);

path_plot = plot(x_s(1), y_s(1), 'k-', 'LineWidth', lw);

legend({'Waypoints', 'Planned Path', 'Initial Point', 'Final Point', 'Bicycle', 'Direction'}, ...
       'Location', 'best', 'Orientation', 'vertical');

% Animation
for i = 1:10:length(x_s)
    set(robot_plot, 'XData', x_s(i), 'YData', y_s(i));
    set(direction_arrow, 'XData', x_s(i), 'YData', y_s(i), ...
        'UData', cos(theta(i))*10, 'VData', sin(theta(i))*10);
    set(path_plot, 'XData', x_s(1:i), 'YData', y_s(1:i));
    
    pause(0.01);
end

%% Plot from simulink

error_x = out.error_x;
error_y = out.error_y;
error_theta = out.error_theta;
error_y1 = out.error_y1;
error_y2 = out.error_y2;
v = out.v;
x_g = out.x_g;
y_g = out.y_g;
y1_d = out.y1_d;
y2_d = out.y2_d;
psi = out.psi;

% e_x
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
subplot(3,1,1);
plot(error_x, 'k', 'LineWidth', lw);
title('e_x(t)');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 20);

% e_y
subplot(3,1,2);
plot(error_y, 'k', 'LineWidth', lw);
title('e_y(t)');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 20);

% e_theta
subplot(3,1,3);
plot(error_theta, 'k', 'LineWidth', lw);
title('e_theta(t)');
xlabel('t [s]');
ylabel('Angle [rad]');
grid on;
box on;
set(gca, 'FontSize', 18);

% e_y1
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
subplot(1,2,1);
plot(error_y1, 'k', 'LineWidth', lw);
title('e_x_CoG(t)');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 20);

% e_y2
subplot(1,2,2);
plot(error_y2, 'k', 'LineWidth', lw);
title('e_y_CoG(t)');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 20);

% v
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
subplot(2,1,1);
plot(v, 'k', 'LineWidth', lw);
title('v(t)');
xlabel('t [s]');
ylabel('Heading Velocity [m/s]');
grid on;
box on;
set(gca, 'FontSize', 20);

% psi
subplot(2,1,2);
plot(psi, 'k', 'LineWidth', lw);
title('psi');
xlabel('t [s]');
ylabel('Steering Angle [rad]');
grid on;
box on;
set(gca, 'FontSize', 20);

% Trajectory
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
plot(y1_d.Data, y2_d.Data, '--r', 'LineWidth', lw+1, 'DisplayName', 'Desired Trajectory');
hold on;
plot(x_g.Data, y_g.Data, '-k', 'LineWidth', lw, 'DisplayName', 'Effective Trajectory');
% Marker for initial and final positions
plot(y1_d.Data(1), y2_d.Data(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(y1_d.Data(end), y2_d.Data(end), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'End');

xlabel('x [m]');
ylabel('y [m]');
title('Trajectory in xy-Plane');
legend('Location', 'best');
grid on;
box on;
axis equal;
set(gca, 'FontSize', 20);
