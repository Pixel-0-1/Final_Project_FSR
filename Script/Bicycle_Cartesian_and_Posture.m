close all
clc

% Apply default settings for all plots
set(0, 'DefaultTextInterpreter', 'latex')
set(0, 'DefaultLegendInterpreter', 'latex')
set(0, 'DefaultAxesTickLabelInterpreter', 'latex')
set(0, 'DefaultFigureColor', 'w')
lw = 2; % linewidth

x = out.x;
y = out.y;
theta = out.theta;
v = out.v;
psi = out.psi;
e_x = out.e_x;
e_y = out.e_y;
e_theta = out.e_theta;

% Plots
% v
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
subplot(2,1,1);
plot(v, 'k', 'LineWidth', lw);
title('v(t)');
legend('v(t)', 'Location', 'best', 'Orientation', 'vertical');
xlabel('t [s]');
ylabel('Heading Velocity [m/s]');
grid on;
box on;
set(gca, 'FontSize', 20);

% psi
subplot(2,1,2);
plot(psi, 'k', 'LineWidth', lw);
title('psi(t)');
legend('psi(t)', 'Location', 'best', 'Orientation', 'vertical');
xlabel('t [s]');
ylabel('Steering Angle [rad]');
grid on;
box on;
set(gca, 'FontSize', 20);

% Trajectory + Animation
figure('Renderer', 'painters', 'Position', [10 10 900 700]);
hold on;
axis equal;
grid on;
box on;
xlabel('x [m]');
ylabel('y [m]');
title('Actual Trajectory');
set(gca, 'FontSize', 20);

% Limits
padding = 0.1;
xlim([min(x.Data)-padding, max(x.Data)+padding]);
ylim([min(y.Data)-padding, max(y.Data)+padding]);

plot(x.Data, y.Data, 'k--', 'LineWidth', lw);  % Path 
plot(x.Data(1), y.Data(1), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');  % Initial point
plot(x.Data(end), y.Data(end), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');  % Final point

% Bicycle
robot_plot = plot(x.Data(1), y.Data(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
direction_arrow = quiver(x.Data(1), y.Data(1), cos(theta.Data(1))*0.2, sin(theta.Data(1))*0.2, ...
                         'b', 'LineWidth', lw, 'MaxHeadSize', 2);

path_plot = plot(x.Data(1), y.Data(1), 'k-', 'LineWidth', lw);

legend({'Path', 'Initial Point', 'Final Point', 'Bicycle'}, 'Location', 'best', 'Orientation', 'horizontal');

% Animation
for i = 1:10:length(x.Data)
    set(robot_plot, 'XData', x.Data(i), 'YData', y.Data(i));
    set(direction_arrow, 'XData', x.Data(i), 'YData', y.Data(i), ...
        'UData', cos(theta.Data(i))*0.2, 'VData', sin(theta.Data(i))*0.2);
    set(path_plot, 'XData', x.Data(1:i), 'YData', y.Data(1:i));
    
    pause(0.01); % Animation time (To see the orientation better, increase up to 0.25)
end

% Position x
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
subplot(3,1,1);
plot(x, 'k', 'LineWidth', lw);
legend('x(t)', 'Location', 'best', 'Orientation', 'vertical');
title('Position x');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 14);

% Position y
subplot(3,1,2);
plot(y, 'k', 'LineWidth', lw);
legend('y(t)', 'Location', 'best', 'Orientation', 'vertical');
title('Position y');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 14);

% Orientation
subplot(3,1,3);
plot(theta, 'k', 'LineWidth', lw);
legend('theta(t)', 'Location', 'best', 'Orientation', 'horizontal');
title('Orientation');
xlabel('t [s]');
ylabel('Orientation [rad]');
grid on;
box on;
set(gca, 'FontSize', 14);

% Position error e_x
figure('Renderer', 'painters', 'Position', [10 10 900 600]);
subplot(3,1,1);
plot(e_x, 'k', 'LineWidth', lw);
legend('e_x(t)', 'Location', 'best', 'Orientation', 'vertical');
title('Position error e_x');
xlabel('t [s]');
ylabel('Position [m]');
grid on;
box on;
set(gca, 'FontSize', 14);

% Position error e_y
subplot(3,1,2);
plot(e_y, 'k', 'LineWidth', lw);
legend('e_y(t)', 'Location', 'best', 'Orientation', 'vertical');
title('Position error y');
xlabel('t [s]');
ylabel('Position error [m]');
grid on;
box on;
set(gca, 'FontSize', 14);

% Orientation error e_theta
subplot(3,1,3);
plot(e_theta, 'k', 'LineWidth', lw);
legend('e_theta(t)', 'Location', 'best', 'Orientation', 'horizontal');
title('Orientation error');
xlabel('t [s]');
ylabel('Orientation [rad]');
grid on;
box on;
set(gca, 'FontSize', 14);