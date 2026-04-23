clc; close all; clear;

g0 = readmatrix('Experiment1/alpha_data.csv');% [t , y]
g1 = readmatrix('Experiment1/u_data.csv');% [t , u]

t = g0(:,1);% time vector
y = deg2rad(g0(:,2));% output (convert deg → rad)
u = g1(:,2);% input (in degrees)

h = mean(diff(t));% sample time

Yraw   = y(2:end);
Phiraw = [ y(1:end-1) , u(1:end-1) ];
theta_raw = Phiraw \ Yraw;

a_raw = theta_raw(1);
b_raw = theta_raw(2);

Gz_raw = tf([0 b_raw], [1 -a_raw], h);
Gs_raw = d2c(Gz_raw, 'zoh');

y_hat_raw = lsim(Gz_raw, u, t);
fit_raw = 100 * (1 - norm(y - y_hat_raw)/norm(y - mean(y)));

fprintf('Manual LS:\n  a = %.4f, b = %.4f, fit = %.2f%%\n', a_raw, b_raw, fit_raw);
disp('Gs_raw (manual LS, continuous):'); disp(Gs_raw);

data_id = iddata(y, u, h);
model_arx = arx(data_id, [1 1 1]);
Gs_arx = d2c(tf(model_arx), 'zoh');

[y_hat_arx, fit_arx, ~] = compare(data_id, model_arx);

fprintf('ARX Model:\n  Fit = %.2f%%\n', fit_arx);
disp('Gs_arx (ARX, continuous):'); disp(Gs_arx);

figure('Name','Continuous-Time Comparison','Position',[100 100 900 600]);
subplot(2,1,1)
plot(t, y, 'k', t, y_hat_raw, 'b--', t, y_hat_arx.OutputData, 'r-.', 'LineWidth', 1.2);
legend('Measured', sprintf('LS (%.1f%%)', fit_raw), sprintf('ARX (%.1f%%)', fit_arx), 'Location','best');
xlabel('Time [s]'); ylabel('\alpha [rad]');
title('Model Output vs Measured Data'); grid on;
