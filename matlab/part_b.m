clc; clear; close all;

alpha = readmatrix('Experiment2/alpha_data_1.csv');%inputs
theta = readmatrix('Experiment2/theta_data_1.csv');%outputs
bad = any(isnan([alpha theta]) | isinf([alpha theta]), 2);
alpha(bad,:) = [];
theta(bad,:) = [];

Ts = 0.05;
t = (0:length(alpha)-1)' * Ts;

u = deg2rad(alpha(:,2));
y = deg2rad(theta(:,2));
u = u - mean(u);
y = y - mean(y);

N = length(y);
Y = zeros(N-2,1);
S = zeros(N-2,1);

for k = 3:N
    Y(k-2) = y(k) - 2*y(k-1) + y(k-2);
    S(k-2) = (Ts^2 / 2) * (u(k-1) + u(k-2));
end
c_est = S \ Y;

s = tf('s');
b = 15.85;
a = 17.39;
Gs = b / (s + a);
Gtheta_s = (c_est / s^2) * Gs;

Gtheta_z = c2d(Gtheta_s, Ts, 'zoh');
rlocus(Gtheta_s);

y_hat = lsim(Gtheta_z, u, t);

figure('Name','Part B: Model vs Measured','Position',[100 100 950 450]);
plot(t(3:end), y(3:end), 'r', ...
     t, y_hat, 'b--', 'LineWidth', 1.4);
legend('Measured \theta', sprintf('Model Prediction (c = %.4f)', c_est), 'Location', 'best');
xlabel('Time (s)'); ylabel('\theta (rad)');
title('Part B: Model Prediction vs. Measured Output');
grid on;

fprintf('Estimated parameter c = %.4f\n', c_est);
disp('Continuous transfer function (theta/u):');
disp(Gtheta_s);
disp('Discrete transfer function (theta/u):');
disp(Gtheta_z);
