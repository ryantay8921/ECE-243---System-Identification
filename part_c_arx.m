clc; clear; close all;

Ts = 0.05;

theta_all = []; u_all = [];
for i = 1:5
    alpha = readmatrix(sprintf('Experiment2/alpha_data_%d.csv', i));
    theta = readmatrix(sprintf('Experiment2/theta_data_%d.csv', i));
    bad = any(isnan([alpha theta]) | isinf([alpha theta]), 2);
    alpha(bad,:) = []; theta(bad,:) = [];
    N = min(size(alpha,1), size(theta,1));
    u_all = [u_all; alpha(1:N,2)];
    theta_all = [theta_all; theta(1:N,2)];
end

u_all = movmean(u_all, 5);
theta_all = movmean(theta_all, 5);

na = 2; nb = 2;
N = length(theta_all);
theta_rls = zeros(na+nb,1);
P = 1e3 * eye(na+nb);
lam = 1;

theta_hist = zeros(na,1); u_hist = zeros(nb,1);
theta_est = zeros(N,1);
param_log = zeros(N, na+nb);

for k = 3:N
    % Regressor vector: [-y(k-1), -y(k-2), u(k-1), u(k-2)]
    phi_k = [-theta_all(k-1); -theta_all(k-2); u_all(k-1); u_all(k-2)];

    y_pred = phi_k' * theta_rls;
    e_k = theta_all(k) - y_pred;

    K = P * phi_k / (lam + phi_k' * P * phi_k);
    theta_rls = theta_rls + K * e_k;
    P = (P - K * phi_k' * P) / lam;
    theta_est(k) = y_pred;
    param_log(k,:) = theta_rls';
    K_log(k,:) = K';
    P_diag_log(k,:) = diag(P)';
end


a1 = theta_rls(1); a2 = theta_rls(2);
b1 = theta_rls(3); b2 = theta_rls(4);
c = b1+b2;
Gq_rls = tf([0 b1 b2], [1 a1 a2], Ts);

figure;
plot((0:N-1)*Ts, theta_all, 'b', ...
     (0:N-1)*Ts, theta_est, 'r--', 'LineWidth', 1.3)
legend('Measured \theta', 'RLS Prediction')
xlabel('Time (s)'), ylabel('\theta (rad)')
title('RLS Model vs Measured Output'), grid on