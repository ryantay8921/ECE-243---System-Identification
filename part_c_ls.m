close all; clear; clc;

T = 0.05;
L = 1;

for i = 1:5
    fprintf('\n=== Dataset %d ===\n', i);

    alpha = readmatrix(sprintf('Experiment2/alpha_data_%d.csv', i));
    theta = readmatrix(sprintf('Experiment2/theta_data_%d.csv', i));
    alpha = deg2rad(alpha(2:end-3,2));% [rad]
    theta = deg2rad(theta(2:end-3,2));% [rad]
    N     = numel(theta);

    c = 36.5616;
    R = 4.2859e-07;

    for k = 3:N
        Y = theta(k) - 2*theta(k-1) + theta(k-2);
        S = (T^2/2) * (alpha(k-1) + alpha(k-2));
        e = Y - S*c;
        R = L*R + S^2;
        c = c + (S/R)*e;
    end

    fprintf('Final RLS estimate of c: %.4f\n', c);
    fprintf('Final RLS value of R:   %.4e\n', R);

    th_hat        = zeros(N,1);
    th_hat(1:2)   = theta(1:2);
    gain          = (c*T^2/2);
    for k = 3:N
        th_hat(k) = gain*(alpha(k-1)+alpha(k-2)) + 2*th_hat(k-1) - th_hat(k-2);
    end
end
