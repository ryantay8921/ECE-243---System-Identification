close all; clear;

%% Load in input/output
outputs  = table2array(readtable("Experiment1/alpha_data.csv")); %input
inputs = table2array(readtable("Experiment1/u_data.csv")); %output

% Preprocess to radians and zero-mean
outputs  = outputs(5:end-3,2) - mean(outputs(5:end-3,2));
inputs = inputs(5:end-3,2)*pi/180 - mean(inputs(5:end-3,2)*pi/180);

%% Least Squares Parameter Estimation
T = 0.05;
Y = [];
S = [];
for k = 3:size(inputs,1)
    Y(k-2,1) = inputs(k,1) - 2*inputs(k-1,1) + inputs(k-2,1);
    S(k-2,1) = (T^2 / 2) * (outputs(k-1,1) + outputs(k-2,1));
end
c_ls = inv(S'*S)*S'*Y;

%% Transfer Function Estimation
s = tf('s');
G2 = c_ls / s^2;                        % outputs -> inputs
G1 = 15.85 / (s + 17.39);              % u -> outputs
G = G1 * G2;                            % u -> inputs
c2d(G, T, 'zoh');
rlocus(G);

%% Recursive Least Squares (RLS)
% Replace data files if needed
outputs  = table2array(readtable("Experiment2/alpha_data_2.csv"));
inputs = table2array(readtable("Experiment2/theta_data_2.csv"));

outputs  = outputs(2:end-3,2)*pi/180;
inputs = inputs(2:end-3,2)*pi/180;

% Initialize RLS
Y = [];
S = [];
e = [];
c = [];
R = [];
c0 = 29.1376;

c(2,1) = 30.6240;         % Initial parameter estimate
R(2,1) = 4.2859e-07;      % Initial R matrix
L = 1;                    % Forgetting factor

for k = 3:size(inputs,1)
    Y(k,1) = inputs(k,1) - 2*inputs(k-1,1) + inputs(k-2,1);
    S(k,1) = (T^2 / 2) * (outputs(k-1,1) + outputs(k-2,1));
    e(k,1) = Y(k,1) - S(k,1)*c(k-1,1);
    R(k,1) = L * R(k-1,1) + S(k,1)^2;
    c(k,1) = c(k-1,1) + inv(R(k,1)) * S(k,1) * e(k,1);
end

fprintf('Final RLS estimate of c: %.4f\n', c(end,1));
fprintf('Final RLS value of R: %.4e\n', R(end));

