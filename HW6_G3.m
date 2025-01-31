clear all
close all
clc
dataTable=load("arxdata.mat");

data=dataTable.arx;

time=data(:,1);
y=data(:,2);

figure(1)
plot(time,y)
xlabel('time')
ylabel('y')

% Generate the input signal u(t)
u = sin(2 * time);

% Define the lagged data points for the ARX model
n = 3; % Number of lagged output terms (phi coefficients)
p = 2; % Number of lagged input terms (gamma coefficients)

% Ensure we have enough data for the lagging
k_start = max(n, p) + 1; % Start index for y_k to ensure all lags exist
N = length(time) - k_start + 1; % Total number of usable data points

% Create H matrix and y vector
H = zeros(N, n + p); % Matrix for lagged terms
y_vec = zeros(N, 1); % Output vector

for k = k_start:length(time)
    % Fill the H matrix with lagged y and u terms
    H(k - k_start + 1, :) = [
        y(k-1), y(k-2), y(k-3), ...  % Lagged y terms
        u(k-1), u(k-2)              % Lagged u terms
    ];
    % Fill the corresponding y value
    y_vec(k - k_start + 1) = y(k); % y_k
end

% Solve the least squares problem for the coefficients
theta = (H' * H) \ (H' * y_vec);

% Extract the coefficients
phi_1 = theta(1);
phi_2 = theta(2);
phi_3 = theta(3);
gamma_1 = theta(4);
gamma_2 = theta(5);

% Display the results
disp('Estimated Coefficients:');
fprintf('phi_1 = %.4f\n', phi_1);
fprintf('phi_2 = %.4f\n', phi_2);
fprintf('phi_3 = %.4f\n', phi_3);
fprintf('gamma_1 = %.4f\n', gamma_1);
fprintf('gamma_2 = %.4f\n', gamma_2);

% Optional: Compare predicted vs actual output
y_pred = H * theta; % Predicted values using the ARX model
figure(2)
plot(time(k_start:end), y_vec, 'b', 'LineWidth', 1.5); hold on;
plot(time(k_start:end), y_pred, 'r--', 'LineWidth', 1.5);
legend('Actual Output', 'Predicted Output');
xlabel('Time');
ylabel('y');
title('ARX Model: Actual vs Predicted Output');
grid on;