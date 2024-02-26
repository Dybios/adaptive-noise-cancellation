% State space model is defined by following
% x_k+1 = x_k + v_k, where v_k is the process noise
% y_k+1 = x_k+1 + w_k+1, where w_k is the measurement noise

%% input the audio files
[input_signal_i, fs_i] = audioread("../data/babble_0dB.opus");
[input_signal1, fs1] = audioread("../data/babble_10dB.opus");

% averaging size for initial state estimate
N = 5;

% fixed buffer size
T = 10000;

% pad data to align with buffer
rem = mod(length(input_signal_i),T);
input_signal_i = paddata(input_signal_i, length(input_signal_i) + (T - rem));
input_signal1 = paddata(input_signal1, length(input_signal1) + (T - rem));
input_estimate = input_signal_i;
data_count = input_signal_i / T; % number of times to iterate over T

% y_k is Tx1 matrix from the input signal (use this iteratively later)
y = input_signal_i(1:T);
y1 = input_signal1(1:T);

%% initial measurement matrix estimate
% find the estimate coeffs from control input signals
Y_i(:,1) = y1;
y_coeff_hat = (transpose(Y_i) * Y_i) \ transpose(Y_i) * y; % eq4

% find the measurement estimate
y_i_hat = Y_i * y_coeff_hat;

% measurement noise vector
w_hat = y - y_i_hat; % eq3

% Measurement noise covariance
R = cov(w_hat);

%% Initialize state estimate and covariance
% do dummy sinewave stuff; get first N sinewave complexes for x_hat
x0 = zeros(T, N);
for k = 1:N
    x0(:,k) = input_signal_i(1+((k-1)*T) : T+((k-1)*T));
end
x_hat = mean(x0, 2); % initial estimated state

P = cov(x_hat); % initial estimated state covariance
Q = 0.01; % initial process covariance
K = zeros(1, T); % Kalman gain
Y_i = zeros(T, 1);

%% iterate through the recorded data
% Adaptive Kalman Filter
for k = 2:data_count
    %% calculate kalman gain, state estimates and its covariances
    K = (P + Q) / (R + P + Q); % eq9
    x_hat = x_hat + K * (y - x_hat); % eq7
    P = P + Q - K * (P + Q); % eq8

    %% update process noise covariance
    % obtain model residual between measured and estimated
    model_residual = y - x_hat; % eq11
    if (mod(k, N) == 0)
        model_residual = mean(model_residual, 2);
    end
    Q = ((1/T) * transpose(model_residual) * model_residual) - P - R; % eq14
    if (Q < 0)
        Q = 0;
    end

    %% estimate measurement noise covariance
    % update y buffer to hold latest audio sample
    y = input_signal_i(1+((k-1)*T) : T+((k-1)*T));
    y1 = input_signal1(1+((k-1)*T) : T+((k-1)*T));

    % find the estimate coeffs from control input signals
    Y_i(:,1) = y1;
    y_coeff_hat = (transpose(Y_i) * Y_i) \ transpose(Y_i) * y; % eq4

    % find the measurement estimate
    y_i_hat = Y_i * y_coeff_hat;

    % measurement noise vector
    w_hat = y - y_i_hat; % eq3

    % Measurement noise covariance
    R = cov(w_hat); % eq4

    input_estimate(1+((k-1)*T) : T+((k-1)*T)) = x_hat;
end

% subtract estimated with recorded input to extract noise
noise = input_signal1 - input_estimate;

% subtract noise from original signal to get clean output
output_signal = input_signal_i - noise; % todo: fix this

% Plot results
t = 1:size(output_signal);
figure;
subplot(2,1,1);
plot(t, input_signal_i, 'b', t, noise, 'r');
xlabel('Time');
ylabel('Signal');
legend('Input Signal', 'Estimated Signal');
title('Input Signal vs Estimated Signal');

sound(output_signal, fs_i);
