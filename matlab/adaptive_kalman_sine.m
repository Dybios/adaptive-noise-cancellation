% State space model is defined by following
% x_k+1 = x_k + v_k, where v_k is the process noise
% y_k+1 = x_k+1 + w_k+1, where w_k is the measurement noise

%% generate a demo set of input noisy sine wave signals
Fs = 1000; % Sampling frequency (Hz) (arduino = 8600Hz approx)
D = 60;     % Duration of the signal (s)
t = 0:1/Fs:D-1/Fs; % Time vector
f = 60;       % Frequency of the sine wave (Hz)
A = 1;        % Amplitude of the sine wave

sine_wave1 = (A * sin(2*pi*f*t)) + 1; % shifted up for all +ve
sine_wave2 = (A * sin(2*pi*f*t)) + 1;
sine_wave3 = (A * sin(2*pi*f*t)) + 1;
sine_wave4 = (A * sin(2*pi*f*t)) + 1;
noise_1 = 0.5 * randn(size(t)); % Generate random noise with amplitude 0.1
noise_2 = 0.15 * randn(size(t)); % Generate random noise with amplitude 0.15
noise_3 = 0.02 * randn(size(t)); % Generate random noise with amplitude 0.2
noise_4 = 0.2 * randn(size(t)); % Generate random noise with amplitude 0.07

% measured input signals
input_signal_i = transpose(sine_wave1 + noise_1);
input_signal1 = transpose(sine_wave2 + noise_2);
input_signal2 = transpose(sine_wave3 + noise_3);
input_signal3 = transpose(sine_wave4 + noise_4);

% Initialize common stuff
N = 6; % Number of heartbeats to average upon

%% find ECG complexes and their lengths by QRS detection

% call pam-tomkins algo on input_signal to find indexes between two peaks

% set T as 120% of mean interval between consecutive heartbeats

T = 50; % dummy length of each sine wave period; assumed equal

%% Initialize state estimate and covariance
% y_k is Tx1 matrix from the input signal (use this iteratively later)
y = input_signal_i(1:T);
y1 = input_signal1(1:T);
y2 = input_signal2(1:T);
y3 = input_signal3(1:T);

% do dummy sinewave stuff; get first N sinewave complexes for x_hat
% x0 = input_signal_i(1:T);
% x0 = x0 + input_signal_i((1+T) : (T+T));
% x0 = x0 + input_signal_i(1+(2*T) : T+(2*T));
% x0 = x0 + input_signal_i(1+(3*T) : T+(3*T));
% x0 = x0 + input_signal_i(1+(4*T) : T+(4*T));
x0 = zeros(T, N);
for k = 1:N
    x0(:,k) = input_signal_i(1+((k-1)*T) : T+((k-1)*T));
end
x_hat = mean(x0, 2); % initial estimated state

P = cov(x_hat); % initial estimated state covariance
Q = 0.01; % initial process covariance
K = zeros(1, T); % Kalman gain
Y_i = zeros(T, 1);

%% intialize measurement noise covariance
% 
% % find the estimate coeffs from control input signals
% Y_i(1,1) = y1;
% Y_i(1,2) = y2;
% Y_i(1,3) = y3;
% 
% y_coeff_hat = (transpose(Y_i) * Y_i) \ transpose(Y_i) * yi; % eq4
% 
% % find the measurement estimate
% y_i_hat = Y_i * y_coeff_hat;
% 
% % measurement noise vector
% w_hat = y - y_i_hat; % eq3
% 
% % Measurement noise covariance; initial estimate
% R = cov(w_hat);

%% iterate through the measured data

% Adaptive Kalman Filter
for k = 1:D
    if ((k*T) > size(input_signal_i))
        break;
    end

    %% estimate measurement noise covariance
    % update y to hold latest recorded ECG complex of length T
    y = input_signal_i(1+((k-1)*T) : T+((k-1)*T));
    y1 = input_signal1(1+((k-1)*T) : T+((k-1)*T));
    y2 = input_signal2(1+((k-1)*T) : T+((k-1)*T));
    y3 = input_signal3(1+((k-1)*T) : T+((k-1)*T));

    % find the estimate coeffs from control input signals
    Y_i(:,1) = y1;
    Y_i(:,2) = y2;
    Y_i(:,3) = y3;

    y_coeff_hat = (transpose(Y_i) * Y_i) \ transpose(Y_i) * y; % eq4

    % find the measurement estimate
    y_i_hat = Y_i * y_coeff_hat;

    % measurement noise vector
    w_hat = y - y_i_hat; % eq3

    % Measurement noise covariance
    R = cov(w_hat);

    % skip the first iteration
    if (k == 1)
        continue;
    end

    %% update process noise covariance
    % obtain model residual between measured and estimated
    model_residual = y - x_hat; % eq11
    Q = ((1/T) * transpose(model_residual) * model_residual) - P - R; % eq14
    if (Q < 0)
        Q = 0;
    end

    %% calculate kalman gain, state estimates and its covariances
    K = (P + Q) / (R + P + Q); % eq9
    x_hat = x_hat + K * (y - x_hat); % eq7
    P = P + Q - K * (P + Q); % eq8
end

% Plot results
t = 1:T;
figure;
subplot(2,1,1);
plot(t, y, 'b', t, x_hat, 'r');
xlabel('Time');
ylabel('Signal');
legend('Input Signal', 'Estimated Signal');
title('Input Signal vs Estimated Signal');

% subplot(2,1,2);
% plot(t, R, 'g');
% xlabel('Time');
% ylabel('Measurement Noise Covariance (R)');
% title('Measurement Noise Covariance (R)');
% 
% figure;
% plot(t, K, 'm');
% xlabel('Time');
% ylabel('Kalman Gain');
% title('Kalman Gain');
