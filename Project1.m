%% Direction of Arrival
% parameters
clear all, close all, clc;
N = 58; % Number of antennas
dx = 1.948e-3; % Antenna spacing in meters
antenna_positions = (0:N-1) * dx;

%target
theta_true = 50; 
theta_rad = deg2rad(theta_true); 


%recieved signal
lambda = 3.896e-3; % Wavelength in meters
k = 2 * pi / lambda; % Wavenumber
phase_shifts = k * antenna_positions * sin(theta_rad); 

received_signal = exp(1j * phase_shifts); 

% DoA
delta_phi = angle(received_signal(2:end) .* conj(received_signal(1:end-1)));
delta_phi_avg = mean(delta_phi);

theta_est_rad = asin((delta_phi_avg * lambda) / (2 * pi * dx));
theta_est = rad2deg(theta_est_rad); 

%output
fprintf('True Angle: %.2f degrees\n', theta_true);
fprintf('Estimated Angle: %.2f degrees\n', theta_est);

%% Two Targets - Beamforming Method

% parameters
N = 58; 
dx = 1.948e-3; 
antenna_positions = (0:N-1) * dx; 

lambda = 3.896e-3;
k = 2 * pi / lambda; 

theta_true = [37 38]; % Angles of two targets in degrees 
theta_rad = deg2rad(theta_true); 
K = length(theta_true); 

% Reshape vectors to align dimensions for broadcasting
antenna_positions = antenna_positions(:); % N x 1 column vector
theta_rad = theta_rad(:)'; % 1 x K row vector

% Compute phase shifts matrix
phase_shifts = -k * antenna_positions * sin(theta_rad); % N x K matrix

% received signals from all targets
% Each column corresponds to a target
received_signals = exp(1j * phase_shifts); % N x K matrix

% Sum the signals across targets for each antenna
received_signal = sum(received_signals, 2); % N x 1 vector

% Beamforming Method

% scanning angles
theta_scan = -90:0.1:90; 
theta_scan_rad = deg2rad(theta_scan); 


beamforming_output = zeros(length(theta_scan), 1);


for idx = 1:length(theta_scan_rad)
    % steering vector for the scanning angle
    steering_vector = exp(-1j * k * antenna_positions * sin(theta_scan_rad(idx)));
    
    % beamformer output (power)
    beamforming_output(idx) = abs(steering_vector' * received_signal)^2;
end

% Normalizing
beamforming_output = beamforming_output / max(beamforming_output);


[pks, locs] = findpeaks(beamforming_output, 'MinPeakHeight', 0.5, 'SortStr', 'descend', 'NPeaks', K);

% Estimated angles
theta_est = theta_scan(locs);

% Sort the estimated angles
theta_est = sort(theta_est);

% Plot 

figure;
plot(theta_scan, beamforming_output, 'LineWidth', 1.5);
xlabel('Angle (degrees)');
ylabel('Normalized Power');
title('Beamforming Output');
grid on;
hold on;

plot(theta_est, pks, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Beamforming Output', 'Estimated Angles');


% Output
fprintf('True Angles: %.2f degrees\n', theta_true);
fprintf('Estimated Angles: %.2f degrees\n', theta_est);


%% part 4: 2D position estimation

clear all;
close all;
clc;

% Parameters
N = 58; 
dx = 1.948e-3; 
antenna_positions = (0:N-1) * dx; 

% Target Parameters
theta_true = 50; 
theta_rad = deg2rad(theta_true); 
r = 50; % Target range in meters

% Received Signal Parameters
lambda = 3.896e-3; 
k = 2 * pi / lambda; 


tau = 2 * r / 3e8; 

B = 1e9; 
Ts = 1 / (10 * B); 
Fs = 1 / Ts; 

t_max = tau + 1e-6; 
t = 0 : Ts : t_max; 

%Transmitted Signal
g_tx = sinc(B * t) .* exp(1j * 2 * pi * 77e9 * t); 

% Received Signal at Each Antenna
received_signal = zeros(N, length(t));


phase_shifts = -k * antenna_positions * sin(theta_rad); % 1 x N vector

% Generate Received Signals
for n = 1:N
   
    t_delayed = t - tau;
    signal = zeros(1, length(t));

    idx_positive = find(t_delayed >= 0);
    
    signal(idx_positive) = sinc(B * t_delayed(idx_positive)) .* ...
                            exp(1j * (2 * pi * 77e9 * t_delayed(idx_positive) + phase_shifts(n)));
    % Received signal at antenna n
    received_signal(n, :) = signal;
end



% Demodulate Received Signals
% Use repmat to match dimensions
received_signal_demod = received_signal .* repmat(exp(-1j * 2 * pi * 77e9 * t), N, 1);

% Perform Matched Filtering for Range Estimation
% Reference signal (matched filter)
ref_signal = sinc(B * t);

range_profile = zeros(N, length(t));

for n = 1:N
    R = ifft(fft(received_signal_demod(n, :)) .* conj(fft(ref_signal)));
    range_profile(n, :) = abs(R); 
end

range_profile_sum = sum(range_profile, 1);


[~, idx_peak] = max(range_profile_sum);
tau_est = t(idx_peak); 

% Estimate range
r_est = (tau_est * 3e8) / 2; 

% Perform DoA Estimation at Estimated Range
% Extract the received signal across antennas at the estimated time delay
signal_at_tau = received_signal_demod(:, idx_peak); % N x 1 vector

% phase differences between adjacent antennas
delta_phi = angle(signal_at_tau(2:end) .* conj(signal_at_tau(1:end-1)));

delta_phi_avg = mean(delta_phi);

theta_est_rad = asin((-delta_phi_avg * lambda) / (2 * pi * dx));
theta_est = rad2deg(theta_est_rad); % Convert to degrees

% Output Results
fprintf('True Range: %.2f meters\n', r);
fprintf('Estimated Range: %.2f meters\n', r_est);
fprintf('True Angle: %.2f degrees\n', theta_true);
fprintf('Estimated Angle: %.2f degrees\n', theta_est);

% Plot Range Profile
figure;
plot(t * 3e8 / 2, range_profile_sum, 'LineWidth', 1.5);
xlabel('Range (meters)');
ylabel('Amplitude');
title('Range Profile');
grid on;
xlim([r - 1, r + 1]); 
hold on;
plot(r_est, range_profile_sum(idx_peak), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Range Profile', 'Estimated Range');

% Plot Phase Across Antennas at Estimated Range
figure;
plot(1:N, angle(signal_at_tau), 'o-', 'LineWidth', 1.5);
xlabel('Antenna Index');
ylabel('Phase (radians)');
title('Phase Across Antennas at Estimated Range');
grid on;

%% part5 MIMO array
clear all;
close all;
clc;

%Parameters

c = 3e8; 

% Radar Specifications
f0 = 77e9; % Carrier frequency (Hz)
lambda = c / f0;
k = 2 * pi / lambda; 

% Antenna Array Parameters
NTx = 2; % Number of transmit antennas
NRx = 29; % Number of receive antennas
N_virtual = NTx * NRx; % Number of virtual elements

% Virtual Array Spacing
dx = 1.948e-3; % Virtual array spacing in meters

% Transmit Antenna Positions
Tx_positions = (0:NTx-1) * NRx * dx; % 1 x NTx

% Receive Antenna Positions
Rx_positions = (0:NRx-1) * dx; % 1 x NRx

% Target Parameters
theta_true = 41; 
theta_rad = deg2rad(theta_true); 
r = 60; 

% Time Delay Calculation
tau = 2 * r / c; 

% Sampling Parameters
B = 1e9; 
Fs = 10 * B; 
Ts = 1 / Fs; 
t_max = tau + 5e-6; 
t = 0 : Ts : t_max; 


% Transmitted Signal: Sinc pulse modulated at carrier frequency
g_tx = sinc(B * t) .* exp(1j * 2 * pi * f0 * t); % 1 x length(t)

% Received Signals at Each Tx-Rx Pair
received_signal = zeros(NTx, NRx, length(t));

[Tx_matrix, Rx_matrix] = meshgrid(Tx_positions, Rx_positions);
Tx_matrix = Tx_matrix'; % NTx x NRx
Rx_matrix = Rx_matrix'; % NTx x NRx

% Phase Shifts for Each (Tx_i, Rx_j) Pair
phi_matrix = -k * (Tx_matrix + Rx_matrix) * sin(theta_rad); % NTx x NRx

% Received Signals
for i = 1:NTx
    % Shift the transmitted signal by tau in samples
    shift_samples = round(tau / Ts); % Number of samples to shift
    if shift_samples < length(t)
        g_tx_shifted = [zeros(1, shift_samples), g_tx(1:end-shift_samples)];
    else
        g_tx_shifted = zeros(1, length(t));
    end
    
    for j = 1:NRx
        % Compute phase shift for Tx_i to Rx_j
        phi_ij = phi_matrix(i, j);
           
        g_rx = g_tx_shifted .* exp(1j * phi_ij);
        
        received_signal(i, j, :) = g_rx;
    end
end
% Demodulate Each Tx-Rx Pair Signal
carrier_conj = conj(exp(1j * 2 * pi * f0 * t)); % 1 x length(t)
received_signal_demod = received_signal .* reshape(carrier_conj, [1, 1, length(t)]); % NTx x NRx x length(t)


% Range Estimation via Matched Filtering

received_signal_demod_sum = squeeze(sum(sum(received_signal_demod, 1), 2)); % length(t) x 1

received_signal_demod_sum = received_signal_demod_sum.'; % 1 x length(t)

% Reference Signal for Matched Filtering (sinc pulse)
ref_signal = sinc(B * t); % 1 x length(t)

%Matched Filtering
R = ifft(fft(received_signal_demod_sum) .* conj(fft(ref_signal)));

%range profile
range_profile = abs(R);

[~, idx_peak] = max(range_profile);
tau_est = t(idx_peak); 

r_est = (tau_est * c) / 2; 

%DoA Estimation Using Virtual Array


V_pos = zeros(N_virtual, 1);
V_signal = zeros(N_virtual, 1);

n = 1;
for i = 1:NTx
    for j = 1:NRx
        % Virtual Array Position
        V_pos(n) = Tx_positions(i) + Rx_positions(j);
        
        % Signal at Estimated Time Index
        V_signal(n) = received_signal_demod(i, j, idx_peak);
        
        n = n + 1;
    end
end

[V_pos_sorted, idx_sort] = sort(V_pos);
V_signal_sorted = V_signal(idx_sort);

% Verify Uniform Spacing
virtual_spacing = diff(V_pos_sorted);
mean_spacing = mean(virtual_spacing);
fprintf('Virtual array mean spacing: %.6e meters\n', mean_spacing);

% Adjust for Mean Spacing 
dx_virtual = mean_spacing; 

% Phase Differences Between Adjacent Virtual Elements
delta_phi = angle(V_signal_sorted(2:end) .* conj(V_signal_sorted(1:end-1)));

% Unwrap Phase Differences
delta_phi_unwrapped = unwrap(delta_phi);


delta_phi_avg = mean(delta_phi_unwrapped);

% Estimate Angle
theta_est_rad = asin((-delta_phi_avg * lambda) / (2 * pi * dx_virtual));
theta_est_deg = rad2deg(theta_est_rad); % Convert to degrees

% Results

fprintf('True Range: %.2f meters\n', r);
fprintf('Estimated Range: %.2f meters\n', r_est);
fprintf('True Angle: %.2f degrees\n', theta_true);
fprintf('Estimated Angle: %.2f degrees\n', theta_est_deg);


% Plot Range Profile
figure;
plot(t * c / 2, range_profile, 'LineWidth', 1.5);
xlabel('Range (meters)');
ylabel('Amplitude');
title('Range Profile');
grid on;
xlim([r - 1, r + 1]); 
hold on;
plot(r_est, range_profile(idx_peak), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Range Profile', 'Estimated Range');

% Plot Phase Across Virtual Array at Estimated Range
figure;
plot(V_pos_sorted, angle(V_signal_sorted), 'o-', 'LineWidth', 1.5);
xlabel('Virtual Array Position (meters)');
ylabel('Phase (radians)');
title('Phase Across Virtual Array at Estimated Range');
grid on;

%% part 6: orthogonal waveform
clear all;
close all;
clc;

%Parameters
c = 3e8; 
f0 = 77e9; 
lambda = c / f0; 
k = 2 * pi / lambda; 

NTx = 2; 
NRx = 4; 
N_virtual = NTx * NRx; 

dx = lambda / 2; 

Tx_positions = (0:NTx-1) * NRx * dx; % 1 x NTx
Rx_positions = (0:NRx-1) * dx; % 1 x NRx

theta_true = 30; 
theta_rad = deg2rad(theta_true); 
r = 100; 

tau = 2 * r / c;

B = 100e6; 
Fs = 2 * B; 
Ts = 1 / Fs; 
t_max = tau + 1e-6; 
t = 0 : Ts : t_max; 

% OFDM Parameters
N_subcarriers = NTx;
subcarrier_spacing = B / N_subcarriers; % Subcarrier spacing
subcarrier_frequencies = f0 + (-floor(N_subcarriers/2):floor((N_subcarriers-1)/2)) * subcarrier_spacing;
% Generate Orthogonal Waveforms (OFDM)

g_tx = zeros(NTx, length(t)); % NTx x length(t)

for i = 1:NTx
    % Generate OFDM subcarrier signal for each transmitter
    g_tx(i, :) = exp(1j * 2 * pi * subcarrier_frequencies(i) * t);
end

received_signal = zeros(NTx, NRx, length(t));

[Tx_matrix, Rx_matrix] = meshgrid(Tx_positions, Rx_positions);
Tx_matrix = Tx_matrix'; % NTx x NRx
Rx_matrix = Rx_matrix'; % NTx x NRx

phi_matrix = -k * (Tx_matrix + Rx_matrix) * sin(theta_rad); % NTx x NRx

%Received Signals
for i = 1:NTx
    shift_samples = round(tau / Ts); 
    if shift_samples < length(t)
        g_tx_shifted = [zeros(1, shift_samples), g_tx(i, 1:end-shift_samples)];
    else
        g_tx_shifted = zeros(1, length(t));
    end
    
    for j = 1:NRx
        
        phi_ij = phi_matrix(i, j);
               
        g_rx = g_tx_shifted .* exp(1j * phi_ij);
        
        received_signal(i, j, :) = g_rx;
    end
end

% Matched Filtering 

matched_filter_output = zeros(NTx, NRx, length(t));


for i = 1:NTx
    
    ref_signal = conj(g_tx(i, :)); % 1 x length(t)
    
    for j = 1:NRx
        % Received Signal at Rx_j from Tx_i
        rx_signal = squeeze(received_signal(i, j, :)).'; % 1 x length(t)
        
        % Perform Matched Filtering
        R = ifft(fft(rx_signal) .* fft(ref_signal));
        
        matched_filter_output(i, j, :) = R;
    end
end
% Range Estimation

range_profile = zeros(NTx, length(t));

for i = 1:NTx
    range_profile(i, :) = sum(abs(matched_filter_output(i, :, :)), 2);
end

% Estimate Range for Each Transmitter
tau_est = zeros(NTx, 1);
r_est = zeros(NTx, 1);
idx_peak = zeros(NTx, 1);

for i = 1:NTx
    [~, idx_peak(i)] = max(range_profile(i, :));
    tau_est(i) = t(idx_peak(i));
    r_est(i) = (tau_est(i) * c) / 2;
end

r_est_avg = mean(r_est);

V_pos = zeros(N_virtual, 1);
V_signal = zeros(N_virtual, 1);

n = 1;
for i = 1:NTx
    for j = 1:NRx
    
        V_pos(n) = Tx_positions(i) + Rx_positions(j);
        
        V_signal(n) = matched_filter_output(i, j, idx_peak(i));
        
        n = n + 1;
    end
end

[V_pos_sorted, idx_sort] = sort(V_pos);
V_signal_sorted = V_signal(idx_sort);

delta_phi = angle(V_signal_sorted(2:end) .* conj(V_signal_sorted(1:end-1)));

delta_phi_unwrapped = unwrap(delta_phi);

delta_phi_avg = mean(delta_phi_unwrapped);

% Estimate Angle
dx_virtual = mean(diff(V_pos_sorted)); % Use actual spacing
theta_est_rad = asin((-delta_phi_avg * lambda) / (2 * pi * dx_virtual));
theta_est_deg = rad2deg(theta_est_rad); % Convert to degrees
% Display Results

fprintf('True Range: %.2f meters\n', r);
fprintf('Estimated Range: %.2f meters\n', r_est_avg);
fprintf('True Angle: %.2f degrees\n', theta_true);
fprintf('Estimated Angle: %.2f degrees\n', theta_est_deg);

% Plot Range Profiles for Each Transmitter

figure;
for i = 1:NTx
    plot(t * c / 2, range_profile(i, :), 'LineWidth', 1.5);
    hold on;
end
xlabel('Range (meters)');
ylabel('Amplitude');
title('Range Profiles for Each Transmitter');
legend('Tx1', 'Tx2');
grid on;

% Plot Phase Across Virtual Array at Estimated Range

figure;
plot(V_pos_sorted, angle(V_signal_sorted), 'o-', 'LineWidth', 1.5);
xlabel('Virtual Array Position (meters)');
ylabel('Phase (radians)');
title('Phase Across Virtual Array at Estimated Range');
grid on;
