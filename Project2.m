clear all, close all, clc;
load("dataset_UAV.mat");
c = 3e8; 
lambda = c/f0; 

x_vec = -200:0.5:200;   
y_vec = 80:0.5:200;     
z_const = 0;          

[X, Y] = meshgrid(x_vec, y_vec);
Z = zeros(size(X)) + z_const;

SAR_image = zeros(size(X),'like',RCData);  

R = sqrt( (X - Sx(1)).^2 + (Y - Sy(1)).^2 + (Z - Sz(1)).^2 );

src_first = interp1(r_ax, RCData(:,1), R, 'linear', 0);

phase_term = exp( 1j * (4*pi/lambda) * R );

SAR_image = SAR_image + src_first .* phase_term;

num_pulses = size(RCData, 2); 

for idx = 2:num_pulses
    
    R = sqrt( (X - Sx(idx)).^2 + (Y - Sy(idx)).^2 + (Z - Sz(idx)).^2 );
    
    src_current = interp1(r_ax, RCData(:,idx), R, 'linear', 0);
    
    phase_term = exp( 1j * (4*pi/lambda) * R );

    SAR_image = SAR_image + src_current .* phase_term;
end

figure;
imagesc(x_vec, y_vec, abs(SAR_image));
colormap jet;
colorbar;
caxis([0 max(abs(SAR_image(:)))/10]); 
axis equal;
axis tight;
xlabel('X (m)');
ylabel('Y (m)');
title('Back Projection SAR Image (Magnitude)');

figure;
imshow('optical_image.jpg');
title('Optical Image of the Scene');

true_resolution = c/(2*B);

[num_ranges, num_pulses] = size(RCData);
freq_response = zeros(num_ranges, 1);

for pulse = 1:num_pulses
    windowed_pulse = RCData(:, pulse) .* hamming(num_ranges);
    
    fft_pulse = abs(fftshift(fft(windowed_pulse)));
    freq_response = freq_response + fft_pulse;
end

freq_response = freq_response / num_pulses;

% Frequency axis
sample_spacing = r_ax(2) - r_ax(1); % Spacing in range axis (m)
fs = c / (2 * sample_spacing);      
freq_axis = linspace(-fs/2, fs/2, num_ranges); 

% Estimatation of bandwidth
power_spectrum = freq_response.^2;
threshold = max(power_spectrum) * 0.5; % Adjusting threshold (e.g., 50% for -3 dB)
indices = find(power_spectrum >= threshold);
bandwidth = freq_axis(indices(end)) - freq_axis(indices(1)); 

% range resolution
range_resolution = c / (2 * bandwidth);

disp(['True Resolution (m): ', num2str(true_resolution)]);
disp(['Estimated Bandwidth (Hz): ', num2str(bandwidth)]);
disp(['Estimated Range Resolution (m): ', num2str(range_resolution)]);
%% part 4
SAR_amp = abs(SAR_image);  % SAR_image obtained from TDBP
N = 3;  % filter size
avg_filter = ones(N, N) / (N*N);
SAR_desp_3x3 = conv2(SAR_amp, avg_filter, 'same');
N = 5;
avg_filter_5x5 = ones(5, 5)/(5*5);
SAR_desp_5x5 = conv2(SAR_amp, avg_filter_5x5, 'same');
figure;
subplot(3,1,1); imagesc(x_vec, y_vec, SAR_amp); title('Original');
axis equal; axis tight; colorbar; colormap jet;

subplot(3,1,2); imagesc(x_vec, y_vec, SAR_desp_3x3); title('3x3 Filter');
axis equal; axis tight; colorbar; colormap jet;

subplot(3,1,3); imagesc(x_vec, y_vec, SAR_desp_5x5); title('5x5 Filter');
axis equal; axis tight; colorbar; colormap jet;
caxis([0 max(SAR_amp(:))/5]);  

%% part 5
x_min = -10; x_max = -2;  
y_min =  98; y_max = 104; 

dx = 0.1;
dy = 0.1;

x_vec = x_min:dx:x_max;
y_vec = y_min:dy:y_max;

[X, Y] = meshgrid(x_vec, y_vec);
Z = zeros(size(X)); 
SAR_patch = zeros(size(X),'like',RCData);
lambda = 3e8/f0;

num_pulses = size(RCData, 2);
for idx = 1:num_pulses
    R = sqrt( (X - Sx(idx)).^2 + (Y - Sy(idx)).^2 + (Z - Sz(idx)).^2 );
    src_current = interp1(r_ax, RCData(:,idx), R, 'linear', 0);
    phase_term = exp(1j*(4*pi/lambda)*R);
    SAR_patch = SAR_patch + src_current .* phase_term;
end
figure;
imagesc(x_vec, y_vec, abs(SAR_patch));
axis equal; axis tight; colormap jet; colorbar;
title('High-Resolution Patch Around Corner Reflector');
xlabel('X (m)'); ylabel('Y (m)');
% peak coordinates
[max_val, max_ind] = max(abs(SAR_patch(:)));
[peak_y, peak_x] = ind2sub(size(SAR_patch), max_ind);

%horizontal and vertical profiles
horizontal_profile = abs(SAR_patch(peak_y,:));
vertical_profile   = abs(SAR_patch(:,peak_x));

figure; 
subplot(2,1,1);
plot(x_vec, horizontal_profile);
title('Horizontal Profile Through Corner Reflector');
xlabel('X (m)'); ylabel('Amplitude');

subplot(2,1,2);
plot(y_vec, vertical_profile);
title('Vertical Profile Through Corner Reflector');
xlabel('Y (m)'); ylabel('Amplitude');
SAR_fft = fftshift(fft2(SAR_patch));
Nx = length(x_vec);
Ny = length(y_vec);

kx_spacing = 2*pi/(Nx*dx);
ky_spacing = 2*pi/(Ny*dy);

kx_vec = (-floor(Nx/2):ceil(Nx/2)-1)*kx_spacing;
ky_vec = (-floor(Ny/2):ceil(Ny/2)-1)*ky_spacing;

figure;
imagesc(kx_vec, ky_vec, abs(SAR_fft));
axis equal; axis tight; colormap jet; colorbar;
title('Spatial Frequency Domain of the Corner Reflector Patch');
xlabel('kx (rad/m)'); ylabel('ky (rad/m)');

% Horizontal Profile FWHM
[max_h, idx_h_peak] = max(horizontal_profile); 
half_max_h = max_h / 2;

% Finding points where profile crosses half-max
left_h = find(horizontal_profile(1:idx_h_peak) < half_max_h, 1, 'last');
right_h = find(horizontal_profile(idx_h_peak:end) < half_max_h, 1, 'first') + idx_h_peak - 1;

x_left_h = interp1(horizontal_profile(left_h:left_h+1), x_vec(left_h:left_h+1), half_max_h);
x_right_h = interp1(horizontal_profile(right_h-1:right_h), x_vec(right_h-1:right_h), half_max_h);

% Horizontal Resolution
horizontal_resolution = x_right_h - x_left_h;

% Vertical Profile FWHM
[max_v, idx_v_peak] = max(vertical_profile); 
half_max_v = max_v / 2;

top_v = find(vertical_profile(1:idx_v_peak) < half_max_v, 1, 'last');
bottom_v = find(vertical_profile(idx_v_peak:end) < half_max_v, 1, 'first') + idx_v_peak - 1;

y_top_v = interp1(vertical_profile(top_v:top_v+1), y_vec(top_v:top_v+1), half_max_v);
y_bottom_v = interp1(vertical_profile(bottom_v-1:bottom_v), y_vec(bottom_v-1:bottom_v), half_max_v);

% Vertical Resolution
vertical_resolution = y_bottom_v - y_top_v;

% Display Results
fprintf('Horizontal Resolution: %.3f m\n', horizontal_resolution);
fprintf('Vertical Resolution: %.3f m\n', vertical_resolution);

%% part 6
c = 3e8;
lambda = c/f0;  
rng(0); 
noise_std_values = [lambda/10, lambda/6, lambda/4, lambda/2];


Sx_orig = Sx;
Sy_orig = Sy;
Sz_orig = Sz;
for ns = 1:length(noise_std_values)
    noise_std = noise_std_values(ns);

    Sx_noisy = Sx_orig + noise_std*randn(size(Sx_orig));
    Sy_noisy = Sy_orig + noise_std*randn(size(Sy_orig));
    Sz_noisy = Sz_orig + noise_std*randn(size(Sz_orig));

    SAR_image_noisy = zeros(size(X),'like',RCData); 
    
    num_pulses = size(RCData, 2);
    for idx = 1:num_pulses
        R = sqrt( (X - Sx_noisy(idx)).^2 + (Y - Sy_noisy(idx)).^2 + (Z - Sz_noisy(idx)).^2 );
        src_current = interp1(r_ax, RCData(:,idx), R, 'linear', 0);
        phase_term = exp(1j*(4*pi/lambda)*R);
        SAR_image_noisy = SAR_image_noisy + src_current.*phase_term;
    end
    
    SAR_images_noisy{ns} = SAR_image_noisy; 
end
figure;
subplot(5,1,1);
imagesc(x_vec, y_vec, abs(SAR_image)); 
title('Original');
axis equal; axis tight; colormap jet; colorbar;
caxis([0 max(abs(SAR_image(:)))/10]);

for ns = 1:length(noise_std_values)
    subplot(5,1,ns+1);
    imagesc(x_vec, y_vec, abs(SAR_images_noisy{ns}));
    title(['Noise std: ', num2str(noise_std_values(ns), '%.2f'), ' m']);
    axis equal; axis tight; colormap jet; colorbar;
    caxis([0 max(abs(SAR_image(:)))/10]);
end
