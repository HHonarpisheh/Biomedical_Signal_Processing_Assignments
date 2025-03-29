close all;
clear all;

% Load the EEG data
eeg_data = load('EEG_128Hz.txt');

% Choose one channel (e.g., channel 1)
eeg_channel = eeg_data(:, 1);

% Define filter specifications
center_freq = 10;    % Center frequency of the bandpass filter (10 Hz)
bandwidth = 5;       % Bandwidth of the bandpass filter (5 Hz)
fs = 128;            % Sampling frequency of the EEG data

% Design an IIR Butterworth filter
order = 4; 
[b, a] = butter(order, [(center_freq - bandwidth/2) (center_freq + bandwidth/2)] / (fs/2), 'bandpass');

% Design a linear-phase FIR filter using fir1
fir_order = 64; 
fir_filter = fir1(fir_order, [(center_freq - bandwidth/2) (center_freq + bandwidth/2)] / (fs/2), 'bandpass');

% Apply the IIR filter to the EEG data
filtered_eeg_iir = filtfilt(b, a, eeg_channel);

% Apply the FIR filter to the EEG data
filtered_eeg_fir = conv(eeg_channel, fir_filter, 'same');

% Calculate power in the 10Hz band for all EEG channels
power_in_band = zeros(size(eeg_data, 2), 1);

for i = 1:size(eeg_data, 2)
    channel = eeg_data(:, i);
    filtered_channel = filtfilt(b, a, channel); %filter(fir_filter, 1, channel);
    power_in_band(i) = sum(abs(filtered_channel).^2) / length(filtered_channel);
end

disp('Power in the 10Hz band across channels:');
disp(power_in_band);


% Plot the original EEG signal and the filtered signals
t = (0:length(eeg_channel)-1) / fs;
figure;
subplot(3,1,1);
plot(t, eeg_channel);
title('Original EEG Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,2);
plot(t, filtered_eeg_iir);
title('Filtered EEG (IIR)');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3,1,3);
plot(t, filtered_eeg_fir);
title('Filtered EEG (FIR)');
xlabel('Time (s)');
ylabel('Amplitude');

% Calculate and plot the magnitude and phase responses of the FIR filter
freq_response = fft(fir_filter, 1024);
freq = linspace(0, fs/2, length(freq_response)/2);
magnitude_response = abs(freq_response(1:length(freq_response)/2));
phase_response = angle(freq_response(1:length(freq_response)/2));

figure;
subplot(3,1,1);
stem(fir_filter);
title('Impulse Response (FIR)');
xlabel('Samples');
ylabel('Amplitude');

subplot(3,1,2);
plot(freq, magnitude_response);
title('Magnitude Response (FIR)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

subplot(3,1,3);
plot(freq, phase_response);
title('Phase Response (FIR)');
xlabel('Frequency (Hz)');
ylabel('Phase (radians)');

% Calculate and plot the magnitude and phase responses of the IIR filter
impulse_response = filter(b, a , [1 zeros(1, 1024-1)]);
% Compute the magnitude and phase responses from the impulse response
H = fft(impulse_response);
magnitude_resp = abs(H);
phase_resp = angle(H);
% Plot the impulse response
figure;
subplot(3, 1, 1);
stem(impulse_response);
title('Impulse Response (IIR)');
xlabel('Samples');
ylabel('Amplitude');
% Plot the magnitude response
subplot(3, 1, 2);
freq = (0:length(magnitude_resp)-1) * (fs / length(magnitude_resp));
plot(freq, magnitude_resp);
title('Magnitude Response (IIR)');
xlabel('Frequency (Hz)');
xlim([0, fs/2]); % Limit x-axis to Nyquist frequency
ylabel('Magnitude');
% Plot the phase response
subplot(3, 1, 3);
plot(freq, phase_resp);
title('Phase Response (IIR)');
xlabel('Frequency (Hz)');
xlim([0, fs/2]); % Limit x-axis to Nyquist frequency
ylabel('Phase (radians)');
