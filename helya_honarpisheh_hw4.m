close all;
clear all;

load('gamma.mat');
load('NotchFilter.mat');

% Extract the signal with 60Hz noise
signal = osc(:,1);

% Apply the notch filter to the signal using sosfilt (IIR filter)
filtered_signal = sosfilt(SOS, signal)*prod(G);
% Calculate the impulse response of the notch filter
% Compute the filter coefficients for an IIR notch filter
b0 = SOS(1,1);
b1 = SOS(1,2);
b2 = SOS(1,3);
a0 = SOS(1,4);
a1 = SOS(1,5);
a2 = SOS(1,6);
impulse_response = filter([b0 b1 b2],[a0 a1 a2], [1 zeros(1, 1000-1)]);
% Compute the magnitude and phase responses from the impulse response
H = fft(impulse_response);
magnitude_response = abs(H);
phase_response = angle(H);

% Compute and display the spectrograms
figure;
subplot(1, 2, 1);
specgram(signal, 512, fs);
title('Spectrogram of Original Signal');
colorbar;
subplot(1, 2, 2);
specgram(filtered_signal, 512, fs);
title('Spectrogram of Filtered Signal');
colorbar;
% Plot the impulse response
figure;
subplot(3, 1, 1);
stem(impulse_response);
title('Impulse Response');
xlabel('Samples');
ylabel('Amplitude');
% Plot the magnitude response
subplot(3, 1, 2);
freq = (0:length(magnitude_response)-1) * (fs / length(magnitude_response));
plot(freq, magnitude_response);
title('Magnitude Response');
xlabel('Frequency (Hz)');
xlim([0, fs/2]); % Limit x-axis to Nyquist frequency
ylabel('Magnitude');
% Plot the phase response
subplot(3, 1, 3);
plot(freq, phase_response);
title('Phase Response');
xlabel('Frequency (Hz)');
xlim([0, fs/2]); % Limit x-axis to Nyquist frequency
ylabel('Phase (radians)');

% Plot the original and filtered signals
figure;
t = (0:length(signal) - 1) / fs;
plot(t, signal, 'b', t, filtered_signal, 'r');
title('Original and Filtered Signals');
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal', 'Filtered Signal');
