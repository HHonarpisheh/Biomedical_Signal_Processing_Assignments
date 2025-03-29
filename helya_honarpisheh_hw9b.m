close all;
clear all;

load hr_rest_stroop.mat; 
fs = mean(hr_rest);
x = hr_rest-median(hr_rest);


% Order of the AR/LPC model
P = 15;

% Autocorrelation matrix
R = xcorr(x, P, 'biased');
R_matrix = toeplitz(R(P+1:end-1));

% Autocorrelation vector
r_vector = R(P+2:end);

% Solve normal equations for AR
a_AR = -inv(R_matrix) * r_vector';

% AR model coefficients
a_AR = [1;a_AR]';

% Estimated power spectrum
[Pxx, f] = pwelch(x*sqrt(pi),128,[],[],fs);

% LPC coefficients using lpc
[a_LPC,e] = lpc(x, P);
h_LPC = filter(1,a_LPC,x);
[PLPC, fLPC] = pwelch(h_LPC*sqrt(pi),[], [], [], fs);

AR_model = filter(1,a_AR(2:end),x);
[PAR, fAR] = pwelch(AR_model*sqrt(pi),[], [], [], fs);


% Estimation AR and its power spectrum for AR model
estimation_AR = filter(a_AR, 1, x);
[P_estimation_AR, f_estimation_AR] = pwelch(estimation_AR*sqrt(pi),[], [], [], fs);
% Estimation LPC and its power spectrum for LPC model
estimation_LPC = filter(a_LPC, 1, x);
[P_estimation_LPC, f_estimation_LPC] = pwelch(estimation_LPC*sqrt(pi),[], [], [], fs);

% Linear prediction estimation error for AR
innovation_AR = filter([0 -a_AR(1:end)], 1, x);
[P_innovation_AR, f_innovation_AR] = pwelch(innovation_AR*sqrt(pi),128, [], [], fs);
% Linear prediction estimation error for LPC
innovation_LPC = filter([0 -a_LPC(1:end)], 1, x);
[P_innovation_LPC, f_innovation_LPC] = pwelch(innovation_LPC*sqrt(pi),128, [], [], fs);


% Plot estimated power spectrum and AR model spectrum and LPC Model Spectrum
figure;
subplot(4, 2, 1);
plot(f, db(Pxx), 'b');
hold on;
plot(fAR, db(PAR), 'r');
hold on;
plot(fLPC, db(PLPC), 'g');
title('Estimated Power Spectrum and AR Model Spectrum and LPC Model Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Estimated Power Spectrum', 'AR Model Spectrum','LPC Model Spectrum');


subplot(4, 2, 3);
plot(f_innovation_AR, (P_innovation_AR), 'g');
hold on;
plot(f_innovation_LPC, (P_innovation_LPC), 'r*');
title('Power Spectrum of AR & LPC Innovation Process');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('AR','LPC')
hold off;

subplot(4, 2, 4);
plot(f_estimation_AR, (P_estimation_AR), 'm');
hold on;
plot(f_estimation_LPC, (P_estimation_LPC), 'b*');
title(' Power Spectrum of AR & LPC estimation');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('AR','LPC')
hold off;

% Comparison plot for AR and LPC coefficients
subplot(4, 2, 2);
stem(a_AR,'b*');
hold on;
stem(a_LPC,'r.');
title('Comparison of AR and LPC Coefficients');
xlabel('Coefficient Index');
ylabel('Coefficient Value');
legend('AR', 'LPC');


% Original signal and estimations
subplot(4, 2, 5);
plot(x, 'b', 'DisplayName', 'Original Signal');
hold on;
plot(estimation_AR, 'r.');
plot(estimation_LPC, 'g--');
title('Original signal and estimations');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Original Signal', 'AR Estimation', 'LPC Estimation');

% Signal-to-noise ratio (SNR) as a function of model order P
P_values = 1:20; % Range of model orders
SNR_values_AR = zeros(size(P_values));
SNR_values_LPC = zeros(size(P_values));

for i = 1:length(P_values)
    current_P = P_values(i);
    
    % AR model SNR
    R = xcorr(x, current_P, 'biased');
    R_matrix = toeplitz(R(current_P+1:end-1));

% Autocorrelation vector
    r_vector = R(current_P+2:end);

% Solve normal equations for AR
    a_AR = -inv(R_matrix) * r_vector';

% AR model coefficients
    a_AR = [1; a_AR]';
    current_a_AR = a_AR;
    current_error_AR = filter(current_a_AR, 1, x);
    SNR_values_AR(i) = 10*log10(mean(abs(x-current_error_AR).^2) / mean(abs(current_error_AR).^2));
    
    % LPC model SNR
    current_a_LPC = lpc(x, current_P);
    current_error_LPC = filter(current_a_LPC, 1, x);
    SNR_values_LPC(i) = 10*log10(mean(abs(x-current_error_LPC).^2) / mean(abs(current_error_LPC).^2));
end

% Plot SNR as a function of model order P for both AR and LPC
subplot(4, 2, 6);
plot(P_values, SNR_values_AR, 'o-');
hold on;
plot(P_values, SNR_values_LPC, 's-');
title('Signal-to-Noise Ratio (SNR) vs. Model Order');
xlabel('Model Order (P)');
ylabel('SNR (dB)');
legend('AR Model', 'LPC Model');
grid on;

subplot(4,2,7);
plot((innovation_AR), 'g*');
hold on;
plot((innovation_LPC), 'm');
title('AR & LPC Innovation Process');
legend('AR Innovation','LPC Innovation');
hold off;

subplot(4,2,8);


%%

fs = mean(hr_stroop);
x = hr_stroop-median(hr_stroop);

% Order of the AR/LPC model
P = 10;

% Autocorrelation matrix
R = xcorr(x, P, 'biased');
R_matrix = toeplitz(R(P+1:end-1));

% Autocorrelation vector
r_vector = R(P+2:end);

% Solve normal equations for AR
a_AR = -inv(R_matrix) * r_vector';

% AR model coefficients
a_AR = [1;a_AR]';

% Estimated power spectrum
[Pxx, f] = pwelch(x*sqrt(pi),128,[],[],fs);

% LPC coefficients using lpc
[a_LPC,e] = lpc(x, P);
h_LPC = filter(1,a_LPC,x);
[PLPC, fLPC] = pwelch(h_LPC*sqrt(pi),[], [], [], fs);

AR_model = filter(1,a_AR(2:end),x);
[PAR, fAR] = pwelch(AR_model*sqrt(pi),[], [], [], fs);


% Estimation AR and its power spectrum for AR model
estimation_AR = filter(a_AR, 1, x);
[P_estimation_AR, f_estimation_AR] = pwelch(estimation_AR*sqrt(pi),[], [], [], fs);
% Estimation LPC and its power spectrum for LPC model
estimation_LPC = filter(a_LPC, 1, x);
[P_estimation_LPC, f_estimation_LPC] = pwelch(estimation_LPC*sqrt(pi),[], [], [], fs);

% Linear prediction estimation error for AR
innovation_AR = filter([0 -a_AR(1:end)], 1, x);
[P_innovation_AR, f_innovation_AR] = pwelch(innovation_AR*sqrt(pi),128, [], [], fs);
% Linear prediction estimation error for LPC
innovation_LPC = filter([0 -a_LPC(1:end)], 1, x);
[P_innovation_LPC, f_innovation_LPC] = pwelch(innovation_LPC*sqrt(pi),128, [], [], fs);


% Plot estimated power spectrum and AR model spectrum and LPC Model Spectrum
figure;
subplot(4, 2, 1);
plot(f, db(Pxx), 'b');
hold on;
plot(fAR, db(PAR), 'r');
hold on;
plot(fLPC, db(PLPC), 'g');
title('Estimated Power Spectrum and AR Model Spectrum and LPC Model Spectrum');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('Estimated Power Spectrum', 'AR Model Spectrum','LPC Model Spectrum');


subplot(4, 2, 3);
plot(f_innovation_AR, (P_innovation_AR), 'g');
hold on;
plot(f_innovation_LPC, (P_innovation_LPC), 'r*');
title('Power Spectrum of AR & LPC Innovation Process');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('AR','LPC')
hold off;

subplot(4, 2, 4);
plot(f_estimation_AR, (P_estimation_AR), 'm');
hold on;
plot(f_estimation_LPC, (P_estimation_LPC), 'b*');
title(' Power Spectrum of AR & LPC estimation');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
legend('AR','LPC')
hold off;

% Comparison plot for AR and LPC coefficients
subplot(4, 2, 2);
stem(a_AR,'b*');
hold on;
stem(a_LPC,'r.');
title('Comparison of AR and LPC Coefficients');
xlabel('Coefficient Index');
ylabel('Coefficient Value');
legend('AR', 'LPC');


% Original signal and estimations
subplot(4, 2, 5);
plot(x, 'b', 'DisplayName', 'Original Signal');
hold on;
plot(estimation_AR, 'r.');
plot(estimation_LPC, 'g--');
title('Original signal and estimations');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Original Signal', 'AR Estimation', 'LPC Estimation');

% Signal-to-noise ratio (SNR) as a function of model order P
P_values = 1:20; % Range of model orders
SNR_values_AR = zeros(size(P_values));
SNR_values_LPC = zeros(size(P_values));

for i = 1:length(P_values)
    current_P = P_values(i);
    
    % AR model SNR
    R = xcorr(x, current_P, 'biased');
    R_matrix = toeplitz(R(current_P+1:end-1));

% Autocorrelation vector
    r_vector = R(current_P+2:end);

% Solve normal equations for AR
    a_AR = -inv(R_matrix) * r_vector';

% AR model coefficients
    a_AR = [1; a_AR]';
    current_a_AR = a_AR;
    current_error_AR = filter(current_a_AR, 1, x);
    SNR_values_AR(i) = 10*log10(mean(abs(x-current_error_AR).^2) / mean(abs(current_error_AR).^2));
    
    % LPC model SNR
    current_a_LPC = lpc(x, current_P);
    current_error_LPC = filter(current_a_LPC, 1, x);
    SNR_values_LPC(i) = 10*log10(mean(abs(x-current_error_LPC).^2) / mean(abs(current_error_LPC).^2));
end

% Plot SNR as a function of model order P for both AR and LPC
subplot(4, 2, 6);
plot(P_values, SNR_values_AR, 'o-');
hold on;
plot(P_values, SNR_values_LPC, 's-');
title('Signal-to-Noise Ratio (SNR) vs. Model Order');
xlabel('Model Order (P)');
ylabel('SNR (dB)');
legend('AR Model', 'LPC Model');
grid on;

subplot(4,2,7);
plot((innovation_AR), 'g*');
hold on;
plot((innovation_LPC), 'm');
title('AR & LPC Innovation Process');
legend('AR Innovation','LPC Innovation');
hold off;

subplot(4,2,8);
