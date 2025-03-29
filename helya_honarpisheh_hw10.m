close all;
clear all;

% Load your speech signal
% For example, if you have a WAV file:
[x, fs] = audioread('speech.wav');
x = x - mean(x);

% Step 1: Select a 16 ms segment
segment_duration = 0.016; % 16 ms
segment_samples = round(segment_duration * fs);
segment = x(0.7 * fs : (0.7 + segment_duration) * fs - 1);

% Step 2: Model as a deterministic harmonic process
% Step 3: Estimate pitch frequency (omega)
% Step 4: Determine the number of harmonics (K)
k_range = 1:15; % Test different values of k
f_range = 1 ./ (0.7 : 0.001 : 4.7);
T = segment_samples;

for k = k_range
    s = zeros(size(f_range));
    
    for i = 1:length(f_range)
        phase = (0.7*fs:0.7*fs+T-1)' * [1:k] * 2 * pi * f_range(i) / fs;
        S = [sin(phase), cos(phase)];
        b = S \ segment;
        s(i) = mean((segment - S * b).^2);
    end
    
    [~, idx] = min(s);
    fundamental_frequency = f_range(idx);
    omega = fundamental_frequency

    phase = (0.7*fs:0.7*fs+T-1)' * [1:k] * 2 * pi * f_range(idx) / fs;
    S = [sin(phase), cos(phase)];
    b = S \ segment;
    

    % Step 5: Measure orthogonality
    orthgnl = xcorr(S * b, segment - S * b);
    
    % Step 6: Model remaining noise as an AR process of order P
    P_range = 1:10; % Test different values of P
    Pnoise = zeros(length(segment), length(P_range));
    
    for p = P_range
        r = xcorr(segment - S * b,p);
        a_normal = [1; -toeplitz(r(p+1:end-1))\r(p+2:end)]';
        %a_ar = aryule(segment, p);
        err = filter(1, a_normal, segment - S * b);
        Pnoise(:, p) = abs(fft(err)).^2 / length(segment);
    end
    
    % Plotting
    figure;
    subplot(3, 2, 1);
    plot((0:length(segment)-1) / fs, [segment, segment - S * b]);
    xlabel('Time (s)');
    legend('Data', 'Linear prediction');
    
    subplot(3, 2, 2);
    fbin = (0:length(segment)-1) / length(segment) * fs;
  plot(fbin, db([abs(fft(segment)).^2 / length(segment), abs(fft(S * b)).^2 / length(segment), sum(Pnoise, 2)]));
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    ylabel('Signal Spectrum');
    legend('Data', 'Harmonic', 'Harmonic + Noise');
    
    subplot(3, 2, 3);
    plot(fbin, db(sum(Pnoise, 2)));
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    ylabel('Noise Spectrum');
    
    subplot(3, 2, 4);
    plot(fbin, db([abs(fft(segment)).^2 / length(segment), abs(fft(S * b)).^2 / length(segment), sum(Pnoise, 2)]));
    ylabel('Spectrum (dB)');
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    legend('Data', 'Harmonic', 'Harmonic + Noise');
    sgtitle(['k = ' num2str(k)]);
    
    % Display orthogonality plot
    subplot(3, 2, 5);
    plot(orthgnl);
    title('Orthogonality');
    xlabel('Sample Index');
    ylabel('Amplitude');
    xlim([1, length(orthgnl)]);
    ylim([-max(abs(orthgnl)), max(abs(orthgnl))]);
    
    % Display power spectra for different P
    subplot(3, 2, 6);
    plot(fbin, db(Pnoise));
    ylabel('Spectrum (dB)');
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    legend(arrayfun(@(p) ['P = ' num2str(p)], P_range, 'UniformOutput', false));
    title(['k = ' num2str(k)]);
end

%% Final figure with best k & P
k=6;
p=2;

    s = zeros(size(f_range));
    
    for i = 1:length(f_range)
        phase = (0.7*fs:0.7*fs+T-1)' * [1:k] * 2 * pi * f_range(i) / fs;
        S = [sin(phase), cos(phase)];
        b = S \ segment;
        s(i) = mean((segment - S * b).^2);
    end
    
    [~, idx] = min(s);
    fundamental_frequency = f_range(idx);
    
    omega = fundamental_frequency
    
    phase = (0.7*fs:0.7*fs+T-1)' * [1:k] * 2 * pi * f_range(idx) / fs;
    S = [sin(phase), cos(phase)];
    b = S \ segment;

    % Step 5: Measure orthogonality
    orthgnl = xcorr(S * b, segment - S * b);
    
    % Step 6: Model remaining noise as an AR process of order P
    Pnoise = zeros(length(segment), length(p));
    
        
        r = xcorr(segment - S * b,p);
        a_normal = [1; -toeplitz(r(p+1:end-1))\r(p+2:end)]';
        %a_ar = aryule(segment, p);
        err = filter(1, a_normal, segment - S * b);
        Pnoise(:, p) = abs(fft(err)).^2 / length(segment);
        
    
    % Plotting
    figure;
    subplot(3, 2, 1);
    plot((0:length(segment)-1) / fs, [segment, segment - S * b]);
    xlabel('Time (s)');
    legend('Data', 'Linear prediction');
    
    subplot(3, 2, 2);
    fbin = (0:length(segment)-1) / length(segment) * fs;
  plot(fbin, db([abs(fft(segment)).^2 / length(segment), abs(fft(S * b)).^2 / length(segment), sum(Pnoise, 2)]));
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    ylabel('Signal Spectrum');
    legend('Data', 'Harmonic', 'Harmonic + Noise');
    
    subplot(3, 2, 3);
    plot(fbin, db(sum(Pnoise, 2)));
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    ylabel('Noise Spectrum');
    
    subplot(3, 2, 4);
    plot(fbin, db([abs(fft(segment)).^2 / length(segment), abs(fft(S * b)).^2 / length(segment), sum(Pnoise, 2)]));
    ylabel('Spectrum (dB)');
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    legend('Data', 'Harmonic', 'Harmonic + Noise');
    sgtitle(['best k = ' num2str(k)]);
    
    % Display orthogonality plot
    subplot(3, 2, 5);
    plot(orthgnl);
    title('Orthogonality');
    xlabel('Sample Index');
    ylabel('Amplitude');
    xlim([1, length(orthgnl)]);
    ylim([-max(abs(orthgnl)), max(abs(orthgnl))]);
    
    % Display power spectra for different P
    subplot(3, 2, 6);
    plot(fbin, db(Pnoise));
    ylabel('Spectrum (dB)');
    xlabel('Frequency (Hz)');xlim([0 fs/2]);
    legend(arrayfun(@(p) ['best P = ' num2str(p)], p, 'UniformOutput', false));
    title(['best k = ' num2str(k)]);
