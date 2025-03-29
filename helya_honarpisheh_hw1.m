clear all;
close all;

load ('eeg.mat','fs','x');

x = x(1:64,:);

T = size(x,2)/fs; % length of signal in seconds
dt = 1/fs; % sampling interval in seconds
t = (1:size(x,2))*dt*1000; % sampling times in ms

% figure
% subplot(3,1,1); plot(t,x(1,:)); xlabel('time (ms)');
% subplot(3,1,2); plot(t,x(2,:)); xlabel('time (ms)');
% subplot(3,1,3); imagesc(x); xlabel('time (ms)'); ylabel('channels');

figure(1)
subplot(2,1,1); 
plot(detrend(x(2,:))/max(x(2,:)),'color',"red"); hold on;
plot(x(8,:)/max(x(8,:))+2,'color',"green"); hold on;
plot(-x(1,:)/max(-x(1,:))+4,'color',"blue"); hold off;
legend('VEOG','F1','F5','Location','northwest');
subplot(2,2,4); 
plot(-x(1,:),x(2,:),'.','color',"blue"); xlabel('VEOG'); 
ylabel('F5'); 
xlim([-1000 2000]); ylim([-500 2000]);
subplot(2,2,3); 
plot(-x(1,:),x(8,:),'.','color',"blue"); 
xlabel('VEOG'); ylabel('F1'); 
xlim([-1000 2000]); ylim([-500 1500]);

figure(2)
t = [-1:1/200:1]; y=sin(2*pi*10*t); %frequency = 10 Hz & fs = 200
subplot(2,3,1);plot(t,y);xlim([0 1]);title('sine');hold on;
subplot(2,3,2);plot(t,sin(t));xlim([-1 1]);title('squashing non-linearity');hold on; 
subplot(2,3,3);plot(t,atan(y));xlim([0 1]);title('squashed sine');hold on; 
subplot(2,3,4);periodogram(y,[],200/5,200); %fs = 200
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Spectrum');
grid on;
hold on;
subplot(2,3,6);periodogram(atan(y),[],200/5,200); %fs = 200
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Spectrum');
grid on;
hold off;

figure(3)
subplot(2,1,1);
mean_value = 0;
standard_deviation = 1; 
sample_size = 2000;      % Number of samples
gaussian_noise = mean_value + standard_deviation * randn(sample_size, 1);
plot(gaussian_noise);
xticks([500,1000,1500,2000]);
ylim([-5 5]);
title('Gaussian Noise');
subplot(2,1,2);
samples = randn(2000,1);
hist(randn(2000,1), 50,'FaceColor','blue');
xlabel('σ, 50 bins'); xlim([-4 4]); xticks([-4:2:4]);
ylabel('count'); ylim([0 140]);
title('Histogram');
mu = mean(samples); sigma = std(samples);
x = linspace(-3, 3, 2000);  % Create a range of x values
pdf_curve = exp(-(x-mu).^2./2./sigma^2)./sqrt(2*pi)./sigma;
% Scale the PDF curve
bin_centers = linspace(-3, 3, 51);
bin_width = bin_centers(2) - bin_centers(1);
scaling_factor = bin_width * length(samples);  % Scale factor
scaled_pdf_curve = pdf_curve * scaling_factor;
hold on;
plot(x, scaled_pdf_curve, "red");
hold off;

figure(4)
load ('tong_emg.mat','tong_emg');
subplot(2,2,1);
plot(tong_emg/std(tong_emg));
xlim([0 length(tong_emg)]); xticks([500,1000,1500,2000]);
ylim([-5 5]);
title('Tongue EMG');
subplot(2,2,2);
for i = 11:length(tong_emg)-10
    s = std(tong_emg(i-10:i+10));
    y1(i) = tong_emg(i)/s;
end
plot(y1);
xlim([0 length(tong_emg)]); xticks([500,1000,1500,2000]);
ylim([-5 5]);
title('Normalized Tongue EMG');
subplot(2,2,3);
tong = tong_emg/std(tong_emg);
hist(tong, 50,'FaceColor','blue');
xlabel('σ, 50 bins'); xlim([-4 6]); xticks([-4:2:6]);
ylabel('count'); 
title('Histogram');
mu1 = mean(tong); sigma1 = std(tong);
x1 = linspace(-3, 3, length(tong));
pdf_curve1 = exp(-(x1-mu1).^2./2./sigma1^2)./sqrt(2*pi)./sigma1;
bin_centers1 = linspace(-3, 3, 51);
bin_width1 = bin_centers(2) - bin_centers(1);
scaling_factor1 = bin_width1 * length(tong);
scaled_pdf_curve1 = pdf_curve1 * scaling_factor1;
hold on;
plot(x1, scaled_pdf_curve1, "red");
hold off;
subplot(2,2,4);
hist(y1, 50,'FaceColor','blue');
xlabel('σ, 50 bins'); xlim([-4 4]); xticks([-4:2:4]);
ylabel('count'); ylim([0 140]); yticks([0:20:140]);
title('Histogram');
mu2 = mean(y1); sigma2 = std(y1);
x2 = linspace(-3, 3, length(y1));
pdf_curve2 = exp(-(x2-mu2).^2./2./sigma2^2)./sqrt(2*pi)./sigma2;;
bin_centers2 = linspace(-3, 3, 51);
bin_width2 = bin_centers(2) - bin_centers(1);
scaling_factor2 = bin_width2 * length(y1);
scaled_pdf_curve2 = pdf_curve2 * scaling_factor2;
hold on;
plot(x2, scaled_pdf_curve2, "red");
hold off;
