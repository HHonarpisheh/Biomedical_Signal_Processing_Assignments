% myhomework #2 part #1
clear all;
close all;
% signal is a stereo audio signal with two channels (microphones)
[signal, Fs] = audioread('bird-stereo.wav');

%xcorr_result = xcorr(signal(:, 1), signal(:, 2));

% Find the peak in the cross-correlation
%[~, lag] = max(xcorr_result);

% MA filter order = Q
Q = 500; %abs(lag);
x = signal(:,1);
y = signal(:,2);
t = 0:1/Fs:(length(signal)-1)/Fs;
[b,e,yest] = myfit(y,x,Q);
% Calculate microphone spacing (assuming speed of sound ~ 343 m/s)
speed_of_sound = 343;  % meters per second
microphone_spacing = (speed_of_sound * Q) / Fs;  % in meters

figure;
subplot(2, 2, 1);
plot(t,y,"blue"); xlim([0 1.5]); xticks([0,0.5,1]); xlabel('time (s)');
title('input - mic 1');
legend('mic 1');
legend boxoff;
subplot(2, 2, 2);
plot(t,x,"green"); xlim([0 1.5]); xticks([0,0.5,1]); xlabel('time (s)');
hold on;
plot(t,e,"red");
hold off;
title('output - mic 2');
legend('mic 2','residual');
legend boxoff;
subplot(2, 1, 2);
tb = [0:0.02:10];
plot(tb,b,"blue"); xlim([0 11]); xticks([0:2:10]); xlabel('time (ms)');
ylim([-0.3 0.4]); yticks([-0.2:0.1:0.3]);
title('input/output impulse response, Q=500');

%fprintf('Estimated Microphone Spacing: %.2f meters\n', microphone_spacing);

function [b,e,yest] = myfit(y,x,Q)
% this function will predict y from x using Q delays
b = toeplitz(x(Q+1:end),x(Q+1:-1:1))\y(Q+1:end);
yest =  toeplitz(x,[x(1) zeros(1,Q)]) * b;
e = y - yest;
end