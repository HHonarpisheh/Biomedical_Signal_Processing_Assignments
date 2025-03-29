clear all;
close all;

%% TWO conjugate poles, one in quadrant I and the conjugate in quadrant IV
% Define the pole locations
p1 = [0.5 + 0.5i, 0.5 - 0.5i];

% Define the zero locations 
z1 = [0.5 0.5]; 

% Generate and plot the impulse response
b1 = 1;  % Coefficient for the zero
figure;
plot(filter(b1 * poly(z1), poly(p1), [1 zeros(1, 100)]));
title('Impulse Response - Quadrant I and Quadrant IV Poles');
xlabel('Time (n)');
ylabel('Amplitude');

%% TWO conjugate poles, one in quadrant II and the conjugate in quadrant III
% Define the pole locations
p2 = [-0.5 + 0.5i, -0.5 - 0.5i];

% Define the zero locations
z2 = [0.5 0.5]; 

% Generate and plot the impulse response
b2 = 1;  % Coefficient for the zero
figure;
plot(filter(b2 * poly(z2), poly(p2), [1 zeros(1, 100)]));
title('Impulse Response - Quadrant II and Quadrant III Poles');
xlabel('Time (n)');
ylabel('Amplitude');

%% ONE pole in quadrant I
% Define the pole location
p3 = 0.5 + 0.5i;

% Define the zero locations
z3 = [0.5]; 

% Generate and plot the impulse response
b3 = 1;  % Coefficient for the zero
figure;
plot(real(filter(b3 * poly(z3), poly(p3), [1 zeros(1, 100)])));
title('Impulse Response - Quadrant I Pole');
xlabel('Time (n)');
ylabel('Amplitude');

%% TWO conjugate poles outside of the unit circle
% Define the pole locations outside of the unit circle
p4 = [1.2 + 0.5i, 1.2 - 0.5i];

% Define the zero locations 
z4 = [0.5 0.5];

% Generate and plot the impulse response
b4 = 1;  % Coefficient for the zero
figure;
plot(filter(b4 * poly(z4), poly(p4), [1 zeros(1, 100)]));
title('Impulse Response - Conjugate Poles Outside of Unit Circle');
xlabel('Time (n)');
ylabel('Amplitude');
