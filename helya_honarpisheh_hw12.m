close all;
clear all;

load LATE_example.mat

True_ACE = 3;

% Estimate ITT
itt_estimate_simulated = mean(Y(Z == 1)) - mean(Y(Z == 0));

% Estimate PP
pp_estimate_simulated = mean(Y(X == 1)) - mean(Y(X == 0));

% Estimate ETT 
ett_estimate_simulated = mean(Y(X == 1 & Z == 1)) - mean(Y(X == 0 & Z == 1));

model = fitlm([X U], Y);
beta_yx = model.Coefficients.Estimate(2); % causal effect of X on Y
beta_yu = model.Coefficients.Estimate(3);

m = fitlm(Z, X);
r_xz = m.Coefficients.Estimate(2);
m = fitlm(Z, Y);
r_yz = m.Coefficients.Estimate(2);


% Estimate LATE
m = fitlm([X U],Y);
r_yx_u = m.Coefficients.Estimate(2);
r_yu_x = m.Coefficients.Estimate(3);
late_estimate = r_yz / r_xz;

m = fitlm([Z U],Y);
r_yz_u = m.Coefficients.Estimate(2);
r_yu_z = m.Coefficients.Estimate(3);

m = fitlm([Z U],X);
r_xz_u = m.Coefficients.Estimate(2);
r_xu_z = m.Coefficients.Estimate(3);

% Estimate ACE Regression using perfect knowledge of U
ace_regression_estimate = r_yz_u / r_xz_u;


% Display results
fprintf('True Average Causal Effect: %.4f\n', True_ACE);
fprintf('Simulated ITT Estimate: %.4f\n', itt_estimate_simulated);
fprintf('Simulated PP Estimate: %.4f\n', pp_estimate_simulated);
fprintf('Simulated ETT Estimate: %.4f\n', ett_estimate_simulated);
fprintf('Simulated LATE Estimate: %.4f\n', late_estimate);
fprintf('ACE Regression Estimate: %.4f\n', ace_regression_estimate);


% Create notBoxPlots for PP:P(Y|X) and ITT:P(Y|Z)
pp_prob_data_x0 = Y(X == 0);
pp_prob_data_x1 = Y(X == 1);
itt_prob_data_z0 = Y(Z == 0);
itt_prob_data_z1 = Y(Z == 1);

figure;
% Per Protocol (PP) Probability Box Plot for X=0
subplot(1, 2, 1);
h = notBoxPlot(pp_prob_data_x0,0);
h(1).data.Marker = '.';
h(2).data.Marker = '.';
hold on;
plot([-0.5, 0.5], [0,0], 'Color', [0.8, 0, 0.4]); 
hold on;
% Per Protocol (PP) Probability Box Plot for X=1
h = notBoxPlot(pp_prob_data_x1,1);
h(1).data.Marker = '.';
h(2).data.Marker = '.';
hold on;
plot([0.5, 1.5], [3,3], 'Color', [0.8, 0, 0.4]);
legend('','','','','True Effect', 'Location', 'SouthEast');
title('PP: P(Y|X)');
xlabel('X');ylabel('Y');ylim([-4 8]);
hold off;
% Intention to Treat (ITT) Probability Box Plot for Z=0
subplot(1, 2, 2);
h = notBoxPlot(itt_prob_data_z0,0);
h(1).data.Marker = '.';
h(2).data.Marker = '.';
hold on;
plot([-0.5, 0.5], [0,0], 'Color', [0.8, 0, 0.4]);
hold on;
% Intention to Treat (ITT) Probability Box Plot for Z=1
h = notBoxPlot(itt_prob_data_z1,1);
h(1).data.Marker = '.';
h(2).data.Marker = '.';
hold on;
plot([0.5, 1.5], [3,3], 'Color', [0.8, 0, 0.4]);
legend('','','','','True Effect', 'Location', 'SouthEast');
title('ITT: P(Y|Z)');
xlabel('Z');ylabel('Y');ylim([-4 8]);
hold off;

