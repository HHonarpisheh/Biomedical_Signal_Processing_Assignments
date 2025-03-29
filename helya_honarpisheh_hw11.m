clear all;
close all;

load cause_or_confound.mat

U_z = Z(:,1); % confounder variation
U_y = Y(:,1); % outcome variation
U_x = X(:,1); % treatment variation

% fit the data with linear model while “adjusting”, “controlling” for Z
m = fitlm([U_x U_z],U_y);
beta_yx = m.Coefficients.Estimate(2)% causal effect
beta_yz = m.Coefficients.Estimate(3);
m = fitlm([U_z],U_x);
beta_xz = m.Coefficients.Estimate(2);

% Fit the naive associations
m = fitlm(U_x, U_y);
r_xy = m.Coefficients.Estimate(2)% total association
m = fitlm(U_x, U_z);
r_xz = m.Coefficients.Estimate(2);% association between X and Z
m = fitlm(U_z, U_y);
r_zy = m.Coefficients.Estimate(2);% association between Z and Y

figure
subplot(1,3,1)
plot(U_x,U_y,'.')
title( ['r_x_y' num2str( r_xy )] )
xlabel('X');ylabel('Y');
xlim([-5 5]);ylim([-4 4]);
subplot(1,3,2)
plot(U_x,U_z,'.')
title( ['r_x_z' num2str( r_xz )] )
xlabel('X');ylabel('Z');
xlim([-5 5]);ylim([-4 2]);
%lsline
subplot(1,3,3)
plot(U_z,U_y,'.')
title( ['r_z_y' num2str( r_zy )] )
xlabel('Z');ylabel('Y');
xlim([-4 2]);ylim([-4 4]);

disp('The estimated causal effect of first dataset is negative,indicating a negative causal relationship. However, the total association is almost zero, suggesting a weak overall association between X and Y. (The negative causal effect implies that an increase in X is associated with a decrease in Y.)')

%% For second dataset
U_z = Z(:,2); % confounder variation
U_y = Y(:,2); % outcome variation
U_x = X(:,2); % treatment variation

% fit the data with linear model while “adjusting”, “controlling” for Z
m = fitlm([U_x U_z],U_y);
beta_yx = m.Coefficients.Estimate(2)% causal effect
beta_yz = m.Coefficients.Estimate(3);
m = fitlm([U_z],U_x);
beta_xz = m.Coefficients.Estimate(2);

% Fit the naive associations
m = fitlm(U_x, U_y);
r_xy = m.Coefficients.Estimate(2)% total association
m = fitlm(U_x, U_z);
r_xz = m.Coefficients.Estimate(2);% association between X and Z
m = fitlm(U_z, U_y);
r_zy = m.Coefficients.Estimate(2);% association between Z and Y

figure
subplot(1,3,1)
plot(U_x,U_y,'.')
title( ['r_x_y' num2str( r_xy )] )
xlabel('X');ylabel('Y');
xlim([-5 5]);ylim([-4 4]);
subplot(1,3,2)
plot(U_x,U_z,'.')
title( ['r_x_z' num2str( r_xz )] )
xlabel('X');ylabel('Z');
xlim([-5 5]);ylim([-3 3]);
%lsline
subplot(1,3,3)
plot(U_z,U_y,'.')
title( ['r_z_y' num2str( r_zy )] )
xlabel('Z');ylabel('Y');
xlim([-2.5 2.5]);ylim([-4 4]);

disp('The estimated causal effect of second dataset is negative,indicating a negative causal relationship. However, the total association is higher, indicating a stronger overall association between X and Y. ')
