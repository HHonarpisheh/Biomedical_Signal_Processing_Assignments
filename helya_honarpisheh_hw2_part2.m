% myhomework #2 part #2
clear all;
close all;

load ('gamma.mat','fs','osc');
N = length(osc);   % Total number of data points
y = osc(:,1);
x = osc(:,2);
t = 0:1/fs:(length(osc)-1)/fs;
Q_range = 1:150; 

% Initialize variables
train_errors = zeros(length(Q_range), 5);
test_errors = zeros(length(Q_range), 5);

% Perform 5-fold cross-validation
num_folds = 5;
fold_size = length(x) / num_folds;

for fold = 1:num_folds
    % Split data 
    test_indices = (fold - 1) * fold_size + 1 : fold * fold_size;
    train_indices = setdiff(1:length(x), test_indices);
    x_train = x(train_indices);
    y_train = y(train_indices);
    x_test = x(test_indices);
    y_test = y(test_indices);

    for i = 1:length(Q_range)
        Q = Q_range(i);
        [b, ~, yest] = myfit(y_train, x_train, Q); % Train the model on the training data
        train_errors(i, fold) = mean((y_train - yest).^2); % Compute training error (mean squared error)
        yest_test = toeplitz(x_test, [x_test(1) zeros(1, Q)]) * b; % Predict on test data
        test_errors(i, fold) = mean((y_test - yest_test).^2); % Compute test error (mean squared error)
    end
end

% Calculate the average training and test errors for each model order
avg_train_errors = mean(train_errors,2);
avg_test_errors = mean(test_errors,2);

% Find the optimal model order (Q) that minimizes the average test error
[min_test_error, optimal_Q_idx] = min(avg_test_errors);
optimal_Q = Q_range(optimal_Q_idx);

figure;
subplot(2, 1, 1);
plot(t,y,"blue"); xlim([0 2.5]); xlabel('Time (s)');
% Train the model on the full dataset with the optimal model order
[b_optimal, ~, yest_optimal] = myfit(y, x, optimal_Q);
hold on;
plot(t,yest_optimal,"green");
hold on;
residual = y - yest_optimal;
plot(t,residual,"red");
hold off;
legend('original','estimate','residual');
legend boxoff;
subplot(2, 2, 3);
plot(Q_range,avg_train_errors, "blue");
hold on;
plot(Q_range,avg_test_errors, "green");
xlabel('Model Order (Q)');
ylabel('Error');
hold on;
plot(optimal_Q,min_test_error, 'o','Color',[0 0 0]);
legend('train', 'test','optimal model');
hold off;
legend boxoff;
xlabel('Model order Q'); ylabel('Error');
subplot(2, 2, 4);
stem([1:optimal_Q]/fs,b(1:optimal_Q),"blue"); yticks([-0.1:0.02:0.04]);
title('Optimal MA filter model');
xlabel('Time (s)');

function [b,e,yest] = myfit(y,x,Q)
% this function will predict y from x using Q delays
b = toeplitz(x(Q+1:end),x(Q+1:-1:1))\y(Q+1:end);
yest =  toeplitz(x,[x(1) zeros(1,Q)]) * b;
e = y - yest;
end