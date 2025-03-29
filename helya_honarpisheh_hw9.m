close all;
clear all;

load lfp
x = lfp;

figure
subplot(2,2,1)
stackedplot(x)
title('Original Data')
% plot(x(:,5), x(:,6), '.')
% title('Original Data: Electrodes 5 vs. 6')

Rxx = cov(x);
x = x - mean(x);
N = length(x);
Rxx = 1/N * x' * x;
[U, D] = eig(Rxx);
z = x * [U(:,1) U(:,2)];
subplot(2,2,2)
plot(z(:,1), z(:,2), '.')
title('Projection on First Two Principal Components')

Rxx = cov(x(1:200000,:));
[U, D] = eig(Rxx);
[~, indx] = sort(diag(D), 'descend');
U = U(:, indx);
D = diag(D);
D = D(indx);

P = eye(6) - U(:, 1:3) * U(:, 1:3)';
z = x * P;
subplot(2,2,3)
stackedplot(z)
title('Cleaned Data after Noise Projection')


% plot(z(:,5), z(:,6), '.')
% title('Cleaned Data: Electrodes 5 vs. 6')

eig_values = eig(cov(z));
% figure
% plot(eig_values)
% title('Eigenvalues of Cleaned Data Covariance Matrix')

Rnoise = cov(x(1:200000,:));
Rsignal = cov(z(400000:end,:));

[U, D] = eig(Rsignal);
eig_values_signal = eig(Rsignal);
[~, indx] = sort(diag(D), 'descend');
U = U(:, indx);
D = diag(D);
D = D(indx);
signal = z*U(:,1);
subplot(2,2,4)
plot(signal)
title('Projection of the cleaned signal on the first principal component')
% figure
% plot(eig_values_signal)
% title('Eigenvalues of Signal Covariance Matrix')


