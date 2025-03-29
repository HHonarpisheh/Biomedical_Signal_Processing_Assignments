close all;
clear all;
%% slide 14
load('eeg.mat');
y=[x(5,:);x(6,:);x(7,:)];
x=x(1,:);
Rxx=x*x'/length(x);
Ryx=y*x'/length(x);
A=Ryx*inv(Rxx);
n=y-A*x;
subplot(2,2,1)
plot(y(2,:)'+1000,'b');hold on;title('y, x');
plot(y(1,:)'+3000,'g');hold on;yticks([-1000 1000 3000 5000]);yticklabels({'x', 'y1', 'y2', 'y3'});
plot(y(3,:)'+5000,'r');hold on;
plot(x'-1000,'c');hold off;
subplot(2,2,2)
plot(n(2,:)'+1000,'b');hold on;title('n = y - Ax');
plot(n(1,:)'+3000,'g');hold on;yticks([1000 3000 5000]);yticklabels({'n1', 'n2', 'n3'});
plot(n(3,:)'+5000,'r');hold off;
subplot(2,2,3)
scatter(x,y(2,:)','b','.');hold on;xlim([-2000 500]);xlabel('x');
scatter(x,y(1,:)','g','.');hold on;ylim([-500 2000]);ylabel('y');
scatter(x,y(3,:)','r','.');hold off;
subplot(2,2,4)
scatter(x,n(2,:)','b','.');hold on;xlim([-2000 500]);xlabel('x');
scatter(x,n(1,:)','g','.');hold on;ylim([-400 800]);ylabel('n = y - Ax');
scatter(x,n(3,:)','r','.');hold off;
%Rnx=n*x';
%% slide 15
mu = [0;0];
Sigma = [1 0.5;0.5 0.4];
x1 = -2:0.1:2;
x2 = -2:0.1:2; 
Dt = det(Sigma);
I = inv(Sigma);

pdf_values = zeros(length(x1), length(x2));

for i = 1:length(x1)
    for j = 1:length(x2)
        x_test = [x1(i); x2(j)];
        x_minus_mu = x_test - mu;
        q = 1/ sqrt (((2*pi) ^2) *Dt);
        b = dot (x_minus_mu'*I,x_minus_mu);
        pdf_values(i, j) = q*exp ((-1/2) *b);
    end
end

figure;
subplot(1,2,1)
mesh(x1, x2, pdf_values);
Z=randn(2,1000);
[U,D]=eig(Sigma);
W=U*sqrt(D);
Y=W*Z;
subplot(1,2,2)
contour(x1, x2, pdf_values, 'LineWidth', 1.5);hold on;xlim([-2 2]);ylim([-2 2]);
plot(Y(1,:),Y(2,:),'.');hold off;

%% slide 18
%Rxx = cov(Z');
Rxx=Sigma;
[U,D] = eig(Rxx);
%[~,indx] = sort(diag(D),'descend'); U=U(:,indx);D=diag(D);D=D(indx);
W = U*sqrt(D);
x = inv(W)*Z;

figure;
scatter(x(1, :), x(2, :),'.'); 
xlim([-4 4]);
ylim([-2 2.5]);
hold on;
plot([0;-U(1,1)],[0;-U(1,2)],'black');hold on;
plot([0;-U(2,1)],[0;-U(2,2)],'black');hold off;
%%