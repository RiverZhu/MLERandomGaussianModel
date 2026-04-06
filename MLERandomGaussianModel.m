% Verify results for solving KKT conditions in multiplicative noise environments
% Ruohai Guo
% e-mail: rhguo@zju.edu.cn
% Date: 2025. 05. 19

clc;clear;close all;
%% uncommented for loading data x, y, H
% load('Example1.mat')
% load('Example2.mat')
% load('Example3.mat')
% load('Example4.mat')
% load('Example5.mat')
% n = length(x);         % length of x
% m = length(y);         % length of y

% rng(1)
%% generation
%% uncommented for randomly generating x, H
n = 2;            % length of x
m = 1;            % length of y
x = rand(n,1);    % creating x
r_man = min(m,n); % creating H with rank(H) = r_man
H_A = randn(m,r_man);
H_B = randn(r_man,n);
H = H_A * H_B;
% x(1) = 0;
%% parameters setting
sigma_e = 0.1;                              % covariance of E
sigma_dB = 15;                              % set SNR
sigma_n = 1/10^(sigma_dB/10);               % calculate sigma_n under the SNR
sigma_x = sigma_e * norm(x) ^ 2 + sigma_n;

%% uncommented for randomly generating E, y
E = sqrt(sigma_e) * randn(m,n);       % creating error matrix E
noise = sqrt(sigma_n) * randn(m,1);   % creating noise
y = (H+E) * x + noise;                % generating measurements y

[x_KKT, NLL_value_KKT] = GRVMLE(y, H, sigma_e, sigma_n);

%% contour plot
% x1 = -1:0.01:1;
% x2 = -1:0.01:1;
% [X1,X2] = meshgrid(x1,x2);
% for i = 1:length(x1)
%     for j = 1:length(x2)
%         x_all = [x1(i);x2(j)];
%         y_NLL(i,j) = Negloglikelihood(x_all, y, H, sigma_e, sigma_n, m);
%     end
% end
% figure;
% [C,h] = contour(X1, X2, y_NLL', 'LevelStepMode','manual','LevelStep',0.1);
% hold on
% plot(x_KKT(1),x_KKT(2),'r*','Markersize',8,'LineWidth', 1.5)
% plot(x(1),x(2),'b*','Markersize',8,'LineWidth', 1.5)
% 
% xlabel('$x_1$','interpreter', 'latex','FontSize',14)
% ylabel('$x_2$','interpreter', 'latex','FontSize',14)
