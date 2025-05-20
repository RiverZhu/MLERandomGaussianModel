% Verify results for solving KKT conditions in multiplicative noise environments
% Ruohai Guo
% e-mail: rhguo@zju.edu.cn
% Date: 2025. 05. 19

clc;clear;close all;
%% uncommented for loading data x, y, H
load('Example1.mat')
% load('Example2.mat')
% load('Example3.mat')
% load('Example4.mat')
% load('Example5.mat')
n = length(x);         % length of x
m = length(y);         % length of y

% rng(1)
%% generation
%% uncommented for randomly generating x, H
% n = 2;            % length of x
% m = 1;            % length of y
% x = rand(n,1);    % creating x
% r_man = min(m,n); % creating H with rank(H) = r_man
% H_A = randn(m,r_man);
% H_B = randn(r_man,n);
% H = H_A * H_B;

%% parameters setting
sigma_e = 0.1;                              % covariance of E
sigma_dB = 15;                              % set SNR
sigma_n = 1/10^(sigma_dB/10);               % calculate sigma_n under the SNR
sigma_x = sigma_e * norm(x) ^ 2 + sigma_n;

%% uncommented for randomly generating E, y
% E = sqrt(sigma_e) * randn(m,n);       % creating error matrix E
% noise = sqrt(sigma_n) * randn(m,1);   % creating noise
% y = (H+E) * x + noise;                % generating measurements y

%% SVD on H
[U,Sigma,V] = svd(H);                   % svd on H
[a,b] = find(Sigma < 1e-10 & Sigma > 0);
Sigma(a,b) = 0;
r = size(nonzeros(Sigma),1);            % get rank(H)
sigma_h = nonzeros(Sigma);
sigma_h_reduced = sigma_h(1:r);
sigma_h_reduced_square = sigma_h_reduced.^2;

%% generate y_tilde
y_tilde = U' * y;               % y_tilde
y1_tilde = y_tilde(1:r);        % get y1_tilde
y2_tilde = y_tilde(r+1:end);    % get y2_tilde

%% calculate S
sum_w_i = (1 - (sigma_e * sum((y1_tilde .^ 2) ./ (sigma_h .^ 2)) + sigma_n) * ...
    m / norm(y2_tilde) ^2) / sigma_e

%% function g(nu) Eq.50
KKT_nu = @(nu)nu_picture(nu,y1_tilde,y2_tilde,sigma_h,sigma_e,sigma_n,m); 

%% main body of the algorithm
if r <= n-1
    if sum_w_i >= 0     % case 1 corresponding to Subsubsection IIIA-1)
        z_star = m / norm(y2_tilde) ^2;
        w_RN(1:n-r) = sum_w_i / (n-r);                              % average for each w_i 
        x_tilde_1R = y1_tilde ./ sigma_h;
        x_tilde_RN = sqrt(w_RN / z_star)' * (randi([0,1]) * 2 - 1); % change the sign if one wants to get anther result
        x_KKT_tilde = [x_tilde_1R; x_tilde_RN];
        x_KKT = V * x_KKT_tilde;
        NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
        NMSE = norm(x_KKT - x)^2 / norm(x)^2
    else                % case 2 corresponding to Subsubsection IIIB-1)
        [nu_root,n_iter] = bisect(KKT_nu, -min(sigma_h) ^ 2 / 2 / (sigma_e+eps), m/2, 1e-3);   
        x_tilde_1R = y1_tilde .* sigma_h ./ (sigma_h .^ 2 + 2 * nu_root *sigma_e);  
        x_KKT_tilde = [x_tilde_1R; zeros(n-r,1)];
        x_KKT = V * x_KKT_tilde;   
        NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
        NMSE = norm(x_KKT - x)^2 / norm(x)^2
    end
elseif sum_w_i == 0     % case 3 corresponding to Subsubsection IIIA-2)
    x_tilde_1N = y1_tilde ./ sigma_h;
    x_KKT = V * x_tilde_1N; 
    NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
    NMSE = norm(x_KKT - x)^2 / norm(x)^2
elseif sum_w_i > 0      % case 4 corresponding to Subsubsection IIIB-2) nu<0
    [nu_root, n_iter] = bisect(KKT_nu, -min(sigma_h) ^ 2 / 2 / (sigma_e+eps), 0, 1e-3);
    x_tilde_1N = y1_tilde .* sigma_h./ (sigma_h .^ 2 + 2 * nu_root *sigma_e);
    x_KKT = V * x_tilde_1N; 
    NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
    NMSE = norm(x_KKT - x)^2 / norm(x)^2
else                    % case 4 corresponding to Subsubsection IIIB-2) nu>0
    [nu_root, n_iter] = bisect(KKT_nu, 0, m/2, 1e-3);
    x_tilde_1N = y1_tilde .* sigma_h./ (sigma_h .^ 2 + 2 * nu_root *sigma_e);
    x_KKT = V * x_tilde_1N; 
    NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
    NMSE = norm(x_KKT - x)^2 / norm(x)^2
end
x_LS = pinv(H) * y;     % LS result
NLL_value_LS = Negloglikelihood(x_LS, y, H, sigma_e, sigma_n, m)



%% contour plot
x1 = -1:0.001:1;
x2 = -1:0.001:1;
[X1,X2] = meshgrid(x1,x2);
for i = 1:length(x1)
    for j = 1:length(x2)
        x_all = [x1(i);x2(j)];
        y_NLL(i,j) = Negloglikelihood(x_all, y, H, sigma_e, sigma_n, m);
    end
end
figure;
[C,h] = contour(X1, X2, y_NLL', 'LevelStepMode','manual','LevelStep',0.01);
hold on
plot(x_KKT(1),x_KKT(2),'r*','Markersize',8,'LineWidth', 1.5)
plot(x(1),x(2),'b*','Markersize',8,'LineWidth', 1.5)

xlabel('$x_1$','interpreter', 'latex','FontSize',14)
ylabel('$x_2$','interpreter', 'latex','FontSize',14)
