function [x_KKT, NLL_value_KKT] = GRVMLE(y, H, sigma_e, sigma_n)
%% SVD on H
[U,Sigma,V] = svd(H);                   % svd on H
[m,n] = size(H);
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


%% function g(nu) Eq.50
KKT_nu = @(nu)nu_picture(nu,y1_tilde,y2_tilde,sigma_h,sigma_e,sigma_n,m); 
%% calculate S
if r < m
    sum_w_i = KKT_nu(0)
else 
    sum_w_i = -inf;
end
% (1 - (sigma_e * sum((y1_tilde .^ 2) ./ (sigma_h .^ 2)) + sigma_n) * ...
%     m / norm(y2_tilde) ^2) / sigma_e

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
            % NMSE = norm(x_KKT - x)^2 / norm(x)^2
        else                % case 2 corresponding to Subsubsection IIIB-1)
            [nu_root,n_iter] = bisect(KKT_nu, -min(sigma_h) ^ 2 / 2 / (sigma_e+eps), m/2, 1e-3);   
            x_tilde_1R = y1_tilde .* sigma_h ./ (sigma_h .^ 2 + 2 * nu_root *sigma_e);  
            x_KKT_tilde = [x_tilde_1R; zeros(n-r,1)];
            x_KKT = V * x_KKT_tilde;   
            NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
            % NMSE = norm(x_KKT - x)^2 / norm(x)^2
        end
    elseif sum_w_i == 0     % case 3 corresponding to Subsubsection IIIA-2)
        x_tilde_1N = y1_tilde ./ sigma_h;
        x_KKT = V * x_tilde_1N; 
        NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
        % NMSE = norm(x_KKT - x)^2 / norm(x)^2
    elseif sum_w_i > 0      % case 4 corresponding to Subsubsection IIIB-2) nu<0
        [nu_root, n_iter] = bisect(KKT_nu, -min(sigma_h) ^ 2 / 2 / sigma_e, 0, 1e-1);
        x_tilde_1N = y1_tilde .* sigma_h./ (sigma_h .^ 2 + 2 * nu_root *sigma_e);
        x_KKT = V * x_tilde_1N; 
        NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
        % NMSE = norm(x_KKT - x)^2 / norm(x)^2
    else                    % case 4 corresponding to Subsubsection IIIB-2) nu>0
        [nu_root, n_iter] = bisect(KKT_nu, 0, m/2, 1e-3);
        x_tilde_1N = y1_tilde .* sigma_h./ (sigma_h .^ 2 + 2 * nu_root *sigma_e);
        x_KKT = V * x_tilde_1N; 
        NLL_value_KKT = Negloglikelihood(x_KKT, y, H, sigma_e, sigma_n, m)
        % NMSE = norm(x_KKT - x)^2 / norm(x)^2
    end
end