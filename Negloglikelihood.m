% Negloglikelihood Eq.5
% Ruohai Guo
% e-mail: rhguo@zju.edu.cn
% Date: 2025. 05. 19

function NLL = Negloglikelihood(x, y, H, sigma_e, sigma_n, m)
    NLL = norm(y - H * x) ^ 2 / (sigma_e * norm(x) ^ 2 + sigma_n)...
        + m * log(sigma_e * norm(x) ^ 2 + sigma_n);
end