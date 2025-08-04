function [pi_ck] = get_shares(w_ck, theta, rho)
    w_ck = w_ck.^(theta/ (1 - rho));
    lambda_k = sum(w_ck, 1);
    pi_ck = w_ck .* (lambda_k .^ (- rho) ./ sum(lambda_k .^ (- rho)));
end