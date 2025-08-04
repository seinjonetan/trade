function [w] = equilibrium(w_init,T_input,alpha,theta,eta,rho,tolfun,tolx)
    % T_input_gpu = gpuArray(T_input);
    % fcf = @(w) gather(model_iteration(gpuArray(w), T_input_gpu, alpha, theta, eta, rho));
    fcf = @(w) model_iteration(w, T_input, alpha, theta, eta, rho);
    
    options = optimset('Display','off','TolFun',tolfun,'TolX',tolx);
    [w] = fsolve(fcf, double(w_init), options);
end