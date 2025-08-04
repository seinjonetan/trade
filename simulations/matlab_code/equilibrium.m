function [w] = equilibrium(w_init,T_input,alpha,theta,eta,rho,tolfun,tolx)
%This function takes a counterfactual matrix of T_ck matrix and an initial guess of wages and solves for
%the new wages associated with this counterfactual.

    %Set optimization options, such as tolerance levels.
    options = optimset('Display','off','TolFun',tolfun,'TolX',tolx);
    
    %Call model_iteration function to solve each iteration of a wage guess.
    fcf = @(w)model_iteration(w,T_input,alpha,theta,eta,rho);
    [w] = fsolve(fcf,w_init,options);

end