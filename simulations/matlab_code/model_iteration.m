function [XLS] = model_iteration(w_input,T_input,alpha,theta,eta,rho)
%This function takes a wage vector and a matrix of technology parameters and generates a measure of excess labour supply,
%which is then fed into a solver to solve for the new wages. The fulle
%model solves when the input to this function delivers an excess labour
%supply matrix that is within the tolerance bounds.

%   Prices at occupation-city level
    Pck = w_input ./ T_input; 
   
    % Prices at occupation level
    Pky = (sum(Pck.^(1-alpha))).^(1/(1-alpha));
    
    % Prices at city level.
    Pc = (sum(Pky.^(1-eta))).^(1/(1-eta));

    LD_wagebill_share = (Pck.^(1-alpha)).*(Pky.^(alpha-eta))*(Pc^(eta-1));

    % Calculate lambda from wage matrix.
    lambda = sum(w_input.^(theta/(1-rho)));

    % Calculate aggregate labour supply price index:
    lambda_sum = sum(lambda.^(1-rho));

    % Calculate labour supply share numerator 
    pi_ls = (w_input.^(theta/(1-rho))).*(lambda.^(-rho))/lambda_sum;
    ss_num = w_input.*(pi_ls.^((1+theta)/theta));
    
    % Calculate labour supply denominator;
    ss_den = sum(sum(ss_num));

    LS_wagebill_share = ss_num / ss_den;

    % Re-calculate new Excess Labour Supply
    XLS = LS_wagebill_share - LD_wagebill_share;
    
end