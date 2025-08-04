function [T] = technology(pick,yck,alpha,theta,eta,rho)
%This function takes as input an initial matrix of choice shares at the
%city-occupation level (pick) and a matrix of initial labour income at the city-occupation level (yck). 
%This function generates a matrix of technology parameters of the same dimensions.
%Note that the largest city-occupation pair, in terms of choice shares, is
%normalized such that T_ck=1.
    
    %Generate choice shares at the occupation level.
    omegak = sum(pick);
    %Generate income at the occupation level.
    yk = sum(yck); 
   
    % Identify numeraire city-occupation pair
    maximum = max(max(pick));
    [x,y]=find(pick==maximum);
    pi_max = pick(x,y);
    y_max = yck(x,y);
    omega_max = omegak(y);
    yk_max = yk(y);

    % Generate the technology matrix;
    T = (((pick/pi_max).^(1-rho).*(omegak/omega_max).^(rho)).^(1/theta)).*(((yck/y_max).*(yk/yk_max).^((alpha-eta)/(eta-1)))).^(1/(alpha-1));

end