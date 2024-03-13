%%% Solve model

function [Cw_i,CL_i,Ctradesh,dCtradesh,CLconverge,Cwconverge,xtic] = cftualLw(param,fund,L,w,tradesh,dist,Cdist,nobs)

global alpha sigma theta epsilon LL;

xtic = tic();

% param=[alpha theta epsilon];
alpha=param(1);
theta=param(2);
epsilon=param(3);

% fund(:,1)=a; fund(:,2)=b; fund(:,3)=H; 
a=fund(:,1); b=fund(:,2); H=fund(:,3); 

% trade costs;
dd=double(dist.^(-theta));

% Counterfactual trade costs;
Cdd=double(Cdist.^(-theta));

ddhat=Cdd./dd;

% gdp;
gdp=w.*L;

% domestic trade share;
dtradesh=diag(tradesh);

% population share;
lambda=L./LL;

% convergence indicator;
Cwconverge=0;
CLconverge=0;

% Initializations;
CL_i=L;
Cw_i=w;

% trade costs;
dd=double(dist.^(-theta));

display('>>>> Start Counterfactual Wage and Population Convergence <<<<');

% *****************************************;
% **** START OUTER LOOP FOR POPULATION ****;
% *****************************************;

%display('>>>> Solving region populations <<<<');
xx=1;
while xx<2000;
    
% ************************************;
% **** START INNER LOOP FOR WAGES ****;
% ************************************;

%display('>>>> Solving region wages <<<<');
x=1;
while x<2000;
    
% Trade share;
pwmat=double(((Cw_i./w).^(-theta))*ones(1,nobs));
nummat=double(tradesh.*ddhat.*pwmat);
denom=double(sum(nummat));
denommat=double(ones(nobs,1)*denom);
Ctradesh=double(nummat./denommat);
% test;
test=sum(Ctradesh);
mntest=mean(test);
% Income equals expenditure;
income=double((Cw_i./w).*(CL_i./L).*gdp);
expend=double(Ctradesh*income);    

% Convergence criterion;
income_r=round(income.*(10.^6));
expend_r=round(expend.*(10.^6));

%display('x');
%[x max(abs(income_r-expend_r))]

% Update loop;
if income_r==expend_r;
    %display('>>>> Wage Convergence Achieved <<<<');
    x=10000;
    Cwconverge=1;
else;  
Cw_e=double(Cw_i.*(expend./income).^(1./theta));    
Cw_i=(0.25.*Cw_e)+(0.75.*Cw_i);
% Normalization;
Cw_i=Cw_i./geomean(Cw_i);
Cwconverge=0;
x=x+1;
end;

end;

% **********************************;
% **** END INNER LOOP FOR WAGES ****;
% **********************************;

% domestic trade share;
dCtradesh=diag(Ctradesh);

% population;
num=((dCtradesh./dtradesh).^(-alpha.*epsilon./theta)).*((CL_i./L).^(-epsilon.*(1-alpha))).*lambda;
CL_e=zeros(nobs,1);
CL_e=num./sum(num);
CL_e=CL_e.*LL;

% Convergence criterion;
CL_i_r=round(CL_i.*(10.^6));
CL_e_r=round(CL_e.*(10.^6));

%display('xx');
%[xx max(abs(CL_e_r-CL_i_r))]

% Update loop;
if CL_i_r==CL_e_r;
    %display('>>>> Population Convergence Achieved <<<<');
    xx=10000;
    CLconverge=1;
else;  
CL_e=double(CL_i.*(CL_e./CL_i).^(1./(epsilon.*(1-alpha))));    
CL_i=(0.25.*CL_e)+(0.75.*CL_i);
CLconverge=0;
xx=xx+1;
end;

end;

xtic=toc(xtic);
xtic

% ****************************************;
% **** END OUTER LOOP FOR POPULATIONS ****;
% ****************************************;
