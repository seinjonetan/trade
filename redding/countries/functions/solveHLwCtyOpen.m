%%% Solve model

function [w_i,L_i,tradesh,dtradesh,Lconverge,wconverge,xtic] = solveHLwCtyOpen(param,fund,dwght,dist,nobs)

global alpha theta epsilon LL LLwest LLeast;

xtic = tic();

% param=[alpha theta epsilon];
alpha=param(1);
theta=param(2);
epsilon=param(3);

% fund(:,1)=a; fund(:,2)=b; fund(:,3)=H; fund(:,4)=Iwest; fund(:,5)=Ieast;
a=fund(:,1); b=fund(:,2); H=fund(:,3); Iwest=fund(:,4); Ieast=fund(:,5); 

% convergence indicator;
wconverge=0;
Lconverge=0;

% Initializations;
L_i=double(ones(nobs,1)).*(LL./nobs);
w_i=double(ones(nobs,1));

% trade costs;
dd=double(dist.^(-theta));
dd=dwght.*dd;

display('>>>> Start Wage and Population Convergence <<<<');

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
pwmat=double((L_i.*(a.^theta).*(w_i.^(-theta)))*ones(1,nobs));
nummat=double(dd.*pwmat);
denom=double(sum(nummat));
denommat=double(ones(nobs,1)*denom);
tradesh=double(nummat./denommat);
% test;
test=sum(tradesh);
mntest=mean(test);
% Income equals expenditure;
income=double(w_i.*L_i);
expend=double(tradesh*income);    

% Convergence criterion;
income_r=round(income.*(10.^6));
expend_r=round(expend.*(10.^6));

[x max(abs(income_r-expend_r))];

% Update loop;
if income_r==expend_r;
    %display('>>>> Wage Convergence Achieved <<<<');
    x=10000;
    wconverge=1;
else;  
w_e=double(w_i.*(expend./income).^(1./theta));    
w_i=(0.25.*w_e)+(0.75.*w_i);
% Normalization;
% Separate countries;
w_i(Iwest==1)=w_i(Iwest==1)./geomean(w_i(Iwest==1));
wconverge=0;
x=x+1;
end;

end;

% **********************************;
% **** END INNER LOOP FOR WAGES ****;
% **********************************;

% domestic trade share;
dtradesh=diag(tradesh);

% population;
num=b.*(a.^(alpha.*epsilon)).*(H.^(epsilon.*(1-alpha))).*(dtradesh.^(-alpha.*epsilon./theta)).*(L_i.^(-((epsilon.*(1-alpha))-(alpha.*epsilon./theta))));
L_e=zeros(nobs,1);
L_e(Iwest==1)=num(Iwest==1)./sum(num(Iwest==1));
L_e(Ieast==1)=num(Ieast==1)./sum(num(Ieast==1));
L_e(Iwest==1)=L_e(Iwest==1).*LLwest;
L_e(Ieast==1)=L_e(Ieast==1).*LLeast;

% Convergence criterion;
L_i_r=round(L_i.*(10.^6));
L_e_r=round(L_e.*(10.^6));

[xx max(abs(L_e_r-L_i_r))];

% Update loop;
if L_i_r==L_e_r;
    %display('>>>> Population Convergence Achieved <<<<');
    xx=10000;
    Lconverge=1;
else;  
L_e=double(L_i.*((L_e./L_i).^(1./(epsilon.*(1-alpha)))));    
L_i=(0.25.*L_e)+(0.75.*L_i);
Lconverge=0;
xx=xx+1;
end;

end;

xtic=toc(xtic);
xtic

% ****************************************;
% **** END OUTER LOOP FOR POPULATIONS ****;
% ****************************************;
