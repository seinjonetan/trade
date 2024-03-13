%%% Solve model

function [w_i,P,r,tradesh,dtradesh,wconverge,xtic] = solveImmobileCtyOpen(param,fund,L,dwght,dist,nobs)

global alpha sigma theta epsilon LL;

xtic = tic();

% param=[alpha theta epsilon];
alpha=param(1);
theta=param(2);
epsilon=param(3);

% fund(:,1)=a; fund(:,2)=b; fund(:,3)=H; fund(:,4)=Iwest; fund(:,5)=Ieast;
a=fund(:,1); b=fund(:,2); H=fund(:,3); Iwest=fund(:,4); Ieast=fund(:,5); 

% convergence indicator;
wconverge=0;

% Initializations;
w_i=double(ones(nobs,1));

% trade costs;
dd=double(dist.^(-theta));
dd=dwght.*dd;

display('>>>> Start Wage and Population Convergence <<<<');

% ************************************;
% **** START INNER LOOP FOR WAGES ****;
% ************************************;

%display('>>>> Solving region wages <<<<');
x=1;
while x<2000;
   
% Trade share;
pwmat=double((a.*(w_i.^(-theta)))*ones(1,nobs));
nummat=double(dd.*pwmat);
denom=double(sum(nummat));
denommat=double(ones(nobs,1)*denom);
tradesh=double(nummat./denommat);
% test;
test=sum(tradesh);
mntest=mean(test);
% Income equals expenditure;
income=double(w_i.*L);
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
w_i=w_i./geomean(w_i(Iwest==1));
wconverge=0;
x=x+1;
end;

end;

% **********************************;
% **** END INNER LOOP FOR WAGES ****;
% **********************************;

% domestic trade share;
dtradesh=diag(tradesh);

% gamma function;
gammaf=gamma((theta+1-sigma)./theta);

% price index;
P=((gammaf.^-theta).*a.*(w_i.^-theta)./dtradesh).^(-1./theta);

% Land price;
r=((1-alpha)./alpha).*((w_i.*L)./H);

xtic=toc(xtic);
xtic

