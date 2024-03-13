%%% Solve model

function [a_i,b_i,tradesh,aconverge,bconverge,xtic] = solveHab(param,observe,dist,nobs)

global alpha theta epsilon LL;

xtic = tic();

% param=[alpha theta epsilon];
alpha=param(1);
theta=param(2);
epsilon=param(3);

% observe(:,1)=L; observe(:,2)=w; observe(:,3)=H; 
L=observe(:,1); w=observe(:,2); H=observe(:,3); 

% convergence indicator;
aconverge=0;
bconverge=0;

% Initializations;
a_i=double(ones(nobs,1));
b_i=double(ones(nobs,1));

% trade costs;
dd=double(dist.^(-theta));

display('>>>> Start productivity and amenities Convergence <<<<');

% ****************************************;
% **** START OUTER LOOP FOR AMENITIES ****;
% ****************************************;

%display('>>>> Solving region amenities <<<<');
xx=1;
while xx<2000;

% *******************************************;
% **** START INNER LOOP FOR PRODUCTIVITY ****;
% *******************************************;

%display('>>>> Solving region productivities <<<<');
x=1;
while x<2000;
   
% Trade share;
pwmat=double((L.*(a_i.^theta).*(w.^(-theta)))*ones(1,nobs));
nummat=double(dd.*pwmat);
denom=double(sum(nummat));
denommat=double(ones(nobs,1)*denom);
tradesh=double(nummat./denommat);
% test;
test=sum(tradesh);
mntest=mean(test);
% Income equals expenditure;
income=double(w.*L);
expend=double(tradesh*income);    

% Convergence criterion;
income_r=round(income.*(10.^6));
expend_r=round(expend.*(10.^6));

[x max(abs(income_r-expend_r))];

% Update loop;
if income_r==expend_r;
    %display('>>>> Productivity Convergence Achieved <<<<');
    x=10000;
    aconverge=1;
else;  
a_e=double(a_i.*((income./expend).^(1./theta)));    
a_i=(0.25.*a_e)+(0.75.*a_i);
% Normalization;
a_i=a_i./geomean(a_i);
aconverge=0;
x=x+1;
end;

end;

% *******************************************;
% **** END INNER LOOP FOR PRODUCTIVITIES ****;
% *******************************************;

% domestic trade share;
dtradesh=diag(tradesh);

% population;
num=b_i.*(a_i.^(alpha.*epsilon)).*(H.^(epsilon.*(1-alpha))).*(dtradesh.^(-alpha.*epsilon./theta)).*(L.^(-((epsilon.*(1-alpha))-(alpha.*epsilon./theta))));
L_e=zeros(nobs,1);
L_e=num./sum(num);
L_e=L_e.*LL;

% Convergence criterion;
L_r=round(L.*(10.^6));
L_e_r=round(L_e.*(10.^6));

[xx max(abs(L_e_r-L_r))];

% Update loop;
if L_r==L_e_r;
    %display('>>>> Population Convergence Achieved <<<<');
    xx=10000;
    bconverge=1;
else;  
b_e=double(b_i.*(L./L_e));    
b_i=(0.25.*b_e)+(0.75.*b_i);
% Normalization;
b_i=b_i./geomean(b_i);
bconverge=0;
xx=xx+1;
end;

end;

xtic=toc(xtic);
xtic

% **************************************;
% **** END OUTER LOOP FOR AMENITIES ****;
% **************************************;
