
function [f] = thetaepsopt(BETA)

global alpha sigma theta epsilon LL nobs W;
global observe Cobserve dist0 dist1;

theta=BETA(1);
epsilon=BETA(2);

estparam=[alpha theta epsilon];

% Solve for baseline region productivities and amenities;
[a_i,b_i,abtradesh,aconverge,bconverge,xtic]=solveab(estparam,observe,dist0,nobs);
display('>>>> Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[aconverge bconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Solve for counterfactual region productivities and amenities  
[Ca_i,Cb_i,Cabtradesh,Caconverge,Cbconverge,Cxtic]=solveab(estparam,Cobserve,dist1,nobs);
display('>>>> Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[aconverge bconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Change in fundamentals;
da_i=log(Ca_i./a_i);
db_i=log(Cb_i./b_i);

% Objective function;

ft=zeros(2,1);
ft(1)=da_i'*da_i;
ft(2)=db_i'*db_i;

f=ft' * W * ft;

