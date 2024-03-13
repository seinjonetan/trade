
% Monte Carlo for Quantitative Spatial Model;
% Constant and increasing returns to scale model;
% Regions specification;

% SJR, November, 2015;

% *********************; 
% **** Choose User ****; 
% *********************;

clear all;
clc;
colormap hot;
close all;

user=1;

% User 1 : steve windows desktop; 
% User 2 : steve mac laptop dropbox;

if user==1;
cd C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\regions;
addpath C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\regions\functions;
end;

format shortG;

global alpha sigma theta epsilon LL nobs;
global observe Cobserve dist0 dist1;

% ************************; 
% **** Initialization ****;
% ************************;

% Set default random number stream;
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% *************************; 
% **** Distance matrix ****;
% *************************;

N=11;
NN=N*N;
% Other latitude-longitude grid;
ltd=linspace(0,4,N)';
lgd=linspace(0,4,N);

% Transport weights;
tt0=7.9; 
tt1=1;

tau0=zeros(N);
tau0(:,:)=tt0;

tau1=zeros(N);
tau1(:,:)=tt0;
tau1(:,6)=tt1;
tau1(6,:)=tt1;

dist0=zeros(NN);
dist1=zeros(NN);

for z=1:NN;
    seed=false(N,N);
    seed(z)=true;
    temp=graydist(tau0,seed,'quasi-euclidean');
    dist0(z,:)=reshape(temp,1,NN); 
    temp=graydist(tau1,seed,'quasi-euclidean');
    dist1(z,:)=reshape(temp,1,NN); 
end;
dist0(eye(size(dist0))==1)=1;
dist1(eye(size(dist0))==1)=1;
rdist=dist1./dist0;

% Idist0=mat2gray(dist0);
% Idist1=mat2gray(dist1);
% Irdist=mat2gray(rdist);
% figure(1);
% imshow(Idist0);
% figure(2);
% imshow(Idist1);
% figure(3);
% imshow(Irdist);

% Define treatment;
treat=zeros(N,N);
treat(:,6)=1;
treat(6,:)=1;
treat=reshape(treat,NN,1);

% Define controls;
X=[ones(size(treat)) treat];

% Trade costs are a power function of effective distance;

dist0=dist0.^0.33;
dist1=dist1.^0.33;

% **************************;
% **** Parameterization ****;
% **************************;

% Share of goods in consumption expenditure (1-housing share);
alpha=0.75; 
% Elasticity of substitution;
sigma=4;

% ***********************;
% **** Random shocks ****;
% ***********************;

a=normrnd(0,1,NN,1);
a=exp(a);
a=a./geomean(a);

b=normrnd(0,1,NN,1);
b=exp(b);
b=b./geomean(b);

display('Summary statistics a');
[mean(a) std(a) max(a) min(a)]

display('Summary statistics a');
[mean(b) std(b) max(b) min(b)]

% **************************;
% **** Other Parameters ****;
% **************************:

% Observations
nobs=NN;
% Land area;
H=100.*ones(nobs,1);
% Aggregate labor Supply;
LL=153889;          % US civilian labor force 2010 (Statistical Abstract, millions);

% ************************************;
% **** Loop over parameter values ****;
% ************************************;

% Parameter grid;
% Sigma is 4, so need theta greater than 3;
epsgrid=linspace(3.1,5.1,10);
thetagrid=linspace(3.1,5.1,10)';
K=size(epsgrid,2);
KK=K.*K;
[thetamat,epsmat]=ndgrid(thetagrid,epsgrid);
epsvec=reshape(epsmat,KK,1);
thetavec=reshape(thetamat,KK,1);

% Matrices to store;

dLmat=zeros(size(epsvec,1),1);
drwmat=zeros(size(epsvec,1),1);
dwmat=zeros(size(epsvec,1),1);
drmat=zeros(size(epsvec,1),1);
dPmat=zeros(size(epsvec,1),1);
dwelfmat=zeros(size(epsvec,1),1);
dacrmat=zeros(size(epsvec,1),1);
dmobmat=zeros(size(epsvec,1),1);

convmat=zeros(size(epsvec,1),4);

% Loop;

for p=1:1:size(epsvec,1);

display('>>>> Value of parameters are <<<<');
epsilon=epsvec(p);
theta=thetavec(p);
display('[alpha theta epsilon]');
estparam=[alpha theta epsilon];

% ****************************************;
% **** Solve for Endogenous Variables ****;
% ****************************************;

fund=zeros(nobs,5);
fund(:,1)=a; fund(:,2)=b; fund(:,3)=H;

% Solve for region populations and wages;
[w,L,tradesh,dtradesh,Lconverge,wconverge,xtic]=solveLw(estparam,fund,dist0,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[wconverge Lconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic
convmat(p,1)=wconverge;
convmat(p,2)=Lconverge;

% Price index;
[P] = pindex(estparam,fund,w,dtradesh,nobs);

% Land price;
[r]=landprice(estparam,fund,L,w,dist0,nobs);

% Welfare;
[welf]=welfare(estparam,fund,L,w,tradesh,dist0,nobs);

% Real wage;
[realwage]=realw(estparam,fund,L,w,tradesh,dist0,nobs);

% % % % *********************************;
% % % % **** Solve for Unobservables ****;
% % % % *********************************;
% % % 
% % % observe=zeros(nobs,5);
% % % observe(:,1)=L; observe(:,2)=w; observe(:,3)=H; 
% % % 
% % % % Solve for region productivities and amenities;
% % % [a_i,b_i,abtradesh,aconverge,bconverge,xtic]=solveab(estparam,observe,dist0,nobs);
% % % display('>>>> Productivity and Amenity System Converged <<<<');
% % % display('>>>> Check Productivity and Amenity Convergence <<<<');
% % % [aconverge bconverge]
% % % display('>>>> Elapsed Time in Seconds <<<<');
% % % xtic

% *******************************;
% **** CHANGE IN TRADE COSTS ****;
% *******************************;

% Solve for region populations and wages;
[Cw,CL,Ctradesh,Cdtradesh,CLconverge,Cwconverge,Cxtic]=solveLw(estparam,fund,dist1,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[Cwconverge CLconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic
convmat(p,3)=Cwconverge;
convmat(p,4)=CLconverge;

% Counterfactual price index;
[CP] = pindex(estparam,fund,Cw,Cdtradesh,nobs);

% Counterfactual land prices;
[Cr]=landprice(estparam,fund,CL,Cw,dist1,nobs);

% Counterfactual welfare;
[Cwelf]=welfare(estparam,fund,CL,Cw,Ctradesh,dist1,nobs);

% Counterfactual real wage;
[Crealwage]=realw(estparam,fund,CL,Cw,Ctradesh,dist1,nobs);

% Welfare gains;
[welfgain]=welfaregains(estparam,fund,Ctradesh,tradesh,CL,L,nobs);
welfgain=round(welfgain.*(10.^4));
welfgain=welfgain./(10.^4);

% Perfectly immobile welfare gains (ACR);
[acrwelfgain]=acrwelfaregains(estparam,fund,Ctradesh,tradesh,CL,L,nobs);

% Perfectly mobile welfare gains;
[mobwelfgain]=mobwelfaregains(estparam,fund,Ctradesh,tradesh,CL,L,nobs);

% % % % *****************************;
% % % % **** Solve Unobservables ****;
% % % % *****************************;
% % % 
% % % Cobserve=zeros(nobs,5);
% % % Cobserve(:,1)=CL; Cobserve(:,2)=Cw; Cobserve(:,3)=H; 
% % % 
% % % % Solve for region productivities and amenities;
% % % [Ca_i,Cb_i,Cabtradesh,Caconverge,Cbconverge,Cxtic]=solveab(estparam,Cobserve,dist1,nobs);
% % % display('>>>> Productivity and Amenity System Converged <<<<');
% % % display('>>>> Check Productivity and Amenity Convergence <<<<');
% % % [Caconverge Cbconverge]
% % % display('>>>> Elapsed Time in Seconds <<<<');
% % % xtic

% **************************;
% **** Relative Changes ****;
% **************************;

dL=CL./L; ldL=log(dL);
dw=Cw./w; ldw=log(dw);
dr=Cr./r; ldr=log(dr);
dP=CP./P; ldP=log(dP);
lacrwelfgain=log(acrwelfgain);
lmobwelfgain=log(mobwelfgain);
drealw=Crealwage./realwage; ldrealw=log(drealw);

% Population treatment;
[bL,bintL,rL,rintL,statsL] = regress(ldL,X);
dLmat(p)=bL(2);
% Wage treatment;
[bw,bintw,rw,rintw,statsw] = regress(ldw,X);
dwmat(p)=bw(2);
% Price index treatment;
[bP,bintP,rP,rintP,statsP] = regress(ldP,X);
dPmat(p)=bP(2);
% Land price treatment;
[br,bintr,rr,rintr,statsr] = regress(ldr,X);
drmat(p)=br(2);
% Real wage treatment;
[brw,bintrw,rrw,rintrw,statsrw] = regress(ldrealw,X);
drwmat(p)=brw(2);
% ACR welfare treatment;
[bacr,bintacr,racr,rintacr,statsacr] = regress(lacrwelfgain,X);
dacrmat(p)=bacr(2);
% ACR welfare treatment;
[bmob,bintmob,rmob,rintmob,statsmob] = regress(lmobwelfgain,X);
dmobmat(p)=bmob(2);

end;

% ****************************************************;
% **** Check all parameter combinations converged ****;
% ****************************************************;

display('>>>> Check all parameter combinations converged <<<<');
min(min(convmat))

% ********************************************************;
% **** Analyze how effects vary with parameter values ****;
% ********************************************************;

% POPULATION;
dLmat=reshape(dLmat,size(thetagrid,1),size(epsgrid,2));
XXL=min(epsgrid):0.001:max(epsgrid); 
YYL=min(thetagrid):0.001:max(thetagrid); YYL=YYL'; 
[XXL,YYL,ZZL]=griddata(epsgrid,thetagrid,dLmat,XXL,YYL);

% REAL WAGE;
drwmat=reshape(drwmat,size(thetagrid,1),size(epsgrid,2));
XXrw=min(epsgrid):0.001:max(epsgrid); 
YYrw=min(thetagrid):0.001:max(thetagrid); YYrw=YYrw';
[XXrw,YYrw,ZZrw]=griddata(epsgrid,thetagrid,drwmat,XXrw,YYrw);

% PRICE INDEX;
dPmat=reshape(dPmat,size(thetagrid,1),size(epsgrid,2));
XXP=min(epsgrid):0.001:max(epsgrid);
YYP=min(thetagrid):0.001:max(thetagrid); YYP=YYP';
[XXP,YYP,ZZP]=griddata(epsgrid,thetagrid,dPmat,XXP,YYP);

% WAGE;
dwmat=reshape(dwmat,size(thetagrid,1),size(epsgrid,2));
XXw=min(epsgrid):0.001:max(epsgrid);
YYw=min(thetagrid):0.001:max(thetagrid); YYw=YYw';
[XXw,YYw,ZZw]=griddata(epsgrid,thetagrid,dwmat,XXw,YYw);

% RELATIVE LAND PRICE;
drmat=reshape(drmat,size(thetagrid,1),size(epsgrid,2));
XXr=min(epsgrid):0.001:max(epsgrid);
YYr=min(thetagrid):0.001:max(thetagrid); YYr=YYr';
[XXr,YYr,ZZr]=griddata(epsgrid,thetagrid,drmat,XXr,YYr);

% ACR WELFARE GAIN;
dacrmat=reshape(dacrmat,size(thetagrid,1),size(epsgrid,2));
XXacrb=min(epsgrid):0.001:max(epsgrid);
YYacrb=min(thetagrid):0.001:max(thetagrid); YYacrb=YYacrb';
[XXacrb,YYacrb,ZZacrb]=griddata(epsgrid,thetagrid,dacrmat,XXacrb,YYacrb);

% MOBILITY WELFARE GAIN;
dmobmat=reshape(dmobmat,size(thetagrid,1),size(epsgrid,2));
XXmobb=min(epsgrid):0.001:max(epsgrid);
YYmobb=min(thetagrid):0.001:max(thetagrid); YYmobb=YYmobb';
[XXmobb,YYmobb,ZZmobb]=griddata(epsgrid,thetagrid,dmobmat,XXmobb,YYmobb);

% MULTI-PANEL FIGURE;
% (Figure 5 in the paper)
figure(1);
% Population;
subplot(3,2,1);
C=contourf(YYL,XXL,ZZL,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('Panel A : Population Treatment','FontSize',8);
set(gca,'fontsize',8);
% Wage;
subplot(3,2,2);
C=contourf(YYw,XXw,ZZw,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('Panel B : Wage Treatment','FontSize',8);
set(gca,'fontsize',8);
% Price index;
subplot(3,2,3);
C=contourf(YYP,XXP,ZZP,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('Panel C : Price Treatment','FontSize',8);
set(gca,'fontsize',8);
% Land price;
subplot(3,2,4);
C=contourf(YYr,XXr,ZZr,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('Panel D : Land Rent Treatment','FontSize',8);
set(gca,'fontsize',8);
% Real wage;
subplot(3,2,5);
C=contourf(YYrw,XXrw,ZZrw,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('Panel E : Real Wage Treatment','FontSize',8);
set(gca,'fontsize',8);
% Welfare
subplot(3,2,6);
C=contourf(YYacrb,XXacrb,ZZacrb,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('Panel F : Incorrect Immobile Welfare','FontSize',8);
set(gca,'fontsize',8);

print -dpdf graphs/compstatic.pdf;

% ACR WELFARE;
figure(2);
C=contourf(YYacrb,XXacrb,ZZacrb,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('ACR Welfare Bias','FontSize',8);
set(gca,'fontsize',8);

% MOBILITY WELFARE;
figure(3);
C=contourf(YYmobb,XXmobb,ZZmobb,10);
xlabel('Theta','FontSize',8);
ylabel('Epsilon','FontSize',8);
title('Mobility Welfare Bias','FontSize',8);
set(gca,'fontsize',8);







