
% Monte Carlo for Quantitative Spatial Model;
% Constant and increasing returns to scale model;
% Regions specification;

% SJR, November, 2015;

% *********************; 
% **** Choose User ****; 
% *********************;

clear all;
clc;
colormap jet;
close all;

user=1;

% User 1 : steve desktop; 
% User 2 : steve mac laptop dropbox;
% User 3 : steve home mac;

if user==1;
cd C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\regions;
addpath C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\regions\functions;
elseif user==2;
cd /Users/reddings/Dropbox/QuantSpatialModel/JIErev1/matlab/regions;
addpath /Users/reddings/Dropbox/QuantSpatialModel/JIErev1/matlab/regions/functions;
elseif user==3;
cd /Users/reddings/Dropbox/QuantSpatialModel/JIErev1/matlab/regions;
addpath /Users/reddings/Dropbox/QuantSpatialModel/JIErev1/matlab/regions/functions;
end;

format shortG;

global alpha sigma Hsigma theta epsilon LL nobs F;
global observe Cobserve dist0 dist1;

% ************************; 
% **** Initialization ****;
% ************************;

% Set default random number stream;
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

% *********************************************; 
% **** Trade cost matrix (called distance) ****;
% *********************************************;

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
Hsigma=5;
% Goods Frechet shape parameter;
theta=4; 
% Worker Frechet shape parameter;
epsilon=3;

param=[alpha theta epsilon];

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
% Fixed production cost;
F=1;

% ****************************************;
% **** Solve for Endogenous Variables ****;
% ****************************************;

fund=zeros(nobs,5);
fund(:,1)=a; fund(:,2)=b; fund(:,3)=H;

% Solve for region populations and wages;
[w,L,tradesh,dtradesh,Lconverge,wconverge,xtic]=solveLw(param,fund,dist0,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[wconverge Lconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Price index;
[P] = pindex(param,fund,w,dtradesh,nobs);

% Land price;
[r]=landprice(param,fund,L,w,dist0,nobs);

% Expected utility;
[EU]=expectut(param,fund,L,w,tradesh,dist0,nobs);
display('>>>> Expected utility <<<<');
[EU]

% Welfare;
[welf]=welfare(param,fund,L,w,tradesh,dist0,nobs);
display('>>>> Welfare <<<<');
welf=round(welf.*(10.^4));
welf=welf./(10.^4);
unique(welf)

% Real wage;
[realwage]=realw(param,fund,L,w,tradesh,dist0,nobs);

% *********************************;
% **** Solve for Unobservables ****;
% *********************************;

observe=zeros(nobs,5);
observe(:,1)=L; observe(:,2)=w; observe(:,3)=H; 

% Solve for region productivities and amenities;
[a_i,b_i,abtradesh,aconverge,bconverge,xtic]=solveab(param,observe,dist0,nobs);
display('>>>> Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[aconverge bconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% *****************************************;
% **** Solve for Helpman Unobservables ****;
% *****************************************;

% Solve for region productivities and amenities;
[Ha_i,Hb_i,Habtradesh,Haconverge,Hbconverge,xtic]=solveHab(param,observe,dist0,nobs);
display('>>>> Helpman Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[Haconverge Hbconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

Hfund=zeros(nobs,5);
Hfund(:,1)=Ha_i; Hfund(:,2)=Hb_i; Hfund(:,3)=H;

% ************************************************;
% **** Solve for Helpman Endogenous Variables ****;
% ************************************************;

% Solve for region populations and wages;
[Hw,HL,Htradesh,Hdtradesh,HLconverge,Hwconverge,Hxtic]=solveHLw(param,Hfund,dist0,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[Hwconverge HLconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Price index;
[HP] = Hpindex(param,Hfund,HL,Hw,Hdtradesh,nobs);

% Land price;
[Hr]=landprice(param,Hfund,HL,Hw,dist0,nobs);

% Expected utility;
[HEU]=Hexpectut(param,Hfund,HL,Hw,HP,Hr,dist0,nobs);
display('>>>> Helpman Expected utility <<<<');
[HEU]

% Welfare;
[Hwelf]=Hwelfare(param,Hfund,HL,Hw,Htradesh,dist0,nobs);
display('>>>> Helpman Welfare <<<<');
Hwelf=round(Hwelf.*(10.^4));
Hwelf=Hwelf./(10.^4);
unique(Hwelf)

% Real wage;
[Hrealwage]=Hrealw(param,Hfund,HL,Hw,Htradesh,dist0,nobs);

% *******************************;
% **** CHANGE IN TRADE COSTS ****;
% *******************************;

% Solve for region populations and wages;
[Cw,CL,Ctradesh,Cdtradesh,CLconverge,Cwconverge,Cxtic]=solveLw(param,fund,dist1,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[Cwconverge CLconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Counterfactual price index;
[CP] = pindex(param,fund,Cw,Cdtradesh,nobs);

% Counterfactual land prices;
[Cr]=landprice(param,fund,CL,Cw,dist1,nobs);

% Counterfactual expected utility;
[CEU]=expectut(param,fund,CL,Cw,Ctradesh,dist1,nobs);
display('>>>> Expected utility <<<<');
CEU

% Counterfactual welfare;
[Cwelf]=welfare(param,fund,CL,Cw,Ctradesh,dist1,nobs);
display('>>>> Welfare <<<<');
Cwelf=round(Cwelf.*(10.^4));
Cwelf=Cwelf./(10.^4);
unique(Cwelf)

% Counterfactual real wage;
[Crealwage]=realw(param,fund,CL,Cw,Ctradesh,dist1,nobs);

% Welfare gains;
[welfgain]=welfaregains(param,fund,Ctradesh,tradesh,CL,L,nobs);
display('>>>> Welfare Gains <<<<');
welfgain=round(welfgain.*(10.^4));
welfgain=welfgain./(10.^4);
unique(welfgain)

% Perfectly immobile welfare gains (ACR);
[acrwelfgain]=acrwelfaregains(param,fund,Ctradesh,tradesh,CL,L,nobs);

% Perfectly mobile no preference heterogeneity welfare gains;
[mobwelfgain]=mobwelfaregains(param,fund,Ctradesh,tradesh,CL,L,nobs);

% ***************************************;
% **** HELPMAN CHANGE IN TRADE COSTS ****;
% ***************************************;

% Solve for region populations and wages;
[CHw,CHL,CHtradesh,CHdtradesh,CHLconverge,CHwconverge,CHxtic]=solveHLw(param,Hfund,dist1,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[CHwconverge CHLconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Counterfactual price index;
[CHP] = Hpindex(param,Hfund,CHL,CHw,CHdtradesh,nobs);

% Counterfactual land prices;
[CHr]=landprice(param,Hfund,CHL,CHw,dist1,nobs);

% Counterfactual expected utility;
[CHEU]=Hexpectut(param,Hfund,CHL,CHw,CHP,CHr,dist0,nobs);
display('>>>> Helpman Expected utility <<<<');
CHEU

% Counterfactual welfare;
[CHwelf]=Hwelfare(param,Hfund,CHL,CHw,CHtradesh,dist0,nobs);
display('>>>> Helpman Welfare <<<<');
CHwelf=round(CHwelf.*(10.^4));
CHwelf=CHwelf./(10.^4);
unique(CHwelf)

% Counterfactual real wage;
[CHrealwage]=Hrealw(param,Hfund,CHL,CHw,CHtradesh,dist0,nobs);

% Welfare gains;
[Hwelfgain]=Hwelfaregains(param,Hfund,CHtradesh,tradesh,CHL,L,nobs);
display('>>>> Helpman Welfare Gains <<<<');
Hwelfgain=round(Hwelfgain.*(10.^4));
Hwelfgain=Hwelfgain./(10.^4);
unique(Hwelfgain)

% Perfectly immobile welfare gains (ACR);
[Hacrwelfgain]=acrwelfaregains(param,Hfund,CHtradesh,tradesh,CHL,L,nobs);

% Perfectly mobile welfare gains;
[Hmobwelfgain]=mobwelfaregains(param,Hfund,CHtradesh,tradesh,CHL,L,nobs);

% *****************************;
% **** Solve Unobservables ****;
% *****************************;

Cobserve=zeros(nobs,5);
Cobserve(:,1)=CL; Cobserve(:,2)=Cw; Cobserve(:,3)=H; 

% Solve for region productivities and amenities;
[Ca_i,Cb_i,Cabtradesh,Caconverge,Cbconverge,Cxtic]=solveab(param,Cobserve,dist1,nobs);
display('>>>> Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[Caconverge Cbconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% *************************************;
% **** Solve Helpman Unobservables ****;
% *************************************;

CHobserve=zeros(nobs,5);
CHobserve(:,1)=CHL; CHobserve(:,2)=CHw; CHobserve(:,3)=H; 

% Solve for region productivities and amenities;
[CHa_i,CHb_i,CHabtradesh,CHaconverge,CHbconverge,CHxtic]=solveHab(param,CHobserve,dist1,nobs);
display('>>>> Helpman Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[CHaconverge CHbconverge]
display('>>>> Elapsed Time in Seconds <<<<');
CHxtic

% ***********************************;
% **** MEAN CHANGE IN TRADE COST ****;
% ***********************************;

mnrdist=mean(rdist,2);

mnrdmat=reshape(mnrdist,N,N);
XXL=min(lgd):0.001:max(lgd);
YYL=min(ltd):0.001:max(ltd); YYL=YYL';
[XXD,YYD,ZZD]=griddata(lgd,ltd,mnrdmat,XXL,YYL);

figure(1);
% Mean change in trade cost;
% (Figure 2 in the paper);
C=contourf(XXD,YYD,ZZD,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca,'fontsize',8);

print -dpdf graphs/mean_rdist.pdf;

% ***********************************************;
% **** Three-Dimensional Initial Equilibrium ****;
% ***********************************************;

% PRODUCTIVITY;
amat=reshape(a,N,N);
XXL=min(lgd):0.001:max(lgd);
YYL=min(ltd):0.001:max(ltd); YYL=YYL';
[XXA,YYA,ZZA]=griddata(lgd,ltd,amat,XXL,YYL);

% AMENITIES;
bmat=reshape(b,N,N);
XXL=min(lgd):0.001:max(lgd);
YYL=min(ltd):0.001:max(ltd); YYL=YYL';
[XXB,YYB,ZZB]=griddata(lgd,ltd,bmat,XXL,YYL);

% POPULATION;
Lmat=reshape(L,N,N);
XXL=min(lgd):0.001:max(lgd);
YYL=min(ltd):0.001:max(ltd); YYL=YYL';
[XXL,YYL,ZZL]=griddata(lgd,ltd,Lmat,XXL,YYL);

% PRICE INDEX;
Pmat=reshape(P,N,N);
XXP=min(lgd):0.001:max(lgd);
YYP=min(ltd):0.001:max(ltd); YYP=YYP';
[XXP,YYP,ZZP]=griddata(lgd,ltd,Pmat,XXP,YYP);

% WAGE;
wmat=reshape(w,N,N);
XXw=min(lgd):0.001:max(lgd);
YYw=min(ltd):0.001:max(ltd); YYw=YYw';
[XXw,YYw,ZZw]=griddata(lgd,ltd,wmat,XXw,YYw);

% RELATIVE LAND PRICE;
rmat=reshape(r,N,N);
XXr=min(lgd):0.001:max(lgd);
YYr=min(ltd):0.001:max(ltd); YYr=YYr';
[XXr,YYr,ZZr]=griddata(lgd,ltd,rmat,XXr,YYr);

% MULTI-PANEL FIGURE;
figure(2);
% Productivity;
subplot(3,2,1);
C=contourf(XXA,YYA,ZZA,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel A : Productivity','FontSize',8);
set(gca,'fontsize',8);
% Amenities;
subplot(3,2,2);
C=contourf(XXB,YYB,ZZB,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel B : Amenities','FontSize',8);
set(gca,'fontsize',8);
% Population;
subplot(3,2,3);
C=contourf(XXL,YYL,ZZL,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel C : Population','FontSize',8);
set(gca,'fontsize',8);
% Price index;
subplot(3,2,4);
C=contourf(XXP,YYP,ZZP,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel D : Price Index','FontSize',8);
set(gca,'fontsize',8);
% Wages;
subplot(3,2,5);
C=contourf(XXw,YYw,ZZw,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel E : Wages','FontSize',8);
set(gca,'fontsize',8);
% Land prices;
subplot(3,2,6);
C=contourf(XXr,YYr,ZZr,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel F : Land Prices','FontSize',8);
set(gca,'fontsize',8);

print -dpdf graphs/initial_equil.pdf;

% **************************************************************************;
% **** Three-Dimensional Impact of Change in Trade Costs Counterfactual ****;
% **************************************************************************;

dL=CL./L; ldL=log(dL);
dw=Cw./w; ldw=log(dw);
dr=Cr./r; ldr=log(dr);
dP=CP./P; ldP=log(dP);
lacrwelfgain=log(acrwelfgain);
drealw=Crealwage./realwage; ldrealw=log(drealw);
dtradesh=diag(tradesh);
Cdtradesh=diag(Ctradesh);
ddtradesh=Cdtradesh./dtradesh;

HdL=CHL./HL; lHdL=log(HdL);
Hdw=CHw./Hw; lHdw=log(Hdw);
Hdr=CHr./Hr; lHdr=log(Hdr);
HdP=CHP./HP; lHdP=log(HdP);
lHacrwelfgain=log(Hacrwelfgain);
Hdrealw=CHrealwage./Hrealwage; lHdrealw=log(Hdrealw);
CHdtradesh=diag(CHtradesh);
Hddtradesh=CHdtradesh./dtradesh;

% GRID DATA;

% POPULATION;
dLmat=reshape(dL,N,N);
XXL=min(lgd):0.001:max(lgd);
YYL=min(ltd):0.001:max(ltd); YYL=YYL';
[XXL,YYL,ZZL]=griddata(lgd,ltd,dLmat,XXL,YYL);

% REAL WAGE;
drealwmat=reshape(drealw,N,N);
XXrw=min(lgd):0.001:max(lgd);
YYrw=min(ltd):0.001:max(ltd); YYrw=YYrw';
[XXrw,YYrw,ZZrw]=griddata(lgd,ltd,drealwmat,XXrw,YYrw);

% PRICE INDEX;
dPmat=reshape(dP,N,N);
XXP=min(lgd):0.001:max(lgd);
YYP=min(ltd):0.001:max(ltd); YYP=YYP';
[XXP,YYP,ZZP]=griddata(lgd,ltd,dPmat,XXP,YYP);

% WAGE;
dwmat=reshape(dw,N,N);
XXw=min(lgd):0.001:max(lgd);
YYw=min(ltd):0.001:max(ltd); YYw=YYw';
[XXw,YYw,ZZw]=griddata(lgd,ltd,dwmat,XXw,YYw);

% RELATIVE LAND PRICE;
drmat=reshape(dr,N,N);
XXr=min(lgd):0.001:max(lgd);
YYr=min(ltd):0.001:max(ltd); YYr=YYr';
[XXr,YYr,ZZr]=griddata(lgd,ltd,drmat,XXr,YYr);

% ACR WELFARE;
acrmat=reshape(acrwelfgain,N,N);
XXa=min(lgd):0.001:max(lgd);
YYa=min(ltd):0.001:max(ltd); YYa=YYa';
[XXa,YYa,ZZa]=griddata(lgd,ltd,acrmat,XXa,YYa);

% WELFARE;
welfmat=reshape(welfgain,N,N);
XXwelf=min(lgd):0.001:max(lgd);
YYwelf=min(ltd):0.001:max(ltd); YYwelf=YYwelf';
[XXwelf,YYwelf,ZZwelf]=griddata(lgd,ltd,welfmat,XXwelf,YYwelf);

% MULTI-PANEL FIGURE;
% Impact of the transport improvement in the constant returns to scale model;
% (Figure 3 in the paper)
figure(3);
% Population;
subplot(3,2,1);
C=contourf(XXL,YYL,ZZL,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel A : Population','FontSize',8);
set(gca,'fontsize',8);
% Wage;
subplot(3,2,2);
C=contourf(XXw,YYw,ZZw,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel B : Wages','FontSize',8);
set(gca,'fontsize',8);
% Price index;
subplot(3,2,3);
C=contourf(XXP,YYP,ZZP,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel C : Price Index','FontSize',8);
set(gca,'fontsize',8);
% Land price;
subplot(3,2,4);
C=contourf(XXr,YYr,ZZr,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel D : Land Rents','FontSize',8);
% Real wage;
subplot(3,2,5);
C=contourf(XXrw,YYrw,ZZrw,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel E : Real Wage','FontSize',8);
set(gca,'fontsize',8);
% Welfare;
subplot(3,2,6);
C=contourf(XXa,YYa,ZZa,5);
xlabel('Longitude','FontSize',8);
ylabel('Latitude','FontSize',8);
set(gca, 'XTickLabel',[2 4 6 8 10],'XTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
set(gca, 'YTickLabel',[2 4 6 8 10],'YTick',[lgd(2) lgd(4) lgd(6) lgd(8) lgd(10)]);
title('Panel F : Incorrect Immobile Welfare','FontSize',8);
set(gca,'fontsize',8);

print -dpdf graphs/transport_impact.pdf;

% HISTOGRAM FIGURE;
% Histogram for the impact of the transport improvement in the constant returns to scale model;
% (Figure 4 in the paper)
figure(4);
% Population;
s1=subplot(3,2,1);
bin=linspace(min(dL),max(dL),20); bin=bin';
h=histc(dL,bin);
h=h./sum(h);
h1=histc(dL(treat==1),bin);
h1=h1./sum(h1);
h0=histc(dL(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel A : Population','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Real wage;
s2=subplot(3,2,2);
bin=linspace(min(drealw),max(drealw),20); bin=bin';
h=histc(drealw,bin);
h=h./sum(h);
h1=histc(drealw(treat==1),bin);
h1=h1./sum(h1);
h0=histc(drealw(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel B : Real Wage','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Price index;
s3=subplot(3,2,3);
bin=linspace(min(dP),max(dP),20); bin=bin';
h=histc(dP,bin);
h=h./sum(h);
h1=histc(dP(treat==1),bin);
h1=h1./sum(h1);
h0=histc(dP(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel C : Price Index','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Wage;
s4=subplot(3,2,4);
bin=linspace(min(dw),max(dw),20); bin=bin';
h=histc(dw,bin);
h=h./sum(h);
h1=histc(dw(treat==1),bin);
h1=h1./sum(h1);
h0=histc(dw(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel D : Wage','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Land price;
s5=subplot(3,2,5);
bin=linspace(min(dr),max(dr),20); bin=bin';
h=histc(dr,bin);
h=h./sum(h);
h1=histc(dr(treat==1),bin);
h1=h1./sum(h1);
h0=histc(dr(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel E : Land Rents','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% ACR welfare;
s6=subplot(3,2,6);
bin=linspace(min(acrwelfgain),max(acrwelfgain),20); bin=bin';
h=histc(acrwelfgain,bin);
h=h./sum(h);
h1=histc(acrwelfgain(treat==1),bin);
h1=h1./sum(h1);
h0=histc(acrwelfgain(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel F : Incorrect Immobile Welfare','FontSize',8);
hold on;
%Hard to display actual welfare gain due to bin definition;
%line([welfgain welfgain],[0 0.6],'Color','red','LineStyle','-','LineWidth',1.5);
hold off
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);

print -dpdf graphs/transport_histogram.pdf;

figure(5);
% Perfect mobility welfare;
% Implied welfare gains under the assumption of perfect mobility with 
% no heterogeneity in worker preferences;
% (no longer in the paper)
bin=linspace(min(acrwelfgain),max(mobwelfgain),20); bin=bin';
h1=histc(acrwelfgain,bin);
h1=h1./sum(h1);
h0=histc(mobwelfgain,bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*100); bbin=bbin/100;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
xlabel('Relative Change in Welfare','FontSize',8);
ylabel('Probability','FontSize',8);
axis([0 20 0 0.6]);
hold on;
line([0 20],[0 0],'Color','black','LineWidth',1);
hold off
set(gca,'fontsize',8);

print -dpdf graphs/perfmobil_histogram.pdf;

% HELPMAN GRID DATA;

% HELPMAN POPULATION;
HdLmat=reshape(HdL,N,N);
XXL=min(lgd):0.001:max(lgd);
YYL=min(ltd):0.001:max(ltd); YYL=YYL';
[XXL,YYL,HZZL]=griddata(lgd,ltd,HdLmat,XXL,YYL);

% HELPMAN REAL WAGE;
Hdrealwmat=reshape(Hdrealw,N,N);
XXrw=min(lgd):0.001:max(lgd);
YYrw=min(ltd):0.001:max(ltd); YYrw=YYrw';
[XXrw,YYrw,HZZrw]=griddata(lgd,ltd,Hdrealwmat,XXrw,YYrw);

% HELPMAN PRICE INDEX;
HdPmat=reshape(HdP,N,N);
XXP=min(lgd):0.001:max(lgd);
YYP=min(ltd):0.001:max(ltd); YYP=YYP';
[XXP,YYP,HZZP]=griddata(lgd,ltd,HdPmat,XXP,YYP);

% HELPMAN WAGE;
Hdwmat=reshape(Hdw,N,N);
XXw=min(lgd):0.001:max(lgd);
YYw=min(ltd):0.001:max(ltd); YYw=YYw';
[XXw,YYw,HZZw]=griddata(lgd,ltd,Hdwmat,XXw,YYw);

% HELPMAN RELATIVE LAND PRICE;
Hdrmat=reshape(Hdr,N,N);
XXr=min(lgd):0.001:max(lgd);
YYr=min(ltd):0.001:max(ltd); YYr=YYr';
[XXr,YYr,HZZr]=griddata(lgd,ltd,Hdrmat,XXr,YYr);

% HELPMAN ACR WELFARE;
Hacrmat=reshape(Hacrwelfgain,N,N);
XXa=min(lgd):0.001:max(lgd);
YYa=min(ltd):0.001:max(ltd); YYa=YYa';
[XXa,YYa,HZZa]=griddata(lgd,ltd,Hacrmat,XXa,YYa);

% HELPMAN WELFARE;
Hwelfmat=reshape(Hwelfgain,N,N);
XXwelf=min(lgd):0.001:max(lgd);
YYwelf=min(ltd):0.001:max(ltd); YYwelf=YYwelf';
[XXwelf,YYwelf,HZZwelf]=griddata(lgd,ltd,Hwelfmat,XXwelf,YYwelf);

% HELPMAN HISTOGRAM FIGURE;
% Histogram for impact of transport improvement in the 
% calibrated increasing returns to scale model
% (Figure 7 in the paper)
figure(7);
% Population;
s1=subplot(3,2,1);
bin=linspace(min(HdL),max(HdL),20); bin=bin';
h=histc(HdL,bin);
h=h./sum(h);
h1=histc(HdL(treat==1),bin);
h1=h1./sum(h1);
h0=histc(HdL(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel A : Population','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Real wage;
s2=subplot(3,2,2);
bin=linspace(min(Hdrealw),max(Hdrealw),20); bin=bin';
h=histc(Hdrealw,bin);
h=h./sum(h);
h1=histc(Hdrealw(treat==1),bin);
h1=h1./sum(h1);
h0=histc(Hdrealw(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel B : Real Wage','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Price index;
s3=subplot(3,2,3);
bin=linspace(min(HdP),max(HdP),20); bin=bin';
h=histc(HdP,bin);
h=h./sum(h);
h1=histc(HdP(treat==1),bin);
h1=h1./sum(h1);
h0=histc(HdP(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel C : Price Index','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Wage;
s4=subplot(3,2,4);
bin=linspace(min(Hdw),max(Hdw),20); bin=bin';
h=histc(Hdw,bin);
h=h./sum(h);
h1=histc(Hdw(treat==1),bin);
h1=h1./sum(h1);
h0=histc(Hdw(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel D : Wage','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% Land price;
s5=subplot(3,2,5);
bin=linspace(min(Hdr),max(Hdr),20); bin=bin';
h=histc(Hdr,bin);
h=h./sum(h);
h1=histc(Hdr(treat==1),bin);
h1=h1./sum(h1);
h0=histc(Hdr(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel E : Land Rents','FontSize',8);
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);
% ACR welfare;
s6=subplot(3,2,6);
bin=linspace(min(Hacrwelfgain),max(Hacrwelfgain),20); bin=bin';
h=histc(Hacrwelfgain,bin);
h=h./sum(h);
h1=histc(Hacrwelfgain(treat==1),bin);
h1=h1./sum(h1);
h0=histc(Hacrwelfgain(treat==0),bin);
h0=h0./sum(h0);
hb=bar([h1 h0]);
hbc = get(hb, 'Children');
set(hbc{1},'FaceColor', 'k','EdgeColor', 'k');
set(hbc{2},'FaceColor', 'c','EdgeColor', 'c');
bbin=round(bin*10); bbin=bbin/10;
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
axis([0 20 0 0.6]);
line([0 20],[0 0],'Color','black','LineWidth',1);
ylabel('Probability','FontSize',8);
title('Panel F : Incorrect Immobile Welfare','FontSize',8);
hold on;
%Hard to display actual welfare gain due to bin definition;
%line([welfgain welfgain],[0 0.6],'Color','red','LineStyle','-','LineWidth',1.5);
hold off
set(gca,'YTick',[0 0.2 0.4 0.6],'fontsize',8);

print -dpdf graphs/H_transport_histogram.pdf;

% HELPMAN PRODUCTIVITIES FIGURE;
% Figure graphing calibrated productivity in the constant and increasing; 
% returns models against one another (Figure 6 in the paper);
figure(8);
s1=subplot(1,2,1);
plot(log(a),log(Ha_i),'bo','MarkerSize',4);
hold on;
plot(log(a),log(a),'k');
ylabel('Increasing Returns','FontSize',8);
xlabel('Constant Returns','FontSize',8);
title('Panel A : Log Productivity','FontSize',8);
set(gca,'fontsize',8);
axis([min(log(a)) max(log(a)) min(log(a)) max(log(a))]);
s2=subplot(1,2,2);
plot(log(b),log(Hb_i),'bo','MarkerSize',4);
hold on;
plot(log(b),log(b),'k');
ylabel('Increasing Returns','FontSize',8);
xlabel('Constant Returns','FontSize',8);
title('Panel B : Log Amenities','FontSize',8);
set(gca,'fontsize',8);
axis([min(log(b)) max(log(b)) min(log(b)) max(log(b))]);

print -dpdf graphs/H_prod_amen.pdf;

% *******************************;
% **** Treatment Regressions ****;
% *******************************;

display('>>>> Treatment Regressions <<<<');

% Define controls;
X=[ones(size(treat)) treat];

% Population treatment;
[bL,bintL,rL,rintL,statsL] = regress(ldL,X);
display('>>>> Population Treatment <<<<');
bL(2)
stats = regstats(ldL,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Real wage treatment;
[brw,bintrw,rrw,rintrw,statsrw] = regress(ldrealw,X);
display('>>>> Real Wage Treatment <<<<');
brw(2)
stats = regstats(ldrealw,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Price index treatment;
[bP,bintP,rP,rintP,statsP] = regress(ldP,X);
display('>>>> Price Index Treatment <<<<');
bP(2)
stats = regstats(ldP,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Wage treatment;
[bw,bintw,rw,rintw,statsw] = regress(ldw,X);
display('>>>> Wage Treatment <<<<');
bw(2)
stats = regstats(ldw,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Land price treatment;
[br,bintr,rr,rintr,statsr] = regress(ldr,X);
display('>>>> Land Price Treatment <<<<');
br(2)
stats = regstats(ldr,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% ACR welfare
[bacr,bintacr,racr,rintacr,statsacr] = regress(lacrwelfgain,X);
display('>>>> ACR Welfare Treatment <<<<');
bacr(2)
stats = regstats(lacrwelfgain,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% ***************************************;
% **** Helpman Treatment Regressions ****;
% ***************************************;

display('>>>> Helpman Treatment Regressions <<<<');

% Define controls;
X=[ones(size(treat)) treat];

% Population treatment;
[bL,bintL,rL,rintL,statsL] = regress(lHdL,X);
display('>>>> Population Treatment <<<<');
bL(2)
stats = regstats(lHdL,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Real wage treatment;
[brw,bintrw,rrw,rintrw,statsrw] = regress(lHdrealw,X);
display('>>>> Real Wage Treatment <<<<');
brw(2)
stats = regstats(lHdrealw,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Price index treatment;
[bP,bintP,rP,rintP,statsP] = regress(lHdP,X);
display('>>>> Price Index Treatment <<<<');
bP(2)
stats = regstats(lHdP,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Wage treatment;
[bw,bintw,rw,rintw,statsw] = regress(lHdw,X);
display('>>>> Wage Treatment <<<<');
bw(2)
stats = regstats(lHdw,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Land price treatment;
[br,bintr,rr,rintr,statsr] = regress(lHdr,X);
display('>>>> Land Price Treatment <<<<');
br(2)
stats = regstats(lHdr,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% ACR welfare
[bacr,bintacr,racr,rintacr,statsacr] = regress(lHacrwelfgain,X);
display('>>>> ACR Welfare Treatment <<<<');
bacr(2)
stats = regstats(lHacrwelfgain,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

