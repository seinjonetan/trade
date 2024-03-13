
% Monte Carlo for Quantitative Spatial Model;
% Increasing returns to scale model;
% Countries and regions;

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
% User 3 : steve home mac;

if user==1;
cd C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\countries;
addpath C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\countries\functions;
end;

format shortG;

global alpha Hsigma theta epsilon LL LLwest LLeast F;
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

% Change in transport infrastructure;

% tau0=zeros(N);
% tau0(:,:)=tt0;
% 
% tau1=zeros(N);
% tau1(:,:)=tt0;
% tau1(:,6)=tt1;
% tau1(6,:)=tt1;

% Always transport infrastructure;

tau0=zeros(N);
tau0(:,:)=tt0;
tau0(:,6)=tt1;
tau0(6,:)=tt1;

tau1=zeros(N);
tau1(:,:)=tt0;
tau1(:,6)=tt1;
tau1(6,:)=tt1;

% Low trade costs everywhere;

% tau0=zeros(N);
% tau0(:,:)=tt1;
% 
% tau1=zeros(N);
% tau1(:,:)=tt1;

% Compute transport costs;

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

% Define east and west;
Iwest=zeros(N,N);
Ieast=zeros(N,N);
Iwest(:,1:7)=1;
Ieast(:,8:11)=1;
Iwest=reshape(Iwest,NN,1);
Ieast=reshape(Ieast,NN,1);

% Define distance weights;
dopen=ones(size(dist0));
dclosed=zeros(size(dist0));
dclosed(Iwest==1,Iwest==1)=1;
dclosed(Ieast==1,Ieast==1)=1;

% Trade costs are a power function of effective distance;

dist0=dist0.^0.33;
dist1=dist1.^0.33;

% **************************;
% **** Parameterization ****;
% **************************;

% Share of goods in consumption expenditure (1-housing share);
alpha=0.75; 
% Elasticity of substitution;
Hsigma=5;
% Goods Frechet shape parameter;
theta=Hsigma-1; 
% Worker Frechet shape parameter;
epsilon=3;

param=[alpha theta epsilon];

% ***********************;
% **** Random shocks ****;
% ***********************;

a=normrnd(0,1,NN,1);
a=exp(a);
a(Iwest==1)=a(Iwest==1)./geomean(a(Iwest==1));
a(Ieast==1)=a(Ieast==1)./geomean(a(Ieast==1));

b=normrnd(0,1,NN,1);
b=exp(b);
b(Iwest==1)=b(Iwest==1)./geomean(b(Iwest==1));
b(Ieast==1)=b(Ieast==1)./geomean(b(Ieast==1));

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
LLwest=(sum(Iwest)./(sum(Iwest)+sum(Ieast))).*LL;
LLeast=(sum(Ieast)./(sum(Iwest)+sum(Ieast))).*LL;
% Fixed production cost;
F=1;

% ****************************************************************;
% **** Closed Economy Solve for Endogenous Neoclassical Model ****;
% ****************************************************************;

fund=zeros(nobs,5);
fund(:,1)=a; fund(:,2)=b; fund(:,3)=H; fund(:,4)=Iwest; fund(:,5)=Ieast;

% Solve for region populations and wages;
[w,L,tradesh,dtradesh,Lconverge,wconverge,xtic]=solveLwCtyClosed(param,fund,dclosed,dist0,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[wconverge Lconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% **************************************************************;
% **** Closed Economy Solve for Unobservables Helpman model ****;
% **************************************************************;

observe=zeros(nobs,5);
observe(:,1)=L; observe(:,2)=w; observe(:,3)=H; observe(:,4)=Iwest; observe(:,5)=Ieast;

% Solve for region productivities and amenities;
[Ha,Hb,abtradesh,aconverge,bconverge,xtic]=solveHabCtyClosed(param,observe,dclosed,dist0,nobs);
display('>>>> Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[aconverge bconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% *********************************************************************;
% **** Closed Economy Solve for Endogenous variables Helpman model ****;
% *********************************************************************;

fund=zeros(nobs,5);
fund(:,1)=Ha; fund(:,2)=Hb; fund(:,3)=H; fund(:,4)=Iwest; fund(:,5)=Ieast;

% Solve for region populations and wages;
[w,L,tradesh,dtradesh,Lconverge,wconverge,xtic]=solveHLwCtyClosed(param,fund,dclosed,dist0,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[wconverge Lconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Price index;
[P] = Hpindex(param,fund,L,w,dtradesh,nobs);

% Land price;
[r]=landprice(param,fund,L,w,dist0,nobs);

% Expected utility;
[EU]=Hexpectut(param,fund,L,w,P,r,dist0,nobs);
display('>>>> Expected utility (West, East) <<<<');
unique(EU)

% Welfare;
[welf]=Hwelfare(param,fund,L,w,tradesh,dist0,nobs);
display('>>>> Welfare <<<<');
welf=round(welf.*(10.^4));
welf=welf./(10.^4);
unique(welf)

% Real wage;
[realwage]=Hrealw(param,fund,L,w,tradesh,dist0,nobs);

% ************************************************;
% **** Closed Economy Solve for Unobservables ****;
% ************************************************;

observe=zeros(nobs,5);
observe(:,1)=L; observe(:,2)=w; observe(:,3)=H; observe(:,4)=Iwest; observe(:,5)=Ieast;

% Solve for region productivities and amenities;
[a_i,b_i,abtradesh,aconverge,bconverge,xtic]=solveHabCtyClosed(param,observe,dclosed,dist0,nobs);
display('>>>> Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[aconverge bconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% ***************************************************;
% **** CLOSED ECONOMY AGGREGATE TO COUNTRY LEVEL ****;
% ***************************************************;

NOBSC=2;

WBILL=zeros(NOBSC,1);
WBILL(1,1)=sum(w(Iwest==1).*L(Iwest==1)); WBILL(2,1)=sum(w(Ieast==1).*L(Ieast==1)); 

OBSERVE=zeros(NOBSC,5); 
OBSERVE(1,1)=sum(L(Iwest==1)); OBSERVE(2,1)=sum(L(Ieast==1));
OBSERVE(1,2)=WBILL(1,1)./OBSERVE(1,1); OBSERVE(2,2)=WBILL(2,1)./OBSERVE(2,1);
OBSERVE(1,3)=sum(H(Iwest==1)); OBSERVE(2,3)=sum(H(Ieast==1));
OBSERVE(1,4)=1; OBSERVE(2,4)=0;
OBSERVE(1,5)=0; OBSERVE(2,5)=1;

income=w.*L;
trade=tradesh.*repmat(income',NN,1);
TRADE=zeros(NOBSC,NOBSC);
TRADE(1,1)=sum(sum(trade(Iwest==1,Iwest==1)));
TRADE(2,2)=sum(sum(trade(Ieast==1,Ieast==1)));
TRADE(1,2)=sum(sum(trade(Iwest==1,Ieast==1)));
TRADE(2,1)=sum(sum(trade(Ieast==1,Iwest==1)));
EXPEND=sum(TRADE,2);
TRADESH=TRADE./repmat(EXPEND,1,2);
DTRADESH=diag(TRADESH);

% *****************************************************;
% **** Open Economy Solve for Endogenous Variables ****;
% *****************************************************;

% Solve for region populations and wages;
[Cw,CL,Ctradesh,Cdtradesh,CLconverge,Cwconverge,Cxtic]=solveHLwCtyOpen(param,fund,dopen,dist1,nobs);
display('>>>> Wage and Population System Converged <<<<');
display('>>>> Check Wage and Population Convergence <<<<');
[Cwconverge CLconverge]
display('>>>> Elapsed Time in Seconds <<<<');
xtic

% Counterfactual price index;
[CP] = Hpindex(param,fund,CL,Cw,Cdtradesh,nobs);

% Counterfactual land prices;
[Cr]=landprice(param,fund,CL,Cw,dist1,nobs);

% Counterfactual expected utility;
[CEU]=Hexpectut(param,fund,CL,Cw,CP,Cr,dist1,nobs);
display('>>>> Expected utility (West, East) <<<<');
unique(CEU)

% Counterfactual welfare;
[Cwelf]=Hwelfare(param,fund,CL,Cw,Ctradesh,dist1,nobs);
display('>>>> Welfare <<<<');
Cwelf=round(Cwelf.*(10.^4));
Cwelf=Cwelf./(10.^4);
unique(Cwelf)

% Counterfactual real wage;
[Crealwage]=Hrealw(param,fund,CL,Cw,Ctradesh,dist1,nobs);

% Welfare gains;
[welfgain]=Hwelfaregains(param,Ctradesh,tradesh,CL,L,nobs);
display('>>>> Welfare Gains <<<<');
welfgain=round(welfgain.*(10.^4));
welfgain=welfgain./(10.^4);
unique(welfgain)

% Perfectly immobile welfare gains (ACR);
[acrwelfgain]=acrwelfaregains(param,Ctradesh,tradesh,CL,L,nobs);

% Perfectly mobile welfare gains;
[mobwelfgain]=Hmobwelfaregains(param,Ctradesh,tradesh,CL,L,nobs);

% ******************************************;
% **** Open Economy Solve Unobservables ****;
% ******************************************;

Cobserve=zeros(nobs,5);
Cobserve(:,1)=CL; Cobserve(:,2)=Cw; Cobserve(:,3)=H; Cobserve(:,4)=Iwest; Cobserve(:,5)=Ieast;

% Solve for region productivities and amenities;
[Ca_i,Cb_i,Cabtradesh,Caconverge,Cbconverge,Cxtic]=solveHabCtyOpen(param,Cobserve,dopen,dist1,nobs);
display('>>>> Productivity and Amenity System Converged <<<<');
display('>>>> Check Productivity and Amenity Convergence <<<<');
[Caconverge Cbconverge]
display('>>>> Elapsed Time in Seconds <<<<');
Cxtic

% ************************************;
% **** AGGREGATE TO COUNTRY LEVEL ****;
% ************************************;

CWBILL=zeros(NOBSC,1);
CWBILL(1,1)=sum(Cw(Iwest==1).*CL(Iwest==1)); CWBILL(2,1)=sum(Cw(Ieast==1).*CL(Ieast==1)); 

COBSERVE=zeros(NOBSC,5); 
COBSERVE(1,1)=sum(CL(Iwest==1)); COBSERVE(2,1)=sum(CL(Ieast==1));
COBSERVE(1,2)=CWBILL(1,1)./COBSERVE(1,1); COBSERVE(2,2)=CWBILL(2,1)./COBSERVE(2,1);
COBSERVE(1,3)=sum(H(Iwest==1)); COBSERVE(2,3)=sum(H(Ieast==1));
COBSERVE(1,4)=1; COBSERVE(2,4)=0;
COBSERVE(1,5)=0; COBSERVE(2,5)=1;

Cincome=Cw.*CL;
Ctrade=Ctradesh.*repmat(Cincome',NN,1);
CTRADE=zeros(NOBSC,NOBSC);
CTRADE(1,1)=sum(sum(Ctrade(Iwest==1,Iwest==1)));
CTRADE(2,2)=sum(sum(Ctrade(Ieast==1,Ieast==1)));
CTRADE(1,2)=sum(sum(Ctrade(Iwest==1,Ieast==1)));
CTRADE(2,1)=sum(sum(Ctrade(Ieast==1,Iwest==1)));
CEXPEND=sum(CTRADE,2);
CTRADESH=CTRADE./repmat(CEXPEND,1,2);
CDTRADESH=diag(CTRADESH);

% WELFARE GAINS;
[WELFGAIN]=Hwelfaregains(param,CTRADESH,TRADESH,COBSERVE(:,1),OBSERVE(:,1),NOBSC);
display('>>>> Aggregate Country Welfare Gains <<<<');
WELFGAIN=round(WELFGAIN.*(10.^4));
WELFGAIN=WELFGAIN./(10.^4);
unique(WELFGAIN)

% PERFECTLY IMMOBILE WELFARE GAINS (ACR);
[ACRWELFGAIN]=acrwelfaregains(param,CTRADESH,TRADESH,COBSERVE(:,1),OBSERVE(:,1),NOBSC);

% PERFECTLY MOBILE WELFARE GAINS;
[MOBWELFGAIN]=Hmobwelfaregains(param,CTRADESH,TRADESH,COBSERVE(:,1),OBSERVE(:,1),NOBSC);

% COMPARE REGION AND COUNTRY WELFARE GAINS;
display('>>>> Region and Country Welfare Gains <<<<');
display('[Region Country]');
temp=[unique(welfgain(Iwest==1)); unique(welfgain(Ieast==1))];
display([temp WELFGAIN])

% COUNTRY PERFECTLY IMMOBILE WELFARE GAINS;
display('>>>> Country Perfectly Immobile Welfare Gains <<<<');
display('[Region Country]');
display([unique(ACRWELFGAIN)])

% COUNTRY PERFECTLY MOBILE WELFARE GAINS;
display('>>>> Region and Country Perfectly Mobile Welfare Gains <<<<');
display('[Region Country]');
display([unique(MOBWELFGAIN)])

% *************************;
% **** Geometric Means ****;
% *************************;

gdtradesh=zeros(NOBSC,1);
gdtradesh(1,1)=geomean(dtradesh(Iwest==1));
gdtradesh(2,1)=geomean(dtradesh(Ieast==1));

Cgdtradesh=zeros(NOBSC,1);
Cgdtradesh(1,1)=geomean(Cdtradesh(Iwest==1));
Cgdtradesh(2,1)=geomean(Cdtradesh(Ieast==1));

gL=zeros(NOBSC,1);
gL(1,1)=geomean(L(Iwest==1));
gL(2,1)=geomean(L(Ieast==1));

gCL=zeros(NOBSC,1);
gCL(1,1)=geomean(CL(Iwest==1));
gCL(2,1)=geomean(CL(Ieast==1));

test=((gdtradesh./Cgdtradesh).^(alpha./theta)).*((gL./gCL).^((1./epsilon)+(1-alpha)-(alpha./theta)));

display('>>> Compare Aggregate and Geometric Mean Domestic Trade Share');
display('[Aggregate Geometric]');
[1./CDTRADESH gdtradesh./Cgdtradesh]

display('>>> Compare Aggregate and Geometric Mean Labor Supply');
display('[Aggregate Geometric]');
[ones(2,1) gL./gCL]

display('>>>> Compare Region Welfare, Test and Aggregate Welfare <<<<');
display('[region test aggregate]');
temp=[unique(welfgain(Iwest==1)); unique(welfgain(Ieast==1))];
[temp test ACRWELFGAIN]

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

print -dpdf graphs/H_cnty_closed.pdf;

% ******************************************************;
% **** Three-Dimensional Impact of Opening to Trade ****;
% ******************************************************;

dL=CL./L; ldL=log(dL);
dw=Cw./w; ldw=log(dw);
dr=Cr./r; ldr=log(dr);
dP=CP./P; ldP=log(dP);
lacrwelfgain=log(acrwelfgain);
drealw=Crealwage./realwage; ldrealw=log(drealw);
dtradesh=diag(tradesh);
Cdtradesh=diag(Ctradesh);
ddtradesh=Cdtradesh./dtradesh;

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
title('Panel F : Immobile Welfare','FontSize',8);
set(gca,'fontsize',8);

print -dpdf graphs/H_cnty_open.pdf;

figure(4);
wtruew=unique(welfgain(Iwest==1));
wtruee=unique(welfgain(Ieast==1));
bin=linspace(min(acrwelfgain),max(acrwelfgain),20); bin=bin';
bbin=round(bin*100); bbin=bbin/100;
wbin=abs(bin-wtruew);
ebin=abs(bin-wtruee);
fwbin=find(wbin==min(wbin));
febin=find(ebin==min(ebin));
h=histc(acrwelfgain,bin);
h=h./sum(h);
h1=histc(acrwelfgain(Iwest==1),bin);
h1=h1./sum(h1);
h0=histc(acrwelfgain(Ieast==1),bin);
h0=h0./sum(h0);
subplot(1,2,1);
bar(h1);
hold on;
ylabel('Probability','FontSize',8);
title('Panel A : West','FontSize',8);
axis([0 20 0 0.6]);
set(gca,'fontsize',8);
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
line([fwbin fwbin],[0 0.6],'LineStyle','-','Color','red','LineWidth',1.5);
hold off
subplot(1,2,2);
bar(h0);
hold on;
ylabel('Probability','FontSize',8);
title('Panel B : East','FontSize',8);
axis([0 20 0 0.6]);
set(gca,'fontsize',8);
set(gca, 'XTickLabel',[bbin(1) bbin(5) bbin(10) bbin(15) bbin(20)],'XTick',[1 5 10 15 20]);
line([febin febin],[0 0.6],'LineStyle','-','Color','red','LineWidth',1.5);
hold off;

print -dpdf graphs/H_cnty_welfare_impact.pdf;
