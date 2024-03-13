
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

% User 1 : steve mac pro dropbox; 
% User 2 : steve mac laptop dropbox;
% User 3 : steve home mac;

if user==1;
cd C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\regions;
addpath C:\Users\reddings\Dropbox\QuantSpatialModel\JIErev1\matlab\regions\functions;
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
epsilon=3.1;

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

% Perfectly mobile welfare gains;
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

% *****************;
% **** Impacts ****;
% *****************;

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

% Wage treatment;
[bw,bintw,rw,rintw,statsw] = regress(ldw,X);
display('>>>> Wage Treatment <<<<');
bw(2)
stats = regstats(ldw,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Price index treatment;
[bP,bintP,rP,rintP,statsP] = regress(ldP,X);
display('>>>> Price Index Treatment <<<<');
bP(2)
stats = regstats(ldP,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Land price treatment;
[br,bintr,rr,rintr,statsr] = regress(ldr,X);
display('>>>> Land Price Treatment <<<<');
br(2)
stats = regstats(ldr,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Real wage treatment;
[brw,bintrw,rrw,rintrw,statsrw] = regress(ldrealw,X);
display('>>>> Real Wage Treatment <<<<');
brw(2)
stats = regstats(ldrealw,treat,'linear',{'beta','covb'});
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

% Wage treatment;
[bw,bintw,rw,rintw,statsw] = regress(lHdw,X);
display('>>>> Wage Treatment <<<<');
bw(2)
stats = regstats(lHdw,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Price index treatment;
[bP,bintP,rP,rintP,statsP] = regress(lHdP,X);
display('>>>> Price Index Treatment <<<<');
bP(2)
stats = regstats(lHdP,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Land price treatment;
[br,bintr,rr,rintr,statsr] = regress(lHdr,X);
display('>>>> Land Price Treatment <<<<');
br(2)
stats = regstats(lHdr,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% Real wage treatment;
[brw,bintrw,rrw,rintrw,statsrw] = regress(lHdrealw,X);
display('>>>> Real Wage Treatment <<<<');
brw(2)
stats = regstats(lHdrealw,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

% ACR welfare
[bacr,bintacr,racr,rintacr,statsacr] = regress(lHacrwelfgain,X);
display('>>>> ACR Welfare Treatment <<<<');
bacr(2)
stats = regstats(lHacrwelfgain,treat,'linear',{'beta','covb'});
stats.beta
sqrt(diag(stats.covb))

