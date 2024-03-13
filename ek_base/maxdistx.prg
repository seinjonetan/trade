/* maxdistx.prg */

/* Program to look at maximum price ratios as a proxy for trade */
/*   frictions.  Based on SIMPARM.PRG  */
/* Sam Kortum 7/20/98 */

/* This program generates:

	normalized import shares, (Xni/Xn)/(Xii/Xi), defined on p.12
	Dni, defined on p.14
	
	method of moments estimate of theta, p.15
	OLS estimates of theta (with & without constant), footnote 26, p.15
	unrestricted pure price regressions, footnote 26, p.16
	OLS estimate of theta, with source & destination dummies, p.25
	2SLS estimate of theta with geography instruments, p. 25
	
	table II
	figure II */

new;
cls;
output file = c:\research\eaton\tgt\maxdistx.out reset;

s = 0;

ncnty = 19;         @ we consider 19 countries @
nyear=1;
n1=6860; nn=7220;   @ 1990 data @

"  ***** ";
" ONLY USES 1990 DATA ";
" INCLUDES TALK, EU, and EFTA ";
" ***** ";

loadm trade1 = c:\research\eaton\tgt\trade1;
trade1 = trade1[n1:nn,.];

impcnty = trade1[.,1];      @ importing country index @
expcnty = trade1[.,2];      @ exporting country index @
depvar = trade1[.,4];
depvarn = trade1[.,5];
dist1 = trade1[.,6];
dist2 = trade1[.,7];
dist3 = trade1[.,8];
dist4 = trade1[.,9];
dist5 = trade1[.,10];
dist6 = trade1[.,11];
border = trade1[.,12];
lrrnd = trade1[.,13];
lrhk = trade1[.,14];
lrwork = trade1[.,16];
lrwage = trade1[.,17];

mybeta = .21221;
diff = -(depvarn-depvar)*mybeta/(1-mybeta);
depvarp = depvar + diff;    @ depvarp is ln(Xni/Xn)-ln(Xii/Xi), as defined on p.12 @


@ Create common language variable @

@          AL AU BE CA DE FI FR GE GR IT JA NE NZ NO PO SP SW UK US   @
english = {1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1};
french =  {0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
german =  {0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
talk = vec(english*english' + french*french' + german*german');
talk = ones(20,1).*.talk;

@ Create country trade blocks @

@          AL AU BE CA DE FI FR GE GR IT JA NE NZ NO PO SP SW UK US   @
eu71 =    {0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0};
eu73 =    {0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0};
eu81 =    {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
eu86 =    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
eu7172 = ones(2,1).*.vec(eu71*eu71');
eu7380 = ones(8,1).*.vec((eu71+eu73)*(eu71+eu73)');
eu8185 = ones(5,1).*.vec((eu71+eu73+eu81)*(eu71+eu73+eu81)');
eu8690 = ones(5,1).*.vec((eu71+eu73+eu81+eu86)*(eu71+eu73+eu81+eu86)');
eu = eu7172|eu7380|eu8185|eu8690;

@          AL AU BE CA DE FI FR GE GR IT JA NE NZ NO PO SP SW UK US   @
ef71 =    {0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0};
ef73 =    {0, 0, 0, 0,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0};
ef86 =    {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1, 0, 0, 0};
ef7172 = ones(2,1).*.vec(ef71*ef71');
 ef7385 = ones(13,1).*.vec((ef71+ef73)*(ef71+ef73)');
ef8690 = ones(5,1).*.vec((ef71+ef73+ef86)*(ef71+ef73+ef86)');
ef = ef7172|ef7385|ef8690;

tblock = eu~ef;

tblock = tblock[n1:nn,.];
talk   = talk[n1:nn,.];

@Create the source and destination dummy variables@
onecnty = ones(ncnty,1);
desdum = ones(nyear,1).*.(eye(ncnty).*.onecnty);
srcdum = ones(nyear,1).*.(onecnty.*.eye(ncnty));
reldum = srcdum-desdum;

desdumr = desdum[.,1:(ncnty-1)]-desdum[.,ncnty];
reldumr = reldum[.,1:(ncnty-1)]-reldum[.,ncnty];

@proc to calculate country pair time averages@
proc order2(datf);
   local avar,i,var,varmat,adatf;
     adatf = ones((ncnty*(ncnty-1)),cols(datf));
     i = 1;
     do while i <= cols(datf);
        var = datf[.,i];
        varmat = reshape(var,nyear,(ncnty*(ncnty-1)));
        adatf[.,i] = meanc(varmat);
        i = i+1;
     endo;
   retp(adatf);
endp;

@ Procedure to do gls @
proc(3)=gls(y,x,icov,dof);
    local betahat,eps,sighat2,varcov,stdb;
    betahat = inv(x'*icov*x)*x'*icov*y;
    eps = y - x*betahat;
    print "Raw SS in this regression " eps'*eps;
    print "Adjusted SS " eps'*icov*eps;
    varcov = inv(x'*icov*x);
    /* Next, we get the variance for the usa in GLS */
    print "source effect: usa=";
    -ones(1,18)*betahat[11:28];
    print "source effect: stdb(usa)=";
    sqrt(ones(1,18)*varcov[11:28,11:28]*ones(1,18)');
    print "destination effect: usa=";
    -ones(1,18)*betahat[29:46];
    print "destination effect: stdb(usa)=";
    sqrt(ones(1,18)*varcov[29:46,29:46]*ones(1,18)');
    stdb = sqrt(diag(varcov));
    retp(betahat,stdb,eps);
endp;

@ Procedure to construct the correlation between dresid and usa @
proc correl(resid,b,axmatf);
   local allresid,imp2,exp2,datloc,residmat,sresid,dresid,beta,
            temp1,temp2,temp3,usa;
   @ average the residuals across sources and destinations @
    allresid = zeros((ncnty^2),1);
    imp2 = seqa(1,1,ncnty).*.ones(ncnty,1);
    exp2 = ones(ncnty,1).*.seqa(1,1,ncnty);
    datloc = indexcat((imp2.==exp2),0);
    allresid[datloc] = resid;
    residmat = reshape(allresid,ncnty,ncnty);
    sresid = sumc(residmat)/(ncnty-1);
    dresid = sumc(residmat')/(ncnty-1);
   @ calculate lnmuhat-thetahat*lnwage @
    beta = b[(cols(axmatf)-2):cols(axmatf),1];
    temp1 = delif(lrrnd~lrhk~lrwage,(impcnty.eq expcnty));
    temp2 = order2(temp1);
    temp3 = temp2*beta;
    usa = zeros(ncnty,1);
    usa[1:ncnty-1,.] = temp3[ncnty*(ncnty-1)-ncnty+2:ncnty*(ncnty-1),1];
   @ calculate the correlation between dresid and usa @
    print "correlation between dresid and usa";
    corrx(dresid~usa);
    print "dresid~usa";
  retp(dresid~usa);
endp;

@ Proc to do TSLS given  @
@ returns coefficients in column 1 and std. errors in column 2 @
@ remember that x must include the constant term if you want it @
@ You supply appropriate degrees of freedom @

 proc(3)= tsls(y,x,xhat,dof);
    local betahat,eps,sighat2,varcov,stdb;
    betahat = inv(xhat'*xhat)*xhat'*y;
    eps = y - x*betahat;
    sighat2 = eps'*eps/dof;
    print "SS in this regression" sighat2*dof;
    varcov = sighat2*inv(xhat'*xhat);
    stdb = sqrt(diag(varcov));
    retp(betahat,stdb,eps);
 endp;


@ ******* Get pure distances ******************************** @

@ Get distances between countries in thousands of miles @

load milevec[361,2] = c:\research\eaton\tgt\getdist.dat;

milevec = milevec[.,2];

@ correct mistakes @
milevec[18] = 10.769;
milevec[324] = 10.769;


@ *********** get the price data **************************** @
load pdat[900,1] = c:\research\eaton\tgt\pppdat1.prn;
pdat = pdat|zeros(50,1);   @ tack on the U.S. @
pdat = reshape(pdat,19,50);  @ countries in rows, goods in columns @

@ *** eliminate office and computing *** @
@ pdat = pdat[.,1:16]~pdat[.,18:cols(pdat)]; @


mpdat = meanc(pdat');
dlpavg = mpdat.*.ones(ncnty,1);
slpavg = ones(ncnty,1).*.mpdat;


dpdat = pdat.*.ones(ncnty,1);
spdat = ones(ncnty,1).*.pdat;
@ create log importer price - log exporter price for each item @
rpdat = dpdat - spdat;

ldni1 = maxc(rpdat');
max1 = maxindc(rpdat');

/*The following loops calculate the second, third, etc. highest
values of r_ni.  We use the second highest value, as described on page 14.*/

i = 1;
 do while i <= rows(rpdat);
    bigcol = max1[i];
    rpdat[i,bigcol] = -100;
   i = i+1;
 endo;

ldni2 = maxc(rpdat');
max2 = maxindc(rpdat');

i = 1;
 do while i <= rows(rpdat);
    bigcol = max2[i];
    rpdat[i,bigcol] = -100;
   i = i+1;
 endo;

ldni3 = maxc(rpdat');
max3 = maxindc(rpdat');

i = 1;
 do while i <= rows(rpdat);
    bigcol = max3[i];
    rpdat[i,bigcol] = -100;
   i = i+1;
 endo;

ldni4 = maxc(rpdat');
max4 = maxindc(rpdat');


@ ************************************************************* @

dist = dist2~dist3~dist4~dist5~dist6~border;
dist = dist~talk~tblock;


ldni = ldni2;   @ choose the order statistic @


xmat = slpavg - dlpavg + ldni;




@ want to save matrix for a spreadsheet @
@ destinations are rows, sources are columns @

@
temp = reshape(xmat,19,19);

temp[.,1:4];
" ";
temp[.,5:8];
" ";
temp[.,9:12];
" ";
temp[.,13:16];
" ";
temp[.,17:19];
" ";
@

eupair = vec((eu71+eu73+eu81+eu86)*(eu71+eu73+eu81+eu86)');



@ remove home countries @

depvarpf = delif(depvarp,(impcnty .eq expcnty));
depvarnf = delif(depvarn,(impcnty .eq expcnty));

xmatf = delif(xmat,(impcnty .eq expcnty));

milef = delif(milevec,(impcnty .eq expcnty));

lmilef = ln(milef);



" mean(trade) ";
meanc(depvarpf);
" ";
" mean(price) ";
meanc(xmatf);
" ";
" method of moments theta (see p.15) ";
meanc(depvarpf)/meanc(xmatf);  @ method of moments theta, as described on p.15 @

" ";
" correlation ";
corrx((depvarpf~xmatf));

" ";
"******************************************************************************";
"Countries are: (1) Australia,   (2) Austria,  (3) Belgium,";
"               (4) Canada,      (5) Denmark,  (6) Finland,";
"               (7) France,      (8) Germany,  (9) Greece,";
"              (10) Italy,      (11) Japan,   (12) Netherlands,";
"              (13) New Zealand (14) Norway   (15) Portugal,";
"              (16) Spain       (17) Sweden,  (18) UK,";
"              (19) USA.";
"******************************************************************************";
" ";

@ make index vectors @
dest = delif(impcnty, impcnty .eq expcnty);
srce = delif(expcnty, impcnty .eq expcnty);

@ convert Dni to exp(Dni) as in table II @
expxmatf=exp(xmatf);

"Table II";
" ";
x=expxmatf~srce~dest;
z=ones(19,8);
i=1;
do while i <= 19;
	y=selif(x,x[.,3] .eq i);
	w=y[.,1];
	z[i,1]=y[minindc(w),2];
	z[i,2]=minc(w);
	z[i,3]=y[maxindc(w),2];
	z[i,4]=maxc(w);
    i=i+1;
endo;

i=1;
do while i <= 19;
	y=selif(x,x[.,2] .eq i);
	w=y[.,1];
	z[i,5]=y[minindc(w),3];
	z[i,6]=minc(w);
	z[i,7]=y[maxindc(w),3];
	z[i,8]=maxc(w);
	i=i+1;
endo;

format 4,3;	



$" "~"srce"~" "~" "~" "~"dest"~" "~" ";
$" "~"min"~" "~"max"~" "~"min"~" "~"max";
z;
" ";
"******************************************************************************";

" ";
" ";
" ";
" ";
" linear regression through the scatter in figure 2, no constant ";
"   (see footnote 26, page 15)";
" ";
" dependent variable is normalized import share: ln(Xni/Xn) - ln(Xii/Xi) ";
" independent variable is price measure: Dni ";
" ";
" ";
@" pure price regression (no constant) ";@
" ";
__con=0;
call ols(0,depvarpf,xmatf);
" ";
"************************************************************************";

" ";
" ";
" ";
" ";
" ";
" ";
" linear regression through the scatter in figure 2, with constant ";
"   (see footnote 26, page 14) ";
" ";
" dependent variable is normalized import share: ln(Xni/Xn) - ln(Xii/Xi) ";
" independent variable is price measure: Dni ";
" ";
@" pure price regression (with constant) ";@
" ";
__con=1;
call ols(0,depvarpf,xmatf);

" ";
"************************************************************************";



" ";
" Instrument for the price variable (using actual distance)";

_olsres=1;

__output=0;

{x1,x2,x3,x4,x5,x6,x7,x8,x9,myres,x11} = ols(0,xmatf,lmilef);
myfit = xmatf-myres;
" ";
" myres ";
" ";
" second stage ";

__output=0;

call ols(0,depvarpf,myfit);

" ";

@ do the standard errors right @
x = ones(rows(myfit),1)~xmatf;
xhat = ones(rows(myfit),1)~myfit;
dof = rows(x) - 1;
{betahat,stdb,eps} = tsls(depvarpf,x,xhat,dof);
(betahat~stdb);

@ include source and destination dummies @;

dvars = desdumr~reldumr;
dvars = delif(dvars,(impcnty .eq expcnty));
dvars = lmilef~dvars;

_olsres=1;
__output=0;

{x1,x2,x3,x4,x5,x6,x7,x8,x9,myres,x11} = ols(0,xmatf,dvars);
myfit = xmatf-myres;
" ";
@ myres @
" ";
@ second stage @

__output=0;

call ols(0,depvarpf,myfit);

" ";
@ do the standard errors right @
x = ones(rows(myfit),1)~xmatf;
xhat = ones(rows(myfit),1)~myfit;
dof = rows(x) - 1;
{betahat,stdb,eps} = tsls(depvarpf,x,xhat,dof);
(betahat~stdb);

__output=1;

" ";
" ";
" ";
" ";
"*******************************************************************";
" ";
" ";
" include source and destination dummies ";
" dependent variable is now depvarnf ";
" ";
" 2SLS regression of the price measure on normalized trade shares "; 
"   with observed geography terms as instruments for price measure, Dni ";
"   (see section 5.3, page 24) ";
" ";

dumvars = desdumr~reldumr;
dumvars = delif(dumvars,(impcnty .eq expcnty));

dvars = dist~desdumr~reldumr;
dvars = delif(dvars,(impcnty .eq expcnty)); 

_olsres=1;

{x1,x2,x3,x4,x5,x6,x7,x8,x9,myres,x11} = ols(0,xmatf,dvars);
myfit = xmatf-myres;
" ";
" myres ";
" "
" second stage ";
"    ";
myfit = myfit~dumvars;
call ols(0,depvarnf,myfit);

" ";
" do the standard errors right ";
x = ones(rows(myfit),1)~xmatf~dumvars;
xhat = ones(rows(myfit),1)~myfit;
dof = rows(x) - 1;
{betahat,stdb,eps} = tsls(depvarnf,x,xhat,dof);
(betahat~stdb);
" ";
" ";
" -----------------------------------------------------------------";
" 2SLS estimate of theta: " betahat[2] "          standard error: " stdb[2];
" *****************************************************************************";
" ";
" ";

" ";
"  Do OLS but include source and destination dummies ";
" ";

dvars = xmatf~dumvars;

call ols(0,depvarnf,dvars);




xmat2 = slpavg~dlpavg~ldni;  @ replaced xmat @

@ remove home countries @

xmat2f = delif(xmat2,(impcnty .eq expcnty));  @ replaced xmatf @

" ";
"********************************";
" ";
" ";
" OLS regression of trade shares on log of exporter price (pi), log of importer price (pn), ";
" and log of geographic barriers (dni), no constant  (see footnote 26, p.16)";
" ";
__con=0;
call ols(0,depvarpf,xmat2f);
" ";
"*******************************"
" ";
" ";
" OLS regression of trade shares on log of exporter price (pi), log of importer price (pn), ";
" and log of geographic barriers (dni), with constant  (see footnote 26, p.16)";
" ";
__con=1;
call ols(0,depvarpf,xmat2f);





end;
