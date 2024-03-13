%%% Solve model

function [acrwelfgain] = acrwelfaregains(param,Ctradesh,tradesh,CL,L,nobs)

global alpha sigma theta epsilon LL;

xtic = tic();

% parameters;
alpha=param(1);
theta=param(2);
epsilon=param(3);

% domestic trade share;
dtradesh=diag(tradesh);
Cdtradesh=diag(Ctradesh);

% welfare gains;
acrwelfgain=(dtradesh./Cdtradesh).^(alpha./theta);
