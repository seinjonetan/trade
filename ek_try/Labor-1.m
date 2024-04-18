clear all;

%%%Inputing Structural parameters
theta=4; %%%Theta as suggested by Simonovska and Waugh (2014).
sigma=2;    %%%Elasticity of substitution in consumption. Unimportant in the estimation. Only used when calculating the price index.
gamma=(gamma((theta+1-sigma)/theta))^(1/(1-sigma)); %%%An exogenous parameter. See equation 9 of EK.
bump=0.05; %%%Controls the speed at which the iterative process moves.
N=40; %%%Number of Countries

%%%For the iterative process set a tolerance level and a starting value of
%%%the loss function making sure that it is greater than the tolerable
%%%loss.
tol=0.000000001;
dev=5;

%%%Inputing exogenous technology levels
T=csvread('T.csv'); %%%%States of Technology from EK (2002).
T=T./T(N);
L=csvread('L.csv'); %%%%Employment Levels from EK (2002).
L=L./L(N);

%%%Distance matrix (rows are sources, columns are destinations)
d=csvread('d.csv');

%%%Inputing a candidate wage vector
w_start=ones(N,1);

while dev>tol
    
    %%%Taking the cost of an input bundle to be equal to the wage.
    c=w_start;

    %%%Setting a NxN matrix of fob ("free-on-board") prices.
    p=((c.^(-theta)).*T)*ones(1,N);
    
    %%%Constructing the matrix of deliverable prices by multiplying the fob
    %%%price matrix element-by-element by iceberg costs.
    p=d.*p;
     
    %%%Constucting the Phi term for each destination country. This
    %%%corresponds to Equation (7) in EK (2002).
    %%%
    Phi=sum(p);
    
    %%%%Setting up the matrix of expenditure shares for each country n for
    %%%%goods from i. The elements of pi corrrespond to equation (8) in EK
    %%%%(2002)
    Phi_mtx=ones(N,1)*Phi;
    pi=p./Phi_mtx;
    
    
    RHS=pi*(w_start.*L);%%%Calculating the RHS of equation (18) in EK (2002). This is labor income.
    LHS=w_start.*L; %%%Calculating the LHS of equation (18) in EK (2002). This is labor expenditure.
    dev=sum(abs(LHS-RHS)); %%%Calculating the loss function deviation between the candidate and implied wages.
    w_start=(w_start-bump*(LHS-RHS));   %%%Adjusting the candidate vector to get closer to the implied vector.
    w_start=w_start/w_start(N); %%%Keeping normalization of one wage to 1. This "anchors" everything.
    
end

w_n=w_start;
P_n=gamma*(Phi.^(-1/theta));
real_w=w_n./P_n;

%%%CHECKING MY WORK
c=w_n;
p=((c.^(-theta)).*T)*ones(1,N);
p=d.*p;
Phi=sum(p);
Phi_mtx=ones(N,1)*Phi;
pi=p./Phi_mtx;
RHS=pi*(w_n.*L);
LHS=w_n.*L;
check_wages=sum(abs(LHS-RHS))







