%This script reads the raw data from "master_subset.csv" and generates the
%technology matrix from the raw data using the function "technology.m".
clear

% Set directory
cd C:\Users\tsjac\Desktop\research\correlated_location

% Set filename
filename = "C:\Users\tsjac\Downloads\master_subset.csv";

% Read the data
data = readtable(filename);
% Keep only 1990 data
data(table2array(data(:,"Year") ~= 1990), :) = [];

% Turn string columns into index columns and drop unecessary columns
data(:,"NAME") = [];
data(:,"Year") = [];

colName = {'c'};
c = array2table(findgroups(data(:,"COMZONE")),'VariableNames',colName);
data(:,"COMZONE") = [];

colName = {'k'};
k = array2table(findgroups(data(:,"Occupation")),'VariableNames',colName);
data(:,"Occupation") = [];

data = [c k data];

%%
% Generate matrix with rows c and columns k
N = table2array(max(data(:,"c")));
K = table2array(max(data(:,"k")));

employment = zeros(N,K);
wagebill = zeros(N,K);

for c = 1:N 
    for k = 1:K
        temp = data(table2array(data(:,"c") == c), :);
        temp = temp(table2array(temp(:,"k") == k), :);
        employment(c,k) = table2array(temp(1,"Employed"));
        wagebill(c,k) = table2array(temp(1,"Wage_Bill"));
    end
end

% Clean employment and wagebill data: replace all zeros or mising values
% with employment of one. Can check sensitivity to these assumptions later.
mask = (employment == 0);
employment(mask) = 1;
mask = isnan(employment);
employment(mask) = 1;

mask = (wagebill == 0);
wagebill(mask) = 1;
mask = isnan(wagebill);
wagebill(mask) = 1;

% Turn raw employment data into choice shares.
L = max(max(employment));
employment = employment / L;

% Set parameters
theta = 2; 
eta = 1.85; 
alpha = 7; 
rho = 0.25;

% Use the technology function to recover the initial matrix of technology
% across city-occupation pairs.
[T_init] = technology(employment,wagebill,alpha,theta,eta,rho);

%%
% Run a counterfactual simulation by shocking all technology for occupation
% k=1 and increase by 10%
tic
% Set up counterfactual technology matrix
T_cf = T_init;
T_cf(:,1) = T_cf(:,1)*1.10;

% Set up initial wage guess (all one) and tolerance levels for solver.
w_guess = ones(N,K)/10;
tolfun = 1e-9;
tolx = 1e-9;

[w_init] = equilibrium(w_guess,T_init,alpha,theta,eta,rho,tolfun,tolx);
[w_cf] = equilibrium(w_init,T_cf,alpha,theta,eta,rho,tolfun,tolx);
toc