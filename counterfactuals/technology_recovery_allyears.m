%This script reads the raw data from "est_master.csv" and generates the
%technology matrix from the raw data using the function "technology.m". I
%do this for all years and for two different parameter assumptions: the
%true theta and rho we estimate, and an estimate under the assumption that
%rho = 0 and therefore theta = theta/(1-rho). 
clear

% Set directory
cd C:\Users\tsjac\Desktop\research\correlated_location

% Set filename
filename = "est_master.csv";

% Set parameters
% From Burstein et al (2022, ECMA):
eta = 1.85; 
alpha = 7; 
% From Table 2 of our draft:
theta = 3.475; 
rho = 0.562;
% Generate no correlation equivalent.
theta_norho = theta / (1-rho);
rho_norho = 0;

for y = [1990,2000,2007,2019]
disp(y)

% Read the data
data = readtable(filename);
% Keep only 1990 data
data(table2array(data(:,"year") ~= y), :) = [];

% Turn string columns into index columns and drop unecessary columns
data(:,"name") = [];
data(:,"year") = [];

colName = {'c'};
c = array2table(findgroups(data(:,"comzone")),'VariableNames',colName);
comzone = (unique(data(:,"comzone")));
comzone = table2array(comzone);

colName = {'k'};
k = array2table(findgroups(data(:,"occupation")),'VariableNames',colName);
occupation = array2table(unique(data(:,"occupation")));

% Generate new data
data = [c k data];

% Generate matrix with rows c and columns k
N = table2array(max(data(:,"c")));
K = table2array(max(data(:,"k")));

employment = zeros(N,K);
wagebill = zeros(N,K);

for c = 1:N 
    for k = 1:K
        temp = data(table2array(data(:,"c") == c), :);
        temp = temp(table2array(temp(:,"k") == k), :);
        employment(c,k) = table2array(temp(1,"employed"));
        wagebill(c,k) = table2array(temp(1,"wage_bill"));
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
L = sum(sum(employment));
employment = employment / L;

% Use the technology function to recover the initial matrix of technology
% across city-occupation pairs.
[T] = technology(employment,wagebill,alpha,theta,eta,rho);
[T_norho] = technology(employment,wagebill,alpha,theta_norho,eta,rho_norho);

year = ones(N,1)*y;
city_id = linspace(1,c,c)';
tech = [city_id comzone year T];
tech_norho = [city_id comzone year T_norho];

if y == 1990 
    Tck = tech;
    Tck_norho = tech_norho;
else
    Tck = [Tck; tech];
    Tck_norho = [Tck_norho; tech_norho];
end
end
