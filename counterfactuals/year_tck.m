% This code provides the basic structure for running counterfactuals. You
% need to have equilibrium.m and model_iteration.m saved in directory
% before begining.

clear

% Set directory
% cd counterfactuals/

% Set filename
filename = "est_master.csv";

% Set parameters
% From Burstein et al (2022, ECMA):
eta = 1.85;
alpha = 7;
% From Table 2 of our draft:
theta = 3.475;
rho = 0.562;

% Read the data
data_master = readtable(filename);

% get unique years
years = unique(table2array(data_master(:,"year")));

for i = 1:length(years)
    fprintf("Processing year %d\n", years(i));
    fprintf("Number of years left: %d\n", length(years) - i);

    % Keep only data for the current year
    data = data_master(table2array(data_master(:,"year")) == years(i), :);

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

    % Calculate initial technology matrix.
    [T] = technology(employment,wagebill,alpha,theta,eta,rho);

    % write the technology matrix to a file
    filename = sprintf('tck/tck_%d.csv', years(i));
    writetable(array2table(T), filename);
end