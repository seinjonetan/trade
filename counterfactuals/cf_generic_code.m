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
% rho = 0.562;
rho = 0;

% Read the data
data = readtable(filename);
% Keep only 1990 data
data(table2array(data(:,"year") ~= 1990), :) = [];

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

% We need an initial guess of wages to then solve for the true wages
% associated with the existing data. Just guess that technology level is
% correlated wages and guess T as initial wage vector.
w_init = T;
[w] = equilibrium(w_init,T,alpha,theta,eta,rho,1e-6,1e-6);
[pick_0] = get_shares(w,theta,rho);

% This is a shock to city number 541, which I think is the largest city,
% with a negative 10% shock to all technology across all occupations for
% this city specifically. Can then solve for wages. Notice that with wages
% we can then construct everything else. So wages are the fundamental
% variable.
city_employment_ranking = sum(employment,2);
[~, city_rank] = sort(city_employment_ranking,'descend');

shockvalue = 0.2;
T_cf = T;
pic_0 = sum(pick_0,2);

% Apply shock to every city, starting with the largest city using a loop and store results
% for i = 1:length(city_rank)
%     T_cf(city_rank(i),:) = T_cf(city_rank(i),:)*(1-shockvalue);
%     [w_cf] = equilibrium(w,T_cf,alpha,theta,eta,rho,1e-6,1e-6);
%     [pi_cf] = get_shares(w_cf,theta,rho);
%     [change_cf] = pi_cf ./ pi_0;
% end

% Initialize storage for all results
num_cities = length(city_rank);
% num_cities = 3;  % For testing, use only first 3 cities
num_locations = N;
num_occupations = K;

% Create storage matrix: [change_cf(:), shocked_city, num_cities_shocked, city_id, occupation_id]
all_results = [];

% Create progress bar
h = waitbar(0, 'Running counterfactual simulations...', 'Name', 'Counterfactual Progress');

% Apply shock to every city, starting with the largest city using a loop and store results
for i = 1:num_cities
    waitbar(i/3, h, sprintf('Processing city %d of %d (City ID: %d)', i, num_cities, city_rank(i)));

    T_cf(city_rank(i),1) = T_cf(city_rank(i),1)*(1-shockvalue);
    [w_cf] = equilibrium(w,T_cf,alpha,theta,eta,rho,1e-6,1e-6);
    [pick_cf] = get_shares(w_cf,theta,rho);
    pic_cf = sum(pick_cf,2);
    [change_cf] = pic_cf ./ pic_0;  % N×1 vector

    % Create metadata for city-level results
    num_obs = length(change_cf);  % N cities

    % Create metadata columns
    shocked_city_indicator = ((1:num_obs)' == city_rank(i));
    shocked_city_id = repmat(city_rank(i), num_obs, 1);
    num_cities_shocked = repmat(i, num_obs, 1);
    city_id_col = (1:num_obs)';

    % Combine results (no occupation dimension)
    iteration_results = [change_cf, shocked_city_indicator, shocked_city_id, ...
        num_cities_shocked, city_id_col];

    all_results = [all_results; iteration_results];
end

% Update column names
result_table = array2table(all_results, 'VariableNames', ...
    {'change_ratio', 'is_shocked_city', 'shocked_city_id', 'num_cities_shocked', 'city_id'});

% Close progress bar
close(h);

writetable(result_table, 'shocks/c_shock_rho0.csv')

% T_kf = T;
% T_kf(:, 1) = T_kf(:, 1)*(1-shockvalue);
% [w_kf] = equilibrium(w,T_kf,alpha,theta,eta,rho,1e-6,1e-6);
% [pi_kf] = get_shares(w_kf,theta,rho);
% [change_kf] = pi_kf ./ pi_0;