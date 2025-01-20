% Jiaqi Wu @ Arizona State University

clc;
clear;
rng(0);

% Load the case
%mpc = matpower_MV_rural();
%case_name = 'mv_rural';
%scaler = 0.05;

%mpc = matpower_LV_urban();
%case_name = 'lv_urban';
%scaler = 0.001;

%mpc = case33mg();
%case_name = 'case33';
%scaler = 0.05;

%mpc = case34sa();
%case_name = 'case34';
%scaler = 0.05;

mpc = Node8();
case_name = 'case8';
scaler = 0.05;


% Define file name
file_name = sprintf('data/%s_results.csv', case_name);
edge_name = sprintf('data/%s_branches.csv', case_name);
g_name = sprintf('data/%s_g.csv', case_name);
b_name = sprintf('data/%s_b.csv', case_name);

result = runpf(mpc);
Ybus_sparse = makeYbus(mpc);
Ybus = full(Ybus_sparse);

% Get the base case data
base_V = result.bus(:, 8)';

% Get number of bus/gen
num_gens = size(mpc.gen, 1);
num_buses = size(mpc.bus, 1);
num_samples = 2976;

branch_set = mpc.branch(:, 1:2);
writematrix(branch_set, edge_name);
writematrix(real(Ybus), g_name);
writematrix(imag(Ybus), b_name);


% Load p data
raw_data = readtable('data/sorted_profile_data.csv');
%profile_data = raw_data(1:213, 4:2979);
profile_data = raw_data;
profile_array = table2array(profile_data);


%data_p = 4 * 0.001 * profile_array(1:num_buses, :);
data_p = 4 * scaler * profile_array(1:num_buses, :);


% Generate q data  
data_pf = 0.1 * rand(num_buses, num_samples) + 0.85;
data_s = data_p ./ data_pf;
data_q = sqrt(data_s.^2 - data_p.^2);



% Preallocate storage for the results
results_data = struct();
results_data.load_P = zeros(num_samples, num_buses);
results_data.load_Q = zeros(num_samples, num_buses);
results_data.gen_P = zeros(num_samples, num_gens);
results_data.gen_Q = zeros(num_samples, num_gens);
results_data.voltage = zeros(num_samples, num_buses);
results_data.theta = zeros(num_samples, num_buses);



% Iterations
for i = 1:num_samples

    % Update the real and reactive power
    for n = 1:num_buses
        mpc.bus(n, 3) = data_p(n, i);
        mpc.bus(n, 4) = data_q(n, i);
    end
    
    % Run the power flow
    result = runpf(mpc);
    
    % Store the results
    results_data.load_P(i, :) = result.bus(:, 3);
    results_data.load_Q(i, :) = result.bus(:, 4);
    results_data.gen_P(i, :) = result.gen(:, 2);
    results_data.gen_Q(i, :) = result.gen(:, 3);
    results_data.voltage(i, :) = result.bus(:, 8);
    results_data.theta(i, :) = result.bus(:, 9);

end



% Flatten the results for the CSV
num_gens = size(mpc.gen, 1);
num_buses = size(mpc.bus, 1);

iteration = 2976;
load_P_flat = reshape(results_data.load_P, [], num_buses);
load_Q_flat = reshape(results_data.load_Q, [], num_buses);
gen_P_flat = reshape(results_data.gen_P, [], num_gens);
gen_Q_flat = reshape(results_data.gen_Q, [], num_gens);
voltage_flat = reshape(results_data.voltage, [], num_buses);
theta_flat = reshape(results_data.theta, [], num_buses);


% Create a table with just the iteration
num_samples = size(load_P_flat, 1);
T = table((1:num_samples)', 'VariableNames', {'Iteration'});

% Add load_P and load_Q data for each bus
for k = 1:num_buses
    bus_idx = mpc.bus(k, 1);
    T.(['load_P_bus_' num2str(bus_idx)]) = load_P_flat(:, k);
    T.(['load_Q_bus_' num2str(bus_idx)]) = load_Q_flat(:, k);
end

% Add generation data for each generator
for k = 1:num_gens
    gen_idx = mpc.gen(k, 1);
    T.(['gen_P_bus_' num2str(gen_idx)]) = gen_P_flat(:, k);
    T.(['gen_Q_bus_' num2str(gen_idx)]) = gen_Q_flat(:, k);
end

% Add voltage and theta data for each bus
for k = 1:num_buses
    bus_idx = mpc.bus(k, 1);
    T.(['voltage_bus_' num2str(bus_idx)]) = voltage_flat(:, k);
    T.(['theta_bus_' num2str(bus_idx)]) = theta_flat(:, k);
end

% Write the table to a CSV file
writetable(T, file_name);


