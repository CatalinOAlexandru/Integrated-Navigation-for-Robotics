% Long format for acurate printing and defining the included cosntants.
format long;
Define_Constants;

% This is the start of the initialisation of the starting postion.
% It is the first scrip which will run.

% We are declaring the number of satellites as seen in the
% Pseudo_ranges.csv and initialise the variables by storing the csv's data.
num_of_satellites = 8;
pseudo_ranges = readmatrix('Pseudo_ranges.csv');
pseudo_range_rates = readmatrix('Pseudo_range_rates.csv');
sattelite_numbers = pseudo_ranges(1,2:9);
satellite_indices = zeros(1,30);

% Reverse mapping for outlier detection
for i=1:num_of_satellites
    satellite_indices(sattelite_numbers(i)) = i;
end

num_rows = size(pseudo_ranges,1);

% The current location set as zeros and previous locations as inifinte just
% as a starting point.
r_ea = zeros(3, 1);
prev_r_ea = [Inf
    Inf
    Inf];

% This is only for initial estimate - hardcode predicted receiver
% clock offset (rho_c)
rho_c = 1;

% The while loop will go through all the data and compute the least square
% position until we converge by having the difference between the previous
% r_ea and current r_ea.
while true
    % More details about least square in its file.
    x_plus = CW1_Least_squares_position(sattelite_numbers,...
        pseudo_ranges(2,2:9), r_ea, rho_c, omega_ie, c);
    % After using the unweighted least square we save the new position in
    % r_ea
    r_ea = x_plus(1:3);
    rho_c = x_plus(4);
    
    % The convergence check. If the difference between the computations is 
    % less then 10cm then we converged.
    if sqrt(sum((prev_r_ea-r_ea).^2)) < 0.1
        break;
    end
    
    % Replacing the previous r_ea for next iteration.
    prev_r_ea = r_ea;
end

% Final run of least squares with outlier detection
outliers = CW1_Detect_outliers(0, sattelite_numbers, pseudo_ranges(2,...
    2:num_of_satellites+1), r_ea, rho_c, omega_ie, c);
num_of_correct_satellites = length(find(outliers == 0));
correct_satellite_numbers = zeros(1, num_of_correct_satellites);
k = 1;
for j=1:num_of_satellites
    if outliers(j) == 0
        correct_satellite_numbers(k) = sattelite_numbers(j);
        k = k + 1;
    end
end

% Pseudo ranges of non-outliers
correct_pseudo_ranges = zeros(1, num_of_correct_satellites);
for j=1:num_of_correct_satellites
    correct_pseudo_ranges(1, j) = pseudo_ranges(2,...
        satellite_indices(correct_satellite_numbers(j))+1);
end

x_plus = CW1_Least_squares_position(correct_satellite_numbers,...
    correct_pseudo_ranges, r_ea, rho_c, omega_ie, c);

% Saving the initial position which will be used in CW1_GNSS
initial_x_plus = x_plus;