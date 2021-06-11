format long;
Define_Constants;

num_of_satellites = 8;
% load `csv` data files
pseudo_ranges = readmatrix('Pseudo_ranges.csv');
pseudo_range_rates = readmatrix('Pseudo_range_rates.csv');
sattelite_numbers = pseudo_ranges(1,2:9);
num_rows = size(pseudo_ranges,1);

% load dead reckoning data from csv file
dead_reckoning_data = readmatrix('Dead_reckoning.csv');

% we use an aggregated speed, 
% which is the mean of the 4 whells
avg_speed_4_wheels = sum(dead_reckoning_data(:,2:5), 2) / 4;
% we use the compass only for vehicle heading
compass_bearing = dead_reckoning_data(:,7);
timestamp = dead_reckoning_data(:,1);
num_epochs = size(dead_reckoning_data, 1);

% the initial position comes from the GNSS initialisatin's
% result, the speed initialisation is 0 (in the dead 
% reckoning csv file, the first row's speed is 0)
[L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(initial_x_plus(1:3), zeros(3,1));

% arrays to store the speed and position data
averageSpeed = zeros(2,num_epochs); % for north and east velocities
latitude = zeros(num_epochs,1);
longitude = zeros(num_epochs,1);

% the initial position is the result of GNSS initialisation
latitude(1) = L_b;
longitude(1) = lambda_b;
h = 24; %h_b; the GNSS is most inaccuarate in the up-down
        % dimension, as the vehicle does not travel a large
        % distance, we use the geographics height of the 
        % starting location (the original GNSS was 38.8 m)
        % This will be used for updating the lattitude and
        % longitude data

% arrays for average speed components
V_N_k = zeros(num_epochs,1);
V_E_k = zeros(num_epochs,1);
% for the damped instantaneous DR velocity, the 0th
% value is just the appropriate component of the speed at time=0
V_N_k(1) = avg_speed_4_wheels(1) * cos(compass_bearing(1) * deg_to_rad);
V_E_k(1) = avg_speed_4_wheels(1) * sin(compass_bearing(1) * deg_to_rad);

% after initialisation we start the loop
for index=2:num_epochs
    % heading in previous epoch
    psi_minus = compass_bearing(index-1) * deg_to_rad;
    % current heading
    psi = compass_bearing(index) * deg_to_rad;
    % current average speed
    v_k = avg_speed_4_wheels(index);
    
    % the average velocity between previous and current epoch
    % by north and east components [eq (1) in week 3 workshop pdf]
    v_north = 1/2 * (cos(psi) + cos(psi_minus)) * v_k;
    v_east = 1/2 * (sin(psi) + sin(psi_minus)) * v_k;
    
    % store them in the pre-defined array
    averageSpeed(:,index) = [v_north, v_east]'; 
    
    % for updating the lattitude and longitude data
    % we add time*velocity to the postions
    % for that we need `time` (delta_t here) 
    delta_t = timestamp(index) - timestamp(index-1);
    % from lattitude we compute the meridian and
    % transverse curvatures, to correct with them the update
    % of lattitude and longitude [eq. (2) in week 3 workshop]
    [R_N,R_E]= Radii_of_curvature(latitude(index-1));
    latitude(index) = latitude(index-1) + (v_north * delta_t) / (R_N + h);
    longitude(index) = longitude(index-1)...
        + (v_east * delta_t) / ((R_E + h) * cos(latitude(index)));
    
    % we calculate the damped instantaneous DR velocity
    % this we'll use for the integration part
    V_N_k(index) = 1.7 * v_north - 0.7 * V_N_k(index-1);
    V_E_k(index) = 1.7 * v_east - 0.7 * V_E_k(index-1);
    
%     disp(latitude(index) * rad_to_deg)
%     disp(longitude(index) * rad_to_deg)
%     disp(V_N_k(index))
%     disp(V_E_k(index))
end

dead_reckoning_state = horzcat(timestamp, latitude*rad_to_deg,...
    longitude*rad_to_deg, V_N_k, V_E_k, compass_bearing);

% writematrix(horzcat(latitude*rad_to_deg, longitude*rad_to_deg),'CW1_Dead_reckoning_coords.csv');