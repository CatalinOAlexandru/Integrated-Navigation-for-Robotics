% we run the GNSS initialisation
% using its result we calculate the GNSS only
% and DR only positions. Here we integrate the two.
CW1_Initial_Position;
CW1_GNSS;
CW1_Dead_Reckoning;

% dead_reckoning_state
% Time in seconds | Lat* | Long* | V_north | V_east | phi
%       1            2      3         4        5       6

% gnss_state
% Time in seconds | Lat* | Long* | Hight/m | V_north | V_east | V_down
%       1            2      3         4         5          6       7

% state vector, initially we don't know the DR errors,
% so we initialise an all zeros vector
x = zeros(4);
% sigma_v and sigma_r are the same as in CW1_GNSS's 
% error covariance matrix (here we have a square roots
% of those)
sigma_v = 0.06;
sigma_r = 10;
tau = 0.5;  % time step, we have new data in every 0.5 sec
S_DR = 0.01;  % Coming from the coursework PDF (wheel speed PSD)

% we construct the state estimation error covariance matrix
% [eq. (5) in week 3 workshop]
P_plus = eye(4);
P_plus(1,1) = sigma_v^2;
P_plus(2,2) = sigma_v^2;
[R_N,R_E]= Radii_of_curvature(gnss_state(1,2) * deg_to_rad);
P_plus(3,3) = (sigma_r^2) / (R_N + gnss_state(1,4))^2;
P_plus(4,4) = (sigma_r^2) / ((R_E + gnss_state(1,4))^2*cos(gnss_state(1,2))^2);

sizeOf = size(gnss_state, 1);

% to save the state vector in each epoch
x_plus = zeros(4,1);

% For step 6
% Assuming the same GNSS position and velocity std everywhere (10m, 0.1m/s)
sigma_Gr = 10;
sigma_Gv = 0.1;

integrated_state = zeros(sizeOf, 6);
integrated_state(1, :) = dead_reckoning_state(1, :);

for index=2:sizeOf
    % the GNSS height of actual time step
    h = gnss_state(index-1,4);
    % lattitude in degrees
    L = gnss_state(index-1,2) * deg_to_rad;
    % meridian and transverse radius from the lattitude
    [R_N,R_E]= Radii_of_curvature(L);

    % Step 1
    % we construct the transition matrix
    phi = eye(4);
    phi(3,1) = tau / (R_N + h);
    phi(4,2) = tau / ((R_E + h) * cos(L));

    % Step 2
    % we construct the system noise covariance matrix
    Q = eye(4);
    Q(1,1) = S_DR * tau;
    Q(1,3) = 1/2 * ((S_DR * tau^2) / (R_N + h));
    Q(2,2) = S_DR * tau;
    Q(2,4) = 1/2 * ((S_DR * tau^2) / ((R_E + h) * cos(L)));
    Q(3,1) = 1/2 * ((S_DR * tau^2) / (R_N + h));
    Q(3,3) = 1/3 * S_DR * tau^3 / ((R_N + h)^2);
    Q(4,2) = 1/2 * ((S_DR * tau^2) / ((R_E + h) * cos(L)));
    Q(4,4) = 1/3 * S_DR * tau^3 / ((R_E + h)^2 * cos(L)^2);
    
    % Step 3 - NOTE: Might be wrong here. Initial x_plus value vs x
    % Propagate the state estimates
    x_minus = phi * x_plus;
    
    % Step 4
    % Propagate the error covariance matrix
    P_minus = phi * P_plus * phi' + Q;
    
    % Step 5
    % Measurement matrix (note: this could be defined once
    % outside of the loop)
    H = zeros(4,4);
    H(1,3) = -1;
    H(2,4) = -1;
    H(3,1) = -1;
    H(4,2) = -1;
    
    % Overrite h for index
    % R_E and R_N might NOT need to be changed to index (keep index-1?)
    h = gnss_state(index,4);
    L = gnss_state(index,2) * deg_to_rad;
    [R_N,R_E] = Radii_of_curvature(L);
    
    % Step 6
    % Compute the measurement noise covariance matrix
    R = zeros(4,4);
    R(1,1) = sigma_Gr^2 / (R_N + h)^2;
    R(2,2) = sigma_Gr^2 / ((R_E + h)^2 * cos(L)^2);
    R(3,3) = sigma_Gv^2;
    R(4,4) = sigma_Gv^2;
    
    % Step 7
    % Compute the Kalman gain matrix
    K = P_minus * H' * inv(H * P_minus * H' + R);
    
    % Step 8
    % Formulate the measurement innovation vector
    L_D = dead_reckoning_state(index, 2) * deg_to_rad; 
    lambda_D = dead_reckoning_state(index, 3) * deg_to_rad; 
    V_N_D = dead_reckoning_state(index, 4); 
    V_E_D = dead_reckoning_state(index, 5); 
    
    delta_z = zeros(4,1);
    delta_z(1,1) = L - L_D;
    delta_z(2,1) = gnss_state(index,3) * deg_to_rad - lambda_D;
    delta_z(3,1) = gnss_state(index,5) - V_N_D;
    delta_z(4,1) = gnss_state(index,6) - V_E_D;
    
    delta_z = delta_z - H * x_minus;
    
    % Step 9
    % Update the state estimates
    x_plus = x_minus + K * delta_z;
    
    % Step 10
    % Update the error covariance matrix
    P_plus = (eye(4) - K * H) * P_minus;
    
    % Use the Kalman filter estimates to correct the
    % DR solution in each epoch
    % Save to integrated state
    integrated_state(index, 1) = dead_reckoning_state(index, 1);
    integrated_state(index, 2) = (L_D - x_plus(3)) * rad_to_deg;
    integrated_state(index, 3) = (lambda_D - x_plus(4)) * rad_to_deg;
    integrated_state(index, 4) = V_N_D - x_plus(1);
    integrated_state(index, 5) = V_E_D - x_plus(2);
    integrated_state(index, 6) = dead_reckoning_state(index, 6);
end

% writematrix(integrated_state(:,2:3),'CW1_integrated_coords.csv');