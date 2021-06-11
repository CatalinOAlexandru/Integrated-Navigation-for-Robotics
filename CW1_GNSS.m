
% From the dead reckoning data we obtain that the robot has no initial
% velocity.

% Initialise state estimates
x_plus = zeros(8,1);
x_plus(1:3) = initial_x_plus(1:3);
x_plus(7) = initial_x_plus(4);

% Initialise error covariance matrix
P_plus =  zeros(8);
P_plus(1,1) = 100;
P_plus(2,2) = 100;
P_plus(3,3) = 100;
P_plus(4,4) = 0.01;
P_plus(5,5) = 0.01;
P_plus(6,6) = 0.01;
P_plus(7,7) = 100000^2; % Clock offset std is 100000m
P_plus(8,8) = 200^2; % Clock drift std is 200m/s

% Compute the transition matrix using
tau_s = 0.5;
phi = eye(8);
phi(1:3,4:6) = tau_s * eye(3);
phi(7,8) = tau_s;

S_a = 1.5; % We assume acceleration PSD is 1.5 m^2s^−3
S_c_phi = 0.01; % We assume clock phase PSD is 0.01 m^2s^−1
S_cf = 0.04; % We assume clock frequency PSD is 0.04 m^2s^−3

% Compute the system noise covariance matrix
Q = zeros(8,8);
Q(1:3,1:3) = 1/3 * S_a * tau_s^3 * eye(3);
Q(1:3,4:6) = 1/2 * S_a * tau_s^2 * eye(3);
Q(4:6,1:3) = 1/2 * S_a * tau_s^2 * eye(3);
Q(4:6,4:6) = S_a * tau_s * eye(3);
Q(7,7) = S_c_phi * tau_s + 1/3 * S_cf * tau_s^3;
Q(7,8) = 1/2 * S_cf * tau_s^2;
Q(8,7) = 1/2 * S_cf * tau_s^2;
Q(8,8) = S_cf * tau_s;

gnss_state = zeros(num_rows-1,7);

for row_idx=2:num_rows
    time = pseudo_ranges(row_idx,1);
    
    outliers = CW1_Detect_outliers(time,sattelite_numbers,...
        pseudo_ranges(row_idx, 2:num_of_satellites+1),...
        x_plus(1:3), x_plus(7), omega_ie, c);
    num_of_correct_satellites = length(find(outliers == 0));
    correct_satellite_numbers = zeros(1, num_of_correct_satellites);
    k = 1;
    for j=1:num_of_satellites
        if outliers(j) == 0
            correct_satellite_numbers(k) = sattelite_numbers(j);
            k = k + 1;
        end
    end
    
    % Propagate the state estimates
    x_minus = phi * x_plus;

    % Propagate the error covariance matrix:
    P_minus = phi * P_plus * phi' + Q;

    % Predict the ranges from approximate user position to each satellite
    r_ej = zeros(3, num_of_correct_satellites);
    v_ej = zeros(3, num_of_correct_satellites);
    for i=1:num_of_correct_satellites
        [sat_r_es_e,sat_v_es_e] = Satellite_position_and_velocity(time,...
            correct_satellite_numbers(i));
        r_ej(:, i) = sat_r_es_e;
        v_ej(:, i) = sat_v_es_e;
    end

    % Initialise the matrices
    r_ea = x_minus(1:3);
    r_aj = zeros(num_of_correct_satellites, 1);
    % 3x3 matrix for each satellite
    C_e = zeros(num_of_correct_satellites, 3, 3);

    for i=1:num_of_correct_satellites
        C_e(i, :, :) = eye(3);
        % First compute r_aj with compensation matrix set to identity
        inner_exp = squeeze(C_e(i, :, :)) * r_ej(:, i) - r_ea;
        r_aj(i, 1) = sqrt(transpose(inner_exp) * inner_exp);
        % Then recompute correct Sagnac effect compensation matrix
        C_e(i, 1, 2) = omega_ie * r_aj(i, 1) / c;
        C_e(i, 2, 1) = -omega_ie * r_aj(i, 1) / c;
        inner_exp = squeeze(C_e(i, :, :)) * r_ej(:, i) - r_ea;
        % And calculate valid range values
        r_aj(i, 1) = sqrt(transpose(inner_exp) * inner_exp);
    end

    % Compute the line-of-sight unit vectors for each satellite
    u_aj = zeros(3, num_of_correct_satellites);

    for i=1:num_of_correct_satellites
        u_aj(:, i) = (squeeze(C_e(i, :, :)) * r_ej(:, i) - r_ea) / r_aj(i);
    end

    % Compute range rates from receiver to each satellite
    r_dot_caret = zeros(1, num_of_correct_satellites);

    for i=1:num_of_correct_satellites
        v_ea = x_minus(4:6);
        C = squeeze(C_e(i, :, :));
        r_dot_caret(1, i) = u_aj(:, i)' * (C * (v_ej(:, i) +...
            Omega_ie * r_ej(:, i)) - (v_ea + Omega_ie * r_ea));
    end

    % Compute measurement matrix H
    H = zeros(num_of_correct_satellites*2,8);
    H(1:num_of_correct_satellites,1:3) = -u_aj';
    H(num_of_correct_satellites+1:2*num_of_correct_satellites,4:6) = -u_aj';
    H(1:num_of_correct_satellites,7) = ones(num_of_correct_satellites,1);
    H(num_of_correct_satellites+1:2*num_of_correct_satellites,8) =...
        ones(num_of_correct_satellites,1);

    % Compute measurement noise covariance matrix R
    % We assume noise standard deviation of 10m on all pseudo-range
    % measurements and 0.05 m/s on all pseudo-range rate measurements
    sigma_rho = 10;
    sigma_r = 0.05;
    R = zeros(2*num_of_correct_satellites,2*num_of_correct_satellites);
    R(1:num_of_correct_satellites,1:num_of_correct_satellites) =...
        eye(num_of_correct_satellites) * sigma_rho^2;
    R(num_of_correct_satellites+1:2*num_of_correct_satellites,...
        num_of_correct_satellites+1:2*num_of_correct_satellites) =...
        eye(num_of_correct_satellites) * sigma_r^2;

    % Calculate Kalman gain
    K = P_minus * H' * inv(H * P_minus * H' + R);

    % Measurement innovation vector delta_z
    delta_z = zeros(2*num_of_correct_satellites,1);
    rho = zeros(1, num_of_correct_satellites);
    rho_dot = zeros(1, num_of_correct_satellites);
    
    for j=1:num_of_correct_satellites
        rho(1, j) = pseudo_ranges(row_idx,...
            satellite_indices(correct_satellite_numbers(j))+1);
        rho_dot(1, j) = pseudo_range_rates(row_idx,...
            satellite_indices(correct_satellite_numbers(j))+1);
    end
    
    delta_z(1:num_of_correct_satellites,1) = rho' - r_aj...
        - ones(num_of_correct_satellites,1) * x_minus(7,1);
    delta_z(num_of_correct_satellites+1:2*num_of_correct_satellites,1)...
        = rho_dot' - r_dot_caret'...
        - ones(num_of_correct_satellites,1)*x_minus(8,1);

    % Update the state vector
    x_plus = x_minus + K * delta_z;

    % Update error covariance matric
    P_plus = (eye(8) - K * H) * P_minus;

    % Convert new state predictions to NED
    [L_b,lambda_b,h_b,v_eb_n] = pv_ECEF_to_NED(x_plus(1:3,1), x_plus(4:6,1));
    
    % Save in GNSS state matrix
    gnss_state(row_idx-1,:) = [time L_b*rad_to_deg lambda_b*rad_to_deg...
        h_b v_eb_n(1) v_eb_n(2) v_eb_n(3)];
end

% writematrix(gnss_state(:, 2:3),'CW1_GNSS_coords.csv');