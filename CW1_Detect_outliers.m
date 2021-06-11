function outliers = CW1_Detect_outliers(time, sattelite_numbers, ranges,...
    est_pos, est_clock_offset, omega_ie, c)

    num_of_satellites = length(sattelite_numbers);
    r_ej = zeros(3, num_of_satellites);
    
    % Compute the cartesian positions of the sattelites
    for j=1:length(sattelite_numbers)
        r_ej(:, j) = Satellite_position_and_velocity(time,...
            sattelite_numbers(j));
    end
    
    r_aj = zeros(num_of_satellites, 1);
    C_e = zeros(num_of_satellites, 3, 3);

    % Predict the ranges from the satellites to the aproximate
    % user position. This is done similarly as in
    % CW1_Least_squares_position.
    for i=1:length(sattelite_numbers)
        C_e(i, :, :) = eye(3);
        inner_exp = squeeze(C_e(i, :, :)) * r_ej(:, i) - est_pos;
        r_aj(i, 1) = sqrt(transpose(inner_exp) * inner_exp);
        C_e(i, 1, 2) = omega_ie * r_aj(i, 1) / c;
        C_e(i, 2, 1) = -omega_ie * r_aj(i, 1) / c;
        inner_exp = squeeze(C_e(i, :, :)) * r_ej(:, i) - est_pos;
        r_aj(i, 1) = sqrt(transpose(inner_exp) * inner_exp);
    end
    
    % Compute the line of sight unit vectors.
    u_aj = zeros(3, num_of_satellites);

    for i=1:length(num_of_satellites)
        u_aj(:, i) =...
            (squeeze(C_e(i, :, :)) * r_ej(:, i) - est_pos) / r_aj(i);
    end
    
    % Compute measurement innovation vector
    delta_z = zeros(num_of_satellites, 1);
    for i=1:num_of_satellites
        delta_z(i, 1) = ranges(i) - r_aj(i, 1) - est_clock_offset;
    end
    
    H_G = -transpose(u_aj);
    H_G(:,4) = ones(8, 1);
    
    % Compute the residuals vector
    v = (H_G * inv(H_G' * H_G) * H_G' - eye(8)) * delta_z;

    % Compute the residuals covariance matrix
    sigma_rho = 10; % <-- measurement error standard deviation
    C_v = (eye(8) - H_G * inv(H_G' * H_G)* H_G') * sigma_rho^2;

    % Compute the normalised residuals
    threshold = 2.23;

    outliers = zeros(1, num_of_satellites);
    for j=1:num_of_satellites
       if abs(v(j)) > (sqrt(C_v(j,j)) * threshold)
           disp('outlier');
           disp(j)
           outliers(1, j) = 1;
       end
    end
end