% The following code is part of our implementation of Workshop 1 Task 1 B
% potentially slightly modified to accomodate the new data.

function x_plus = CW1_Least_squares_position(sattelite_numbers,...
    pseudo_ranges, r_ea, rho_c, omega_ie, c)
    % Tast 1A b  --------------------------------------------------------------
    % The initial step 1A_a has been achived at the beginning CW1_Initial_Position
    % In this step we compute the cartesian positions of the sattelites using
    % the pre-defined function Satellite_position_and_velocity which returns
    % the position and velocity given the time and satellite.
    r_ej = zeros(3, 8);
    for i=1:length(sattelite_numbers)
        [position, velocity] = Satellite_position_and_velocity(0,...
            sattelite_numbers(i));
        % r_ej will hold the all the position values for each satellite.
        r_ej(:, i) = position;
    end


    % Task 1A c  --------------------------------------------------------------
    % I this step we predict the ranges from the satellites to the aproximate
    % user position.
    % r_aj will hold these prediction and C_e is the Sagnac effect compensation
    % matrix as based on the equation (1) from [1]
    r_aj = zeros(8, 1);
    C_e = zeros(8, 3, 3);

    for i=1:length(sattelite_numbers)
        % We start by creating and computing C_e
        C_e(i, :, :) = eye(3);
        % inner_exp is based on the r_aj formula from [1] Task 1Ac and
        % initialised C_e
        inner_exp = squeeze(C_e(i, :, :)) * r_ej(:, i) - r_ea;
        % We now use the temporary inner_exp to calculate r_aj
        r_aj(i, 1) = sqrt(transpose(inner_exp) * inner_exp);
        % As we have an r_aj now, we add the missing elements to C_e as in the
        % formula.
        C_e(i, 1, 2) = omega_ie * r_aj(i, 1) / c;
        C_e(i, 2, 1) = -omega_ie * r_aj(i, 1) / c;
        % We re-compute the inner_exp with the new C_e
        inner_exp = squeeze(C_e(i, :, :)) * r_ej(:, i) - r_ea;
        % Now with an accurate inner_exp we can get our final r_aj
        r_aj(i, 1) = sqrt(transpose(inner_exp) * inner_exp);
        % now we repeat for all satellites and all data will be stored in r_aj
    end



    % Task 1A d  --------------------------------------------------------------
    % Compute the line of sight unit vector using the aproximate user position
    % for each sattelite using the equation (2) as seen in [1] and previously
    % calculated values.
    u_aj = zeros(3, 8);

    for i=1:length(sattelite_numbers)
        u_aj(:, i) = (squeeze(C_e(i, :, :)) * r_ej(:, i) - r_ea) / r_aj(i);
    end


    % Task 1A e  --------------------------------------------------------------
    % Compute the predicted state vector x_minus and the measurement innovation
    % vectore delta_z using the measurement matrix H_G

    % Equation 1
    % We initialise and set the x_minus values from r_ea and rho_c (defined in 
    % CW1_Initial_Position.m)
    x_minus = zeros(4, 1);
    x_minus(1:3, 1) = r_ea;
    x_minus(4, 1) = rho_c;

    % Equation 2
    % Compute delta_z (measurement innovation vector) for all satellites,
    % using rho_a (measured pseudo-ranges from satellite to user) and rho_c.
    % We start by initiallising rho_a from the pseudo_ranges csv and delta_z
    % with zeros.
    rho_a = pseudo_ranges;
    delta_z = zeros(8, 1);
    for i=1:length(sattelite_numbers)
        % We iterate through all the sattelites and compute delta_z using the
        % equation 4 in [1].
        delta_z(i, 1) = rho_a(i) - r_aj(i, 1) - rho_c;
    end

    % Equation 3
    % We compute the measurement matrix as shown in equations (4) from [1] and 
    % add ones as the 4th column. 
    H_G = -transpose(u_aj);
    H_G(:,4) = ones(8, 1);


    % Task 1A f  --------------------------------------------------------------
    % We compute the position and reciver clock offset using the unweighted
    % least squares based on the equation (5) from [1].
    x_plus = x_minus + inv(transpose(H_G) * H_G) * transpose(H_G) * delta_z;
end