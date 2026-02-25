function transformed_voltage = applyRotorToVoltage(data, rotor)
    % APPLYROTORTOVOLTAGE Applies a rotor to three-phase voltage data
    %
    % Inputs:
    %   data  - 3xN matrix with the three voltage phases [Va; Vb; Vc]
    %   rotor - Multivector object representing the rotor
    %
    % Output:
    %   transformed_voltage - 3xN matrix with the transformed data

    % Create the rotor reverse (needed for the sandwich transformation)
    rotor_reverse = rotor.reverse();

    % Get dimensions
    [n_phases, n_samples] = size(data);

    % Initialize output matrix
    transformed_voltage = zeros(size(data));

    % Process each voltage sample
    for i = 1:n_samples
        % Convert current sample to a GA vector
        v = ga3.EncodeVector([data(1,i), data(2,i), data(3,i)]);

        % Apply the rotor transformation: v' = R v R̃
        v_transformed = rotor.gp(v).gp(rotor_reverse);

        % Extract components of the transformed vector
        v_components = v_transformed.getDataArray();

        % Store components in the output matrix
        % Vector components are at positions 2-4
        transformed_voltage(1,i) = v_components(2); % x component
        transformed_voltage(2,i) = v_components(3); % y component
        transformed_voltage(3,i) = v_components(4); % z component
    end
end
