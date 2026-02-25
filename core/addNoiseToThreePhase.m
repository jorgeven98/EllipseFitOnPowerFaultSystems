function [signal_with_noise] = addNoiseToThreePhase(signal, noise_level, noise_type)
%ADDNOISETOTHREEPHASE Adds Gaussian noise to a three-phase signal
%   [SIGNAL_WITH_NOISE] = ADDNOISETOTHREEPHASE(SIGNAL, NOISE_LEVEL, NOISE_TYPE)
%
%   Parameters:
%   SIGNAL - Three-phase signal matrix of dimensions [n_phases, n_samples]
%            where n_phases must be 3
%   NOISE_LEVEL - Noise level (in dB or as standard deviation)
%   NOISE_TYPE - 'dB' for decibels or 'std' for standard deviation
%
%   Returns:
%   SIGNAL_WITH_NOISE - Three-phase signal with added Gaussian noise

    % Validate inputs
    [n_phases, n_samples] = size(signal);
    if n_phases ~= 3
        error('Signal must have exactly 3 phases (rows)');
    end

    if ~ischar(noise_type) && ~isstring(noise_type)
        error('noise_type must be "dB" or "std"');
    end

    % Calculate signal power
    signal_power = mean(sum(signal.^2, 1)) / 3;

    % Determine noise standard deviation
    if strcmpi(noise_type, 'dB')
        % Convert dB to linear factor
        SNR_linear = 10^(noise_level/10);
        noise_power = signal_power / SNR_linear;
        noise_std = sqrt(noise_power);
    elseif strcmpi(noise_type, 'std')
        noise_std = noise_level;
    else
        error('noise_type must be "dB" or "std"');
    end

    % Generate Gaussian noise
    noise = noise_std * randn(n_phases, n_samples);

    % Add noise to signal
    signal_with_noise = signal + noise;

    % If specified in dB, verify the resulting SNR
    if strcmpi(noise_type, 'dB')
        % Calculate actual SNR
        noise_power_actual = mean(sum(noise.^2, 1)) / 3;
        SNR_actual = 10 * log10(signal_power / noise_power_actual);

        fprintf('Target SNR: %.2f dB, Actual SNR: %.2f dB\n', noise_level, SNR_actual);
    end
end
