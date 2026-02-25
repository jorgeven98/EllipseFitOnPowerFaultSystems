function [va, vb, vc, t, faultType] = generateFault(faultType, faultStartTime, faultDuration, samplingFreq, faultMagnitude, faultImpedance, groundImpedance)
% GENERATEFAULT Generates three-phase voltage waveforms for different fault types
%
% Inputs:
%   faultType      - Fault type ('A-G', 'B-G', 'C-G', 'A-B', 'B-C', 'C-A',
%                    'A-B-G', 'B-C-G', 'C-A-G', 'ABC')
%   faultStartTime - Fault start time [s]
%   faultDuration  - Fault duration [s]
%   samplingFreq   - Sampling frequency [Hz]
%   faultMagnitude - Fault magnitude (0 to 1), where 0 is full fault and 1 is no fault
%   faultImpedance - Fault impedance [p.u.]
%   groundImpedance- Ground impedance [p.u.]
%
% Outputs:
%   va, vb, vc - Normalized three-phase voltages
%   t          - Time vector [s]
%   faultType  - Generated fault type string

% Default parameters
if nargin < 7, groundImpedance = 0; end
if nargin < 6, faultImpedance = 0.001; end
if nargin < 5, faultMagnitude = 0.2; end
if nargin < 4, samplingFreq = 10000; end
if nargin < 3, faultDuration = 0.1; end
if nargin < 2, faultStartTime = 0.05; end
if nargin < 1, faultType = 'A-G'; end

% System parameters
f = 50; % System frequency [Hz]
totalTime = 0.3; % Total simulation time [s]

% Time vector
t = 0:1/samplingFreq:totalTime;
omega = 2*pi*f;

% Initialize voltages under normal conditions (p.u. values)
Va_mag = 1.0;
Vb_mag = 1.0;
Vc_mag = 1.0;
Va_angle = 0;
Vb_angle = -2*pi/3;
Vc_angle = 2*pi/3;

% Generate normal three-phase voltages
va = Va_mag * cos(omega*t + Va_angle);
vb = Vb_mag * cos(omega*t + Vb_angle);
vc = Vc_mag * cos(omega*t + Vc_angle);

% Indices for the fault duration
faultStartIndex = find(t >= faultStartTime, 1);
faultEndIndex = find(t >= (faultStartTime + faultDuration), 1);

% Apply fault according to type
switch upper(faultType)
    case 'A-G' % Single-phase fault: phase A to ground
        % Calculate voltage during fault considering fault and ground impedance
        Z_ratio = faultImpedance / (faultImpedance + groundImpedance);
        va(faultStartIndex:faultEndIndex) = faultMagnitude * va(faultStartIndex:faultEndIndex);

        % Effect on other phases due to grounding
        if groundImpedance > 0
            vb(faultStartIndex:faultEndIndex) = vb(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * va(faultStartIndex:faultEndIndex);
            vc(faultStartIndex:faultEndIndex) = vc(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * va(faultStartIndex:faultEndIndex);
        end

    case 'B-G' % Single-phase fault: phase B to ground
        Z_ratio = faultImpedance / (faultImpedance + groundImpedance);
        vb(faultStartIndex:faultEndIndex) = faultMagnitude * vb(faultStartIndex:faultEndIndex);

        if groundImpedance > 0
            va(faultStartIndex:faultEndIndex) = va(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * vb(faultStartIndex:faultEndIndex);
            vc(faultStartIndex:faultEndIndex) = vc(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * vb(faultStartIndex:faultEndIndex);
        end

    case 'C-G' % Single-phase fault: phase C to ground
        Z_ratio = faultImpedance / (faultImpedance + groundImpedance);
        vc(faultStartIndex:faultEndIndex) = faultMagnitude * vc(faultStartIndex:faultEndIndex);

        if groundImpedance > 0
            va(faultStartIndex:faultEndIndex) = va(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * vc(faultStartIndex:faultEndIndex);
            vb(faultStartIndex:faultEndIndex) = vb(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * vc(faultStartIndex:faultEndIndex);
        end

    case 'A-B' % Phase-to-phase fault between A and B
        % During a phase-to-phase fault, voltages are adjusted to maintain zero sum
        faultVoltage = faultMagnitude * mean([va(faultStartIndex:faultEndIndex); vb(faultStartIndex:faultEndIndex)]);
        va(faultStartIndex:faultEndIndex) = va(faultStartIndex:faultEndIndex).*(1-faultMagnitude) + faultVoltage;
        vb(faultStartIndex:faultEndIndex) = vb(faultStartIndex:faultEndIndex).*(1-faultMagnitude) + faultVoltage;

    case 'B-C' % Phase-to-phase fault between B and C
        faultVoltage = faultMagnitude * (vb(faultStartIndex:faultEndIndex) - vc(faultStartIndex:faultEndIndex))/2;
        vb(faultStartIndex:faultEndIndex) = vb(faultStartIndex:faultEndIndex).*(1-faultMagnitude) + faultVoltage;
        vc(faultStartIndex:faultEndIndex) = vc(faultStartIndex:faultEndIndex).*(1-faultMagnitude) + faultVoltage;

    case 'C-A' % Phase-to-phase fault between C and A
        faultVoltage = faultMagnitude * (vc(faultStartIndex:faultEndIndex) - va(faultStartIndex:faultEndIndex))/2;
        vc(faultStartIndex:faultEndIndex) = vc(faultStartIndex:faultEndIndex).*(1-faultMagnitude) + faultVoltage;
        va(faultStartIndex:faultEndIndex) = va(faultStartIndex:faultEndIndex).*(1-faultMagnitude) + faultVoltage;

    case 'A-B-G' % Double phase-to-ground fault: phases A and B
        Z_ratio = faultImpedance / (faultImpedance + groundImpedance);

        % Faulted phase voltages are reduced
        va(faultStartIndex:faultEndIndex) = faultMagnitude * va(faultStartIndex:faultEndIndex);
        vb(faultStartIndex:faultEndIndex) = faultMagnitude * vb(faultStartIndex:faultEndIndex);

        % Phase C may be affected due to ground impedance
        if groundImpedance > 0
            vc(faultStartIndex:faultEndIndex) = vc(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * (va(faultStartIndex:faultEndIndex) + vb(faultStartIndex:faultEndIndex));
        end

    case 'B-C-G' % Double phase-to-ground fault: phases B and C
        Z_ratio = faultImpedance / (faultImpedance + groundImpedance);

        vb(faultStartIndex:faultEndIndex) = faultMagnitude * vb(faultStartIndex:faultEndIndex);
        vc(faultStartIndex:faultEndIndex) = faultMagnitude * vc(faultStartIndex:faultEndIndex);

        if groundImpedance > 0
            va(faultStartIndex:faultEndIndex) = va(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * (vb(faultStartIndex:faultEndIndex) + vc(faultStartIndex:faultEndIndex));
        end

    case 'C-A-G' % Double phase-to-ground fault: phases C and A
        Z_ratio = faultImpedance / (faultImpedance + groundImpedance);

        vc(faultStartIndex:faultEndIndex) = faultMagnitude * vc(faultStartIndex:faultEndIndex);
        va(faultStartIndex:faultEndIndex) = faultMagnitude * va(faultStartIndex:faultEndIndex);

        if groundImpedance > 0
            vb(faultStartIndex:faultEndIndex) = vb(faultStartIndex:faultEndIndex) - ...
                (1-faultMagnitude) * Z_ratio * (vc(faultStartIndex:faultEndIndex) + va(faultStartIndex:faultEndIndex));
        end

    case 'ABC' % Three-phase fault
        % In a balanced three-phase fault, all three voltages are reduced equally
        va(faultStartIndex:faultEndIndex) = faultMagnitude * va(faultStartIndex:faultEndIndex);
        vb(faultStartIndex:faultEndIndex) = faultMagnitude * vb(faultStartIndex:faultEndIndex);
        vc(faultStartIndex:faultEndIndex) = faultMagnitude * vc(faultStartIndex:faultEndIndex);

    otherwise
        error('Unrecognized fault type');
end

end
