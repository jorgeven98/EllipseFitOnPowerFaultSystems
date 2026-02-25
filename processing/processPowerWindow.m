function [vec_plane, bivector_components, center_3D, semiaxes, angle] = processPowerWindow(data, t, freq)

    % PROCESSPOWERWINDOW Processes a signal window to extract ellipse parameters

    % Parameters
    cycle_duration = 1/freq; % Duration of one cycle in seconds
    n_samples = length(t);

    % Compute samples per cycle
    windows_size_percentage = 0.75;
    dt = t(2) - t(1);
    samples_per_cycle = round(cycle_duration/ dt);
    windows_size = round(samples_per_cycle*windows_size_percentage);

    % Initialize output arrays
    d2 = round(windows_size*0.375);
    d4 = round(windows_size*0.1875);
    data_size = length(1:round(windows_size/6):n_samples - windows_size);
    bivector_components = zeros(3,data_size);
    semiaxes = zeros(2,data_size);
    center_3D = zeros(3,data_size);
    angle = zeros(1,data_size);
    vec_plane = zeros(3,data_size);

    k = 1;

    for i = 1:round(windows_size/8):n_samples - windows_size

    % Extract window
    data_window = data(:,i:i + windows_size);
    v1 = ga3.EncodeVector([data_window(1,1) - data_window(1,d2), data_window(2,1) - data_window(2,d2), data_window(3,1) - data_window(3,d2)]);
    v2 = ga3.EncodeVector([data_window(1,d4) - data_window(1,end), data_window(2,d4) - data_window(2,end), data_window(3,d4) - data_window(3,end)]);

    % Compute the normalized bivector
    bivector = v2.op(v1);
    bivector_norm = bivector.norm();
    %bivector_norm = 1;
    bivector_normalized = bivector/bivector_norm;
    bivector_components_pre = bivector_normalized.getDataArray;

    bivector_components(1, k) = bivector_components_pre(1);
    bivector_components(2, k) = bivector_components_pre(3);
    bivector_components(3, k) = -bivector_components_pre(2);

    n = bivector_normalized.dual;
    p1 = ga3.EncodeVector([data_window(1,1),data_window(2,1),data_window(3,1)]);
    d = p1.sp(n);
    N = n*d;

    vec_plane(:,k) = N.Data;
    data_0 = data_window - vec_plane(:,k);

    % Get the rotor to align the bivector with the XY plane
    % NOTE: alignBivectorWithPlaneXY is missing from the project and must be implemented
    R = alignBivectorWithPlaneXY(bivector_normalized);

    % Apply the rotor to the window data
    data_XY = applyRotorToVoltage(data_0, R);

    % Extract X and Y components for ellipse fitting
    X = data_XY(1,:);
    Y = data_XY(2,:);

    % Fit an ellipse to the transformed data
    try
        [center, semiaxes(:,k), angle(k)] = fitEllipseGAC([X',Y']);
    catch
        disp("Ellipse not found")
        semiaxes(:,k) = semiaxes(:,k-1);
        angle(k) = angle(k-1);
        center(k) = center(k-1);
    end

    if semiaxes(:,k) == [0;0]

        semiaxes(:,k) = semiaxes(:,k-1);
        angle(k) = angle(k-1);

    end

    center(3)= 0;
    center_3D(:,k) = applyRotorToVoltage(center, R.inverse);

    k = k + 1;

    end

end
