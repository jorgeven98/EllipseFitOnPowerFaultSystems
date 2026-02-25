function [center, axes, angle] = fitEllipseGAC(points)
    % FITELLIPSEGAC Fits a 2D ellipse using the GAC (Geometric Algebra for Conics) method
    % points: Nx2 matrix where each row is a point [x, y]
    % Returns center, semi-axes, and angle of the fitted ellipse

    % Validate input
    arguments
        points (:,2) double
    end
    points = points';

    % Number of points
    n = size(points, 2);

    % Matrix denoting CGA space
    I3 =[0 0 1; 0 1 0; 1 0 0];
    B(1:3, 6:8) = -I3;
    B(4:5, 4:5) = eye(2);
    B(6:8, 1:3) = -I3;
    Bc = B(3:6, 3:6);

    % Points in CGA notation [0 0 1 x y 0.5(x^2 + y^2) 0.5(x^2 - y^2) xy]
    D = zeros([8, n]);
    D(3,:) = ones(1,n);
    D(4:5,:) = points;
    D(6,:) = 0.5 * (points(1,:).^2 + points(2,:) .^ 2);
    D(7,:) = 0.5 * (points(1,:).^2 - points(2,:) .^ 2);
    D(8,:) = points(1,:) .* points(2,:);

    P = 1/n * B * (D * D') * B;
    Pc = P( 3 : 6 , 3 : 6 );
    P0 = P( 1 : 2 , 1 : 2 );
    P1 = P( 1 : 2 , 3 : 6 );

    Pcon = Bc * (Pc - P1' / P0 * P1);

    % Constraint matrix C to ensure the result is an ellipse
    % Implements the constraint (v̄⁺)² - (v̄⁻)² - (v̄ˣ)² = 1
    % Equivalent to 4AC - B² = 1 in standard form
    C = zeros(4, 4);
    C(1, 1) = -0.5;  % Coefficient for (v̄ˣ)²
    C(2, 2) = -0.5;  % Coefficient for (v̄⁻)²
    C(3, 3) = 0.5;   % Coefficient for (v̄⁺)²

    [EV, ED] = eig(Pcon,C);
    EW = diag(ED);

    % Find the smallest positive eigenvalue
    valid_indices = find(isfinite(EW));
    [~, min_idx] = min(abs(EW));

    k_opt = valid_indices(min_idx);
    v_opt = EV(:,k_opt);
    w = -P0 \ P1 * v_opt;

    Q = zeros(8, 1);
    Q(1) = w(1);      % v̄ˣ
    Q(2) = w(2);      % v̄⁻
    Q(3) = v_opt(1);  % v̄⁺
    Q(4) = v_opt(2);  % v¹
    Q(5) = v_opt(3);  % v²
    Q(6) = v_opt(4);  % v⁺

    Q = Q/Q(3);

    [center, axes, angle] = extractEllipseParameters(Q);

end
