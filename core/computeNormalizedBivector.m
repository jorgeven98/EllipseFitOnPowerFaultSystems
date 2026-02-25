function bivector = computeNormalizedBivector(point1, point2)
    % COMPUTENORMALIZEDBIVECTOR Computes the normalized bivector from two 3D points
    %
    % Inputs:
    %   point1: Vector [x1, y1, z1] of the first point
    %   point2: Vector [x2, y2, z2] of the second point
    %
    % Output:
    %   bivector: Normalized bivector object (outer product of the two vectors)

    % Validate inputs
    if length(point1) ~= 3 || length(point2) ~= 3
        error('Points must have 3 coordinates (x,y,z)');
    end

    % Encode points as vectors in geometric algebra
    % If points are relative to the origin, vectors are directly the point coordinates
    v1 = ga3.EncodeVector(point1);
    v2 = ga3.EncodeVector(point2);

    % Compute the outer product (wedge product)
    bv = v1.op(v2);

    % Compute the bivector norm
    bv_norm = bv.norm();

    % Normalize the bivector if its norm is not too small
    if bv_norm > 1e-10  % Threshold to avoid division by zero
        bivector = bv * (1/bv_norm);
    else
        % If vectors are (nearly) linearly dependent, norm will be close to zero
        % warning('Vectors are nearly collinear, bivector is not well-defined');
        bivector = bv; % Return unnormalized bivector
    end
    %
    % % Display bivector information (optional)
    % fprintf('Normalized bivector components:\n');
    %
    % % Get bivector components
    % components = bivector.getDataArray();
    %
    % % Convention: [scalar, e1, e2, e3, e12, e23, e31, e123]
    % fprintf('  σ12 (xy): %.6f\n', components(1));
    % fprintf('  σ23 (yz): %.6f\n', components(3));
    % fprintf('  σ31 (zx): %.6f\n', -components(2));
end
