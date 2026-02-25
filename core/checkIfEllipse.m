function isellipse = checkIfEllipse(Q)

    % Extract scalar components of GAC vector
    vx = Q(1);
    vn = Q(2);
    v1 = Q(4);
    v2 = Q(5);
    vp = Q(6);

    if isnan(vx) || isnan(vn) || isnan(v1) || isnan(v2) || isnan(vp)
       isellipse = false;
       return
    end

    % Calculate ellipse parameters
    alpha = sqrt(vn^2+vx^2);

    if alpha == -1 || alpha > 1
       isellipse = false;
       return
    end

    % Center of ellipse; should be near the origin
    center = inv([1+vn, vx; vx, 1-vn])*[v1; v2];

    if isnan(center(1)) || isnan(center(2))
        isellipse = false;
        return
    end

    centerX = center(1);
    centerY = center(2);
    beta = centerX^2 + centerY^2 + (centerX^2 - centerY^2) * vn + 2 * centerX * centerY * vx - 2 * vp;

    if beta <= 0
       isellipse = false;
       return 
    else
        isellipse = true;
    end

    % Major axis length of ellipse
    rM = sqrt(beta / (1 - alpha));
    % Minor axis length of ellipse
    rm = sqrt(beta / (1 + alpha));

    if rm/rM < 1e-3
        isellipse = false;
        return 
    else
        isellipse = true;
    end

end