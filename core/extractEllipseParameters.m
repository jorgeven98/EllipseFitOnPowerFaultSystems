function [center, axes, angle] = extractEllipseParameters(Q)

    % Extract scalar components of GAC vector
    vx = real(Q(1));
    vn = real(Q(2));
    v1 = real(Q(4));
    v2 = real(Q(5));
    vp = real(Q(6));

    % Calculate ellipse parameters
    alpha = sqrt(vn^2 + vx^2);
    cos2Theta = -vn / alpha;
    sin2Theta = -vx / alpha;
    
    % Major axis rotation angle of ellipse
    angle = real(atan2(sin2Theta, cos2Theta) / 2);
    %majorTheta * 180 / pi;

    % Center of ellipse; 
    center = inv([1+vn, vx; vx, 1-vn])*[v1; v2];

    centerX = center(1);
    centerY = center(2);

    % Signo de los términos de centro corregido (rev. TPWRD 2026): con el signo
    % negativo los ejes salían ~3x en elipses descentradas (inofensivo con centro=0).
    beta = centerX^2 + centerY^2 + (centerX^2 - centerY^2) * vn + 2 * centerX * centerY * vx - 2 * vp;
    % if beta < 0
    %     beta= -beta;
    % end

    % Major axis length of ellipse
    rM = sqrt(beta / (1 - alpha));

    % Minor axis length of ellipse
    rm = sqrt(beta / (1 + alpha));

    axes = real([rM, rm]);

    % if (rM - rm)/rM < 1e-3
    %     angle = 0;
    % end
    if angle < 15*pi/180
        angle = angle + pi;
    end
end