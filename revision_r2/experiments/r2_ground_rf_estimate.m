function E = r2_ground_rf_estimate(Vw, fs, f, metodo)
% R2_GROUND_RF_ESTIMATE  u0 = V0/V1 y u2 = V2/V1 de una ventana (T23b).
% En los dos metodos se estima el vector complejo U de la elipse del lugar
% geometrico, x(t) = Re{U e^{j tau}}: Re U y -Im U son dos semidiametros
% conjugados de la elipse, y sus componentes simetricas dan V0, V1 y V2 salvo
% una fase comun que no afecta a los cocientes.
% metodo 'conj' (por defecto): semidiametros conjugados separados un cuarto de
%   ciclo, que son los dos vectores x1 y x2 cuyo producto exterior da el
%   bivector del paper, estimados por minimos cuadrados con todas las muestras
%   de la ventana (x(t) = x1 cos wt + x2 sin wt). Usa la frecuencia nominal,
%   igual que la separacion de un cuarto de ciclo del clasificador.
% metodo 'geom': solo geometria, sin usar la frecuencia, con el plano del
%   lugar geometrico y la elipse ajustada en ese plano:
%   1. Plano por el origen que mejor ajusta la ventana (minimos cuadrados; su
%      normal es el dual del bivector del plano), orientado con el sentido de
%      recorrido.
%   2. Elipse centrada en ese plano con fit_ellipse_gac_improved (el ajuste
%      del paper): semiejes a >= b y orientacion theta.
%   3. Los semiejes como vectores de R^3, A (mayor) y B (menor, a +90 grados
%      de A en el sentido de recorrido), dan el vector complejo U = A - jB,
%      con x(t) = Re{U e^{j tau}}. Sus componentes simetricas dan V0, V1, V2
%      salvo una fase comun, que no afecta a los cocientes.
% Devuelve ademas rK = redondez del lugar proyectado en el plano de Kirchhoff,
% (1 - |u2|)/(1 + |u2|).
if nargin < 4 || isempty(metodo), metodo = 'conj'; end
a = exp(2j*pi/3);
E = struct('u0', NaN, 'u2', NaN, 'rK', NaN, 'U', nan(1,3), 'a', NaN, 'b', NaN);
switch lower(metodo)
    case 'geom'
        [~, ~, Vv] = svd(Vw, 'econ');
        n = Vv(:,3)'; e1 = Vv(:,1)';
        c = cross(Vw(1:end-1,:), Vw(2:end,:), 2);
        if sum(c*n') < 0, n = -n; end
        e2 = cross(n, e1);
        X = Vw*[e1', e2'];
        ws = warning('off', 'all');
        try
            [~, ax, th] = fit_ellipse_gac_improved(X);
        catch
            ax = [NaN NaN]; th = NaN;
        end
        warning(ws);
        aE = ax(1); bE = ax(2);
        if ~(isfinite(aE) && isfinite(bE) && isfinite(th)) || aE <= 0
            th = 0; aE = max(abs(X(:,1))); bE = 0;      % recta: direccion principal
        end
        if bE > aE, [aE, bE] = deal(bE, aE); th = th + pi/2; end
        A = aE*(cos(th)*e1 + sin(th)*e2); B = bE*(-sin(th)*e1 + cos(th)*e2);
        U = A - 1j*B;
        E.a = aE; E.b = bE;
    case {'conj','fasor'}
        nn = size(Vw, 1); t = (0:nn-1)'/fs;
        M = [cos(2*pi*f*t), -sin(2*pi*f*t)];
        Xc = M \ Vw; U = Xc(1,:) + 1j*Xc(2,:);
    otherwise
        error('metodo desconocido');
end
V0 = sum(U)/3; V1 = (U(1) + a*U(2) + a^2*U(3))/3; V2 = (U(1) + a^2*U(2) + a*U(3))/3;
E.U = U; E.u0 = V0/V1; E.u2 = V2/V1;
E.rK = (1 - abs(E.u2))/(1 + abs(E.u2));
end
