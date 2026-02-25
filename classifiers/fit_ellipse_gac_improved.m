function [center, axes, angle] = fit_ellipse_gac_improved(points)
    % FITELLIPSE00GAC_IMPROVED Ajuste robusto de elipses centradas en el origen
    % Versión mejorada que maneja correctamente casos degenerados:
    % - Círculos pequeños
    % - Elipses muy aplastadas (líneas)
    % - Datos con ruido numérico
    %
    % Input:
    %   points: Matriz Nx2 donde cada fila es un punto [x, y]
    %
    % Output:
    %   center: Centro de la elipse [cx, cy] (siempre [0,0] en esta versión)
    %   axes:   Semiejes [a, b] donde a >= b
    %   angle:  Ángulo de orientación en radianes

    % Verificar entrada
    arguments
        points (:,2) double
    end

    % Transponer para trabajar con columnas
    points = points';
    n = size(points, 2);

    % VALIDACIÓN 1: Verificar que hay suficientes puntos
    if n < 5
        warning('fitEllipse:InsufficientPoints', ...
            'Se necesitan al menos 5 puntos. Usando fallback para %d puntos.', n);
        [center, axes, angle] = fallback_fit(points);
        return;
    end

    % VALIDACIÓN 2: Detectar caso degenerado (todos los puntos son iguales o colineales)
    point_spread = max(std(points, 0, 2));
    if point_spread < 1e-10
        warning('fitEllipse:DegenerateData', 'Todos los puntos son prácticamente idénticos.');
        center = [0, 0];
        axes = [point_spread, 0];
        angle = 0;
        return;
    end

    % MEJORA 1: Normalización de datos para estabilidad numérica
    % Escalar los puntos al rango [-1, 1]
    scale_factor = max(max(abs(points)));

    if scale_factor < 1e-12
        warning('fitEllipse:TooSmall', 'Puntos demasiado pequeños (escala < 1e-12).');
        [center, axes, angle] = fallback_fit(points);
        return;
    end

    points_normalized = points / scale_factor;

    % Matriz para el espacio GAC inner product
    B = zeros(8,8);
    I3 = [0 0 1; 0 1 0; 1 0 0];
    B(1:3, 6:8) = -I3;
    B(4:5, 4:5) = eye(2);
    B(6:8, 1:3) = -I3;

    % Puntos en notación GAC [0 0 1 x y 0.5(x^2 + y^2) 0.5(x^2 - y^2) xy]
    D = zeros(8, n);
    D(3,:) = ones(1,n);
    D(4:5,:) = points_normalized;
    D(6,:) = 0.5 * (points_normalized(1,:).^2 + points_normalized(2,:).^2);
    D(7,:) = 0.5 * (points_normalized(1,:).^2 - points_normalized(2,:).^2);
    D(8,:) = points_normalized(1,:) .* points_normalized(2,:);

    % Matriz de datos P
    P = 1/n * B * (D * D') * B;

    % Para una elipse centrada en el origen, v¹ = 0 y v² = 0
    % Extraemos solo las filas/columnas correspondientes a [v̄ˣ, v̄⁻, v̄⁺, v⁺]
    indices = [1, 2, 3, 6];
    P_reduced = P(indices, indices);

    % MEJORA 2: Regularización para evitar matrices singulares
    % Agregar pequeña perturbación a la diagonal
    epsilon = 1e-12 * norm(P_reduced, 'fro');
    P_reduced = P_reduced + epsilon * eye(size(P_reduced));

    % Matriz de restricción C para garantizar que el resultado sea una elipse
    % Implementa la restricción (v̄⁺)² - (v̄⁻)² - (v̄ˣ)² = 1
    C = zeros(4, 4);
    C(1, 1) = -1;  % Coeficiente para (v̄ˣ)²
    C(2, 2) = -1;  % Coeficiente para (v̄⁻)²
    C(3, 3) = 1;   % Coeficiente para (v̄⁺)²

    % MEJORA 3: Resolver problema de eigenvalores con manejo de errores
    try
        [EV, ED] = eig(P_reduced, C);
        EW = diag(ED);
    catch ME
        warning('fitEllipse:EigenError', 'Error en cálculo de eigenvalores: %s', ME.message);
        [center, axes, angle] = fallback_fit(points);
        return;
    end

    % MEJORA 4: Selección robusta de eigenvalor
    % Filtrar eigenvalores válidos con tolerancia ajustada
    tolerance = 1e-6;
    valid_indices = find(isfinite(EW) & real(EW) > -tolerance);

    if isempty(valid_indices)
        warning('fitEllipse:NoValidEigenvalues', 'No se encontraron eigenvalores válidos.');
        [center, axes, angle] = fallback_fit(points);
        return;
    end

    % Buscar el eigenvalor positivo más pequeño
    [~, min_idx] = min(real(EW(valid_indices)));
    k_opt = valid_indices(min_idx);
    v_opt = EV(:, k_opt);

    % Reconstruir el vector Q completo en formato GAC
    Q = zeros(8, 1);
    Q(1) = v_opt(1);  % v̄ˣ
    Q(2) = v_opt(2);  % v̄⁻
    Q(3) = v_opt(3);  % v̄⁺
    Q(6) = v_opt(4);  % v⁺

    % MEJORA 5: Normalización segura
    if abs(Q(3)) < 1e-15
        warning('fitEllipse:NormalizationFailed', 'Q(3) muy pequeño para normalizar.');
        [center, axes, angle] = fallback_fit(points);
        return;
    end

    Q = Q / Q(3);

    % Extraer componentes escalares del vector GAC
    vx = real(Q(1));
    vn = real(Q(2));
    vp = real(Q(6));

    % MEJORA 6: Calculo robusto de parametros de elipse
    % Calcular alpha con proteccion contra valores muy pequenos
    alpha = sqrt(vn^2 + vx^2);

    % Centro (siempre origen en esta version)
    centerX = 0;
    centerY = 0;
    center = [centerX, centerY];

    % CASO ESPECIAL: Circulo casi perfecto (alpha muy pequeno)
    % Cuando alpha -> 0, la elipse es un circulo
    CIRCLE_THRESHOLD = 1e-8;

    if alpha < CIRCLE_THRESHOLD
        % Para un circulo centrado en origen: beta = -2*vp
        % y radius = sqrt(beta) = sqrt(-2*vp)
        beta = -2 * vp;

        if beta > 0
            radius = sqrt(beta) * scale_factor;
            axes = real([radius, radius]);
            angle = 0;
            return;
        else
            % Datos problematicos, usar fallback
            [center, axes, angle] = fallback_fit(points);
            return;
        end
    end

    % Caso general: elipse no circular
    cos2Theta = -vn / alpha;
    sin2Theta = -vx / alpha;
    angle = real(atan2(sin2Theta, cos2Theta) / 2);

    % MEJORA 7: Calculo protegido de beta
    % Para centro en origen: beta = -2*vp
    beta = -2 * vp;

    % Verificar beta positivo
    if beta <= 0
        % Intentar recuperar con valor absoluto si es pequeno
        if abs(beta) < 1e-10
            [center, axes, angle] = fallback_fit(points);
            return;
        end
        beta = abs(beta);
    end

    % MEJORA 8: Calculo seguro de semiejes
    % Para elipse: rM^2 = beta/(1-alpha), rm^2 = beta/(1+alpha)
    denominator_major = 1 - alpha;
    denominator_minor = 1 + alpha;

    % Proteger contra division por valores muy pequenos
    MIN_DENOM = 1e-10;
    if denominator_major < MIN_DENOM
        % Alpha muy cercano a 1: elipse muy elongada
        denominator_major = MIN_DENOM;
    end

    rM = sqrt(beta / denominator_major);
    rm = sqrt(beta / denominator_minor);

    % MEJORA 9: Des-normalizar los semiejes
    rM = rM * scale_factor;
    rm = rm * scale_factor;

    % Asegurar que rM >= rm (eje mayor >= eje menor)
    if rm > rM
        [rM, rm] = deal(rm, rM);
        angle = angle + pi/2;
    end

    axes = real([rM, rm]);

    % MEJORA 10: Validación final de resultados
    if any(~isfinite(axes)) || any(axes < 0)
        warning('fitEllipse:InvalidResult', 'Resultado inválido. Usando fallback.');
        [center, axes, angle] = fallback_fit(points);
        return;
    end

    % Normalizar ángulo al rango [0, pi)
    angle = mod(real(angle), pi);

    % Ajuste de ángulo si es muy pequeño
    if angle < 15*pi/180
        angle = angle + pi;
    end
    angle = mod(angle, pi);
end

function [center, axes, angle] = fallback_fit(points)
    % FALLBACK_FIT Ajuste alternativo usando PCA para casos problemáticos
    %
    % Usa análisis de componentes principales para estimar la orientación
    % y los semiejes de la elipse

    % Centrar los datos (aunque deberían estar centrados en origen)
    points_mean = mean(points, 2);
    points_centered = points - points_mean;

    % Calcular matriz de covarianza
    C = cov(points_centered');

    % Eigenvalores y eigenvectores
    [V, D] = eig(C);
    eigenvalues = diag(D);

    % Ordenar por magnitud (mayor primero)
    [eigenvalues, idx] = sort(eigenvalues, 'descend');
    V = V(:, idx);

    % Semiejes basados en desviaciones estándar
    % Factor 2 para aproximar el semieje (contiene ~95% de los datos)
    axes = 2 * sqrt(max(eigenvalues, 0));

    % Ángulo del primer componente principal
    angle = atan2(V(2,1), V(1,1));

    % Normalizar ángulo
    angle = mod(angle, pi);

    % Centro (origen)
    center = [0, 0];

    % Validar resultado
    if any(~isfinite(axes))
        axes = [1, 0.5];  % Default razonable
    end

    if ~isfinite(angle)
        angle = 0;
    end
end
