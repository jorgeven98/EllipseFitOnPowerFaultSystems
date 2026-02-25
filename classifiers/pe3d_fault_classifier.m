function [fault_type, params] = pe3d_fault_classifier(V_abc)
% PE3D_FAULT_CLASSIFIER Clasificación de faltas mediante parámetros de elipse 3D
% Basado en: Alam et al., "Classification and Localization of Fault-Initiated
% Voltage Sags Using 3-D Polarization Ellipse Parameters", IEEE Trans. Power Del., 2020
%
% Input:
%   V_abc: Matriz [N x 3] de voltajes trifásicos (Va, Vb, Vc) en p.u.
%
% Output:
%   fault_type: String con tipo de falta
%   params: Estructura con parámetros de la elipse de polarización

%% 1. EXTRAER FASORES (1 ciclo de datos)
N = size(V_abc, 1);

% Calcular fasores mediante análisis de Fourier
Va_phasor = (2/N) * sum(V_abc(:,1) .* exp(-1j*2*pi*(0:N-1)'/N));
Vb_phasor = (2/N) * sum(V_abc(:,2) .* exp(-1j*2*pi*(0:N-1)'/N));
Vc_phasor = (2/N) * sum(V_abc(:,3) .* exp(-1j*2*pi*(0:N-1)'/N));

% Magnitudes y ángulos
V_mag = [abs(Va_phasor), abs(Vb_phasor), abs(Vc_phasor)];
V_ang = [angle(Va_phasor), angle(Vb_phasor), angle(Vc_phasor)];

%% 2. CALCULAR VECTOR NORMAL AL PLANO DE POLARIZACIÓN
Wx = -2 * V_mag(2) * V_mag(3) * sin(V_ang(2) - V_ang(3));
Wy = -2 * V_mag(3) * V_mag(1) * sin(V_ang(3) - V_ang(1));
Wz = -2 * V_mag(1) * V_mag(2) * sin(V_ang(1) - V_ang(2));

%% 3. CALCULAR ÁNGULOS DE AZIMUTH Y ELEVACIÓN
theta = atan2(sqrt(Wx^2 + Wy^2), Wz);
phi = atan2(Wy, Wx);

theta_deg = rad2deg(theta);
phi_deg = rad2deg(phi);

%% 4. CALCULAR COMPONENTES TRANSVERSALES
Va = V_mag(1) * exp(1j*V_ang(1));
Vb = V_mag(2) * exp(1j*V_ang(2));
Vc = V_mag(3) * exp(1j*V_ang(3));

S_theta = Va*cos(theta)*cos(phi) + Vb*cos(theta)*sin(phi) - Vc*sin(theta);
S_phi = -Va*sin(phi) + Vb*cos(phi);

%% 5. CALCULAR PARÁMETROS DE STOKES Y ELIPSE
beta = angle(S_phi) - angle(S_theta);

numerator = 2 * abs(S_theta) * abs(S_phi) * cos(beta);
denominator = abs(S_theta)^2 - abs(S_phi)^2;

if abs(denominator) < 1e-10
    psi = 0;
else
    psi = 0.5 * atan2(numerator, denominator);
end
psi_abs = abs(psi);
psi_abs_deg = rad2deg(psi_abs);

% Semiejes mayor y menor
g0 = abs(S_theta)^2 + abs(S_phi)^2;
sin_arg = 2 * abs(S_theta) * abs(S_phi) * sin(beta) / (g0 + 1e-10);
sin_arg = max(-1, min(1, sin_arg));

ellipticity_arg = 0.5 * asin(sin_arg);

Ax = sqrt(g0 + 1e-10) / (sin(ellipticity_arg) + 1e-10);
Ay = sqrt(g0 + 1e-10) / (cos(ellipticity_arg) + 1e-10);

if Ay > Ax
    temp = Ax;
    Ax = Ay;
    Ay = temp;
end

%% 6. CLASIFICACIÓN
Ay_Ax_ratio = Ay / (Ax + 1e-10);

% MODIFICACIÓN: Umbrales ajustados para coherencia con GA
% Criterio de circularidad: Ay_Ax_ratio ≈ 0.99 (equivalente a roundness > 0.99 en GA)
TH_NORMAL = 0.99;      % MODIFICADO: 0.933 → 0.99 (coherente con TH_BALANCED de GA)
TH_AY_MIN_RATIO = 0.90; % MODIFICADO: Nuevo umbral para magnitud reducida (coherente con GA)
TH_AY_MAX = 1.2;

% MODIFICACIÓN: Lógica mejorada para distinguir Normal vs ABCG usando ratio de magnitud
if Ay_Ax_ratio >= TH_NORMAL
    % Sistema circular (alta circularidad)
    % Distinguir Normal vs ABCG por magnitud del semieje mayor
    % Valor nominal esperado Ay ≈ 1.1 - 1.2 para sistema normal
    magnitude_ratio = Ay / 1.15;  % 1.15 como valor nominal promedio

    if magnitude_ratio >= TH_AY_MIN_RATIO
        % Magnitud normal (≥90% del nominal) → Sistema balanceado sin falta
        fault_type = 'Normal';
    else
        % Magnitud reducida (<90% del nominal) → Falta trifásica simétrica ABCG
        fault_type = 'ABCG';
    end

elseif Ay_Ax_ratio < TH_NORMAL
    % Sistema elíptico o lineal → Faltas asimétricas (L-G, L-L, L-L-G)
    fault_type = classify_asymmetric_faults_paper_method(psi_abs_deg, phi_deg, Ax, theta_deg, V_abc);
else
    fault_type = 'Normal';
end

%% 7. EMPAQUETAR PARÁMETROS
params.Ax = Ax;
params.Ay = Ay;
params.theta = theta_deg;
params.phi = phi_deg;
params.psi = rad2deg(psi);
params.psi_abs = psi_abs_deg;
params.Ay_Ax_ratio = Ay_Ax_ratio;
params.V_mag = V_mag;
params.fault_type = fault_type;

end

function fault_type = classify_asymmetric_faults_paper_method(psi_abs, phi, Ax, theta, V_abc)
% Clasificación EXACTA según Fig. 6 del paper
% Usa fronteras de decisión de Tabla II + análisis de fases para desambiguar

phi = mod(phi + 180, 360) - 180;

%% Análisis de fases afectadas para desambiguación
Va_rms = rms(V_abc(:,1));
Vb_rms = rms(V_abc(:,2));
Vc_rms = rms(V_abc(:,3));
V_rms = [Va_rms, Vb_rms, Vc_rms];

% Detectar cuántas fases tienen sag (< 85% del valor máximo)
V_max = max(V_rms);
threshold = 0.85;
phase_with_sag = V_rms < (threshold * V_max);
num_phases_with_sag = sum(phase_with_sag);

%% Fronteras de decisión (Tabla II del paper)
DB_CG_CAG = 0.0003885*phi^2 - 0.2182*phi + 38.14;
DB_CAG_AG = -0.003477*phi^2 - 1.617*phi - 110.9;
DB_AG_ABG = -0.003866*phi^2 + 1.399*phi - 104;
DB_CG_BCG = 0.0003885*phi^2 + 0.4284*phi + 125.4;
DB_BCG_BG = -0.003477*phi^2 - 0.2605*phi + 72.27;
DB_BG_ABG = -0.003866*phi^2 - 0.6885*phi + 8.13;

%% CLASIFICACIÓN POR ZONAS DE PHI (Fig. 6)

if phi >= -180 && phi < -135
    % Zona: AG, CAG, CG, AB/ABG
    if psi_abs > DB_CG_CAG
        fault_type = 'CG';
    elseif psi_abs > DB_CAG_AG
        fault_type = 'CAG';
    elseif psi_abs > DB_AG_ABG
        fault_type = 'AG';
    else
        % Zona AB/ABG: CLAVE - verificar número de fases con sag
        if num_phases_with_sag >= 2
            fault_type = 'ABG'; % 2+ fases con sag → con tierra
        else
            fault_type = classify_AB_ABG_paper(Ax, theta);
        end
    end
    
elseif phi >= -135 && phi <= -90
    % Zona: BG, BCG, CG, AB/ABG, BC/BCG
    
    if abs(phi + 135) < 1 % φ ≈ -135°
        % Punto especial en -135°
        if psi_abs > 80
            fault_type = 'CG';
        elseif psi_abs < 10
            if num_phases_with_sag >= 2
                fault_type = 'ABG';
            else
                fault_type = classify_AB_ABG_paper(Ax, theta);
            end
        else
            fault_type = 'Unclassified';
        end
    elseif psi_abs > DB_CG_BCG
        fault_type = 'CG';
    elseif psi_abs > DB_BCG_BG
        % Zona BC/BCG: verificar número de fases con sag
        if num_phases_with_sag >= 2 && phase_with_sag(2) && phase_with_sag(3)
            fault_type = 'BCG'; % B y C con sag → con tierra
        else
            fault_type = classify_BC_BCG_paper(Ax, theta, V_abc);
        end
    elseif psi_abs > DB_BG_ABG
        fault_type = 'BG';
    else
        if num_phases_with_sag >= 2
            fault_type = 'ABG';
        else
            fault_type = classify_AB_ABG_paper(Ax, theta);
        end
    end
    
elseif phi > -90 && phi < -45
    % Zona: BC/BCG
    if num_phases_with_sag >= 2 && phase_with_sag(2) && phase_with_sag(3)
        fault_type = 'BCG';
    else
        fault_type = classify_BC_BCG_paper(Ax, theta, V_abc);
    end
    
elseif phi >= -45 && phi < 45
    % Zona: AB o CA (dependiendo del signo de φ)
    if phi < 0 && phi >= -45
        % Cerca de -45° → AB
        fault_type = classify_AB_ABG_paper(Ax, theta);
    elseif phi >= 0 && phi < 45
        % Cerca de 0° → CA
        if num_phases_with_sag >= 2 && phase_with_sag(3) && phase_with_sag(1)
            fault_type = 'CAG';
        else
            fault_type = classify_CA_CAG_paper(Ax, theta, V_abc);
        end
    else
        fault_type = classify_AB_ABG_paper(Ax, theta);
    end
    
elseif phi >= 45 && phi < 135
    % Zona: CA/CAG
    if num_phases_with_sag >= 2 && phase_with_sag(3) && phase_with_sag(1)
        fault_type = 'CAG';
    else
        fault_type = classify_CA_CAG_paper(Ax, theta, V_abc);
    end
    
else
    % φ >= 135 o φ < -180: zonas extremas
    if phi >= 135
        fault_type = classify_CA_CAG_paper(Ax, theta, V_abc);
    else
        fault_type = 'Unclassified';
    end
end

end

function fault_type = classify_AB_ABG_paper(Ax, theta)
% Clasificar AB vs ABG según Fig. 7 del paper

DB_Ax = 0.0001698*theta^2 - 0.0327*theta + 2.62;

if Ax >= DB_Ax
    fault_type = 'AB';
else
    fault_type = 'ABG';
end
end

function fault_type = classify_BC_BCG_paper(Ax, theta, V_abc)
% Clasificar BC vs BCG con remapeo

N = size(V_abc, 1);
V_remapped = [V_abc(:,2), V_abc(:,3), V_abc(:,1)];

[Ax_r, theta_r] = recalculate_params(V_remapped);

DB_Ax = 0.0001698*theta_r^2 - 0.0327*theta_r + 2.62;

if Ax_r >= DB_Ax
    fault_type = 'BC';
else
    fault_type = 'BCG';
end
end

function fault_type = classify_CA_CAG_paper(Ax, theta, V_abc)
% Clasificar CA vs CAG con remapeo

N = size(V_abc, 1);
V_remapped = [V_abc(:,3), V_abc(:,1), V_abc(:,2)];

[Ax_r, theta_r] = recalculate_params(V_remapped);

DB_Ax = 0.0001698*theta_r^2 - 0.0327*theta_r + 2.62;

if Ax_r >= DB_Ax
    fault_type = 'CA';
else
    fault_type = 'CAG';
end
end

function [Ax, theta_deg] = recalculate_params(V_abc)
% Recalcular Ax y theta con voltajes remapeados

N = size(V_abc, 1);

Va = (2/N) * sum(V_abc(:,1) .* exp(-1j*2*pi*(0:N-1)'/N));
Vb = (2/N) * sum(V_abc(:,2) .* exp(-1j*2*pi*(0:N-1)'/N));
Vc = (2/N) * sum(V_abc(:,3) .* exp(-1j*2*pi*(0:N-1)'/N));

V_mag = [abs(Va), abs(Vb), abs(Vc)];
V_ang = [angle(Va), angle(Vb), angle(Vc)];

Wx = -2 * V_mag(2) * V_mag(3) * sin(V_ang(2) - V_ang(3));
Wy = -2 * V_mag(3) * V_mag(1) * sin(V_ang(3) - V_ang(1));
Wz = -2 * V_mag(1) * V_mag(2) * sin(V_ang(1) - V_ang(2));

theta = atan2(sqrt(Wx^2 + Wy^2), Wz);
phi = atan2(Wy, Wx);

Va_c = V_mag(1) * exp(1j*V_ang(1));
Vb_c = V_mag(2) * exp(1j*V_ang(2));
Vc_c = V_mag(3) * exp(1j*V_ang(3));

S_theta = Va_c*cos(theta)*cos(phi) + Vb_c*cos(theta)*sin(phi) - Vc_c*sin(theta);
S_phi = -Va_c*sin(phi) + Vb_c*cos(phi);

beta = angle(S_phi) - angle(S_theta);

g0 = abs(S_theta)^2 + abs(S_phi)^2;
sin_arg = 2 * abs(S_theta) * abs(S_phi) * sin(beta) / (g0 + 1e-10);
sin_arg = max(-1, min(1, sin_arg));

ellipticity_arg = 0.5 * asin(sin_arg);
Ax = sqrt(g0 + 1e-10) / (sin(ellipticity_arg) + 1e-10);

theta_deg = rad2deg(theta);
end