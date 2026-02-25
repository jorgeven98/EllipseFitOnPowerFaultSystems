function metrics = compute_detection_metrics(classifications, time_vector, expected_fault, fault_start, fault_end)
% COMPUTE_DETECTION_METRICS Calcula métricas de rendimiento de detección
%
% Input:
%   classifications: Cell array de strings con clasificaciones detectadas
%   time_vector: Vector de tiempos centrales de cada ventana
%   expected_fault: String con tipo de falta esperado
%   fault_start: Tiempo de inicio de falta (s)
%   fault_end: Tiempo de fin de falta (s)
%
% Output:
%   metrics: Estructura con todas las métricas de rendimiento
%
% Autor: Experiment Script
% Fecha: 2026-02-10

%% Inicializar estructura de métricas
metrics = struct();

%% Definir regiones temporales
% Pre-falta: 0 - fault_start
% Durante falta: fault_start - fault_end
% Post-falta: fault_end - fin

idx_pre = time_vector < fault_start;
idx_during = (time_vector >= fault_start) & (time_vector <= fault_end);
idx_post = time_vector > fault_end;

%% 1. PRECISIÓN GENERAL
n_total = length(classifications);
n_correct = sum(strcmp(classifications, expected_fault));
metrics.overall_accuracy = n_correct / n_total;

%% 2. TIEMPO DE DETECCIÓN
% Primera ventana que detecta correctamente DESPUÉS del inicio de falta
idx_correct_during = find(idx_during & strcmp(classifications, expected_fault), 1, 'first');

if ~isempty(idx_correct_during)
    metrics.detection_time = time_vector(idx_correct_during) - fault_start;
else
    % No se detectó la falta
    metrics.detection_time = NaN;
end

%% 3. TASA DE FALSOS POSITIVOS (Pre-falta)
% Definir lo que es "normal" antes de la falta
normal_states = {'A-B-C'};  % Sistema balanceado

n_pre = sum(idx_pre);
if n_pre > 0
    classifications_pre = classifications(idx_pre);
    n_false_positives = sum(~ismember(classifications_pre, normal_states));
    metrics.false_positive_rate = n_false_positives / n_pre;
else
    metrics.false_positive_rate = NaN;
end

%% 4. PRECISIÓN DURANTE FALTA
n_during = sum(idx_during);
if n_during > 0
    classifications_during = classifications(idx_during);
    n_correct_during = sum(strcmp(classifications_during, expected_fault));
    metrics.fault_accuracy = n_correct_during / n_during;
else
    metrics.fault_accuracy = NaN;
end

%% 5. TIEMPO DE RECUPERACIÓN (Post-falta)
% Primera ventana que vuelve a estado normal DESPUÉS de fin de falta
idx_normal_post = find(idx_post & ismember(classifications, normal_states), 1, 'first');

if ~isempty(idx_normal_post)
    metrics.recovery_time = time_vector(idx_normal_post) - fault_end;
else
    % No se recuperó a estado normal
    metrics.recovery_time = NaN;
end

%% 6. MÉTRICAS POR REGIÓN
% Pre-falta
if n_pre > 0
    classifications_pre = classifications(idx_pre);
    metrics.pre_fault_stats = compute_region_stats(classifications_pre, normal_states);
else
    metrics.pre_fault_stats = struct('n_windows', 0, 'accuracy', NaN);
end

% Durante falta
if n_during > 0
    classifications_during = classifications(idx_during);
    metrics.during_fault_stats = compute_region_stats(classifications_during, {expected_fault});
else
    metrics.during_fault_stats = struct('n_windows', 0, 'accuracy', NaN);
end

% Post-falta
if sum(idx_post) > 0
    classifications_post = classifications(idx_post);
    metrics.post_fault_stats = compute_region_stats(classifications_post, normal_states);
else
    metrics.post_fault_stats = struct('n_windows', 0, 'accuracy', NaN);
end

%% 7. MATRIZ DE CONFUSIÓN SIMPLIFICADA
% Solo durante la falta
if n_during > 0
    unique_classes = unique(classifications(idx_during));
    conf_matrix = zeros(length(unique_classes), 1);

    for i = 1:length(unique_classes)
        conf_matrix(i) = sum(strcmp(classifications(idx_during), unique_classes{i}));
    end

    metrics.confusion_matrix.classes = unique_classes;
    metrics.confusion_matrix.counts = conf_matrix;
    metrics.confusion_matrix.expected = expected_fault;
else
    metrics.confusion_matrix = struct('classes', {{}}, 'counts', [], 'expected', expected_fault);
end

%% 8. ESTADÍSTICAS ADICIONALES
metrics.n_windows_total = n_total;
metrics.n_windows_pre = n_pre;
metrics.n_windows_during = n_during;
metrics.n_windows_post = sum(idx_post);

end

function stats = compute_region_stats(classifications, expected_classes)
% Calcular estadísticas para una región específica

n = length(classifications);
n_correct = sum(ismember(classifications, expected_classes));

stats = struct();
stats.n_windows = n;
stats.accuracy = n_correct / n;
stats.classifications = classifications;

end
