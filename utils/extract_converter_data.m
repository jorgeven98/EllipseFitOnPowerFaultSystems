function [V_abc, time_vector] = extract_converter_data(data, converter_name)
% EXTRACT_CONVERTER_DATA Extrae datos de voltaje de un convertidor específico
%
% Input:
%   data: Estructura cargada desde archivo .mat (campo GridConvs)
%   converter_name: 'DG1', 'DG2', o 'Pex'
%
% Output:
%   V_abc: Matriz [N x 3] de voltajes trifásicos
%   time_vector: Vector [N x 1] de tiempos correspondientes
%
% Autor: Experiment Script
% Fecha: 2026-02-10

% Validar entrada
if ~isfield(data, 'GridConvs')
    error('extract_converter_data:NoGridConvs', ...
        'El archivo .mat no contiene el campo GridConvs');
end

GridConvs = data.GridConvs;

% Validar que el convertidor existe
if ~isfield(GridConvs, converter_name)
    error('extract_converter_data:ConverterNotFound', ...
        'Convertidor %s no encontrado en GridConvs', converter_name);
end

converter_data = GridConvs.(converter_name);

% Extraer voltajes según el convertidor
switch converter_name
    case {'DG1', 'DG2'}
        % DG1 y DG2 usan campo Vac
        if ~isfield(converter_data, 'Vac')
            error('extract_converter_data:NoVac', ...
                'Convertidor %s no tiene campo Vac', converter_name);
        end

        V_timeseries = converter_data.Vac;
        V_abc = V_timeseries.Data;
        time_vector = V_timeseries.Time;

    case 'Pex'
        % Pex usa campo Vabc
        if ~isfield(converter_data, 'Vabc')
            error('extract_converter_data:NoVabc', ...
                'Convertidor Pex no tiene campo Vabc');
        end

        V_timeseries = converter_data.Vabc;
        V_abc = V_timeseries.Data;
        time_vector = V_timeseries.Time;

    otherwise
        error('extract_converter_data:InvalidConverter', ...
            'Convertidor desconocido: %s. Use DG1, DG2 o Pex', converter_name);
end

% Validar dimensiones
if size(V_abc, 2) ~= 3
    error('extract_converter_data:InvalidDimensions', ...
        'V_abc debe tener 3 columnas (fases A, B, C), tiene %d', size(V_abc, 2));
end

% Asegurar que time_vector es columna
if size(time_vector, 2) > size(time_vector, 1)
    time_vector = time_vector';
end

end
