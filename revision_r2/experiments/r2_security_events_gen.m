function E = r2_security_events_gen(nPerType)
% R2_SECURITY_EVENTS_GEN  Generador de la suite de eventos de seguridad (T41, ronda 2).
%
% Copia del generador synth_event del proyecto DGC
%   code/fault_detection/matlab/experiments/exp_security_events.m
% (solo lectura alli). Se conserva rng(7) y el ORDEN EXACTO de las llamadas a
% rand/randn, incluidas las del canal de corriente que aqui no se usa, para
% obtener los mismos 320 eventos (8 familias x 40). La unica diferencia es que
% se devuelven los parametros sorteados de cada evento.
%
% Ningun evento es una falta: escalon de carga, inrush, condensadores,
% arranque de motor, hueco remoto equilibrado, escalon de desequilibrio,
% deriva de frecuencia y salto de fase. Todos llevan 1 % de desequilibrio
% permanente de modulo, armonicos 5 y 7 (3 y 2 %) y ruido de 60 dB.
% Tension en pu (pico de fase 1), fs = 4 kHz, 2 s, evento en t = 1.0 s.
%
%   E = r2_security_events_gen()      % 40 por familia
%   E(k).kind, .q, .t, .v [N x 3], .tE, .par (struct con lo sorteado)
if nargin < 1 || isempty(nPerType), nPerType = 40; end
rng(7);
types = {'load_step','inrush','cap_switch','motor_start','remote_sag', ...
         'unbal_step','freq_drift','phase_jump'};
dt = 250e-6; E = struct('kind', {}, 'q', {}, 't', {}, 'v', {}, 'tE', {}, 'par', {});
for ty = 1:numel(types)
    for q = 1:nPerType
        [v, ~, t, tE, par] = synth_event(types{ty}, dt);
        E(end+1) = struct('kind', types{ty}, 'q', q, 't', t, 'v', v, 'tE', tE, 'par', par); %#ok<AGROW>
    end
end
end

% ------------------------------------------------------------- generador ----
function [v,i,t,tE,par] = synth_event(kind, dt)
    f0=50; Tend=2.0; t=(0:dt:Tend).'; N=numel(t); tE=1.0;
    ph0=[0 -2*pi/3 2*pi/3];
    Vmag = 1 + 0.01*[1 -0.4 -0.6];
    Imag = 0.9 + 0.01*[0.5 -0.5 0];
    phiI = 0.2;
    h5 = 0.03; h7 = 0.02;
    after = t>=tE;
    par = struct('p1', NaN, 'p2', NaN, 'p3', NaN, 'desc', '');

    f = f0*ones(N,1);
    if strcmp(kind,'freq_drift')
        fEnd = 49.3 + 0.4*rand;
        f(after) = f0 + (fEnd-f0) * min(1,(t(after)-tE)/0.8);
        par.p1 = fEnd; par.desc = 'p1 = frecuencia final (Hz)';
    end
    phase = 2*pi*cumsum(f)*dt;

    v=zeros(N,3); i=zeros(N,3);
    for k=1:3
        vk = Vmag(k)*( sin(phase+ph0(k)) + h5*sin(5*(phase+ph0(k))) + h7*sin(7*(phase+ph0(k))) );
        ik = Imag(k)*( sin(phase+ph0(k)-phiI) + 2*h5*sin(5*(phase+ph0(k))) );
        v(:,k)=vk; i(:,k)=ik;
    end

    switch kind
        case 'load_step'
            stp = (0.2+0.4*rand) * sign(rand-0.3);
            g = 1 + stp*after;
            i = i.*g;
            v = v.*(1 - 0.02*abs(stp)*after);
            par.p1 = stp; par.p2 = 0.02*abs(stp); par.desc = 'p1 = escalon de carga (pu); p2 = caida de tension (pu)';
        case 'inrush'
            tau = 0.1+0.3*rand; A = 1+2*rand;
            for k=1:3
                phk = ph0(k)+2*pi*rand*0.2;
                burst = A*(0.5+0.5*rand)*exp(-(t-tE)/tau) .* max(0,sin(phase+phk)).^3;
                i(:,k) = i(:,k) + after.*burst;
            end
            par.p1 = A; par.p2 = tau; par.desc = 'p1 = pico de corriente (pu); p2 = tau (s); la tension no cambia';
        case 'cap_switch'
            fc = 400+500*rand; tau = 0.005+0.015*rand; A = 0.05+0.10*rand;
            osc = exp(-(t-tE)/tau).*sin(2*pi*fc*(t-tE));
            for k=1:3
                v(:,k) = v(:,k) + after.*A.*osc*cos(ph0(k));
                i(:,k) = i(:,k) + after.*(A*fc/f0/8).*exp(-(t-tE)/tau).*cos(2*pi*fc*(t-tE))*cos(ph0(k));
            end
            par.p1 = A; par.p2 = fc; par.p3 = tau; par.desc = 'p1 = amplitud (pu); p2 = frecuencia (Hz); p3 = tau (s)';
        case 'motor_start'
            A = 2.5+1.5*rand; tau = 0.3+0.5*rand;
            g = 1 + after.*( (A-1).*exp(-(t-tE)/tau) + 0.2 );
            i = i.*g;
            sag = 0.05+0.07*rand;
            v = v.*(1 - after.*sag.*exp(-(t-tE)/(0.8*tau)));
            par.p1 = sag; par.p2 = tau; par.p3 = A; par.desc = 'p1 = hueco inicial (pu); p2 = tau (s); p3 = corriente (pu)';
        case 'remote_sag'
            depth = 0.2+0.2*rand; dur = 0.1+0.2*rand;
            ramp = 0.005;
            prof = 1 - depth*0.5*(tanh((t-tE)/ramp) - tanh((t-tE-dur)/ramp));
            v = v.*prof;
            par.p1 = depth; par.p2 = dur; par.desc = 'p1 = profundidad (pu); p2 = duracion (s)';
        case 'unbal_step'
            drop = 0.02+0.03*rand; kph = randi(3);
            v(:,kph) = v(:,kph).*(1 - drop*after);
            i(:,kph) = i(:,kph).*(1 - 0.5*drop*after);
            par.p1 = drop; par.p2 = kph; par.desc = 'p1 = caida de una fase (pu); p2 = fase (1 = A)';
        case 'freq_drift'
            % ya aplicado en el perfil de fase
        case 'phase_jump'
            jmp = (5+10*rand)*pi/180 * sign(rand-0.5);
            phase2 = phase + jmp*after;
            for k=1:3
                v(:,k) = Vmag(k)*( sin(phase2+ph0(k)) + h5*sin(5*(phase2+ph0(k))) + h7*sin(7*(phase2+ph0(k))) );
                i(:,k) = Imag(k)*( sin(phase2+ph0(k)-phiI) + 2*h5*sin(5*(phase2+ph0(k))) );
            end
            par.p1 = jmp*180/pi; par.desc = 'p1 = salto de fase (grados)';
    end
    % ruido de medida de 60 dB (mismas llamadas que el original)
    v = v + 10^(-60/20)*std(v(:))*randn(N,3);
    i = i + 10^(-60/20)*std(i(:))*randn(N,3); %#ok<NASGU>
end
