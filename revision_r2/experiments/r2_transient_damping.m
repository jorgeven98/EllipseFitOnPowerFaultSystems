function T = r2_transient_damping(casedir, outcsv, figdir)
% R2_TRANSIENT_DAMPING  Characterize the fault-inception transient of the EMT records
% of the resistance sweeps (S2a, S2d: resistance per phase; S2c: ground resistance)
% and of the inception-angle sweep (S2b): 85 records, two measuring points each
% (Pex bus = fault point, DG1 terminals), 170 events.
%
% The voltages logged at 4 kHz (what the classifier sees) are split, over the first
% 140 ms of the fault (it lasts 160 ms), into
%   v(t) = v_slow(t) + r(t).
% v_slow is a least-squares fit of
%   - a fundamental whose amplitude and phase may change with time: cos/sin times
%     1, t, t^2, exp(-t/10 ms) and exp(-t/30 ms);
%   - an offset with a linear trend and the same two exponentials;
%   - harmonics 2 to 7 of constant amplitude.
% It holds the fault steady state and the response of the converter controls, which
% modulates the fundamental over tens of milliseconds. r(t) is the fast inception
% transient that remains (oscillation of the network, inception spike).
%
% Two results per event (amplitudes in p.u. of the pre-fault phase peak):
% (1) fast transient r(t), observed during 100 ms
%   rpk_pu     largest |r| in the first 20 ms
%   fdom_Hz    dominant frequency of r in the first 40 ms, searched above 100 Hz
%              (NaN if rpk < 5 %)
%   nhalf      half-cycles of r (worst phase) that exceed +-5 %
%   t5_ms      end of the last 2.5 ms block whose peak |r| exceeds 5 %  (0 if none)
%   regime     'oscillatory'  : at least 4 half-cycles above 5 % (two periods)
%              'damped'       : 1 to 3 half-cycles above 5 %
%              'none'         : r never exceeds 5 %
% (2) fundamental-frequency signature (model-free): phasors fitted on a sliding
%     window of one cycle (step 2.5 ms) lying inside the fault, reduced to the
%     quantities that do not depend on a common rotation of the three phasors,
%     g = [ |V1|/Vb, V2/V1, V0/V1 ]  (the two ratios are complex),
%     and compared with g in the window 120 to 140 ms (fault steady state)
%   step_pu    largest change of a phase phasor from pre-fault to the first window
%   dev0_pu    max |g - g_steady| in the first window (0 to 20 ms)
%   ts_ms      end of the last window with max |g - g_steady| above 2 %
%              (0 if the first window is already within 2 %; cannot be below 20)
%   df_Hz      frequency deviation during the fault, from the slope of arg(V1)
%              between the windows that start at 20 and 120 ms
% The 5 % level is above the switching ripple of the converter models at the Pex bus.
%
% Usage:
%   r2_transient_damping('<folder with the EMT cases>', 'transient_damping.csv', '<folder for the figure>')
% No toolbox is needed.

if nargin < 3, figdir = ''; end
files = [dir(fullfile(casedir,'S2a_*.mat')); dir(fullfile(casedir,'S2d_*.mat')); ...
         dir(fullfile(casedir,'S2b_*.mat')); dir(fullfile(casedir,'S2c_*.mat'))];
rows = {};
for k = 1:numel(files)
    d = load(fullfile(files(k).folder, files(k).name));
    c = d.caso;
    ang = NaN; if isfield(c,'meta') && isfield(c.meta,'ang_inc_pedido'), ang = c.meta.ang_inc_pedido; end
    rg = NaN; if isfield(c,'Rtierra'), rg = c.Rtierra; end
    buses = {'Pex', d.GridConvs.Pex.Vabc; 'DG1', d.GridConvs.DG1.Vac};
    for b = 1:2
        m = one_bus(buses{b,2}, c.t_falta, c.red_f);
        rows(end+1,:) = {c.id, c.meta.conjunto, c.tipo, buses{b,1}, c.Rf, rg, ang, ...
            m.rpk, m.fdom, m.nhalf, m.t5, m.regime, m.step, m.dev0, m.ts, m.df}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'id','set','type','bus','Rf_ohm','Rg_ohm','ang_deg', ...
    'rpk_pu','fdom_Hz','nhalf','t5_ms','regime','step_pu','dev0_pu','ts_ms','df_Hz'});
T = sortrows(T, {'set','type','bus','Rf_ohm','Rg_ohm','ang_deg'});
writetable(T, outcsv);
fprintf('%d events from %d records -> %s\n', height(T), numel(files), outcsv);
if ~isempty(figdir), plot_examples(casedir, figdir); end
end

function m = one_bus(ts, tf, f0)
[tt, r, Vb, t, V] = split_transient(ts, tf, f0);
% ---- (1) fast transient
m.rpk = max(max(abs(r(tt <= 0.020,:))));
[~, ph] = max(sum(r(tt <= 0.040,:).^2, 1)); x = r(:,ph);
s = 0; n = 0;
for i = 1:numel(x)
    if x(i) > 0.05, sn = 1; elseif x(i) < -0.05, sn = -1; else, sn = 0; end
    if sn ~= 0 && sn ~= s, n = n + 1; s = sn; end
end
m.nhalf = n;
a = max(abs(r), [], 2); blk = 0.0025; nb = round(0.100/blk); m.t5 = 0;
for q = 1:nb
    pk = max(a(tt > (q-1)*blk & tt <= q*blk));
    if pk > 0.05, m.t5 = q*blk*1e3; end
end
m.fdom = NaN;
if m.rpk >= 0.05
    j = tt <= 0.040; y = x(j); y = y - mean(y); N = 8192; fs = 1/median(diff(tt));
    wn = 0.5*(1 - cos(2*pi*(0:numel(y)-1)'/(numel(y)-1)));   % Hann window, no toolbox needed
    Y = abs(fft(y .* wn, N)); f = (0:N-1)'*fs/N; kk = f > 100 & f < fs/2;
    [~, im] = max(Y(kk)); fk = f(kk); m.fdom = fk(im);
end
if n >= 4, m.regime = 'oscillatory'; elseif n >= 1, m.regime = 'damped'; else, m.regime = 'none'; end
% ---- (2) fundamental-frequency signature
w = 2*pi*f0;
P0 = phasor(t, V, tf-0.040, tf, tf, w, 0);            % pre-fault
gf = signature(phasor(t, V, tf+0.120, tf+0.140, tf, w, 0), Vb);
t0 = 0:0.0025:0.120; m.ts = 0; a1 = zeros(size(t0));
for q = 1:numel(t0)
    P = phasor(t, V, tf+t0(q), tf+t0(q)+0.020, tf, w, 0);
    dev = max(abs(signature(P, Vb) - gf));
    a1(q) = angle(P*[1; exp(2j*pi/3); exp(-2j*pi/3)]);
    if q == 1, m.dev0 = dev; m.step = max(abs(P - P0)) / Vb; end
    if dev > 0.02 && t0(q) <= 0.100, m.ts = (t0(q) + 0.020)*1e3; end
end
k = t0 >= 0.020; pf = polyfit(t0(k), unwrap(a1(k)), 1); m.df = pf(1)/(2*pi);
end

function g = signature(P, Vb)
% P: row of phase phasors [Va Vb Vc]; g = [|V1|/Vb, V2/V1, V0/V1]
a = exp(2j*pi/3);
V1 = P*[1; a; a^2]/3; V2 = P*[1; a^2; a]/3; V0 = sum(P)/3;
g = [abs(V1)/Vb, V2/V1, V0/V1];
end

function P = phasor(t, V, ta, tb, tf, w, nh)
% fundamental phasor of each column of V on (ta, tb], time origin tf; nh harmonics fitted too
i = t > ta + 1e-9 & t <= tb + 1e-9; x = t(i) - tf;
B = [ones(nnz(i),1), cos(w*x), sin(w*x)];
for h = 2:nh, B = [B, cos(h*w*x), sin(h*w*x)]; end %#ok<AGROW>
c = B \ V(i,:); P = c(2,:) - 1j*c(3,:);
end

function [tt, r, Vb, t, V] = split_transient(ts, tf, f0)
t = ts.Time(:); V = squeeze(ts.Data); if size(V,1) ~= numel(t), V = V.'; end
w = 2*pi*f0;
ipre = t >= tf-0.04 & t < tf-1e-9;
Bp = [ones(nnz(ipre),1), cos(w*(t(ipre)-tf)), sin(w*(t(ipre)-tf))];
cp = Bp \ V(ipre,:); Vb = mean(hypot(cp(2,:), cp(3,:)));
i = t > tf+1e-9 & t <= tf+0.140; tt = t(i) - tf;
B = slow_basis(tt, w);
r = (V(i,:) - B*lsqminnorm(B, V(i,:))) / Vb;
j = tt <= 0.100; tt = tt(j); r = r(j,:);
end

function B = slow_basis(tt, w)
tt = tt(:); s = tt/0.140; c = cos(w*tt); q = sin(w*tt);
B = [ones(size(tt)), s, c, q, s.*c, s.*q, s.^2.*c, s.^2.*q];
for tau = [0.010 0.030], e = exp(-tt/tau); B = [B, e, e.*c, e.*q]; end %#ok<AGROW>
for h = 2:7, B = [B, cos(h*w*tt), sin(h*w*tt)]; end %#ok<AGROW>
end

function plot_examples(casedir, figdir)
% bolted A-G; A-G with 1, 3 and 30 ohm; A-B with 3 ohm; A-B-G with 3 ohm
ids = {'S2a_AG_001','S2a_AG_004','S2a_AG_005','S2a_AG_007','S2a_AB_015','S2a_ABG_025'};
f = figure('Visible','off','Position',[0 0 1500 1300]);
for k = 1:numel(ids)
    fn = fullfile(casedir,[ids{k} '.mat']); if ~isfile(fn), continue; end
    d = load(fn); c = d.caso;
    [tt, r, Vb, t, V] = split_transient(d.GridConvs.Pex.Vabc, c.t_falta, c.red_f);
    i = t > c.t_falta-0.02 & t <= c.t_falta+0.10;
    subplot(numel(ids),2,2*k-1); plot((t(i)-c.t_falta)*1e3, V(i,:)/Vb); grid on; xlim([-20 100]);
    title(sprintf('%s  %s  R_f = %g \\Omega: v(t) at the fault bus (p.u.)', strrep(ids{k},'_','\_'), c.tipo, c.Rf));
    subplot(numel(ids),2,2*k); plot(tt*1e3, r); grid on; ylim([-0.4 0.4]); xlim([-20 100]);
    title('fast transient r(t) (p.u.)');
end
exportgraphics(f, fullfile(figdir,'transient_residual_examples.png'), 'Resolution', 150);
end
