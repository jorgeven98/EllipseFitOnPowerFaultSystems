function T = r2_lab_phase_order(labdir, toncsv, outcsv)
% R2_LAB_PHASE_ORDER  Order of the phases in the voltage channels of the laboratory
% records, obtained without the classifier.
%
% Two measured facts fix the order of the three channels of the grid-side measurement
% box (rows 183 to 185 of the records, phase voltages referred to a floating star):
%  (1) phase sequence: positive- and negative-sequence voltage of the channels in the
%      recorded order, over one cycle before the fault;
%  (2) faulted pair: the pair of channels whose line-to-line voltage collapses during
%      the fault (one cycle starting 40 ms after inception, over the cycle before it).
% The pair of (2) is compared with the pair in the name of the record (ab, bc, ac).
% The inception instants are read from the per-record results of the laboratory
% (column t_on_s of lab_por_registro.csv, written by r2_final_lab.m).
%
%   T = r2_lab_phase_order('<folder with Fault_emulator_Outkontaktor and
%        Fault_emulator_Autokontaktor>', 'lab_por_registro.csv', 'lab_phase_order.csv')
L = readtable(toncsv, 'TextType', 'string', 'Delimiter', ',');
L = L(L.senal == "box2" & L.puerta == "adaptive" & L.tipo ~= "abc" & ~isnan(L.t_on_s), :);
a = exp(2j*pi/3); rows = {};
for k = 1:height(L)
    fn = fullfile(labdir, "Fault_emulator_" + L.conjunto(k), L.registro(k) + ".mat");
    d = load(fn); R = d.results; t = R(1,:)' - R(1,1); V = R(183:185,:)';
    assert(contains(string(d.legendText{183}), 'From_meas_box2') && contains(string(d.legendText{183}), 'V_abc_meas_SI'));
    w = 2*pi*50; ton = L.t_on_s(k);
    P0 = phasor(t, V, ton - 0.025, ton - 0.005, w); P1 = phasor(t, V, ton + 0.040, ton + 0.060, w);
    pos = abs(P0*[1; a; a^2])/3; neg = abs(P0*[1; a^2; a])/3;                 % sequence in the recorded order
    pr = [1 2; 2 3; 3 1]; ll0 = abs(P0(pr(:,1)) - P0(pr(:,2))); ll1 = abs(P1(pr(:,1)) - P1(pr(:,2)));
    ratio = ll1./ll0; [rs, is] = sort(ratio); q = pr(is(1),:);
    rows(end+1,:) = {L.registro(k), L.conjunto(k), L.tipo(k), pos, neg, ratio(1), ratio(2), ratio(3), sprintf('%d-%d', min(q), max(q)), rs(1), rs(2)}; %#ok<AGROW>
end
T = cell2table(rows, 'VariableNames', {'record','set','name_pair','V_pos_recorded_order_V','V_neg_recorded_order_V', ...
    'LL12_fault_over_pre','LL23_fault_over_pre','LL31_fault_over_pre','collapsed_channels','collapsed_ratio','next_ratio'});
writetable(T, outcsv);
fprintf('%d line-to-line records\n', height(T));
fprintf('sequence in the recorded order: positive %.1f to %.1f V, negative %.1f to %.1f V (peak)\n', min(T.V_pos_recorded_order_V), max(T.V_pos_recorded_order_V), min(T.V_neg_recorded_order_V), max(T.V_neg_recorded_order_V));
fprintf('collapsed line-to-line voltage: %.3f to %.3f of pre-fault; the other two: %.3f at least\n', min(T.collapsed_ratio), max(T.collapsed_ratio), min(T.next_ratio));
G = groupsummary(T, {'set','name_pair','collapsed_channels'}); disp(G(:, {'set','name_pair','collapsed_channels','GroupCount'}));
end

function P = phasor(t, V, ta, tb, w)
i = t > ta & t <= tb; x = t(i); B = [ones(nnz(i),1), cos(w*x), sin(w*x)]; c = B \ V(i,:); P = c(2,:) - 1j*c(3,:);
end
