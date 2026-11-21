function [lab, n] = r2_ground_rf_legacy(u0, u2)
% R2_GROUND_RF_LEGACY  Etiqueta que da la rama de tierra del clasificador
% publicado (classify_by_bivector) a partir de u0 = V0/V1 y u2 = V2/V1 (T23b).
% La normal del plano del lugar geometrico es, salvo un factor positivo,
%   n = N0 (1,1,1) + |T| (cos(d), cos(d - 120), cos(d + 120)),
%   N0 = 1 - |u2|^2,  T = |T| e^{j d} = 2 (u2* u0 - u0*)
% (n_a = Im(Vb* Vc), etc.). El clasificador ordena |n| y decide:
%   L-L-G con la fase de menor |n| como fase sana si dif1 < dif2 o dif2 > 0.3;
%   si no, L-G en la fase de mayor |n|.
N0 = 1 - abs(u2)^2; T = 2*(conj(u2)*u0 - conj(u0)); d = angle(T);
n = N0 + abs(T)*[cos(d), cos(d - 2*pi/3), cos(d + 2*pi/3)];   % (n_a, n_b, n_c)
n = n/norm(n);
biv = n([3 1 2]);                                              % [B12 B23 B31] = [n_c n_a n_b]
[bs, idx] = sort(abs(biv), 'descend');
dif1 = abs(bs(1) - bs(2)); dif2 = abs(bs(2) - bs(3));
if dif1 < dif2 || dif2 > 0.3
    L = {'A-B-G','B-C-G','C-A-G'}; lab = L{idx(3)};
else
    L = {'C-G','A-G','B-G'}; lab = L{idx(1)};
end
end
