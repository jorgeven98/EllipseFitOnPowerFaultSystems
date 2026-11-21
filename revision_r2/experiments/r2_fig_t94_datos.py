"""Datos de las figuras nuevas o recuperadas de la segunda revision (T94).

uv run --with pandas python experiments/r2_fig_t94_datos.py

Escribe en paper/data/ (nombres nuevos; no sobrescribe CSV existentes del paper):
  r2_heat_q.csv            mapa de calor de precision por tipo, ventana de cuarto de ciclo,
                           desde los ids HEAT_* de resultados/T90/cifras_finales.csv
  r2_val_curves.csv        curvas de giro 10/20/30 grados y limite de deteccion (copia de T21/crossing_curves.csv)
  r2_val_region_wrong.csv  franja de par adyacente (copia de T21/fig_region_wrong.csv)
  r2_val_region_undet.csv  region sin deteccion (copia de T21/fig_region_undetected.csv)
  r2_val_dataset.csv       registros originales R (copia de T21/fig_markers_dataset.csv)
  r2_val_emt.csv           eventos L-L de la campana EMT en el bus de la falta (Pex), resistivos,
                           situados en (arg Z2, Rf/|Z2|); Z2 de cada escenario de red estimado con (14)
                           a partir de V2/V1 medido en sus eventos resistivos (mediana de Re e Im)
  r2_rf_out.csv            resultado por evento frente a Rf por fase (barrido S2a + S2d, Pex y DG1)
Resumen y comprobaciones en _revision_TPWRD/ronda2/resultados/T94/r2_fig_t94_datos_log.txt
"""
import os, shutil, io, contextlib
import numpy as np
import pandas as pd

RAIZ = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
RES = os.path.join(RAIZ, '_revision_TPWRD', 'ronda2', 'resultados')
OUT = os.path.join(RAIZ, 'paper', 'data')
LOG = os.path.join(RES, 'T94', 'r2_fig_t94_datos_log.txt')
os.makedirs(os.path.dirname(LOG), exist_ok=True)
buf = io.StringIO()


def p(*a):
    print(*a, file=buf)


# ---------- 1. mapa de calor (cuarto de ciclo) ----------
cf = pd.read_csv(os.path.join(RES, 'T90', 'cifras_finales.csv'))
cols = ['AB', 'CA', 'BC', 'AG', 'BG', 'CG', 'ABG', 'CAG', 'BCG']
rows = [('GA', 'DG1'), ('CLARKE', 'DG1'), ('PE3D', 'DG1'), ('GA', 'Pex'), ('CLARKE', 'Pex'), ('PE3D', 'Pex')]
H = []
for i, (m, c) in enumerate(rows):
    for j, k in enumerate(cols):
        v = float(cf.loc[cf.id == f'HEAT_{m}_{c}_{k}', 'valor'].iloc[0])
        H.append({'col': j, 'row': i, 'val': v, 'dark': int(v >= 70)})
H = pd.DataFrame(H)
H.to_csv(os.path.join(OUT, 'r2_heat_q.csv'), index=False)
for i, (m, c) in enumerate(rows):
    p(f'heat {m:6s} {c}: media {H[H.row == i].val.mean():.2f}')
# comprobacion contra la Tabla IV (columna 1/4 de ciclo)
nom = {('GA', 'DG1'): 'NOM_GA_DG1_W20', ('GA', 'Pex'): 'NOM_GA_Pex_W20'}
for k, idn in nom.items():
    q = cf.loc[cf.id == idn, 'valor']
    if len(q):
        p(f'  Tabla IV {idn} = {q.iloc[0]}')

# ---------- 2. mapa de validez ----------
T21 = os.path.join(RES, 'T21')
for src, dst in [('crossing_curves.csv', 'r2_val_curves.csv'), ('fig_region_wrong.csv', 'r2_val_region_wrong.csv'),
                 ('fig_region_undetected.csv', 'r2_val_region_undet.csv'), ('fig_markers_dataset.csv', 'r2_val_dataset.csv')]:
    shutil.copyfile(os.path.join(T21, src), os.path.join(OUT, dst))

s1 = pd.read_csv(os.path.join(RES, 'T90', 'emt', 'S1_final.csv'))
e = s1[(s1.terminal == 'Pex') & (s1.valida == 1)].copy()
p(f'S1 Pex validos: {len(e)}; resultado: {e.resultado4.value_counts().to_dict()}')
ind = e[e.Lf_H > 0]
p(f'excluidos por rama de falta inductiva (Lf > 0): {len(ind)}; resultado: {ind.resultado4.value_counts().to_dict()}')
e = e[e.Lf_H == 0].copy()
phi = {'B-C': 0.0, 'A-B': -120.0, 'C-A': 120.0}
q = e.V2V1_falta_pct.values / 100 * np.exp(1j * np.deg2rad(e.argV2V1_falta_deg.values - e.tipo.map(phi).values))
e['Z2ev'] = q * e.R_bucle_ohm.values / (1 - q)
e['cond'] = np.abs(1 - q)
key = ['red_Z_ohm', 'red_XR', 'red_f_Hz', 'Id_ref_pu', 'dg1_on', 'dg2_on', 'dg1_dP_W']
ok = e[e.cond > 0.05]
gz = ok.groupby(key).Z2ev.apply(lambda s: np.median(s.values.real) + 1j * np.median(s.values.imag)).rename('Z2')
e = e.join(gz, on=key)
assert e.Z2.notna().all()
e['argZ2_deg'] = np.rad2deg(np.angle(e.Z2.values.astype(complex)))
e['Rf_over_absZ2'] = e.R_bucle_ohm / np.abs(e.Z2.values.astype(complex))
a = e.argZ2_deg.values
ycr = np.where(a > 60, np.sin(np.pi / 3) / np.sin(np.deg2rad(np.clip(a - 60, 1e-9, None))), np.inf)
e['mas_alla_30'] = e.Rf_over_absZ2.values > ycr
p(f'eventos dibujados: {len(e)}; escenarios de red: {e.groupby(key).ngroups}')
p('resultado frente a la curva de 30 grados (columnas: mas alla de la curva):')
p(pd.crosstab(e.resultado4, e.mas_alla_30).to_string())
nz = e[e.R_bucle_ohm.isin([2.12])]
red = e[(e.red_Z_ohm.round(4) == 5.2847) & (e.dg1_on == 1) & (e.dg2_on == 1)]
p(f'arg Z2 de la red del paper (escenarios con los dos GFM): {red.argZ2_deg.min():.2f} a {red.argZ2_deg.max():.2f} deg')
p('sin declaracion:')
p(e[e.resultado4 == 'sin_declaracion'][['id', 'factor', 'R_bucle_ohm', 'argZ2_deg', 'Rf_over_absZ2', 'mas_alla_30']].to_string())
out = e[['argZ2_deg', 'Rf_over_absZ2', 'resultado4', 'id', 'factor', 'R_bucle_ohm']].copy()
out['res'] = out.resultado4.map({'correcta': 0, 'par_adyacente': 1, 'sin_declaracion': 2})
out.to_csv(os.path.join(OUT, 'r2_val_emt.csv'), index=False)
# Comprobacion no circular: Z2 de cada evento estimada solo con los OTROS eventos resistivos de su
# escenario que tienen otra Rf (deja uno fuera). Salida en resultados/T94/r2_val_emt_loo.csv
e = e.reset_index(drop=True)
e['g'] = e.groupby(key).ngroup()
lx, ly = [], []
for i, r in e.iterrows():
    o = e[(e.g == r.g) & (e.index != i) & (e.cond > 0.05) & (e.R_bucle_ohm != r.R_bucle_ohm)]
    if len(o) == 0:
        lx.append(np.nan); ly.append(np.nan); continue
    z = np.median(o.Z2ev.values.real) + 1j * np.median(o.Z2ev.values.imag)
    lx.append(np.rad2deg(np.angle(z))); ly.append(r.R_bucle_ohm / abs(z))
e['argZ2_loo_deg'] = lx; e['Rf_over_absZ2_loo'] = ly
a = e.argZ2_loo_deg.values
e['mas_alla_30_loo'] = np.where(np.isnan(a), np.nan,
    (e.Rf_over_absZ2_loo.values > np.where(a > 60, np.sin(np.pi / 3) / np.sin(np.deg2rad(np.clip(a - 60, 1e-9, None))), np.inf)).astype(float))
e[['id', 'tipo', 'factor', 'R_bucle_ohm', 'resultado4', 'argZ2_loo_deg', 'Rf_over_absZ2_loo', 'mas_alla_30_loo']].to_csv(
    os.path.join(RES, 'T94', 'r2_val_emt_loo.csv'), index=False)
lo = e.dropna(subset=['argZ2_loo_deg'])
p(f'deja uno fuera: {len(lo)} de {len(e)} eventos con Z2 de otros eventos de su escenario con otra Rf')
p(pd.crosstab(lo.resultado4, lo.mas_alla_30_loo).to_string())

# ---------- 3. resultado frente a Rf, con latencia ----------
s2 = pd.read_csv(os.path.join(RES, 'T90', 'emt', 'S2_final.csv'))
w = s2[s2.conjunto.isin(['S2a', 'S2d']) & (s2.valida == 1)].copy()
NOM = [('Pex', 'A-G'), ('Pex', 'A-B-G'), ('Pex', 'A-B'), ('DG1', 'A-G'), ('DG1', 'A-B-G'), ('DG1', 'A-B')]
w['fila'] = [NOM.index((t, k)) + 1 for t, k in zip(w.terminal, w.tipo)]
w['res'] = w.resultado4.map({'correcta': 0, 'par_adyacente': 1, 'erronea_otra': 1, 'sin_declaracion': 2})
w['ms'] = np.where(w.res == 0, w.correct_decl_ms, np.where(w.res == 1, w.decl_ms, np.nan))
w['serie'] = [{'A-G': 0, 'A-B-G': 1, 'A-B': 2}[k] + (0 if t == 'Pex' else 3) for t, k in zip(w.terminal, w.tipo)]
r = w[['Rf_ohm', 'fila', 'res', 'ms', 'serie', 'terminal', 'tipo', 'decl_class']].sort_values(['fila', 'Rf_ohm'])
r.to_csv(os.path.join(OUT, 'r2_rf_out.csv'), index=False)
p(f'S2a + S2d: {len(w)} eventos; por resultado {w.resultado4.value_counts().to_dict()}')
p('clases vecinas declaradas: ' + str(w[w.res == 1].groupby(['tipo', 'decl_class']).size().to_dict()))
p(f'latencia de las declaraciones (correctas y vecinas): {np.nanmin(r.ms):.2f} a {np.nanmax(r.ms):.2f} ms')
for (t, k), g in w.groupby(['terminal', 'tipo']):
    okk = g[g.res == 0]
    p(f'  {t} {k}: correcta hasta {okk.Rf_ohm.max() if len(okk) else None} ohm; vecina en {sorted(g[g.res == 1].Rf_ohm)}; '
      f'sin declaracion desde {g[g.res == 2].Rf_ohm.min() if (g.res == 2).any() else None}; '
      f'lat. correctas {okk.correct_decl_ms.min():.2f} a {okk.correct_decl_ms.max():.2f} ms' if len(okk) else '')
# comprobacion contra la figura de T23/T90 (fig_rf_*.csv)
for res_, f in [(0, 'fig_rf_ok.csv'), (1, 'fig_rf_wrong.csv'), (2, 'fig_rf_none.csv')]:
    ref = pd.read_csv(os.path.join(RES, 'T90', 'emt', f)).sort_values(['fila', 'Rf_ohm']).reset_index(drop=True)
    mine = r[r.res == res_][['Rf_ohm', 'fila']].sort_values(['fila', 'Rf_ohm']).reset_index(drop=True)
    p(f'  {f}: {len(ref)} filas; coincide con r2_rf_out: {ref.equals(mine.astype(ref.dtypes.to_dict()))}')

with open(LOG, 'w', encoding='utf-8') as fh:
    fh.write(buf.getvalue())
print(buf.getvalue())
