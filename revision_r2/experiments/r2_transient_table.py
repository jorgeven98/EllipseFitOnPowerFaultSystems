"""Builds the body of the inception-transient table (Table S13 of the supplementary document) from
the CSV written by r2_transient_damping.m.
Usage: python3 r2_transient_table.py [transient_damping.csv] [transient_table.tex]
Fault bus (Pex), resistance sweep at the base inception instant (sets S2a and S2d).
Per fault type: n, peak and t5 of the fast transient, and ts of the fundamental-frequency signature.
A cell reads n/a while the record of that case has not been processed."""
import csv
import sys
src = sys.argv[1] if len(sys.argv) > 1 else 'results/transient_damping.csv'
dst = sys.argv[2] if len(sys.argv) > 2 else 'results/transient_table.tex'
rows = [r for r in csv.DictReader(open(src, encoding='utf-8')) if r['set'] in ('S2a', 'S2d') and r['bus'] == 'Pex']
types = ['A-G', 'A-B-G', 'A-B']
rfs = sorted({float(r['Rf_ohm']) for r in rows})
cell = {(r['type'], float(r['Rf_ohm'])): r for r in rows}
def fmt(x):
    return '%g' % x
def t5(x):   # observed for 100 ms
    return r'$\geq$100' if x >= 100 else fmt(x)
def ts(x):   # 0 means that the first one-cycle window is already within 2 % of the steady state
    return '--' if x == 0 else fmt(x)
out, missing = [], 0
for rf in rfs:
    line = [fmt(rf)]
    for t in types:
        r = cell.get((t, rf))
        if r is None:
            line += ['n/a', '', '', '']; missing += 1
        else:
            line += [r['nhalf'], '%.2f' % float(r['rpk_pu']), t5(float(r['t5_ms'])), ts(float(r['ts_ms']))]
    out.append(' & '.join(line) + r' \\')
open(dst, 'w', encoding='utf-8').write('\n'.join(out) + '\n\\bottomrule\n')
print('%d rows, %d cells without record -> %s' % (len(rfs), missing, dst))
