# Second-round revision (IEEE TPWRD, TPWRD-00232-2026.R1)

Scripts that regenerate the numbers, tables and figures added or changed in the second revision of the manuscript. The final classifier is the GA classifier with the zero-sequence gate (G5, threshold 0.05) and persistence M = 10, fixed in `experiments/r2_final_opts.m`.

Comments and console messages in these scripts are in Spanish.

## Requirements

- MATLAB R2024b. The EMT campaign uses Simulink (Specialized Power Systems), and the model must be run in R2024b.
- No extra toolbox for the classifiers (`src/fault_classifiers/ga_fault_classifier_pure.m`). The optional GA-FuL toolbox is only used to re-run the original implementation (`compare_classifier_accuracy`, timing): set the environment variable `GAFUL_PATH` to its folder.
- Python 3 with pandas for `experiments/r2_fig_t94_datos.py`.

## Quick start

```matlab
run('revision_r2/setup_data.m')            % copies ../data/*.mat to revision_r2/data/simulation_curated
addpath('revision_r2/experiments')
r2_final_run_all('bloques', {'nominal', 'robustez', 'instrumentacion', 'degenerada', 'tiempos'})
```

Outputs are written to `revision_r2/_revision_TPWRD/ronda2/resultados/` (ignored by git). `r2_final_run_all` without arguments runs every block; `emt`, `lab`, `lab_fig`, `sintef` and `puerta` need non-public data (see below) and stop with an error without it.

## Where each result comes from

| Manuscript item | Script | Data |
|---|---|---|
| Table IV (accuracy, GA/Clarke/PE3D, four windows) and Fig. 8 (accuracy per fault type) | `r2_final_nominal.m`; `compare_classifier_accuracy.m` for Clarke and PE3D | public simulation records |
| Table V (robustness suite) | `r2_final_robustez.m`, `r2_final_instrumentacion.m` | public simulation records and synthetic events |
| Security of the degenerate branch and of the zero-sequence gate | `r2_degenerate_security.m`, `r2_degenerate_guard_test.m`, `r2_gate_sensitivity.m` | synthetic events (`r2_security_events_gen.m`) |
| Computation time per window | `r2_final_tiempos.m` | public simulation records |
| Fig. 5 (A-B trajectory, R and RL cases) | `r2_fig_rrl.m` | public simulation records |
| Line-to-line sector law and validity region, Fig. 9(a) curves | `r2_sector_analytic.m`, `r2_sector_seqsolve.m`, `r2_ground_rf_laws.m` | analytical and synthetic |
| EMT campaign (Table III, Fig. 9) | `r2_emt_campaign.m`, `r2_emt_run_case.m`, `r2_final_emt.m` | Simulink model, available from the authors on request |
| Laboratory results and Fig. 7 | `r2_final_lab.m`, `r2_final_lab_fig.m` | laboratory records, available from the authors on request |
| SINTEF data set | `r2_t20_run_all.m` (set `SINTEF_DATA`) | available from the authors on request |
| Cross-country illustration | `r2_t73_crosscountry_final.m` | synthetic |
| Figure data files of Fig. 8 and Fig. 9 | `r2_fig_t94_datos.py` | outputs of the scripts above |
| Data set inventory (Table III) | `r2_dataset_inventory.m` | public simulation records |

## Not included

- Laboratory records (about 2 GB), the 294 EMT cases (about 1.1 GB), the SINTEF data set and the Simulink model. They are available from the authors on reasonable request.
- The scripts that write LaTeX fragments for the manuscript. The blocks `macros` and `cifras` of `r2_final_run_all` depend on them and are not in its default list; the `tex` step of `r2_final_emt` prints a warning and continues.
