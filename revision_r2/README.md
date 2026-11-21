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

## Supplementary document

`supplementary.pdf` collects the complete tables of the test campaigns that the manuscript summarizes in Section IV (Tables S1 to S13 and Fig. S1).

| Item of the supplementary document | Script | Data |
|---|---|---|
| Table S1, factors of the line-to-line sector law | `r2_sector_scan.m`, `r2_sector_analytic.m`, `r2_sector_seqsolve.m`, `r2_sector_wave.m`, `r2_sector_event.m`; `r2_sector_paper_network.m` must be run before `r2_sector_analytic.m`, in the same folder | analytical and synthetic; public simulation records |
| Table S2, sector law in the EMT campaign | `r2_final_emt.m`, `r2_final_emt_resumen.m` | EMT cases |
| Table S3, phase-to-ground and floating-star voltages | `r2_t15_star_final_gate.m`; `r2_floating_star_check.m` and `r2_dy_transformer_check.m` for the derivation and the Dy transformer | synthetic |
| Tables S4 and S5, times of the decision chain | `r2_final_nominal.m`, `r2_final_lab.m`, `r2_final_emt.m` | public simulation records; laboratory records; EMT cases |
| Table S6, degenerate branch | `r2_degenerate_security.m`, `r2_degenerate_guard_test.m` | synthetic events (`r2_security_events_gen.m`) |
| Table S7, earlier and revised zero-sequence gate | campaign that selected the gate, not included (see below); `r2_final_instrumentacion.m` repeats part of its grid with the final settings | public simulation records |
| Table S8, adaptive roundness gate | `r2_gate_sensitivity.m` | public simulation records; laboratory records |
| Eq. (S1) and Table S9, cross-country faults | `r2_t73_crosscountry_final.m`; `r2_t73_xc.m` for the two EMT cases | synthetic; EMT cases |
| Table S10, data sets not used for the design | `r2_final_lab.m`, `r2_t20_run_all.m` | laboratory records; SINTEF data set |
| Tables S11 and S12, fault resistance and inception angle | `r2_final_emt.m`, `r2_final_emt_resumen.m` | EMT cases |
| Table S13 and Fig. S1, inception transient | `r2_transient_damping.m`, `r2_transient_table.py` | EMT cases (85 cases of the sweeps) |

`reference_results/S2_final.csv` holds the 170 events of the fault-resistance and inception-angle sweeps as reported in the paper (outcome, declared class and latency per event), and `reference_results/transient_damping.csv` the transient figures of the same events.

## Check of reproducibility

`reproduce_sweep_events.m` finds the files that the decision chain needs (eight files), copies them to an empty folder, resets the MATLAB path and recomputes the 170 events of the sweeps from the EMT cases:

```matlab
cd revision_r2
reproduce_sweep_events('<folder with the EMT cases>')
```

On MATLAB R2024a (macOS, Intel) the outcome, the declared class and the latency of the 170 events are identical to those of `reference_results/S2_final.csv`, which were obtained on MATLAB R2024b (Windows). The largest latency difference is 4e-13 ms.

## Not included

- Laboratory records (about 2 GB), the 294 EMT cases (about 1.1 GB), the SINTEF data set and the Simulink model. They are available from the authors on reasonable request.
- The scripts of the campaign that selected the zero-sequence gate, whose grid of 6528 events is the one of Table S7 of the supplementary document.
- The scripts that write LaTeX fragments for the manuscript. The blocks `macros` and `cifras` of `r2_final_run_all` depend on them and are not in its default list; the `tex` step of `r2_final_emt` prints a warning and continues.
