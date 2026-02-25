# Ellipse Fitting for Fault Classification in Power Systems with Inverter-Based Resources

MATLAB implementation of the methods presented in:

> J. Ventura, J. Hrdina, A. Návrat, M. Stodola, A. Eid, S. Sanchez-Acevedo, F. G. Montoya,
> "Understanding the Geometry of Faulted Power Systems under High Penetration of Inverter-Based Resources via Ellipse Fitting and Geometric Algebra,"
> *arXiv:2509.10044 [eess.SY]*, 2025.
> https://arxiv.org/abs/2509.10044

## Overview

This repository contains three voltage-based fault classification methods for three-phase power systems with high penetration of inverter-based resources (IBR). The approach addresses the failure of traditional distance protection under asymmetric conditions by exploiting the geometry of the voltage polarization ellipse.

Each period of a three-phase voltage waveform traces an ellipse in 3D space. The shape, size, and orientation of this ellipse encode information about the type of fault, enabling detection within a **quarter-cycle** without relying on Clarke transform assumptions.

**Methods implemented:**

1. **GA (Geometric Algebra)** — Extracts the instantaneous voltage plane via bivectors, projects onto that plane, and fits an ellipse using Geometric Algebra for Conics (GAC). Classifies faults from roundness, magnitude, and orientation angle.
2. **Clarke** — Transforms to α-β-0 coordinates and applies PCA to extract ellipse parameters for classification.
3. **PE3D** — Characterizes the 3D polarization ellipse via Stokes parameters (following Alam et al., 2020).

All three methods share coherent thresholds for direct comparison.

## Supported Fault Types

| Code | Description |
|------|-------------|
| Normal | Balanced three-phase system |
| ABCG | Three-phase symmetric fault |
| AG, BG, CG | Single-phase to ground |
| AB, BC, CA | Two-phase line-to-line |
| ABG, BCG, CAG | Two-phase with ground |

## Requirements

- MATLAB R2020b or higher
- No additional toolboxes required

## Repository Structure

```
classifiers/       - Three fault classification methods (GA, Clarke, PE3D)
core/              - Core geometric algebra and ellipse fitting functions
processing/        - Sliding window voltage processing pipeline
signal_generation/ - Synthetic three-phase fault waveform generator
experiments/       - Main experiment scripts for method comparison
utils/             - Shared utility functions (metrics, plotting, export)
tests/             - Validation and threshold coherence verification scripts
gui/               - Interactive MATLAB GUI for real-time analysis
data/              - Simulation data (.mat files) for all 11 fault types
```

### Classifiers (`classifiers/`)

| File | Method | Description |
|------|--------|-------------|
| `ga_fault_classifier.m` | GA | Geometric Algebra with bivectors and GAC ellipse fitting |
| `clarke_fault_classifier.m` | Clarke | Clarke transform + PCA ellipse extraction |
| `pe3d_fault_classifier.m` | PE3D | 3D Polarization Ellipse via Stokes parameters |
| `fit_ellipse_gac.m` | — | GAC ellipse fitting wrapper |
| `fit_ellipse_gac_improved.m` | — | Enhanced fitting with degenerate case handling |

### Core Functions (`core/`)

| File | Description |
|------|-------------|
| `fitEllipseGAC.m` | Fits a 2D ellipse using Geometric Algebra for Conics |
| `extractEllipseParameters.m` | Extracts center, semi-axes, and rotation angle |
| `checkIfEllipse.m` | Validates that the fitted conic is an ellipse |
| `computeNormalizedBivector.m` | Computes bivector of the instantaneous voltage plane |
| `applyRotorToVoltage.m` | Applies GA rotor to align bivector to XY plane |
| `addNoiseToThreePhase.m` | Adds synthetic noise to three-phase signals |

## Data (`data/`)

Simulink simulation files for all supported fault types under a grid-connected system with distributed generation:

| File pattern | Description |
|---|---|
| `DataGrid_*.mat` | 11 basic fault types (AG, BG, CG, AB, BC, CA, ABG, BCG, CAG, ABCG, ABC) |
| `DataGridRL_*.mat` | Faults with resistive-inductive load impedance |
| `DataGridRF*.mat` | Faults with resistive fault impedance |

Each file contains three-phase voltage waveforms sampled at **4000 Hz**.

> **Note:** Laboratory hardware measurement files (~2 GB) are not included in this repository. Contact the authors for access.

## Quick Start

### Run the main comparison experiment

```matlab
addpath(genpath('.'))
run_sliding_window_detection
```

This script runs all three classifiers (GA, Clarke, PE3D) with four window sizes (1/4, 1/2, 3/4, 1 cycle) on all simulation data files and generates comparison metrics.

### Classify a single fault window

```matlab
addpath(genpath('.'))

% Load a data file
load('data/DataGrid_AB.mat')
% Extract V_abc (Nx3 matrix) from the loaded data
V_abc = extract_converter_data(data, 'voltage');

% Classify using the GA method
result = ga_fault_classifier(V_abc, 4000, 50);
disp(result.fault_type)
```

### Generate and classify a synthetic fault

```matlab
addpath(genpath('.'))

% Generate a synthetic A-B fault starting at t=3s, duration 0.16s
[va, vb, vc, t] = generateFault('A-B', 3.0, 0.16, 4000, 0.5);
V_abc = [va, vb, vc];

result = ga_fault_classifier(V_abc, 4000, 50);
disp(result.fault_type)   % Expected: 'AB'
```

### Launch the interactive GUI

```matlab
addpath(genpath('.'))
sliding_window_analyzer_gui
```

## Classifier Interface

All three classifiers share the same function signature:

```matlab
result = ga_fault_classifier(V_abc, fs, f)
result = clarke_fault_classifier(V_abc, fs, f)
result = pe3d_fault_classifier(V_abc, fs, f)
```

**Inputs:**
- `V_abc` — `[N×3]` matrix of three-phase voltages (columns: Va, Vb, Vc)
- `fs` — Sampling frequency in Hz (e.g., `4000`)
- `f` — System frequency in Hz (e.g., `50`)

**Output:** Structure with fields:
- `result.fault_type` — Detected fault type string (`'Normal'`, `'AB'`, `'AG'`, etc.)
- `result.params` — Diagnostic parameters (roundness, magnitude, angle, etc.)

## System Parameters

| Parameter | Value |
|-----------|-------|
| Sampling frequency | 4000 Hz |
| System frequency | 50 Hz |
| Window sizes tested | 1/4, 1/2, 3/4, 1 cycle |
| Window overlap | 75% |

## Classification Thresholds (coherent across all methods)

| Parameter | Threshold |
|-----------|-----------|
| Circularity — Normal/ABCG boundary | 0.99 |
| Magnitude — fault detection | 0.90 pu |
| Degeneracy — line-to-line detection | 0.05 |
| L-L angle margin | ±30° |

## Validation

Run the threshold coherence check:

```matlab
addpath(genpath('.'))
verify_threshold_coherence
```

Run the GA classifier unit tests:

```matlab
addpath(genpath('.'))
test_ga_classifier
```

## Citation

If you use this code, please cite:

```bibtex
@misc{ventura2025ellipse,
  title={Understanding the Geometry of Faulted Power Systems under High Penetration
         of Inverter-Based Resources via Ellipse Fitting and Geometric Algebra},
  author={Ventura, Jorge and Hrdina, Jaroslav and N{\'a}vrat, Ale{\v{s}} and
          Stodola, Marek and Eid, Ahmad and Sanchez-Acevedo, Santiago and
          Montoya, Francisco G.},
  year={2025},
  eprint={2509.10044},
  archivePrefix={arXiv},
  primaryClass={eess.SY},
  url={https://arxiv.org/abs/2509.10044}
}
```

## License

This code is provided for research and reproducibility purposes. Please cite the paper above if you use it in your work.
