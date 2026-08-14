# Phase-delays shape multistability and basin sizes in Kuramoto networks: analytical estimates from network structure

This repository contains the MATLAB code used to generate the analyses and figures for the paper:

**Phase-delays shape multistability and basin sizes in Kuramoto networks: analytical estimates from network structure**

Kalel L. Rossi, Antonio Mihara, Lyle E. Muller, Rene O. Medrano-T, and Roberto C. Budzinski

## Repository structure

Every script adds a single folder to the path:

```
./functions   % ring_graph, heterogeneous_phase_lag, random_circulant_graph_sparse,
              % circulant_eigensystem, linear_stability, order_parameter,
              % simulate_KM, simulate_KM_final_state,
              % plot_phase_scatter, plot_abs_mu_stem, plot_order_parameter_time, imagesc_phase
```

Scripts are organized by figure and split into **simulation** scripts, which run the numerical simulations and save the resulting data to `files_scans/`, and **figure** scripts, which load that data and produce the plots.

The `files_scans/` data is already included in this repository, so the figures can be reproduced directly from the figure scripts without rerunning the simulations. Rerun the simulation scripts only if you want to regenerate the data from scratch.

### Figure 1 — stability and basin size of q-states

- **`figure_1_stability_analysis.m`** — self-contained. Reproduces panels **b–e**: linear stability of q-states with and without phase-delays, the corresponding spatiotemporal patterns, and the stability diagram as a function of coupling range `k`.
- **`simulations_basin_size.m`** — generates the basin-size data for a single coupling range `k`. The script expects `ii_k` to already be defined in the workspace (it uses `clearvars -except ii_k`), so it should be called from a loop, e.g.:
  ```matlab
  for ii_k = 1:25
      simulations_basin_size
  end
  ```
  Each run saves `files_scans/basin_size_data/N_51_k_<k>.mat`.
- **`figure_1_basin_size_plot.m`** — reproduces panels **f–i**: basin-size stem plots for representative `k` values and the mean basin-size estimation error as a function of `k`. Loads data from `files_scans/basin_size_data/`.
- **`figure_1_basin_geometry.m`** — self-contained. Reproduces the basin-geometry insets in panel **i** (`k = 15` and `k = 19`) by directly perturbing the two first oscillators and integrating the resulting initial conditions. This is the most computationally demanding script in the repository.

### Figure 2 — multi-q states

- **`figure_2_ratio_q_states.m`** — reproduces panel **a**: the fraction of simulations converging to single q-states as a function of `k`. Reuses the basin-size data from Figure 1 (`files_scans/basin_size_data/`).
- **`simulations_multi_q_states.m`** — generates the q- and multi-q-state data for the two example cases in the paper (`k = 4` and `k = 6`), saving `files_scans/multi_q_states_data/N_51_k_<k>.mat`. To keep the saved files small, only the per-realization order-parameter mean/std over the post-transient window are stored, not the full 1000 × 50001 time-series trace.
- **`figure_2_multi_q_states.m`** — reproduces panels **b–g**: order-parameter histograms, nearest-neighbour phase differences, eigenmode projections, and the spatiotemporal evolution and eigenmode contributions over time for the `k = 6` example. Loads data from `files_scans/multi_q_states_data/`.

### Figure 3 — chiral dynamics with random phase-lags

- **`simulations_random_phase_lags_chiral_dynamics.m`** — generates the stability, basin-size, and simulation data for the two random circulant phase-lag matrices, saving `files_scans/random_phase_lags_data/figure_3.mat`.
- **`figure_3_random_phase_lags_chiral_dynamics.m`** — reproduces panels **a–g**: the two phase-delay matrices, their stability spectra, basin-size predictions vs. numerical estimates, and the resulting spatiotemporal patterns. Loads `files_scans/random_phase_lags_data/figure_3.mat`.

## Requirements

MATLAB, with the Statistics and Machine Learning Toolbox (for `nanmean`, `nanstd`, `nansum`). No other toolboxes are required.

Runtimes vary substantially across scripts. For example, `simulations_basin_size.m` is run once per `k` and involves 20,000 simulations each time; `figure_1_basin_geometry.m` integrates a 1,000 × 1,000 grid of initial conditions. These script take a long time to run.

## Reproducing the figures

The `files_scans/` folder ships with the data used in the paper, so the figure scripts can be run directly. To regenerate the data instead:

1. Run the simulation script(s) for the figure of interest — they save to `files_scans/` automatically.
2. Run the corresponding figure script.
