# Analyze
Analyze replica exchange simulations.

## Python scripts
- `0_cinnabar_plots_10ns_replicate_1.ipynb` - Calculated vs experiment and Forward vs reverse $\Delta\Delta G$ plots for barnase:barstar (10 ns/replica AREX). Figure 3B, 6A-B
- `0_cinnabar_plots_50ns.ipynb` - Calculated vs experiment and Forward vs reverse $\Delta\Delta G$ plots for barnase:barstar (50 ns/replica AREX). Figure 6C-D
- `0_cinnabar_plots_rest_50ns.ipynb` - Calculated vs experiment and Forward vs reverse $\Delta\Delta G$ plots for barnase:barstar (50 ns/replica AREST). Figure 6E-F
- `0_cinnabar_plots_terminally_blocked.ipynb` - Forward vs reverse $\Delta\Delta G$ plots for terminally-blocked amino acids (5 ns/replica AREX). Figure 3A
- `0_protonation_state_correction.ipynb` - Protonation state $\Delta G_{phase}$ s used for computing $\Delta\Delta G_{binding}$ s for ASH35A and LYN27A. Supplemental Information
- `1_free_energy_discrepancy_timeseries_plots.ipynb` - Figure 6A, 6D, Supplementary Figure 11
- `1_free_energy_rmse_mue_timeseries_plots.ipynb` - Figure 6C,F
- `1_free_energy_timeseries_plots.ipynb` - $\Delta G_{phase}$ time series for barnase:barstar and terminally-blocked amino acids. Figure 3C-D, 3E-F, 4A, 4D
- `1_free_energy_timeseries_summary_plot.ipynb` - Summary plots of $\Delta G_{phase}$ time series slopes for barnase:barstar and terminally-blocked amino acids. Figure 3G-H, Supplementary Figure 3, 9
- `2_phi_psi_angle_timeseries_plots.ipynb` - Supplementary Figure 4
- `3_replica_mixing_matrix_plots.ipynb` - Supplementary Figure 2, 5, 10
- `3_replica_mixing_summary_plots.ipynb` - Supplementary Figure 2, 5, 10
- `4_rest_parameter_combo_plots.ipynb` - Supplementary Figure 8
- `5_alchemical_rest_protocol_plot.ipynb` - Supplementary Figure 1
- `6_outlier_residue_distance_plot.ipynb` - Supplementary Figure 7
- `7_ddg_per_mutation_10ns_arex.ipynb` - Supplementary Figure 6A
- `7_ddg_per_mutation_50ns_arest.ipynb` - Supplementary Figure 6C
- `7_ddg_per_mutation_50ns_arex.ipynb` - Supplementary Figure 6B
- `8_compute_cis_for_pccs.ipynb` -  Figure 5
- `8_correlate_du_dlambda_with_features_per_replica_50ns.ipynb` - Figure 5
- `8_generate_heatmap_per_replica_final_50ns.ipynb` - Figure 5
- `analysis_tools.py` - 
- `analyze_dg.py` - 
- `analyze_dg_timeseries.py` - 
- `generate_heatmap_data_50ns.py` - 
- `get_residue_dihedrals_per_replica.py` - 
- `get_residue_distances_per_replica.py` - 
- `get_water_counts_per_replica.py` - 
- `run_du_dlambda_analysis_per_replica.py` - 
- `run_make_traj_per_replica.py` - 

## Bash scripts
Bash scripts for analyzing experiment are located in `perses-barnase-barstar-paper/data/`.
