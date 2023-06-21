# Analyze
Analyze replica exchange simulations.

## Python scripts
- `0_cinnabar_plots_10ns_replicate_1.ipynb` - Calculated vs experiment and Forward vs reverse $\Delta\Delta G$ plots for barnase:barstar (10 ns/replica AREX for both apo and complex phases). Figure 3B, 6A-B
- `0_cinnabar_plots_50ns.ipynb` - Calculated vs experiment and Forward vs reverse $\Delta\Delta G$ plots for barnase:barstar (50 ns/replica AREX for complex phase, 10 ns/replica AREX for apo phase). Figure 6C-D
- `0_cinnabar_plots_rest_10ns.ipynb` - Calculated vs experiment and Forward vs reverse $\Delta\Delta G$ plots for barnase:barstar (10 ns/replica AREST for complex phase, 10 ns/replica AREX for apo phase). Supplementary Figure 10A-B
- `0_cinnabar_plots_rest_50ns.ipynb` - Calculated vs experiment and Forward vs reverse $\Delta\Delta G$ plots for barnase:barstar (50 ns/replica AREST, 10 ns/replica AREX for apo phase). Figure 6E-F
- `0_cinnabar_plots_terminally_blocked.ipynb` - Forward vs reverse $\Delta\Delta G$ plots for terminally-blocked amino acids (5 ns/replica AREX for both phases). Figure 3A
- `0_protonation_state_correction.ipynb` - Protonation state $\Delta G_{phase}$ s used for computing $\Delta\Delta G_{binding}$ s for ASH35A and LYN27A. Supplemental Information Section C, Supplementary Table 2
- `10_arex_arest_correlation_plot.ipynb` - AREST vs AREX $\Delta\Delta G$ comparison plot for barnase:barstar (50 ns/replica for complex phase, 10 ns/replica AREX for apo phase). Supplementary Figure 14
- `11_arex_arest_statistical_inefficiency_plot.ipynb` - AREST vs AREX statistical inefficiency comparison plot for barnase:barstar (50 ns/replica for complex phase). Supplementary Figure 13
- `1_free_energy_discrepancy_timeseries_plots.ipynb` - $\Delta\Delta G$ discrepancy time series plots for barnase:barstar (50 ns/replica AREX and AREST for complex phase, 10 ns/replica AREX for apo phase). Figure 7A, 7D, Supplementary Figure 15
- `1_free_energy_rmse_mue_timeseries_plots.ipynb` - RMSE and MUE time series plots for barnase:barstar (50 ns/replica AREX and AREST for complex phase, 10 ns/replica AREX for apo phase). Figure 7C,F
- `1_free_energy_timeseries_plots.ipynb` - $\Delta G_{phase}$ time series plots for barnase:barstar and terminally-blocked amino acids. Figure 3C-D, 3E-F, 4A, 4D
- `1_free_energy_timeseries_summary_plot_reordered.ipynb` - Summary plots of $\Delta G_{phase}$ time series slopes for barnase:barstar and terminally-blocked amino acids. Figure 3G-H, Supplementary Figure 3, 11
- `2_phi_psi_angle_timeseries_plots.ipynb` - Phi and psi time series plots for terminally-blocked amino acids (5 ns/replica AREX). Supplementary Figure 4
- `3_replica_mixing_matrix_plots.ipynb` - Replica mixing matrix plots for barnase:barstar and terminally-blocked amino acids. Supplementary Figure 2, 5, 12
- `3_replica_mixing_summary_plots.ipynb` - Summary plots for replica mixing for barnase:barstar and terminally-blocked amino acids. Supplementary Figure 2, 5, 12
- `4_rest_parameter_combo_plots.ipynb` - REST parameter combination comparison plots for barnase:barstar. Supplementary Figure 9
- `5_alchemical_rest_protocol_plot.ipynb` - Functions for alchemical and REST protocols. Supplementary Figure 1
- `6_outlier_residue_distance_plot.ipynb` - A29Y residue pair distance plots. Supplementary Figure 8
- `7_ddg_per_mutation_10ns_arest.ipynb` - $\Delta\Delta G$ per mutation for barnase:barstar (10 ns/replica AREST for complex phase, 10 ns/replica AREX for apo phase). Supplementary Figure 6B
- `7_ddg_per_mutation_10ns_arex.ipynb` - $\Delta\Delta G$ per mutation for barnase:barstar (10 ns/replica AREX for both apo and complex phases). Supplementary Figure 6A
- `7_ddg_per_mutation_50ns_arest.ipynb` - $\Delta\Delta G$ per mutation for barnase:barstar (50 ns/replica AREST for complex phase, 10 ns/replica AREX for apo phase).Supplementary Figure 6D
- `7_ddg_per_mutation_50ns_arex.ipynb` - $\Delta\Delta G$ per mutation for barnase:barstar (50 ns/replica AREX for complex phase, 10 ns/replica AREX for apo phase).Supplementary Figure 6C
- `8_compute_cis_for_pccs.ipynb` - 95% confidence intervals for the correlations in the $\partial U$ / $\partial \lambda$ heatmap (with 50 ns/replica complex phase simulation data). Figure 5
- `8_compute_cis_for_pccs_10ns.ipynb` - 95% confidence intervals for the correlations in the $\partial U$ / $\partial \lambda$ heatmap (with 10 ns/replica complex phase simulation data). Supplementary Figure 7
- `8_correlate_du_dlambda_with_features_per_replica_50ns.ipynb` - Time series plots for specific degrees of freedom and $\partial U$ / $\partial \lambda$ with 50 ns/replica complex phase simulation data. Figure 4B-C, 4E-F, 7B, 7E
- `8_generate_heatmap_per_replica_final_10ns.ipynb` -  Time series plots for specific degrees of freedom and $\partial U$ / $\partial \lambda$ with 10 ns/replica complex phase simulation data. Supplementary Figure 7
- `8_generate_heatmap_per_replica_final_50ns.ipynb` - $\partial U$ / $\partial \lambda$ heatmap. Figure 5
- `8_table_for_heatmap_10ns.ipynb` - Table for the names of the degrees of freedom corresponding to each PCC in the 10 ns/replica complex phase heatmap. Supplementary Figure 7
- `8_table_for_heatmap_50ns.ipynb` - Table for the names of the degrees of freedom corresponding to each PCC in the 50 ns/replica complex phase heatmap. Figure 5
- `analysis_tools.py` - Contains `DataAnalyzer` class for computing free energies with MBAR
- `analyze_dg.py` - Compute (single) free energies with MBAR
- `analyze_dg_timeseries.py` - Generate free energy time series with MBAR
- `generate_heatmap_data_50ns.py` - Generate pearson correlation coefficients for $\partial U$ / $\partial \lambda$ heatmap
- `get_residue_dihedrals_per_replica.py` - Generate dihedral angle time series for each interface residue for $\partial U$ / $\partial \lambda$ heatmap
- `get_residue_distances_per_replica.py` - Generate distance time series for all pairs of residues involving interface residues for $\partial U$ / $\partial \lambda$ heatmap
- `get_water_counts_per_replica.py` - Generate time series of number of waters in the neighborhood of the mutating residue for $\partial U$ / $\partial \lambda$ heatmap
- `run_du_dlambda_analysis_per_replica.py` - Generate du/dlambda time series for each replica for $\partial U$ / $\partial \lambda$ heatmap
- `run_make_traj_per_replica.py` - Generate DCD trajectories for each replica (to use for generating time series for dihedral angles, distances, and water counts for $\partial U$ / $\partial \lambda$ heatmap)

## Bash scripts
Bash scripts for analyzing experiment are located in `perses-barnase-barstar-paper/data/`.
