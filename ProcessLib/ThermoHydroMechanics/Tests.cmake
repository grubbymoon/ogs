# ThermoHydroMechanics; Small deformation, linear poroelastic, homogeneous
AddTest(
    NAME ThermoHydroMechanics_square_1e0
    PATH ThermoHydroMechanics/Linear/Square_sealed_homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu epsilon_xx epsilon_xx 1e-8 1e-8 
    expected_square_1e0_pcs_0_ts_10_t_1000.000000.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu epsilon_yy epsilon_yy 1e-8 1e-8
)

AddTest(
    NAME ThermoHydroMechanics_square_1e0_analyitcal
    PATH ThermoHydroMechanics/Linear/Square_sealed_homogeneous
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e0.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    analytical.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu pressure_ana pressure 1e-4 1e-4
    analytical.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu sigma_xx_ana sigma_xx 1e-2 1e-2
    analytical.vtu square_1e0_pcs_0_ts_10_t_1000.000000.vtu sigma_yy_ana sigma_yy 1e-2 1e-2
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, sealed, bimaterial
AddTest(
    NAME ThermoHydroMechanics_square_1e2_sealed_bimaterial
    PATH ThermoHydroMechanics/Linear/Beam_sealed_bimaterial
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu sigma_xx sigma_xx 1e-8 1e-8 
    expected_square_1e2_pcs_0_ts_10_t_100.000000.vtu square_1e2_pcs_0_ts_10_t_100.000000.vtu sigma_yy sigma_yy 1e-8 1e-8
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, unsealed, bimaterial
AddTest(
    NAME ThermoHydroMechanics_square_1e2_unsealed_bimaterial
    PATH ThermoHydroMechanics/Linear/Beam_unsealed_bimaterial
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu displacement displacement 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu pressure pressure 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu temperature temperature 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu sigma_xx sigma_xx 1e-8 1e-8
    expected_square_1e2_pcs_0_ts_10_t_1000.000000.vtu square_1e2_pcs_0_ts_10_t_1000.000000.vtu sigma_yy sigma_yy 1e-8 1e-8
)

# ThermoHydroMechanics; Small deformation, linear poroelastic, point heat source consolidation
AddTest(
    NAME ThermoHydroMechanics_square_1e2_point_heat_injection
    PATH ThermoHydroMechanics/Linear/Point_injection
    EXECUTABLE ogs
    EXECUTABLE_ARGS square_1e2.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu displacement displacement 1e-6 1e-6
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu pressure pressure 1e-6 1e-6
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu temperature temperature 1e-6 1e-6
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu sigma_xx sigma_xx 1e-6 1e-6
    expected_square_1e0_pcs_0_ts_10_t_10000.000000.vtu square_1e0_pcs_0_ts_10_t_10000.000000.vtu sigma_yy sigma_yy 1e-6 1e-6
)