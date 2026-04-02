!>
!! @file
!! @brief Contains constant values used throughout the code(s).

!> @brief Compile-time constant parameters: default values, tolerances, and physical constants
module m_constants

    use m_precision_select

    character, parameter :: dflt_char = ' '             !< Default string value
    real(wp), parameter  :: dflt_real = -1.e6_wp        !< Default real value
    real(wp), parameter  :: sgm_eps = 1.e-16_wp         !< Segmentation tolerance
    real(wp), parameter  :: small_alf = 1.e-11_wp       !< Small alf tolerance
    real(wp), parameter  :: pi = 3.141592653589793_wp   !< Pi
    real(wp), parameter  :: verysmall = 1.e-12_wp       !< Very small number
    !> Radius cutoff to avoid division by zero for 3D spherical harmonic patch (geometry 14)
    real(wp), parameter :: small_radius = 1.e-32_wp
    integer, parameter  :: num_stcls_min = 5  !< Minimum # of stencils
    integer, parameter  :: path_len = 400  !< Maximum path length
    integer, parameter  :: name_len = 50  !< Maximum name length
    integer, parameter  :: dflt_int = -100  !< Default integer value
    integer, parameter  :: num_fluids_max = 10  !< Maximum number of fluids in the simulation
    integer, parameter  :: num_patches_max = 1000  !< Maximum number of IC patches
    integer, parameter  :: num_bc_patches_max = 10  !< Maximum number of boundary condition patches
    integer, parameter  :: max_sph_harm_degree = 5  !< Max degree L for 3D spherical harmonic patch (geometry 14)
    integer, parameter  :: pathlen_max = 400  !< Maximum path length
    integer, parameter  :: dflt_num_igr_iters = 2  !< number of iterations for IGR elliptic solve
    integer, parameter  :: dflt_num_igr_warm_start_iters = 50  !< default number of iterations for IGR elliptic solve
    real(wp), parameter :: dflt_alf_factor = 10._wp  !< scaling factor for IGR alpha
    real(wp), parameter :: dflt_vcfl_dt = 100._wp  !< value of vcfl_dt when viscosity is off for computing adaptive timestep size
    ! Reconstruction Types
    integer, parameter :: WENO_TYPE = 1   !< Using WENO for reconstruction type
    integer, parameter :: MUSCL_TYPE = 2  !< Using MUSCL for reconstruction type
    ! System constants
    integer, parameter :: CASE_FILE_ERROR_CODE = 22  !< Exit code for case file validation errors

    ! Boundary condition enumeration
    integer, parameter :: BC_PERIODIC = -1
    integer, parameter :: BC_REFLECTIVE = -2
    integer, parameter :: BC_GHOST_EXTRAP = -3
    integer, parameter :: BC_RIEMANN_EXTRAP = -4
    integer, parameter :: BC_CHAR_SLIP_WALL = -5
    integer, parameter :: BC_CHAR_NR_SUB_BUFFER = -6
    integer, parameter :: BC_CHAR_NR_SUB_INFLOW = -7
    integer, parameter :: BC_CHAR_NR_SUB_OUTFLOW = -8
    integer, parameter :: BC_CHAR_FF_SUB_OUTFLOW = -9
    integer, parameter :: BC_CHAR_CP_SUB_OUTFLOW = -10
    integer, parameter :: BC_CHAR_SUP_INFLOW = -11
    integer, parameter :: BC_CHAR_SUP_OUTFLOW = -12
    integer, parameter :: BC_NULL = -13
    integer, parameter :: BC_AXIS = -14
    integer, parameter :: BC_SLIP_WALL = -15
    integer, parameter :: BC_NO_SLIP_WALL = -16
    integer, parameter :: BC_DIRICHLET = -17
end module m_constants
