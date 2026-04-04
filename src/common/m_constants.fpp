
module m_constants

    use m_precision_select

    real(wp), parameter  :: dflt_real = -1.e6_wp        !< Default real value
    real(wp), parameter  :: sgm_eps = 1.e-16_wp         !< Segmentation tolerance
    real(wp), parameter  :: pi = 3.141592653589793_wp   !< Pi
    real(wp), parameter :: small_radius = 1.e-32_wp     !< Radius cutoff for spherical harmonic patch
    integer, parameter  :: num_stcls_min = 5  !< Minimum # of stencils
    integer, parameter  :: path_len = 400  !< Maximum path length
    integer, parameter  :: name_len = 50  !< Maximum name length
    integer, parameter  :: dflt_int = -100  !< Default integer value
    integer, parameter  :: num_fluids_max = 10  !< Maximum number of fluids in the simulation
    integer, parameter  :: num_patches_max = 1000  !< Maximum number of IC patches
    integer, parameter  :: num_bc_patches_max = 10  !< Maximum number of boundary condition patches
    integer, parameter  :: max_sph_harm_degree = 5  !< Max degree L for 3D spherical harmonic patch
    integer, parameter  :: dflt_num_igr_iters = 2  !< number of iterations for IGR elliptic solve
    integer, parameter  :: dflt_num_igr_warm_start_iters = 50  !< default number of iterations for IGR elliptic solve
    real(wp), parameter :: dflt_alf_factor = 10._wp  !< scaling factor for IGR alpha
    integer, parameter :: WENO_TYPE = 1   !< Using WENO for reconstruction type
    integer, parameter :: MUSCL_TYPE = 2  !< Using MUSCL for reconstruction type
    integer, parameter :: CASE_FILE_ERROR_CODE = 22  !< Exit code for case file validation errors

    ! Boundary condition enumeration
    integer, parameter :: BC_PERIODIC = -1
    integer, parameter :: BC_REFLECTIVE = -2
    integer, parameter :: BC_GHOST_EXTRAP = -3
    integer, parameter :: BC_CHAR_SUP_OUTFLOW = -12
    integer, parameter :: BC_AXIS = -14
    integer, parameter :: BC_SLIP_WALL = -15
    integer, parameter :: BC_NO_SLIP_WALL = -16
    integer, parameter :: BC_DIRICHLET = -17
end module m_constants
