!>
!! @file
!! @brief Contains module m_derived_types

#:include "macros.fpp"

!> @brief Shared derived types for field data, patch geometry, and MPI I/O structures
module m_derived_types

    use m_constants
    use m_precision_select

    implicit none

    !> Derived type adding the field position (fp) as an attribute
    type field_position
        real(stp), allocatable, dimension(:,:,:) :: fp  !< Field position
    end type field_position

    !> Derived type annexing a scalar field (SF)
    type scalar_field
        real(stp), pointer, dimension(:,:,:) :: sf => null()
    end type scalar_field

    !> Derived type for pressure/mass-velocity fields (bubbles)
    type pres_field
        real(stp), pointer, dimension(:,:,:,:,:) :: sf => null()
    end type pres_field

    !> Derived type annexing an integer scalar field (SF)
    type integer_field
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), pointer, dimension(:,:,:) :: sf => null()
#else
        integer, pointer, dimension(:,:,:) :: sf => null()
#endif
    end type integer_field

    type mpi_io_var
        integer, allocatable, dimension(:)            :: view
        type(scalar_field), allocatable, dimension(:) :: var
    end type mpi_io_var

    !> Derived type annexing a vector field (VF)
    type vector_field
        type(scalar_field), allocatable, dimension(:) :: vf  !< Vector field
    end type vector_field

    !> Generic 3-component vector (e.g., spatial coordinates or field components) Named _dt (derived types: x,y,z) to differentiate
    !! from t_vec3 (3-component vector)
    type vec3_dt  ! dt for derived types
        real(wp) :: x
        real(wp) :: y
        real(wp) :: z
    end type vec3_dt

    !> Integer bounds for variables
    type int_bounds_info
        integer                             :: beg
        integer                             :: end
        real(wp)                            :: vb1
        real(wp)                            :: vb2
        real(wp)                            :: vb3
        real(wp)                            :: ve1
        real(wp)                            :: ve2
        real(wp)                            :: ve3
        real(wp)                            :: pres_in, pres_out
        real(wp), dimension(3)              :: vel_in, vel_out
        real(wp), dimension(num_fluids_max) :: alpha_rho_in, alpha_in
        logical                             :: grcbc_in, grcbc_out, grcbc_vel_out
    end type int_bounds_info

    type bc_patch_parameters
        integer                :: geometry
        integer                :: type
        integer                :: dir
        integer                :: loc
        real(wp), dimension(3) :: centroid
        real(wp), dimension(3) :: length
        real(wp)               :: radius
    end type bc_patch_parameters

    !> Derived type adding beginning (beg) and end bounds info as attributes
    type bounds_info
        real(wp) :: beg
        real(wp) :: end
    end type bounds_info

    !> Derived type adding initial condition (ic) patch parameters as attributes NOTE: The requirements for the specification of the
    !! above parameters are strongly dependent on both the choice of the multicomponent flow model as well as the choice of the
    !! patch geometry.
    type ic_patch_parameters

        integer :: geometry  !< Type of geometry for the patch
        real(wp) :: x_centroid, y_centroid, z_centroid  !< Geometric center coordinates of the patch
        real(wp) :: length_x, length_y, length_z  !< Dimensions of the patch. x,y,z Lengths.
        real(wp) :: radius  !< Dimensions of the patch. radius.
        real(wp), dimension(3) :: radii  !< Elliptical/ellipsoidal patch radii in x, y, z
        real(wp) :: epsilon, beta  !< The isentropic vortex parameters for the amplitude of the disturbance and domain of influence.
        real(wp), dimension(2:9) :: a  !< Used by hardcoded IC and as temporary variables.
        logical :: non_axis_sym

        ! Geometry 13 (2D modal Fourier): fourier_cos(n), fourier_sin(n) for mode n
        real(wp), dimension(1:max_2d_fourier_modes) :: fourier_cos, fourier_sin
        logical :: modal_clip_r_to_min  !< When true, clip boundary radius: R(theta) = max(R(theta), modal_r_min) (Non-exp form only)
        real(wp) :: modal_r_min  !< Minimum boundary radius when modal_clip_r_to_min is true (Non-exp form only)
        logical :: modal_use_exp_form  !< When true, boundary = radius*exp(Fourier series)
        ! Geometry 14 (3D spherical harmonic): sph_har_coeff(l,m) for real Y_lm
        real(wp), dimension(0:max_sph_harm_degree,-max_sph_harm_degree:max_sph_harm_degree) :: sph_har_coeff
        real(wp), dimension(3) :: normal  !< Patch orientation normal vector (x, y, z)
        logical, dimension(0:num_patches_max - 1) :: alter_patch  !< Overwrite permissions for preceding patches
        logical :: smoothen  !< Whether patch boundaries are smoothed across cells
        integer :: smooth_patch_id  !< Identity (id) of the patch with which current patch is to get smoothed
        real(wp) :: smooth_coeff  !< Smoothing stencil size coefficient
        real(wp), dimension(num_fluids_max) :: alpha_rho
        real(wp) :: rho
        real(wp), dimension(3) :: vel
        real(wp) :: pres
        real(wp), dimension(num_fluids_max) :: alpha
        real(wp) :: gamma
        real(wp) :: pi_inf
        real(wp) :: cv
        real(wp) :: qv
        real(wp) :: qvp  !< Reference entropy per unit mass (SGEOS)
        integer :: hcid  !< Hardcoded initial condition ID
        real(wp) :: cf_val  !< Color function value
        real(wp) :: Y(1:1)  !< Species mass fractions

    end type ic_patch_parameters

    !> Derived type annexing the physical parameters (PP) of the fluids. These include the specific heat ratio function and liquid
    !! stiffness function.
    type physical_parameters
        real(wp)               :: gamma   !< Sp. heat ratio
        real(wp)               :: pi_inf  !< Liquid stiffness
        real(wp), dimension(2) :: Re      !< Reynolds number
        real(wp)               :: cv      !< heat capacity
        real(wp)               :: qv      !< reference energy per unit mass for SGEOS, q (see Le Metayer (2004))
        real(wp)               :: qvp     !< reference entropy per unit mass for SGEOS, q' (see Le Metayer (2004))
        real(wp)               :: G
    end type physical_parameters


    !> Species parameters
    type species_parameters
        character(LEN=name_len) :: name  !< Name of species
    end type species_parameters

    !> Max and min number of cells in a direction of each combination of x-,y-, and z-
    type cell_num_bounds
        integer :: mn_max, np_max, mp_max, mnp_max
        integer :: mn_min, np_min, mp_min, mnp_min
    end type cell_num_bounds

    type simplex_noise_params
        logical, dimension(3)                   :: perturb_vel
        real(wp), dimension(3)                  :: perturb_vel_freq
        real(wp), dimension(3)                  :: perturb_vel_scale
        real(wp), dimension(3, 3)               :: perturb_vel_offset
        logical, dimension(1:num_fluids_max)    :: perturb_dens
        real(wp), dimension(1:num_fluids_max)   :: perturb_dens_freq
        real(wp), dimension(1:num_fluids_max)   :: perturb_dens_scale
        real(wp), dimension(1:num_fluids_max,3) :: perturb_dens_offset
    end type simplex_noise_params
end module m_derived_types
