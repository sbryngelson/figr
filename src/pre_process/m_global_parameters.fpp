!>
!! @file
!! @brief Contains module m_global_parameters

#:include 'case.fpp'

!> @brief Defines global parameters for the computational domain, simulation algorithm, and initial conditions
module m_global_parameters

#ifdef MFC_MPI
    use mpi  ! Message passing interface (MPI) module
#endif

    use m_derived_types  ! Definitions of the derived types
    use m_helper_basic  ! Functions to compare floating point numbers

    implicit none

    ! Logistics
    integer                 :: num_procs                 !< Number of processors
    character(LEN=path_len) :: case_dir                  !< Case folder location
    logical                 :: old_grid                  !< Use existing grid data
    logical                 :: old_ic, non_axis_sym      !< Use existing IC data
    integer                 :: t_step_old, t_step_start  !< Existing IC/grid folder
    logical                 :: cfl_adap_dt, cfl_const_dt, cfl_dt
    integer                 :: n_start, n_start_old

    ! Computational Domain Parameters

    integer :: proc_rank  !< Rank of the local processor Number of cells in the x-, y- and z-coordinate directions
    integer :: m
    integer :: n
    integer :: p

    !> @name Max and min number of cells in a direction of each combination of x-,y-, and z-
    type(cell_num_bounds) :: cells_bounds
    integer(kind=8)       :: nGlobal              !< Global number of cells in the domain
    integer               :: m_glb, n_glb, p_glb  !< Global number of cells in each direction
    integer               :: num_dims             !< Number of spatial dimensions
    integer               :: num_vels             !< Number of velocity components
    integer               :: grid_geometry        !< Grid geometry (always Cartesian in IGR-only build)
    !> Locations of cell-centers (cc) in x-, y- and z-directions, respectively
    real(wp), allocatable, dimension(:) :: x_cc, y_cc, z_cc
    !> Locations of cell-boundaries (cb) in x-, y- and z-directions, respectively
    real(wp), allocatable, dimension(:) :: x_cb, y_cb, z_cb
    real(wp) :: dx, dy, dz                             !< Minimum cell-widths in the x-, y- and z-coordinate directions
    type(bounds_info) :: x_domain, y_domain, z_domain  !< Locations of the domain bounds in the x-, y- and z-coordinate directions
    logical :: stretch_x, stretch_y, stretch_z         !< Grid stretching flags for the x-, y- and z-coordinate directions
    ! Grid stretching: a_x/a_y/a_z = rate, x_a/y_a/z_a = location
    real(wp) :: a_x, a_y, a_z
    integer  :: loops_x, loops_y, loops_z
    real(wp) :: x_a, y_a, z_a
    real(wp) :: x_b, y_b, z_b

    ! Simulation Algorithm Parameters
    integer :: model_eqns   !< Multicomponent flow model
    integer :: num_fluids   !< Number of different fluids present in the flow
    integer :: sys_size     !< Number of unknowns in the system of equations
    integer :: recon_type   !< Reconstruction Type
    integer :: weno_polyn   !< Degree of the WENO polynomials (polyn)
    integer :: muscl_polyn  !< Degree of the MUSCL polynomials (polyn)
    integer :: weno_order   !< Order of accuracy for the WENO reconstruction
    integer :: muscl_order  !< Order of accuracy for the MUSCL reconstruction
    integer :: igr_order    !< IGR reconstruction order
    ! Annotations of the structure, i.e. the organization, of the state vectors
    type(int_bounds_info) :: cont_idx              !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx               !< Indexes of first & last momentum eqns.
    integer               :: E_idx                 !< Index of total energy equation
    integer               :: alf_idx               !< Index of void fraction
    type(int_bounds_info) :: adv_idx               !< Indexes of first & last advection eqns.
    type(int_bounds_info) :: internalEnergies_idx  !< Indexes of first & last internal energy eqns.
    integer               :: gamma_idx             !< Index of specific heat ratio func. eqn.
    integer               :: pi_inf_idx            !< Index of liquid stiffness func. eqn.
    type(int_bounds_info) :: stress_idx            !< Indexes of elastic shear stress eqns.
    type(int_bounds_info) :: xi_idx                !< Indexes of first and last reference map eqns.
    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p). Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwint(1:3)

    ! Cell indices (InDices With BUFFer): includes buffer except in pre_process
    type(int_bounds_info)      :: idwbuff(1:3)
    type(int_bounds_info)      :: bc_x, bc_y, bc_z      !< Boundary conditions in the x-, y- and z-coordinate directions
    logical                    :: parallel_io           !< Format of the data files
    logical                    :: file_per_process      !< type of data output
    integer                    :: precision             !< Precision of output files
    logical                    :: down_sample           !< Down-sample the output data
    logical                    :: mixlayer_vel_profile  !< Set hyperbolic tangent streamwise velocity profile
    real(wp)                   :: mixlayer_vel_coef     !< Coefficient for the hyperbolic tangent streamwise velocity profile
    logical                    :: mixlayer_perturb      !< Superimpose instability waves to surrounding fluid flow
    integer                    :: mixlayer_perturb_nk   !< Number of Fourier modes for perturbation with mixlayer_perturb flag
    real(wp)                   :: mixlayer_perturb_k0   !< Peak wavenumber for mixlayer perturbation (default: most unstable mode)
    logical                    :: simplex_perturb
    type(simplex_noise_params) :: simplex_params
    real(wp)                   :: pi_fac                !< Factor for artificial pi_inf
    logical                    :: viscous

    ! Perturb density of surrounding air so as to break symmetry of grid
    logical                             :: perturb_flow
    integer                             :: perturb_flow_fluid  !< Fluid to be perturbed with perturb_flow flag
    real(wp)                            :: perturb_flow_mag    !< Magnitude of perturbation with perturb_flow flag
    logical                             :: perturb_sph
    integer                             :: perturb_sph_fluid   !< Fluid to be perturbed with perturb_sph flag
    real(wp), dimension(num_fluids_max) :: fluid_rho
    logical                             :: elliptic_smoothing
    integer                             :: elliptic_smoothing_iters
    integer, allocatable, dimension(:)  :: proc_coords         !< Processor coordinates in MPI_CART_COMM
    integer, allocatable, dimension(:)  :: start_idx           !< Starting cell-center index of local processor in global grid
#ifdef MFC_MPI
    type(mpi_io_var), public :: MPI_IO_DATA
    character(LEN=name_len)  :: mpiiofs
    integer                  :: mpi_info_int  !< MPI info for parallel IO with Lustre file systems
#endif

    ! Initial Condition Parameters
    integer                                                  :: num_patches     !< Number of patches composing initial condition
    type(ic_patch_parameters), dimension(num_patches_max)    :: patch_icpp      !< IC patch parameters (max: num_patches_max)
    integer                                                  :: num_bc_patches  !< Number of boundary condition patches
    logical                                                  :: bc_io           !< whether or not to save BC data
    type(bc_patch_parameters), dimension(num_bc_patches_max) :: patch_bc        !< Boundary condition patch parameters

    ! Fluids Physical Parameters
    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp  !< Stiffened gas EOS parameters and Reynolds numbers per fluid
    real(wp)                                             :: rhoref, pref  !< Reference parameters for Tait EOS
    !> @name Index variables used for m_variables_conversion
    !> @{
    integer :: momxb, momxe
    integer :: advxb, advxe
    integer :: contxb, contxe
    integer :: intxb, intxe
    !> @}

    integer, allocatable, dimension(:,:,:) :: logic_grid
    integer                                :: buff_size  !< Number of ghost cells for boundary condition storage
    integer, dimension(2)                  :: Re_size
    integer                                :: Re_size_max
    integer, allocatable, dimension(:,:)   :: Re_idx
    logical                                :: dummy      !< AMDFlang workaround for case-optimization + GPU-kernel bug

contains

    !> Assigns default values to user inputs prior to reading them in. This allows for an easier consistency check of these
    !! parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i  !< Generic loop operator
        ! Logistics

        case_dir = '.'
        old_grid = .false.
        old_ic = .false.
        t_step_old = dflt_int
        t_step_start = dflt_int

        cfl_adap_dt = .false.
        cfl_const_dt = .false.
        cfl_dt = .false.
        n_start = dflt_int

        ! Computational domain parameters
        m = dflt_int; n = 0; p = 0

        call s_update_cell_bounds(cells_bounds, m, n, p)

        x_domain%beg = dflt_real
        x_domain%end = dflt_real
        y_domain%beg = dflt_real
        y_domain%end = dflt_real
        z_domain%beg = dflt_real
        z_domain%end = dflt_real

        stretch_x = .false.
        stretch_y = .false.
        stretch_z = .false.

        a_x = dflt_real
        a_y = dflt_real
        a_z = dflt_real
        loops_x = 1
        loops_y = 1
        loops_z = 1
        x_a = dflt_real
        x_b = dflt_real
        y_a = dflt_real
        y_b = dflt_real
        z_a = dflt_real
        z_b = dflt_real

        ! Simulation algorithm parameters
        model_eqns = dflt_int
        num_fluids = dflt_int
        recon_type = WENO_TYPE
        weno_order = dflt_int
        igr_order = dflt_int
        muscl_order = dflt_int

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int

        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                bc_${DIM}$%vb${DIR}$ = 0._wp
                bc_${DIM}$%ve${DIR}$ = 0._wp
            #:endfor
        #:endfor

        parallel_io = .false.
        file_per_process = .false.
        precision = 2
        down_sample = .false.
        viscous = .false.
        mixlayer_vel_profile = .false.
        mixlayer_vel_coef = 1._wp
        mixlayer_perturb = .false.
        mixlayer_perturb_nk = 100
        mixlayer_perturb_k0 = 0.4446_wp
        perturb_flow = .false.
        perturb_flow_fluid = dflt_int
        perturb_flow_mag = dflt_real
        perturb_sph = .false.
        perturb_sph_fluid = dflt_int
        fluid_rho = dflt_real
        elliptic_smoothing_iters = dflt_int
        elliptic_smoothing = .false.

        dummy = .false.

        simplex_perturb = .false.
        simplex_params%perturb_vel(:) = .false.
        simplex_params%perturb_vel_freq(:) = dflt_real
        simplex_params%perturb_vel_scale(:) = dflt_real
        simplex_params%perturb_vel_offset(:,:) = dflt_real
        simplex_params%perturb_dens(:) = .false.
        simplex_params%perturb_dens_freq(:) = dflt_real
        simplex_params%perturb_dens_scale(:) = dflt_real
        simplex_params%perturb_dens_offset(:,:) = dflt_real

        ! Initial condition parameters
        num_patches = dflt_int

        do i = 1, num_patches_max
            patch_icpp(i)%geometry = dflt_int
            patch_icpp(i)%x_centroid = dflt_real
            patch_icpp(i)%y_centroid = dflt_real
            patch_icpp(i)%z_centroid = dflt_real
            patch_icpp(i)%length_x = dflt_real
            patch_icpp(i)%length_y = dflt_real
            patch_icpp(i)%length_z = dflt_real
            patch_icpp(i)%radius = dflt_real
            patch_icpp(i)%epsilon = dflt_real
            patch_icpp(i)%beta = dflt_real
            patch_icpp(i)%normal = dflt_real
            patch_icpp(i)%radii = dflt_real
            patch_icpp(i)%alter_patch = .false.
            patch_icpp(i)%alter_patch(0) = .true.
            patch_icpp(i)%smoothen = .false.
            patch_icpp(i)%smooth_patch_id = i
            patch_icpp(i)%smooth_coeff = dflt_real
            patch_icpp(i)%alpha_rho = dflt_real
            patch_icpp(i)%rho = dflt_real
            patch_icpp(i)%vel = dflt_real
            patch_icpp(i)%pres = dflt_real
            patch_icpp(i)%alpha = dflt_real
            patch_icpp(i)%gamma = dflt_real
            patch_icpp(i)%pi_inf = dflt_real
            patch_icpp(i)%cv = 0._wp
            patch_icpp(i)%qv = 0._wp
            patch_icpp(i)%qvp = 0._wp
            patch_icpp(i)%a(2) = dflt_real
            patch_icpp(i)%a(3) = dflt_real
            patch_icpp(i)%a(4) = dflt_real
            patch_icpp(i)%a(5) = dflt_real
            patch_icpp(i)%a(6) = dflt_real
            patch_icpp(i)%a(7) = dflt_real
            patch_icpp(i)%a(8) = dflt_real
            patch_icpp(i)%a(9) = dflt_real
            patch_icpp(i)%non_axis_sym = .false.
            patch_icpp(i)%fourier_cos(:) = 0._wp
            patch_icpp(i)%fourier_sin(:) = 0._wp
            patch_icpp(i)%modal_clip_r_to_min = .false.
            patch_icpp(i)%modal_r_min = 1.e-12_wp
            patch_icpp(i)%modal_use_exp_form = .false.
            patch_icpp(i)%sph_har_coeff(:,:) = 0._wp

            patch_icpp(i)%hcid = dflt_int
        end do

        num_bc_patches = 0
        bc_io = .false.

        do i = 1, num_bc_patches_max
            patch_bc(i)%geometry = dflt_int
            patch_bc(i)%type = dflt_int
            patch_bc(i)%dir = dflt_int
            patch_bc(i)%loc = dflt_int
            patch_bc(i)%centroid(:) = dflt_real
            patch_bc(i)%length(:) = dflt_real
            patch_bc(i)%radius = dflt_real
        end do

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        pi_fac = 1._wp

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
            fluid_pp(i)%G = 0._wp
        end do

    end subroutine s_assign_default_values_to_user_inputs

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_global_parameters_module

        integer :: i, j, k, fac

        if (recon_type == WENO_TYPE) then
            weno_polyn = (weno_order - 1)/2
        else if (recon_type == MUSCL_TYPE) then
            muscl_polyn = muscl_order
        end if

        ! Volume Fraction Model (model_eqns == 2)
        cont_idx%beg = 1
        cont_idx%end = num_fluids
        mom_idx%beg = cont_idx%end + 1
        mom_idx%end = cont_idx%end + num_dims
        E_idx = mom_idx%end + 1

        adv_idx%beg = E_idx + 1
        adv_idx%end = E_idx + num_fluids - 1

        sys_size = adv_idx%end

        alf_idx = 1

        momxb = mom_idx%beg
        momxe = mom_idx%end
        advxb = adv_idx%beg
        advxe = adv_idx%end
        contxb = cont_idx%beg
        contxe = cont_idx%end

        Re_size = 0
        Re_size_max = 0
        do i = 1, num_fluids
            if (fluid_pp(i)%Re(1) > 0) Re_size(1) = Re_size(1) + 1
            if (fluid_pp(i)%Re(2) > 0) Re_size(2) = Re_size(2) + 1
        end do
        Re_size_max = maxval(Re_size)
        if (Re_size_max > 0) then
            allocate (Re_idx(1:2,1:Re_size_max))
            k = 0
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(1) > 0) then; k = k + 1; Re_idx(1, k) = i; end if
            end do
            k = 0
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(2) > 0) then; k = k + 1; Re_idx(2, k) = i; end if
            end do
        end if

        call s_configure_coordinate_bounds(igr_order, buff_size, idwint, idwbuff, viscous, m, n, p, num_dims)

#ifdef MFC_MPI
        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))

        if (.not. down_sample) then
            do i = 1, sys_size
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
        end if
#endif

        ! Allocating grid variables for the x-direction
        allocate (x_cc(0:m), x_cb(-1:m))
        ! Allocating grid variables for the y- and z-directions
        if (n > 0) then
            allocate (y_cc(0:n), y_cb(-1:n))
            if (p > 0) then
                allocate (z_cc(0:p), z_cb(-1:p))
            end if
        end if

        grid_geometry = 1  ! Always Cartesian in IGR-only build

    end subroutine s_initialize_global_parameters_module

    !> Configure MPI parallel I/O settings and allocate processor coordinate arrays.
    impure subroutine s_initialize_parallel_io

#ifdef MFC_MPI
        integer :: ierr  !< Generic flag used to identify and report MPI errors
#endif

        num_dims = 1 + min(1, n) + min(1, p)

        num_vels = num_dims

        allocate (proc_coords(1:num_dims))

        if (parallel_io .neqv. .true.) return

#ifdef MFC_MPI
        ! Option for Lustre file system (Darter/Comet/Stampede)
        write (mpiiofs, '(A)') '/lustre_'
        mpiiofs = trim(mpiiofs)
        call MPI_INFO_CREATE(mpi_info_int, ierr)
        call MPI_INFO_SET(mpi_info_int, 'romio_ds_write', 'disable', ierr)

        ! Option for UNIX file system (Hooke/Thomson) WRITE(mpiiofs, '(A)') '/ufs_' mpiiofs = TRIM(mpiiofs) mpi_info_int =
        ! MPI_INFO_NULL

        allocate (start_idx(1:num_dims))
#endif

    end subroutine s_initialize_parallel_io

    !> Deallocate all global grid, index, and equation-of-state parameter arrays.
    impure subroutine s_finalize_global_parameters_module

        integer :: i

        ! Deallocating grid variables for the x-direction

        deallocate (x_cc, x_cb)
        ! Deallocating grid variables for the y- and z-directions
        if (n > 0) then
            deallocate (y_cc, y_cb)
            if (p > 0) then
                deallocate (z_cc, z_cb)
            end if
        end if

        deallocate (proc_coords)

#ifdef MFC_MPI
        if (parallel_io) then
            deallocate (start_idx)
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if
#endif

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
