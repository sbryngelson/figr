!>
!! @file
!! @brief Contains module m_global_parameters

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Global parameters for the computational domain, fluid properties, and simulation algorithm configuration
module m_global_parameters

    use mpi  !< Message passing interface (MPI) module

    use m_derived_types
    use m_helper_basic
    ! $:USE_GPU_MODULE()

    implicit none

    real(wp) :: wall_time = 0
    real(wp) :: wall_time_avg = 0

    ! Logistics
    integer                 :: num_procs      !< Number of processors
    character(LEN=path_len) :: case_dir       !< Case folder location
    logical                 :: run_time_info  !< Run-time output flag
    integer                 :: t_step_old     !< Existing IC/grid folder
    ! Computational Domain Parameters
    integer :: proc_rank  !< Rank of the local processor
    !> @name Number of cells in the x-, y- and z-directions, respectively
    !> @{
    integer :: m, n, p
    !> @}

    !> @name Max and min number of cells in a direction of each combination of x-,y-, and z-
    type(cell_num_bounds) :: cells_bounds

    !> @name Global number of cells in each direction
    !> @{
    integer :: m_glb, n_glb, p_glb
    !> @}

    integer, parameter :: grid_geometry = 1  !< Grid geometry (always Cartesian in IGR-only build)

    !> @name Cell-boundary (CB) locations in the x-, y- and z-directions, respectively
    !> @{
    real(wp), target, allocatable, dimension(:) :: x_cb, y_cb, z_cb
    !> @}

    !> @name Cell-center (CC) locations in the x-, y- and z-directions, respectively
    !> @{
    real(wp), target, allocatable, dimension(:) :: x_cc, y_cc, z_cc
    !> @}
    ! type(bounds_info) :: x_domain, y_domain, z_domain !< Locations of the domain bounds in the x-, y- and z-coordinate directions
    !> @name Cell-width distributions in the x-, y- and z-directions, respectively
    !> @{
    real(wp), target, allocatable, dimension(:) :: dx, dy, dz
    !> @}

    real(wp) :: dt  !< Size of the time-step
    $:GPU_DECLARE(create='[x_cb, y_cb, z_cb, x_cc, y_cc, z_cc, dx, dy, dz, dt, m, n, p]')

    !> @name Starting time-step iteration, stopping time-step iteration and the number of time-step iterations between successive
    !! solution backups, respectively
    !> @{
    integer :: t_step_start, t_step_stop, t_step_save
    !> @}

    !> @name Starting time, stopping time, and time between backups, simulation time, and prescribed cfl respectively
    !> @{
    real(wp) :: t_stop, t_save, cfl_target
    integer  :: n_start
    !> @}
    $:GPU_DECLARE(create='[cfl_target]')

    logical :: cfl_adap_dt, cfl_const_dt, cfl_dt
    integer :: t_step_print  !< Number of time-steps between printouts
    ! Simulation Algorithm Parameters
    integer :: model_eqns  !< Multicomponent flow model
    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: num_dims = ${num_dims}$  !< Number of spatial dimensions
    #:else
        integer :: num_dims  !< Number of spatial dimensions
    #:endif
    ! num_vels removed: always equals num_dims in IGR-only build. Use num_dims directly. Legacy alias kept for compatibility with
    ! GPU routines that declare dimension(num_vels):
    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: num_vels = num_dims
    #:else
        integer :: num_vels
    #:endif
    integer :: time_stepper  !< Time-stepper algorithm
    logical :: prim_vars_wrt

    #:if MFC_CASE_OPTIMIZATION
        integer, parameter :: recon_type = ${recon_type}$    !< Reconstruction type
        integer, parameter :: weno_polyn = ${weno_polyn}$    !< Degree of the WENO polynomials (polyn)
        integer, parameter :: muscl_polyn = ${muscl_polyn}$  !< Degree of the MUSCL polynomials (polyn)
        integer, parameter :: weno_order = ${weno_order}$    !< Order of the WENO reconstruction
        integer, parameter :: muscl_order = ${muscl_order}$  !< Order of the MUSCL order
        !> Number of stencils for WENO reconstruction (only different from weno_polyn for TENO(>5))
        integer, parameter  :: weno_num_stencils = ${weno_num_stencils}$
        integer, parameter  :: muscl_lim = ${muscl_lim}$               !< MUSCL Limiter
        integer, parameter  :: num_fluids = ${num_fluids}$             !< number of fluids in the simulation
        logical, parameter  :: wenojs = (${wenojs}$ /= 0)              !< WENO-JS (default)
        logical, parameter  :: mapped_weno = (${mapped_weno}$ /= 0)    !< WENO-M (WENO with mapping of nonlinear weights)
        logical, parameter  :: wenoz = (${wenoz}$ /= 0)                !< WENO-Z
        logical, parameter  :: teno = (${teno}$ /= 0)                  !< TENO (Targeted ENO)
        real(wp), parameter :: wenoz_q = ${wenoz_q}$                   !< Power constant for WENO-Z
        integer, parameter  :: igr_iter_solver = ${igr_iter_solver}$   !< IGR elliptic solver
        integer, parameter  :: igr_order = ${igr_order}$               !< Reconstruction order for IGR
        logical, parameter  :: igr_pres_lim = (${igr_pres_lim}$ /= 0)  !< Limit to positive pressures for IGR
        logical, parameter  :: viscous = (${viscous}$ /= 0)            !< Viscous effects
    #:else
        integer  :: recon_type         !< Reconstruction Type
        integer  :: weno_polyn         !< Degree of the WENO polynomials (polyn)
        integer  :: muscl_polyn        !< Degree of the MUSCL polynomials (polyn)i
        integer  :: weno_order         !< Order of the WENO reconstruction
        integer  :: muscl_order        !< Order of the MUSCL reconstruction
        integer  :: weno_num_stencils  !< Number of stencils for WENO reconstruction (only different from weno_polyn for TENO(>5))
        integer  :: muscl_lim          !< MUSCL Limiter
        integer  :: num_fluids         !< number of fluids in the simulation
        logical  :: wenojs             !< WENO-JS (default)
        logical  :: mapped_weno        !< WENO-M (WENO with mapping of nonlinear weights)
        logical  :: wenoz              !< WENO-Z
        logical  :: teno               !< TENO (Targeted ENO)
        real(wp) :: wenoz_q            !< Power constant for WENO-Z
        integer  :: igr_iter_solver    !< IGR elliptic solver
        integer  :: igr_order          !< Reconstruction order for IGR
        logical  :: igr_pres_lim       !< Limit to positive pressures for IGR
        logical  :: viscous            !< Viscous effects
    #:endif

    !> @name Variables for our of core IGR computation on NVIDIA
    !> @{
    logical :: nv_uvm_out_of_core       !< Enable out-of-core storage of q_cons_ts(2) in timestepping (default FALSE)
    integer :: nv_uvm_igr_temps_on_gpu  !< 0 => jac, jac_rhs, and jac_old on CPU
    ! 1 => jac on GPU, jac_rhs and jac_old on CPU 2 => jac and jac_rhs on GPU, jac_old on CPU 3 => jac, jac_rhs, and jac_old on GPU
    ! (default)
    logical :: nv_uvm_pref_gpu  !< Enable explicit gpu memory hints (default FALSE)
    !> @}

    real(wp) :: weno_eps                  !< Binding for the WENO nonlinear weights
    real(wp) :: teno_CT                   !< Smoothness threshold for TENO
    logical  :: mp_weno                   !< Monotonicity preserving (MP) WENO
    logical  :: weno_avg                  !< Average left/right cell-boundary states
    logical  :: weno_Re_flux              !< WENO reconstruct velocity gradients for viscous stress tensor
    logical  :: null_weights              !< Null undesired WENO weights
    logical  :: mixture_err               !< Mixture properties correction
    logical  :: shear_stress              !< Shear stresses (computed from Re)
    logical  :: bulk_stress               !< Bulk stresses (computed from Re)
    integer  :: num_igr_iters             !< number of iterations for elliptic solve
    integer  :: num_igr_warm_start_iters  !< number of warm start iterations for elliptic solve
    real(wp) :: alf_factor                !< alpha factor for IGR
    integer  :: cpu_start, cpu_end, cpu_rate

    #:if not MFC_CASE_OPTIMIZATION
        $:GPU_DECLARE(create='[num_dims, weno_polyn, weno_order]')
        $:GPU_DECLARE(create='[weno_num_stencils, num_fluids, wenojs]')
        $:GPU_DECLARE(create='[mapped_weno, wenoz, teno, wenoz_q]')
        $:GPU_DECLARE(create='[igr_iter_solver, igr_order, viscous, igr_pres_lim]')
        $:GPU_DECLARE(create='[recon_type, muscl_order, muscl_polyn, muscl_lim]')
    #:endif

    $:GPU_DECLARE(create='[model_eqns, mixture_err]')
    $:GPU_DECLARE(create='[mp_weno, weno_eps, teno_CT]')
    $:GPU_DECLARE(create='[shear_stress, bulk_stress]')

    integer :: num_bc_patches
    logical :: bc_io
    !> @name Boundary conditions (BC) in the x-, y- and z-directions, respectively
    !> @{
    type(int_bounds_info) :: bc_x, bc_y, bc_z
    !> @}
#if defined(MFC_OpenACC)
    $:GPU_DECLARE(create='[bc_x%vb1, bc_x%vb2, bc_x%vb3, bc_x%ve1, bc_x%ve2, bc_x%ve3]')
    $:GPU_DECLARE(create='[bc_y%vb1, bc_y%vb2, bc_y%vb3, bc_y%ve1, bc_y%ve2, bc_y%ve3]')
    $:GPU_DECLARE(create='[bc_z%vb1, bc_z%vb2, bc_z%vb3, bc_z%ve1, bc_z%ve2, bc_z%ve3]')
#elif defined(MFC_OpenMP)
    $:GPU_DECLARE(create='[bc_x, bc_y, bc_z]')
#endif
    type(bounds_info) :: x_domain, y_domain, z_domain
    $:GPU_DECLARE(create='[x_domain, y_domain, z_domain]')
    real(wp) :: x_a, y_a, z_a
    real(wp) :: x_b, y_b, z_b
    logical  :: parallel_io       !< Format of the data files
    logical  :: file_per_process  !< shared file or not when using parallel io
    integer  :: precision         !< Precision of output files
    logical  :: down_sample       !< down sample the output files
    $:GPU_DECLARE(create='[down_sample]')

    integer, allocatable, dimension(:) :: proc_coords  !< Processor coordinates in MPI_CART_COMM
    integer, allocatable, dimension(:) :: start_idx    !< Starting cell-center index of local processor in global grid
    type(mpi_io_var), public           :: MPI_IO_DATA

    !> @name MPI info for parallel IO with Lustre file systems
    !> @{
    character(LEN=name_len) :: mpiiofs
    integer                 :: mpi_info_int
    !> @}

    !> @name Annotations of the structure of the state and flux vectors in terms of the size and the configuration of the system of
    !! equations to which they belong
    !> @{
    integer               :: sys_size  !< Number of unknowns in system of eqns.
    type(int_bounds_info) :: cont_idx  !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx   !< Indexes of first & last momentum eqns.
    integer               :: E_idx     !< Index of energy equation
    type(int_bounds_info) :: adv_idx   !< Indexes of first & last advection eqns.
    integer               :: alf_idx   !< Index of void fraction
    !> @}
    $:GPU_DECLARE(create='[sys_size, E_idx, alf_idx]')

    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p). Stands for "InDices With INTerior".
    type(int_bounds_info) :: idwint(1:3)
    $:GPU_DECLARE(create='[idwint]')

    ! Cell Indices for the entire (local) domain. In simulation and post_process, this includes the buffer region. idwbuff and
    ! idwint are the same otherwise. Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwbuff(1:3)
    $:GPU_DECLARE(create='[idwbuff]')

    !> @name The number of fluids, along with their identifying indexes, respectively, for which viscous effects, e.g. the shear
    !! and/or the volume Reynolds (Re) numbers, will be non-negligible.
    !> @{
    integer, dimension(2)                :: Re_size
    integer                              :: Re_size_max
    integer, allocatable, dimension(:,:) :: Re_idx
    !> @}

    $:GPU_DECLARE(create='[Re_size, Re_size_max, Re_idx]')

    ! WENO averaging flag: use arithmetic mean or unaltered WENO-reconstructed cell-boundary values
    !> @{
    real(wp) :: wa_flg
    !> @}

    $:GPU_DECLARE(create='[wa_flg]')

    !> @name The coordinate direction indexes and flags (flg), respectively, for which the configurations will be determined with
    !! respect to a working direction and that will be used to isolate the contributions, in that direction, in the dimensionally
    !! split system of equations.
    !> @{
    integer, dimension(3)  :: dir_idx
    real(wp), dimension(3) :: dir_flg
    !> @}

    $:GPU_DECLARE(create='[dir_idx, dir_flg]')

    integer :: buff_size  !< Number of ghost cells for boundary condition storage
    $:GPU_DECLARE(create='[buff_size]')

    ! END: Simulation Algorithm Parameters

    ! Fluids Physical Parameters

    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp  !< Stiffened gas EOS parameters and Reynolds numbers per fluid

    !> @name Reference density and pressure for Tait EOS
    !> @{
    real(wp) :: rhoref, pref
    !> @}
    $:GPU_DECLARE(create='[rhoref, pref]')

    integer :: momxb, momxe
    integer :: advxb, advxe
    integer :: contxb, contxe
    $:GPU_DECLARE(create='[momxb, momxe, advxb, advxe, contxb, contxe]')

    real(wp), allocatable, dimension(:) :: gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps
    $:GPU_DECLARE(create='[gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps]')

    real(wp)                                    :: mytime     !< Current simulation time
    real(wp)                                    :: finaltime  !< Final simulation time
    logical                                     :: rdma_mpi
    type(pres_field), allocatable, dimension(:) :: pb_ts
    type(pres_field), allocatable, dimension(:) :: mv_ts
    $:GPU_DECLARE(create='[pb_ts, mv_ts]')

    logical :: dummy  !< AMDFlang workaround for case-optimization + GPU-kernel bug

contains

    !> Assigns default values to the user inputs before reading them in. This enables for an easier consistency check of these
    !! parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i, j  !< Generic loop iterator
        ! Logistics

        case_dir = '.'
        run_time_info = .false.
        t_step_old = dflt_int

        ! Computational domain parameters
        m = dflt_int; n = 0; p = 0
        call s_update_cell_bounds(cells_bounds, m, n, p)

        dt = dflt_real

        cfl_adap_dt = .false.
        cfl_const_dt = .false.
        cfl_dt = .false.
        cfl_target = dflt_real

        t_step_start = dflt_int
        t_step_stop = dflt_int
        t_step_save = dflt_int
        t_step_print = 1

        n_start = dflt_int
        t_stop = dflt_real
        t_save = dflt_real

        ! NVIDIA UVM options
        nv_uvm_out_of_core = .false.
        nv_uvm_igr_temps_on_gpu = 3  ! => jac, jac_rhs, and jac_old on GPU (default)
        nv_uvm_pref_gpu = .false.

        ! Simulation algorithm parameters
        model_eqns = dflt_int
        time_stepper = dflt_int
        weno_eps = dflt_real
        teno_CT = dflt_real
        mp_weno = .false.
        weno_avg = .false.
        weno_Re_flux = .false.
        null_weights = .false.
        mixture_err = .false.
        parallel_io = .false.
        file_per_process = .false.
        precision = 2
        down_sample = .false.
        rdma_mpi = .false.
        shear_stress = .false.
        bulk_stress = .false.
        num_igr_iters = dflt_num_igr_iters
        num_igr_warm_start_iters = dflt_num_igr_warm_start_iters
        alf_factor = dflt_alf_factor

        #:if not MFC_CASE_OPTIMIZATION
            mapped_weno = .false.
            wenoz = .false.
            teno = .false.
            wenoz_q = dflt_real
            igr_order = dflt_int
            igr_pres_lim = .false.
            viscous = .false.
            igr_iter_solver = 1
        #:endif

        num_bc_patches = 0
        bc_io = .false.

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int

        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                bc_${DIM}$%vb${DIR}$ = 0._wp
                bc_${DIM}$%ve${DIR}$ = 0._wp
            #:endfor
        #:endfor

        x_domain%beg = dflt_real; x_domain%end = dflt_real
        y_domain%beg = dflt_real; y_domain%end = dflt_real
        z_domain%beg = dflt_real; z_domain%end = dflt_real

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
            fluid_pp(i)%Re(:) = dflt_real
        end do

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        #:if not MFC_CASE_OPTIMIZATION
            recon_type = WENO_TYPE
            weno_order = dflt_int
            muscl_order = dflt_int
            muscl_lim = dflt_int
            num_fluids = dflt_int
        #:endif

        dummy = .false.

        ! GRCBC flags
        #:for dir in {'x', 'y', 'z'}
            bc_${dir}$%grcbc_in = .false.
            bc_${dir}$%grcbc_out = .false.
            bc_${dir}$%grcbc_vel_out = .false.
        #:endfor

    end subroutine s_assign_default_values_to_user_inputs

    !> Initialize the global parameters module
    impure subroutine s_initialize_global_parameters_module

        integer :: i, j, k
        integer :: fac

        #:if not MFC_CASE_OPTIMIZATION
            ! Determining the degree of the WENO polynomials

            if (recon_type == WENO_TYPE) then
                weno_polyn = (weno_order - 1)/2
                if (teno) then
                    weno_num_stencils = weno_order - 3
                else
                    weno_num_stencils = weno_polyn
                end if
            else if (recon_type == MUSCL_TYPE) then
                muscl_polyn = muscl_order
            end if
            $:GPU_UPDATE(device='[weno_polyn, muscl_polyn]')
            $:GPU_UPDATE(device='[weno_num_stencils]')
            $:GPU_UPDATE(device='[num_dims, num_fluids]')
            $:GPU_UPDATE(device='[igr_order, igr_iter_solver]')
        #:endif

        ! Initialize counts: viscous fluids, surface-tension interfaces, curvature interfaces
        Re_size = 0
        Re_size_max = 0

        ! Volume Fraction Model (model_eqns == 2)
        cont_idx%beg = 1
        cont_idx%end = num_fluids
        mom_idx%beg = cont_idx%end + 1
        mom_idx%end = cont_idx%end + num_dims
        E_idx = mom_idx%end + 1

        ! IGR: volume fractions after energy (N-1 for N fluids; skipped when num_fluids=1)
        adv_idx%beg = E_idx + 1  ! Alpha for fluid 1
        adv_idx%end = E_idx + num_fluids - 1

        sys_size = adv_idx%end

        alf_idx = 1

        ! Count fluids with non-negligible viscous effects (Re > 0)
        do i = 1, num_fluids
            if (fluid_pp(i)%Re(1) > 0) Re_size(1) = Re_size(1) + 1
            if (fluid_pp(i)%Re(2) > 0) Re_size(2) = Re_size(2) + 1
        end do

        if (Re_size(1) > 0._wp) shear_stress = .true.
        if (Re_size(2) > 0._wp) bulk_stress = .true.

        Re_size_max = maxval(Re_size)

        $:GPU_UPDATE(device='[Re_size, Re_size_max, shear_stress, bulk_stress]')

        if (viscous) then
            @:ALLOCATE(Re_idx(1:2, 1:Re_size_max))

            k = 0
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(1) > 0) then
                    k = k + 1; Re_idx(1, k) = i
                end if
            end do

            k = 0
            do i = 1, num_fluids
                if (fluid_pp(i)%Re(2) > 0) then
                    k = k + 1; Re_idx(2, k) = i
                end if
            end do
        end if

        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))

        if (.not. down_sample) then
            do i = 1, sys_size
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
                MPI_IO_DATA%var(i)%sf => null()
            end do
        end if

        ! Configure WENO averaging flag (arithmetic mean vs. unaltered values)
        wa_flg = 0._wp; if (weno_avg) wa_flg = 1._wp
        $:GPU_UPDATE(device='[wa_flg]')

        ! Resort to default WENO-JS if no other WENO scheme is selected
        #:if not MFC_CASE_OPTIMIZATION
            wenojs = .not. (mapped_weno .or. wenoz .or. teno)
        #:endif

        call s_configure_coordinate_bounds(igr_order, buff_size, idwint, idwbuff, viscous, m, n, p, num_dims)
        $:GPU_UPDATE(device='[idwint, idwbuff]')

        momxb = mom_idx%beg
        momxe = mom_idx%end
        advxb = adv_idx%beg
        advxe = adv_idx%end
        contxb = cont_idx%beg
        contxe = cont_idx%end
        $:GPU_UPDATE(device='[momxb, momxe, advxb, advxe, contxb, contxe, sys_size, buff_size, E_idx, alf_idx]')

        $:GPU_UPDATE(device='[cfl_target, m, n, p]')

        $:GPU_UPDATE(device='[dt, sys_size, buff_size, pref, rhoref, E_idx, alf_idx, model_eqns, mixture_err, mp_weno, weno_eps, teno_CT]')

        #:if not MFC_CASE_OPTIMIZATION
            $:GPU_UPDATE(device='[wenojs, mapped_weno, wenoz, teno]')
            $:GPU_UPDATE(device='[wenoz_q]')
            $:GPU_UPDATE(device='[muscl_order, muscl_lim]')
            $:GPU_UPDATE(device='[igr_order]')
            $:GPU_UPDATE(device='[num_fluids, num_dims, viscous, muscl_lim]')
        #:endif

        $:GPU_UPDATE(device='[dir_idx, dir_flg]')

        ! Allocating grid variables for the x-, y- and z-directions
        @:ALLOCATE(x_cb(-1 - buff_size:m + buff_size))
        @:ALLOCATE(x_cc(-buff_size:m + buff_size))
        @:ALLOCATE(dx(-buff_size:m + buff_size))
        @:PREFER_GPU(x_cb)
        @:PREFER_GPU(x_cc)
        @:PREFER_GPU(dx)

        if (n == 0) return
        @:ALLOCATE(y_cb(-1 - buff_size:n + buff_size))
        @:ALLOCATE(y_cc(-buff_size:n + buff_size))
        @:ALLOCATE(dy(-buff_size:n + buff_size))
        @:PREFER_GPU(y_cb)
        @:PREFER_GPU(y_cc)
        @:PREFER_GPU(dy)

        if (p == 0) return
        @:ALLOCATE(z_cb(-1 - buff_size:p + buff_size))
        @:ALLOCATE(z_cc(-buff_size:p + buff_size))
        @:ALLOCATE(dz(-buff_size:p + buff_size))
        @:PREFER_GPU(z_cb)
        @:PREFER_GPU(z_cc)
        @:PREFER_GPU(dz)

    end subroutine s_initialize_global_parameters_module

    !> Initializes parallel infrastructure
    impure subroutine s_initialize_parallel_io

        integer :: ierr  !< Generic flag used to identify and report MPI errors

        #:if not MFC_CASE_OPTIMIZATION
            num_dims = 1 + min(1, n) + min(1, p)
            num_vels = num_dims
        #:endif

        allocate (proc_coords(1:num_dims))

        if (parallel_io .neqv. .true.) return

        ! Option for Lustre file system (Darter/Comet/Stampede)
        write (mpiiofs, '(A)') '/lustre_'
        mpiiofs = trim(mpiiofs)

        call MPI_INFO_CREATE(mpi_info_int, ierr)
        call MPI_INFO_SET(mpi_info_int, 'romio_ds_write', 'disable', ierr)

        ! Option for UNIX file system (Hooke/Thomson) WRITE(mpiiofs, '(A)') '/ufs_' mpiiofs = TRIM(mpiiofs) mpi_info_int =
        ! MPI_INFO_NULL

        allocate (start_idx(1:num_dims))

    end subroutine s_initialize_parallel_io

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_global_parameters_module

        integer :: i

        ! Deallocating the variables bookkeeping the indexes of any viscous fluids and any pairs of fluids whose interfaces
        ! supported effects of surface tension

        if (viscous) then
            @:DEALLOCATE(Re_idx)
        end if

        deallocate (proc_coords)
        if (parallel_io) then
            deallocate (start_idx)

            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

        ! Deallocating grid variables for the x-, y- and z-directions
        @:DEALLOCATE(x_cb, x_cc, dx)

        if (n == 0) return
        @:DEALLOCATE(y_cb, y_cc, dy)

        if (p == 0) return
        @:DEALLOCATE(z_cb, z_cc, dz)

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
