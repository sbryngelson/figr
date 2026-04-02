
#:include 'case.fpp'

module m_global_parameters

    use mpi  !< Message passing interface (MPI) module

    use m_derived_types
    use m_helper_basic

    implicit none

    integer                 :: num_procs  !< Number of processors
    character(LEN=path_len) :: case_dir   !< Case folder location

    ! Computational Domain Parameters

    integer :: proc_rank  !< Rank of the local processor
    integer :: m, m_root
    integer :: n
    integer :: p

    type(cell_num_bounds) :: cells_bounds
    integer(kind=8)       :: nGlobal  !< Total number of cells in global domain
    integer               :: grid_geometry

    integer :: m_glb, n_glb, p_glb

    integer :: num_dims  !< Number of spatial dimensions
    integer :: num_vels  !< Number of velocity components
    real(wp), allocatable, dimension(:) :: x_cb, x_root_cb, y_cb, z_cb

    real(wp), allocatable, dimension(:) :: x_cc, x_root_cc, y_cc, z_cc
    real(sp), allocatable, dimension(:) :: x_root_cc_s, x_cc_s

    !> Cell-width distributions in the x-, y- and z-coordinate directions
    real(wp), allocatable, dimension(:) :: dx, dy, dz

    integer                              :: buff_size     !< Number of ghost cells for boundary condition storage
    integer, dimension(2)                :: Re_size
    integer                              :: Re_size_max
    integer, allocatable, dimension(:,:) :: Re_idx
    integer                              :: t_step_start  !< First time-step directory
    integer                              :: t_step_stop   !< Last time-step directory
    integer                              :: t_step_save   !< Interval between consecutive time-step directory
    logical  :: cfl_adap_dt, cfl_const_dt, cfl_dt
    real(wp) :: t_save
    real(wp) :: t_stop
    real(wp) :: cfl_target
    integer  :: n_save
    integer  :: n_start

    ! NOTE: m_root, x_root_cb, x_root_cc = defragmented grid (1D only; equals m, x_cb, x_cc in serial)

    integer :: model_eqns  !< Multicomponent flow model
    integer :: num_fluids  !< Number of different fluids present in the flow
    integer :: sys_size    !< Number of unknowns in the system of equations
    integer :: igr_order   !< IGR reconstruction order
    type(int_bounds_info) :: cont_idx  !< Indexes of first & last continuity eqns.
    type(int_bounds_info) :: mom_idx   !< Indexes of first & last momentum eqns.
    integer               :: E_idx     !< Index of energy equation
    type(int_bounds_info) :: adv_idx   !< Indexes of first & last advection eqns.
    integer               :: alf_idx   !< Index of void fraction

    ! Cell Indices for the (local) interior points (O-m, O-n, 0-p). Stands for "InDices With BUFFer".
    type(int_bounds_info) :: idwint(1:3)

    ! Cell indices (InDices With BUFFer): includes buffer in simulation only
    type(int_bounds_info) :: idwbuff(1:3)
    integer               :: num_bc_patches
    logical               :: bc_io
    type(int_bounds_info) :: bc_x, bc_y, bc_z

    logical                            :: parallel_io       !< Format of the data files
    logical                            :: sim_data
    logical                            :: file_per_process  !< output format
    integer, allocatable, dimension(:) :: proc_coords       !< Processor coordinates in MPI_CART_COMM
    integer, allocatable, dimension(:) :: start_idx         !< Starting cell-center index of local processor in global grid
    type(mpi_io_var), public :: MPI_IO_DATA

    character(LEN=name_len) :: mpiiofs
    integer                 :: mpi_info_int

    type(physical_parameters), dimension(num_fluids_max) :: fluid_pp  !< Stiffened gas EOS parameters and Reynolds numbers per fluid
    real(wp), allocatable, dimension(:)                  :: adv       !< Advection variables
    ! Formatted Database File(s) Structure Parameters

    integer               :: format                                    !< Format of the database file(s)
    integer               :: precision                                 !< Floating point precision of the database file(s)
    logical               :: viscous                                   !< Viscous effects
    logical               :: down_sample                               !< down sampling of the database file(s)
    logical               :: output_partial_domain                     !< Specify portion of domain to output for post-processing
    type(bounds_info)     :: x_output, y_output, z_output              !< Portion of domain to output for post-processing
    type(int_bounds_info) :: x_output_idx, y_output_idx, z_output_idx  !< Indices of domain to output for post-processing
    !! necessary when using the Silo database file format in multidimensions. These zones provide VisIt with the subdomain
    !! connectivity information that it requires in order to produce smooth plots.
    type(int_bounds_info) :: offset_x, offset_y, offset_z

    !! momentum, velocity, energy, pressure, volume fraction(s), specific heat ratio function, specific heat ratio, liquid stiffness
    !! function, liquid stiffness, primitive variables, conservative variables, speed of sound, the vorticity, and the numerical
    !! Schlieren function.
    logical, dimension(num_fluids_max) :: alpha_rho_wrt
    logical                            :: rho_wrt
    logical, dimension(3)              :: mom_wrt
    logical, dimension(3)              :: vel_wrt
    integer                            :: flux_lim
    logical, dimension(3)              :: flux_wrt
    logical                            :: E_wrt
    logical, dimension(num_fluids_max) :: alpha_rho_e_wrt
    logical                            :: dummy  !< AMDFlang workaround for case-optimization + GPU-kernel bug
    logical                            :: pres_wrt
    logical, dimension(num_fluids_max) :: alpha_wrt
    logical                            :: gamma_wrt
    logical                            :: heat_ratio_wrt
    logical                            :: pi_inf_wrt
    logical                            :: pres_inf_wrt
    logical                            :: prim_vars_wrt
    logical                            :: cons_vars_wrt
    logical                            :: c_wrt
    logical, dimension(3)              :: omega_wrt
    logical                            :: qm_wrt
    logical                            :: schlieren_wrt

    real(wp), dimension(num_fluids_max) :: schlieren_alpha  !< Per-fluid Schlieren intensity amplitude coefficients
    integer                             :: fd_order         !< Finite-difference order for vorticity and Schlieren derivatives
    integer                             :: fd_number        !< Finite-difference half-stencil size: MAX(1, fd_order/2)
    real(wp) :: rhoref, pref

    integer :: momxb, momxe
    integer :: advxb, advxe
    integer :: contxb, contxe

    real(wp) :: wall_time, wall_time_avg  !< Wall time measurements

contains

    !> Assigns default values to user inputs prior to reading them in. This allows for an easier consistency check of these
    !! parameters once they are read from the input file.
    impure subroutine s_assign_default_values_to_user_inputs

        integer :: i  !< Generic loop iterator
        ! Logistics

        case_dir = '.'

        ! Computational domain parameters
        m = dflt_int; n = 0; p = 0
        call s_update_cell_bounds(cells_bounds, m, n, p)

        m_root = dflt_int

        t_step_start = dflt_int
        t_step_stop = dflt_int
        t_step_save = dflt_int

        cfl_adap_dt = .false.
        cfl_const_dt = .false.
        cfl_dt = .false.
        cfl_target = dflt_real
        t_save = dflt_real
        n_start = dflt_int
        t_stop = dflt_real

        ! Simulation algorithm parameters
        model_eqns = dflt_int
        num_fluids = dflt_int

        bc_x%beg = dflt_int; bc_x%end = dflt_int
        bc_y%beg = dflt_int; bc_y%end = dflt_int
        bc_z%beg = dflt_int; bc_z%end = dflt_int
        bc_io = .false.
        num_bc_patches = dflt_int

        #:for DIM in ['x', 'y', 'z']
            #:for DIR in [1, 2, 3]
                bc_${DIM}$%vb${DIR}$ = 0._wp
                bc_${DIM}$%ve${DIR}$ = 0._wp
            #:endfor
        #:endfor

        ! Fluids physical parameters
        do i = 1, num_fluids_max
            fluid_pp(i)%gamma = dflt_real
            fluid_pp(i)%pi_inf = dflt_real
            fluid_pp(i)%cv = 0._wp
            fluid_pp(i)%qv = 0._wp
            fluid_pp(i)%qvp = 0._wp
        end do

        ! Formatted database file(s) structure parameters
        format = dflt_int

        precision = dflt_int
        down_sample = .false.

        alpha_rho_wrt = .false.
        alpha_rho_e_wrt = .false.
        rho_wrt = .false.
        mom_wrt = .false.
        vel_wrt = .false.
        flux_lim = dflt_int
        flux_wrt = .false.
        parallel_io = .false.
        file_per_process = .false.
        E_wrt = .false.
        dummy = .false.
        pres_wrt = .false.
        alpha_wrt = .false.
        gamma_wrt = .false.
        heat_ratio_wrt = .false.
        pi_inf_wrt = .false.
        pres_inf_wrt = .false.
        prim_vars_wrt = .false.
        cons_vars_wrt = .false.
        c_wrt = .false.
        omega_wrt = .false.
        qm_wrt = .false.
        schlieren_wrt = .false.
        sim_data = .false.

        schlieren_alpha = dflt_real

        fd_order = dflt_int

        ! Tait EOS
        rhoref = dflt_real
        pref = dflt_real

        ! Output partial domain
        output_partial_domain = .false.
        x_output%beg = dflt_real
        x_output%end = dflt_real
        y_output%beg = dflt_real
        y_output%end = dflt_real
        z_output%beg = dflt_real
        z_output%end = dflt_real

    end subroutine s_assign_default_values_to_user_inputs

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_global_parameters_module

        integer :: i

        ! Setting m_root equal to m in the case of a 1D serial simulation

        if (n == 0) m_root = m_glb

        ! Volume Fraction Model (5-equation model) Annotating structure of the state and flux vectors belonging to the system of
        ! equations defined by the selected number of spatial dimensions and the volume fraction model
        cont_idx%beg = 1
        cont_idx%end = num_fluids
        mom_idx%beg = cont_idx%end + 1
        mom_idx%end = cont_idx%end + num_vels
        E_idx = mom_idx%end + 1

        ! Volume fractions are stored in the indices immediately following the energy equation. IGR tracks a total of (N-1) volume
        ! fractions for N fluids, hence the "-1" in adv_idx%end. If num_fluids = 1 then adv_idx%end < adv_idx%beg, which skips all
        ! loops over the volume fractions since there is no volume fraction to track
        adv_idx%beg = E_idx + 1  ! Alpha for fluid 1
        adv_idx%end = E_idx + num_fluids - 1

        sys_size = adv_idx%end

        alf_idx = 1

        ! Chemistry removed (IGR-only build)

        if (output_partial_domain) then
            x_output_idx%beg = 0
            x_output_idx%end = 0
            y_output_idx%beg = 0
            y_output_idx%end = 0
            z_output_idx%beg = 0
            z_output_idx%end = 0
        end if

        momxb = mom_idx%beg
        momxe = mom_idx%end
        advxb = adv_idx%beg
        advxe = adv_idx%end
        contxb = cont_idx%beg
        contxe = cont_idx%end

        allocate (MPI_IO_DATA%view(1:sys_size))
        allocate (MPI_IO_DATA%var(1:sys_size))

        do i = 1, sys_size
            if (down_sample) then
                allocate (MPI_IO_DATA%var(i)%sf(-1:m + 1,-1:n + 1,-1:p + 1))
            else
                allocate (MPI_IO_DATA%var(i)%sf(0:m,0:n,0:p))
            end if
            MPI_IO_DATA%var(i)%sf => null()
        end do

        ! Size of the ghost zone layer is non-zero only when post-processing the raw simulation data of a parallel multidimensional
        ! computation in the Silo-HDF5 format. If this is the case, one must also verify whether the raw simulation data is 2D or
        ! 3D. In the 2D case, size of the z-coordinate direction ghost zone layer must be zeroed out.
        if (num_procs == 1 .or. format /= 1) then
            offset_x%beg = 0
            offset_x%end = 0
            offset_y%beg = 0
            offset_y%end = 0
            offset_z%beg = 0
            offset_z%end = 0
        else if (n == 0) then
            offset_y%beg = 0
            offset_y%end = 0
            offset_z%beg = 0
            offset_z%end = 0
        else if (p == 0) then
            offset_z%beg = 0
            offset_z%end = 0
        end if

        ! Determining the finite-difference number and the buffer size. Note that the size of the buffer is unrelated to the order
        ! of the WENO scheme. Rather, it is directly dependent on maximum size of ghost zone layers and possibly the order of the
        ! finite difference scheme used for the computation of vorticity and/or numerical Schlieren function.
        buff_size = max(offset_x%beg, offset_x%end, offset_y%beg, offset_y%end, offset_z%beg, offset_z%end)

        if (any(omega_wrt) .or. schlieren_wrt .or. qm_wrt) then
            fd_number = max(1, fd_order/2)
            buff_size = buff_size + fd_number
        end if

        ! Configuring Coordinate Direction Indexes
        idwint(1)%beg = 0; idwint(2)%beg = 0; idwint(3)%beg = 0
        idwint(1)%end = m; idwint(2)%end = n; idwint(3)%end = p

        idwbuff(1)%beg = -buff_size
        if (num_dims > 1) then; idwbuff(2)%beg = -buff_size; else; idwbuff(2)%beg = 0; end if
        if (num_dims > 2) then; idwbuff(3)%beg = -buff_size; else; idwbuff(3)%beg = 0; end if

        idwbuff(1)%end = idwint(1)%end - idwbuff(1)%beg
        idwbuff(2)%end = idwint(2)%end - idwbuff(2)%beg
        idwbuff(3)%end = idwint(3)%end - idwbuff(3)%beg

        ! Allocating single precision grid variables if needed
        allocate (x_cc_s(-buff_size:m + buff_size))

        ! Allocating the grid variables in the x-coordinate direction
        allocate (x_cb(-1 - offset_x%beg:m + offset_x%end))
        allocate (x_cc(-buff_size:m + buff_size))
        allocate (dx(-buff_size:m + buff_size))

        ! Allocating grid variables in the y- and z-coordinate directions
        if (n > 0) then
            allocate (y_cb(-1 - offset_y%beg:n + offset_y%end))
            allocate (y_cc(-buff_size:n + buff_size))
            allocate (dy(-buff_size:n + buff_size))

            if (p > 0) then
                allocate (z_cb(-1 - offset_z%beg:p + offset_z%end))
                allocate (z_cc(-buff_size:p + buff_size))
                allocate (dz(-buff_size:p + buff_size))
            end if

            ! Allocating the grid variables, only used for the 1D simulations, and containing the defragmented computational domain
            ! grid data
        else
            allocate (x_root_cb(-1:m_root))
            allocate (x_root_cc(0:m_root))

            if (precision == 1) then
                allocate (x_root_cc_s(0:m_root))
            end if
        end if

        allocate (adv(num_fluids))

        grid_geometry = 1

    end subroutine s_initialize_global_parameters_module

    !> Subroutine to initialize parallel infrastructure
    impure subroutine s_initialize_parallel_io

        integer :: ierr  !< Generic flag used to identify and report MPI errors

        num_dims = 1 + min(1, n) + min(1, p)

        num_vels = num_dims

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

    !> Deallocation procedures for the module
    impure subroutine s_finalize_global_parameters_module

        integer :: i

        ! Deallocating the grid variables for the x-coordinate direction

        deallocate (x_cc, x_cb, dx)

        ! Deallocating grid variables for the y- and z-coordinate directions
        if (n > 0) then
            deallocate (y_cc, y_cb, dy)
            if (p > 0) then
                deallocate (z_cc, z_cb, dz)
            end if
        else
            ! Deallocating the grid variables, only used for the 1D simulations, and containing the defragmented computational
            ! domain grid data
            deallocate (x_root_cb, x_root_cc)
        end if

        deallocate (proc_coords)

        deallocate (adv)

        if (parallel_io) then
            deallocate (start_idx)
            do i = 1, sys_size
                MPI_IO_DATA%var(i)%sf => null()
            end do

            deallocate (MPI_IO_DATA%var)
            deallocate (MPI_IO_DATA%view)
        end if

    end subroutine s_finalize_global_parameters_module

end module m_global_parameters
