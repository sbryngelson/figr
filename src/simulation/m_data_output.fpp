
#:include 'macros.fpp'
#:include 'case.fpp'

module m_data_output

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_compile_specific
    use m_helper
    use m_helper_basic
    use m_sim_helpers
    use m_boundary_common

    implicit none

    private
    public :: s_initialize_data_output_module, s_open_run_time_information_file, s_write_run_time_information, &
        & s_write_data_files, s_write_serial_data_files, s_write_parallel_data_files, s_close_run_time_information_file, &
        & s_finalize_data_output_module

    real(wp), allocatable, dimension(:,:,:) :: icfl_sf  !< ICFL stability criterion
    real(wp), allocatable, dimension(:,:,:) :: vcfl_sf  !< VCFL stability criterion
    real(wp), allocatable, dimension(:,:,:) :: ccfl_sf  !< CCFL stability criterion
    real(wp), allocatable, dimension(:,:,:) :: Rc_sf    !< Rc stability criterion
    $:GPU_DECLARE(create='[icfl_sf, vcfl_sf, ccfl_sf, Rc_sf]')

    real(wp) :: icfl_max_loc, icfl_max_glb  !< ICFL stability extrema on local and global grids
    real(wp) :: vcfl_max_loc, vcfl_max_glb  !< VCFL stability extrema on local and global grids
    real(wp) :: ccfl_max_loc, ccfl_max_glb  !< CCFL stability extrema on local and global grids
    real(wp) :: Rc_min_loc, Rc_min_glb      !< Rc stability extrema on local and global grids
    $:GPU_DECLARE(create='[icfl_max_loc, icfl_max_glb, vcfl_max_loc, vcfl_max_glb]')
    $:GPU_DECLARE(create='[ccfl_max_loc, ccfl_max_glb, Rc_min_loc, Rc_min_glb]')

    real(wp) :: icfl_max  !< ICFL criterion maximum
    real(wp) :: vcfl_max  !< VCFL criterion maximum
    real(wp) :: ccfl_max  !< CCFL criterion maximum
    real(wp) :: Rc_min    !< Rc criterion maximum

    type(scalar_field), allocatable, dimension(:) :: q_cons_temp_ds

contains

    !> Write data files. Dispatch subroutine that replaces procedure pointer.
    impure subroutine s_write_data_files(q_cons_vf, q_prim_vf, t_step, bc_type, beta)

        type(scalar_field), dimension(sys_size), intent(inout)      :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout)      :: q_prim_vf
        integer, intent(in)                                         :: t_step
        type(scalar_field), intent(inout), optional                 :: beta
        type(integer_field), dimension(1:num_dims,-1:1), intent(in) :: bc_type

        if (.not. parallel_io) then
            call s_write_serial_data_files(q_cons_vf, q_prim_vf, t_step, bc_type, beta)
        else
            call s_write_parallel_data_files(q_cons_vf, t_step, bc_type, beta)
        end if

    end subroutine s_write_data_files

    !> Open the run-time information file and write the stability criteria table header
    impure subroutine s_open_run_time_information_file

        character(LEN=name_len), parameter :: file_name = 'run_time.inf'  !< Name of the run-time information file
        character(LEN=path_len + name_len) :: file_path                   !< Relative path to a file in the case directory
        character(LEN=8)                   :: file_date                   !< Creation date of the run-time information file

        file_path = trim(case_dir) // '/' // trim(file_name)

        open (3, FILE=trim(file_path), form='formatted', STATUS='replace')

        write (3, '(A)') 'Description: Stability information at ' // 'each time-step of the simulation. This'
        write (3, '(13X,A)') 'data is composed of the inviscid ' // 'Courant-Friedrichs-Lewy (ICFL)'
        write (3, '(13X,A)') 'number, the viscous CFL (VCFL) number, ' // 'the capillary CFL (CCFL)'
        write (3, '(13X,A)') 'number and the cell Reynolds (Rc) ' // 'number. Please note that only'
        write (3, '(13X,A)') 'those stability conditions pertinent ' // 'to the physics included in'
        write (3, '(13X,A)') 'the current computation are displayed.'

        call date_and_time(DATE=file_date)

        write (3, '(A)') 'Date: ' // file_date(5:6) // '/' // file_date(7:8) // '/' // file_date(3:4)

        write (3, '(A)') ''; write (3, '(A)') ''

        write (3, '(13X,A9,13X,A10,13X,A10,13X,A10)', advance="no") trim('Time-step'), trim('dt'), trim('Time'), trim('ICFL Max')

        if (viscous) then
            write (3, '(13X,A10,13X,A16)', advance="no") trim('VCFL Max'), trim('Rc Min')
        end if

        write (3, *)  ! new line

    end subroutine s_open_run_time_information_file

    !> Write stability criteria extrema to the run-time information file at the given time step
    impure subroutine s_write_run_time_information(q_prim_vf, t_step)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer, intent(in)                                 :: t_step
        real(wp)                                            :: rho  !< Cell-avg. density

            real(wp), dimension(num_fluids) :: alpha  !< Cell-avg. volume fraction
            real(wp), dimension(num_vels)   :: vel    !< Cell-avg. velocity
        real(wp)               :: vel_sum  !< Cell-avg. velocity sum
        real(wp)               :: pres     !< Cell-avg. pressure
        real(wp)               :: gamma    !< Cell-avg. sp. heat ratio
        real(wp)               :: pi_inf   !< Cell-avg. liquid stiffness function
        real(wp)               :: qv       !< Cell-avg. internal energy reference value
        real(wp)               :: c        !< Cell-avg. sound speed
        real(wp)               :: H        !< Cell-avg. enthalpy
        real(wp), dimension(2) :: Re       !< Cell-avg. Reynolds numbers
        integer                :: j, k, l

        ! Computing Stability Criteria at Current Time-step

        $:GPU_PARALLEL_LOOP(collapse=3, private='[j, k, l, vel, alpha, Re, rho, vel_sum, pres, gamma, pi_inf, c, H, qv]')
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    call s_compute_enthalpy(q_prim_vf, pres, rho, gamma, pi_inf, Re, H, alpha, vel, vel_sum, qv, j, k, l)

                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, alpha, vel_sum, 0._wp, c, qv)

                    if (viscous) then
                        call s_compute_stability_from_dt(vel, c, rho, Re, j, k, l, icfl_sf, vcfl_sf, Rc_sf)
                    else
                        call s_compute_stability_from_dt(vel, c, rho, Re, j, k, l, icfl_sf)
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

#ifdef _CRAYFTN
        $:GPU_UPDATE(host='[icfl_sf]')

        if (viscous) then
            $:GPU_UPDATE(host='[vcfl_sf, Rc_sf]')
        end if

        icfl_max_loc = maxval(icfl_sf)

        if (viscous) then
            vcfl_max_loc = maxval(vcfl_sf)
            Rc_min_loc = minval(Rc_sf)
        end if
#else
        #:call GPU_PARALLEL(copyout='[icfl_max_loc]', copyin='[icfl_sf]')
            icfl_max_loc = maxval(icfl_sf)
        #:endcall GPU_PARALLEL
        if (viscous .or. dummy) then
            #:call GPU_PARALLEL(copyout='[vcfl_max_loc, Rc_min_loc]', copyin='[vcfl_sf,Rc_sf]')
                vcfl_max_loc = maxval(vcfl_sf)
                Rc_min_loc = minval(Rc_sf)
            #:endcall GPU_PARALLEL
        end if
#endif

        if (num_procs > 1) then
            call s_mpi_reduce_stability_criteria_extrema(icfl_max_loc, vcfl_max_loc, Rc_min_loc, icfl_max_glb, vcfl_max_glb, &
                & Rc_min_glb)
        else
            icfl_max_glb = icfl_max_loc
            if (viscous) vcfl_max_glb = vcfl_max_loc
            if (viscous) Rc_min_glb = Rc_min_loc
        end if

        if (icfl_max_glb > icfl_max) icfl_max = icfl_max_glb

        if (viscous) then
            if (vcfl_max_glb > vcfl_max) vcfl_max = vcfl_max_glb
            if (Rc_min_glb < Rc_min) Rc_min = Rc_min_glb
        end if

        if (proc_rank == 0) then
            write (3, '(13X,I9,13X,F10.6,13X,F10.6,13X,F10.6)', advance="no") t_step, dt, mytime, icfl_max_glb

            if (viscous) then
                write (3, '(13X,F10.6,13X,ES16.6)', advance="no") vcfl_max_glb, Rc_min_glb
            end if

            write (3, *)  ! new line

            if (.not. f_approx_equal(icfl_max_glb, icfl_max_glb)) then
                call s_mpi_abort('ICFL is NaN. Exiting.')
            else if (icfl_max_glb > 1._wp) then
                print *, 'icfl', icfl_max_glb
                call s_mpi_abort('ICFL is greater than 1.0. Exiting.')
            end if

            if (viscous) then
                if (.not. f_approx_equal(vcfl_max_glb, vcfl_max_glb)) then
                    call s_mpi_abort('VCFL is NaN. Exiting.')
                else if (vcfl_max_glb > 1._wp) then
                    print *, 'vcfl', vcfl_max_glb
                    call s_mpi_abort('VCFL is greater than 1.0. Exiting.')
                end if
            end if
        end if

        call s_mpi_barrier()

    end subroutine s_write_run_time_information

    !> Write grid and conservative variable data files in serial format
    impure subroutine s_write_serial_data_files(q_cons_vf, q_prim_vf, t_step, bc_type, beta)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_prim_vf
        integer, intent(in) :: t_step
        type(scalar_field), intent(inout), optional :: beta
        type(integer_field), dimension(1:num_dims,-1:1), intent(in) :: bc_type
        character(LEN=path_len + 2*name_len) :: t_step_dir  !< Relative path to the current time-step directory
        character(LEN=path_len + 3*name_len) :: file_path   !< Relative path to the grid and conservative variables data files
        logical :: file_exist                               !< Logical used to check existence of current time-step directory
        character(LEN=15) :: FMT
        integer :: i, j, k, l, r
        real(wp) :: gamma, lit_gamma, pi_inf, qv            !< Temporary EOS params

        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir) // '/p_all'
        write (t_step_dir, '(a,i0,a,i0)') trim(case_dir) // '/p_all/p', proc_rank, '/', t_step

        file_path = trim(t_step_dir) // '/.'
        call my_inquire(file_path, file_exist)
        if (file_exist) call s_delete_directory(trim(t_step_dir))
        call s_create_directory(trim(t_step_dir))

        file_path = trim(t_step_dir) // '/x_cb.dat'

        open (2, FILE=trim(file_path), form='unformatted', STATUS='new')
        write (2) x_cb(-1:m); close (2)

        if (n > 0) then
            file_path = trim(t_step_dir) // '/y_cb.dat'

            open (2, FILE=trim(file_path), form='unformatted', STATUS='new')
            write (2) y_cb(-1:n); close (2)

            if (p > 0) then
                file_path = trim(t_step_dir) // '/z_cb.dat'

                open (2, FILE=trim(file_path), form='unformatted', STATUS='new')
                write (2) z_cb(-1:p); close (2)
            end if
        end if

        do i = 1, sys_size
            write (file_path, '(A,I0,A)') trim(t_step_dir) // '/q_cons_vf', i, '.dat'

            open (2, FILE=trim(file_path), form='unformatted', STATUS='new')

            write (2) q_cons_vf(i)%sf(0:m,0:n,0:p); close (2)
        end do

        gamma = gammas(1)
        lit_gamma = gs_min(1)
        pi_inf = pi_infs(1)
        qv = qvs(1)

        if (precision == 1) then
            FMT = "(2F30.3)"
        else
            FMT = "(2F40.14)"
        end if

        write (t_step_dir, '(A,I0,A,I0)') trim(case_dir) // '/D'
        file_path = trim(t_step_dir) // '/.'

        inquire (FILE=trim(file_path), EXIST=file_exist)

        if (.not. file_exist) call s_create_directory(trim(t_step_dir))

        if (n == 0 .and. p == 0) then
            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step, '.dat'

                open (2, FILE=trim(file_path))
                do j = 0, m
                    write (2, FMT) x_cb(j), q_cons_vf(i)%sf(j, 0, 0)
                end do
                close (2)
            end do
        end if

        if (precision == 1) then
            FMT = "(3F30.7)"
        else
            FMT = "(3F40.14)"
        end if

        if ((n > 0) .and. (p == 0)) then
            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_path))
                do j = 0, m
                    do k = 0, n
                        write (2, FMT) x_cb(j), y_cb(k), q_cons_vf(i)%sf(j, k, 0)
                    end do
                    write (2, *)
                end do
                close (2)
            end do
        end if

        if (precision == 1) then
            FMT = "(4F30.7)"
        else
            FMT = "(4F40.14)"
        end if

        if (p > 0) then
            do i = 1, sys_size
                write (file_path, '(A,I0,A,I2.2,A,I6.6,A)') trim(t_step_dir) // '/cons.', i, '.', proc_rank, '.', t_step, '.dat'
                open (2, FILE=trim(file_path))
                do j = 0, m
                    do k = 0, n
                        do l = 0, p
                            write (2, FMT) x_cb(j), y_cb(k), z_cb(l), q_cons_vf(i)%sf(j, k, l)
                        end do
                        write (2, *)
                    end do
                    write (2, *)
                end do
                close (2)
            end do
        end if

    end subroutine s_write_serial_data_files

    !> Write grid and conservative variable data files in parallel via MPI I/O
    impure subroutine s_write_parallel_data_files(q_cons_vf, t_step, bc_type, beta)

        type(scalar_field), dimension(sys_size), intent(inout)      :: q_cons_vf
        integer, intent(in)                                         :: t_step
        type(scalar_field), intent(inout), optional                 :: beta
        type(integer_field), dimension(1:num_dims,-1:1), intent(in) :: bc_type

        integer                              :: ifile, ierr, data_size
        integer, dimension(MPI_STATUS_SIZE)  :: status
        integer(kind=MPI_OFFSET_kind)        :: disp
        integer(kind=MPI_OFFSET_kind)        :: m_MOK, n_MOK, p_MOK
        integer(kind=MPI_OFFSET_kind)        :: WP_MOK, var_MOK, str_MOK
        integer(kind=MPI_OFFSET_kind)        :: NVARS_MOK
        integer(kind=MPI_OFFSET_kind)        :: MOK
        character(LEN=path_len + 2*name_len) :: file_loc
        logical                              :: file_exist, dir_check
        character(len=10)                    :: t_step_string
        integer                              :: i        !< Generic loop iterator
        integer                              :: alt_sys  !< Altered system size for the lagrangian subgrid bubble model
        ! Down sampling variables
        integer :: m_ds, n_ds, p_ds
        integer :: m_glb_ds, n_glb_ds, p_glb_ds
        integer :: m_glb_save, n_glb_save, p_glb_save  !< Global save size

        if (down_sample) then
            call s_downsample_data(q_cons_vf, q_cons_temp_ds, m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds)
        end if

        if (present(beta)) then
            alt_sys = sys_size + 1
        else
            alt_sys = sys_size
        end if

        if (file_per_process) then
            call s_int_to_str(t_step, t_step_string)

            if (down_sample) then
                call s_initialize_mpi_data_ds(q_cons_temp_ds)
            else
                call s_initialize_mpi_data(q_cons_vf)
            end if

            if (proc_rank == 0) then
                file_loc = trim(case_dir) // '/restart_data/lustre_' // trim(t_step_string)
                call my_inquire(file_loc, dir_check)
                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(file_loc))
                end if
                call s_create_directory(trim(file_loc))
            end if
            call s_mpi_barrier()

            call s_initialize_mpi_data(q_cons_vf)

            write (file_loc, '(I0,A,i7.7,A)') t_step, '_', proc_rank, '.dat'
            file_loc = trim(case_dir) // '/restart_data/lustre_' // trim(t_step_string) // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_SELF, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)

            if (down_sample) then
                data_size = (m_ds + 3)*(n_ds + 3)*(p_ds + 3)
                m_glb_save = m_glb_ds + 1
                n_glb_save = n_glb_ds + 1
                p_glb_save = p_glb_ds + 1
            else
                data_size = (m + 1)*(n + 1)*(p + 1)
                m_glb_save = m_glb + 1
                n_glb_save = n_glb + 1
                p_glb_save = p_glb + 1
            end if

            m_MOK = int(m_glb_save + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb_save + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb_save + 1, MPI_OFFSET_KIND)
            WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(sys_size, MPI_OFFSET_KIND)

            if (down_sample) then
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    call MPI_FILE_WRITE_ALL(ifile, q_cons_temp_ds(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                end do
            else
                do i = 1, sys_size
                    var_MOK = int(i, MPI_OFFSET_KIND)

                    call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
                end do
            end if

            call MPI_FILE_CLOSE(ifile, ierr)
        else
            call s_initialize_mpi_data(q_cons_vf)

            write (file_loc, '(I0,A)') t_step, '.dat'
            file_loc = trim(case_dir) // '/restart_data' // trim(mpiiofs) // trim(file_loc)
            inquire (FILE=trim(file_loc), EXIST=file_exist)
            if (file_exist .and. proc_rank == 0) then
                call MPI_FILE_DELETE(file_loc, mpi_info_int, ierr)
            end if
            call MPI_FILE_OPEN(MPI_COMM_WORLD, file_loc, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpi_info_int, ifile, ierr)

            data_size = (m + 1)*(n + 1)*(p + 1)

            m_MOK = int(m_glb + 1, MPI_OFFSET_KIND)
            n_MOK = int(n_glb + 1, MPI_OFFSET_KIND)
            p_MOK = int(p_glb + 1, MPI_OFFSET_KIND)
            WP_MOK = int(storage_size(0._stp)/8, MPI_OFFSET_KIND)
            MOK = int(1._wp, MPI_OFFSET_KIND)
            str_MOK = int(name_len, MPI_OFFSET_KIND)
            NVARS_MOK = int(alt_sys, MPI_OFFSET_KIND)

            do i = 1, sys_size
                var_MOK = int(i, MPI_OFFSET_KIND)

                disp = m_MOK*max(MOK, n_MOK)*max(MOK, p_MOK)*WP_MOK*(var_MOK - 1)

                call MPI_FILE_SET_VIEW(ifile, disp, mpi_p, MPI_IO_DATA%view(i), 'native', mpi_info_int, ierr)
                call MPI_FILE_WRITE_ALL(ifile, MPI_IO_DATA%var(i)%sf, data_size*mpi_io_type, mpi_io_p, status, ierr)
            end do

            call MPI_FILE_CLOSE(ifile, ierr)
        end if

    end subroutine s_write_parallel_data_files

    !> Write footer with stability criteria extrema and run-time to the information file, then close it
    impure subroutine s_close_run_time_information_file

        real(wp) :: run_time  !< Run-time of the simulation

        write (3, '(A)') '    '
        write (3, '(A)') ''

        write (3, '(A,F9.6)') 'ICFL Max: ', icfl_max
        if (viscous) write (3, '(A,F9.6)') 'VCFL Max: ', vcfl_max
        if (viscous) write (3, '(A,F10.6)') 'Rc Min: ', Rc_min

        call cpu_time(run_time)

        write (3, '(A)') ''
        write (3, '(A,I0,A)') 'Run-time: ', int(anint(run_time)), 's'
        write (3, '(A)') '    '
        close (3)

    end subroutine s_close_run_time_information_file

    !> Initialize the data output module
    impure subroutine s_initialize_data_output_module

        integer :: i, m_ds, n_ds, p_ds

        if (run_time_info) then
            @:ALLOCATE(icfl_sf(0:m, 0:n, 0:p))
            icfl_max = 0._wp

            if (viscous) then
                @:ALLOCATE(vcfl_sf(0:m, 0:n, 0:p))
                @:ALLOCATE(Rc_sf  (0:m, 0:n, 0:p))

                vcfl_max = 0._wp
                Rc_min = 1.e3_wp
            end if
        end if

        if (down_sample) then
            m_ds = int((m + 1)/3) - 1
            n_ds = int((n + 1)/3) - 1
            p_ds = int((p + 1)/3) - 1

            allocate (q_cons_temp_ds(1:sys_size))
            do i = 1, sys_size
                allocate (q_cons_temp_ds(i)%sf(-1:m_ds + 1,-1:n_ds + 1,-1:p_ds + 1))
            end do
        end if

    end subroutine s_initialize_data_output_module

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_data_output_module

        integer :: i

        if (run_time_info) then
            @:DEALLOCATE(icfl_sf)
            if (viscous) then
                @:DEALLOCATE(vcfl_sf, Rc_sf)
            end if
        end if

        if (down_sample) then
            do i = 1, sys_size
                deallocate (q_cons_temp_ds(i)%sf)
            end do
            deallocate (q_cons_temp_ds)
        end if

    end subroutine s_finalize_data_output_module

end module m_data_output
