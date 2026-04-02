
module m_data_output

    use m_derived_types
    use m_global_parameters
    use m_derived_variables
    use m_mpi_proxy
    use m_compile_specific
    use m_helper
    use m_variables_conversion

    implicit none

    private; public :: s_initialize_data_output_module, s_define_output_region, s_open_formatted_database_file, &
        & s_open_intf_data_file, s_open_energy_data_file, s_write_grid_to_formatted_database_file, &
        & s_write_variable_to_formatted_database_file, s_write_intf_data_file, s_write_energy_data_file, &
        & s_close_formatted_database_file, s_close_intf_data_file, s_close_energy_data_file, s_finalize_data_output_module

    ! Include Silo-HDF5 interface library
    include 'silo_f9x.inc'

    ! Flow variable storage; q_root_sf gathers to rank 0 in 1D parallel runs
    real(wp), allocatable, dimension(:,:,:), public :: q_sf
    real(wp), allocatable, dimension(:,:,:)         :: q_root_sf

    ! Single precision storage for flow variables
    real(sp), allocatable, dimension(:,:,:), public :: q_sf_s
    real(sp), allocatable, dimension(:,:,:)         :: q_root_sf_s

    ! Spatial and data extents for VisIt visualization
    real(wp), allocatable, dimension(:,:) :: spatial_extents
    real(wp), allocatable, dimension(:,:) :: data_extents

    ! Ghost zone layer sizes (lo/hi) for subdomain connectivity in VisIt
    integer, allocatable, dimension(:) :: lo_offset
    integer, allocatable, dimension(:) :: hi_offset

    ! Track cell-boundary count per active coordinate direction
    integer, allocatable, dimension(:) :: dims

    ! Locations of various folders in the case's directory tree, associated with the choice of the formatted database format. These
    ! include, in order, the location of the folder named after the selected formatted database format, and the locations of two
    ! sub-directories of the latter, the first of which is named after the local processor rank, while the second is named 'root'.
    ! The folder associated with the local processor rank contains only the data pertaining to the part of the domain taken care of
    ! by the local processor. The root directory, on the other hand, will contain either the information about the connectivity
    ! required to put the entire domain back together, or the actual data associated with the entire computational domain. This all
    ! depends on dimensionality and the choice of the formatted database format.
    character(LEN=path_len + name_len)   :: dbdir
    character(LEN=path_len + 2*name_len) :: proc_rank_dir
    character(LEN=path_len + 2*name_len) :: rootdir

    ! Handles of the formatted database master/root file, slave/local processor file and options list. The list of options is
    ! explicitly used in the Silo- HDF5 database format to provide additional details about the contents of a formatted database
    ! file, such as the previously described spatial and data extents.
    integer :: dbroot
    integer :: dbfile
    integer :: optlist

    ! The total number of flow variable(s) to be stored in a formatted database file. Note that this is only needed when using the
    ! Binary format.
    integer :: dbvars

    ! Generic error flags utilized in the handling, checking and the reporting of the input and output operations errors with a
    ! formatted database file
    integer, private :: err

contains

    !> Allocate storage arrays, configure output directories, and count flow variables for formatted database output.
    impure subroutine s_initialize_data_output_module()

        character(LEN=len_trim(case_dir) + 2*name_len) :: file_loc
        logical                                        :: dir_check
        integer                                        :: i

        allocate (q_sf(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end))

        if (precision == 1) then
            allocate (q_sf_s(-offset_x%beg:m + offset_x%end,-offset_y%beg:n + offset_y%end,-offset_z%beg:p + offset_z%end))
        end if

        if (n == 0) then
            allocate (q_root_sf(0:m_root,0:0,0:0))
            if (precision == 1) then
                allocate (q_root_sf_s(0:m_root,0:0,0:0))
            end if
        end if

        ! Allocating the spatial and data extents and also the variables for the offsets and the one bookkeeping the number of
        ! cell-boundaries in each active coordinate direction. Note that all these variables are only needed by the Silo-HDF5 format
        ! for multidimensional data.
        if (format == 1) then
            allocate (data_extents(1:2,0:num_procs - 1))

            if (p > 0) then
                allocate (spatial_extents(1:6,0:num_procs - 1))
                allocate (lo_offset(1:3))
                allocate (hi_offset(1:3))
                allocate (dims(1:3))
            else if (n > 0) then
                allocate (spatial_extents(1:4,0:num_procs - 1))
                allocate (lo_offset(1:2))
                allocate (hi_offset(1:2))
                allocate (dims(1:2))
            else
                allocate (spatial_extents(1:2,0:num_procs - 1))
                allocate (lo_offset(1:1))
                allocate (hi_offset(1:1))
                allocate (dims(1:1))
            end if
        end if

        ! The size of the ghost zone layer in each of the active coordinate directions was set in the module m_mpi_proxy.f90. The
        ! results are now transferred to the local variables of this module when they are required by the Silo-HDF5 format, for
        ! multidimensional data sets. With the same, latter, requirements, the variables bookkeeping the number of cell-boundaries
        ! in each active coordinate direction are also set here.
        if (format == 1) then
            if (p > 0) then
                lo_offset(:) = (/offset_x%beg, offset_y%beg, offset_z%beg/)
                hi_offset(:) = (/offset_x%end, offset_y%end, offset_z%end/)

                dims(:) = (/m + offset_x%beg + offset_x%end + 2, n + offset_y%beg + offset_y%end + 2, &
                     & p + offset_z%beg + offset_z%end + 2/)
            else if (n > 0) then
                lo_offset(:) = (/offset_x%beg, offset_y%beg/)
                hi_offset(:) = (/offset_x%end, offset_y%end/)

                dims(:) = (/m + offset_x%beg + offset_x%end + 2, n + offset_y%beg + offset_y%end + 2/)
            else
                lo_offset(:) = (/offset_x%beg/)
                hi_offset(:) = (/offset_x%end/)
                dims(:) = (/m + offset_x%beg + offset_x%end + 2/)
            end if
        end if

        if (format == 1) then
            dbdir = trim(case_dir) // '/silo_hdf5'

            write (proc_rank_dir, '(A,I0)') '/p', proc_rank

            proc_rank_dir = trim(dbdir) // trim(proc_rank_dir)

            file_loc = trim(proc_rank_dir) // '/.'

            call my_inquire(file_loc, dir_check)
            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(proc_rank_dir))
            end if

            if (proc_rank == 0) then
                rootdir = trim(dbdir) // '/root'

                file_loc = trim(rootdir) // '/.'

                call my_inquire(file_loc, dir_check)
                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(rootdir))
                end if
            end if
        else
            dbdir = trim(case_dir) // '/binary'

            write (proc_rank_dir, '(A,I0)') '/p', proc_rank

            proc_rank_dir = trim(dbdir) // trim(proc_rank_dir)

            file_loc = trim(proc_rank_dir) // '/.'

            call my_inquire(file_loc, dir_check)

            if (dir_check .neqv. .true.) then
                call s_create_directory(trim(proc_rank_dir))
            end if

            if (n == 0 .and. proc_rank == 0) then
                rootdir = trim(dbdir) // '/root'

                file_loc = trim(rootdir) // '/.'

                call my_inquire(file_loc, dir_check)

                if (dir_check .neqv. .true.) then
                    call s_create_directory(trim(rootdir))
                end if
            end if
        end if

        ! Contrary to the Silo-HDF5 database format, handles of the Binary database master/root and slave/local process files are
        ! perfectly static throughout post-process. Hence, they are set here so that they do not have to be repetitively computed in
        ! later procedures.
        if (format == 2) then
            if (n == 0 .and. proc_rank == 0) dbroot = 2
            dbfile = 1
        end if

        if (format == 2) then
            dbvars = 0

            do i = 1, num_fluids
                if (alpha_rho_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    dbvars = dbvars + 1
                end if
            end do

            if (rho_wrt) then
                dbvars = dbvars + 1
            end if

            do i = 1, E_idx - mom_idx%beg
                if (mom_wrt(i) .or. cons_vars_wrt) dbvars = dbvars + 1
            end do

            do i = 1, E_idx - mom_idx%beg
                if (vel_wrt(i) .or. prim_vars_wrt) dbvars = dbvars + 1
            end do

            do i = 1, E_idx - mom_idx%beg
                if (flux_wrt(i)) dbvars = dbvars + 1
            end do

            if (E_wrt .or. cons_vars_wrt) dbvars = dbvars + 1
            if (pres_wrt .or. prim_vars_wrt) dbvars = dbvars + 1

            if (model_eqns == 2) then
                do i = 1, num_fluids - 1
                    if (alpha_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                        dbvars = dbvars + 1
                    end if
                end do

                if (alpha_wrt(num_fluids) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    dbvars = dbvars + 1
                end if
            end if

            if (gamma_wrt) then
                dbvars = dbvars + 1
            end if

            if (heat_ratio_wrt) dbvars = dbvars + 1

            if (pi_inf_wrt) then
                dbvars = dbvars + 1
            end if

            if (pres_inf_wrt) dbvars = dbvars + 1
            if (c_wrt) dbvars = dbvars + 1

            if (p > 0) then
                do i = 1, num_vels
                    if (omega_wrt(i)) dbvars = dbvars + 1
                end do
            else if (n > 0) then
                do i = 1, num_vels
                    if (omega_wrt(i)) dbvars = dbvars + 1
                end do
            end if

            if (schlieren_wrt) dbvars = dbvars + 1
        end if

    end subroutine s_initialize_data_output_module

    !> Compute the cell-index bounds for the user-specified partial output domain in each coordinate direction.
    impure subroutine s_define_output_region

        integer :: i
        integer :: lower_bound, upper_bound

        #:for X, M in [('x', 'm'), ('y', 'n'), ('z', 'p')]
            if (${M}$ == 0) return  ! Early return for y or z if simulation is 1D or 2D

            lower_bound = -offset_${X}$%beg
            upper_bound = ${M}$ + offset_${X}$%end

            do i = lower_bound, upper_bound
                if (${X}$_cc(i) > ${X}$_output%beg) then
                    ${X}$_output_idx%beg = i + offset_${X}$%beg
                    exit
                end if
            end do

            do i = upper_bound, lower_bound, -1
                if (${X}$_cc(i) < ${X}$_output%end) then
                    ${X}$_output_idx%end = i + offset_${X}$%beg
                    exit
                end if
            end do

            ! If no grid points are within the output region
            if ((${X}$_cc(lower_bound) > ${X}$_output%end) .or. (${X}$_cc(upper_bound) < ${X}$_output%beg)) then
                ${X}$_output_idx%beg = 0
                ${X}$_output_idx%end = 0
            end if
        #:endfor

    end subroutine s_define_output_region

    !> Open (or create) the Silo-HDF5 or Binary formatted database slave and master files for a given time step.
    impure subroutine s_open_formatted_database_file(t_step)

        integer, intent(in)                            :: t_step
        character(LEN=len_trim(case_dir) + 3*name_len) :: file_loc
        integer                                        :: ierr

        if (format == 1) then
            write (file_loc, '(A,I0,A)') '/', t_step, '.silo'
            file_loc = trim(proc_rank_dir) // trim(file_loc)

            ierr = DBCREATE(trim(file_loc), len_trim(file_loc), DB_CLOBBER, DB_LOCAL, 'figr', 8, DB_HDF5, dbfile)

            if (dbfile == -1) then
                call s_mpi_abort('Unable to create Silo-HDF5 database ' // 'slave file ' // trim(file_loc) // '. ' // 'Exiting.')
            end if

            if (proc_rank == 0) then
                write (file_loc, '(A,I0,A)') '/collection_', t_step, '.silo'
                file_loc = trim(rootdir) // trim(file_loc)

                ierr = DBCREATE(trim(file_loc), len_trim(file_loc), DB_CLOBBER, DB_LOCAL, 'figr', 8, DB_HDF5, dbroot)

                if (dbroot == -1) then
                    call s_mpi_abort('Unable to create Silo-HDF5 database ' // 'master file ' // trim(file_loc) // '. ' &
                                     & // 'Exiting.')
                end if
            end if
        else
            write (file_loc, '(A,I0,A)') '/', t_step, '.dat'
            file_loc = trim(proc_rank_dir) // trim(file_loc)

            open (dbfile, IOSTAT=err, FILE=trim(file_loc), form='unformatted', STATUS='replace')

            if (err /= 0) then
                call s_mpi_abort('Unable to create Binary database slave ' // 'file ' // trim(file_loc) // '. Exiting.')
            end if

            if (output_partial_domain) then
                write (dbfile) x_output_idx%end - x_output_idx%beg, y_output_idx%end - y_output_idx%beg, &
                       & z_output_idx%end - z_output_idx%beg, dbvars
            else
                write (dbfile) m, n, p, dbvars
            end if

            if (n == 0 .and. proc_rank == 0) then
                write (file_loc, '(A,I0,A)') '/', t_step, '.dat'
                file_loc = trim(rootdir) // trim(file_loc)

                open (dbroot, IOSTAT=err, FILE=trim(file_loc), form='unformatted', STATUS='replace')

                if (err /= 0) then
                    call s_mpi_abort('Unable to create Binary database ' // 'master file ' // trim(file_loc) // '. Exiting.')
                end if

                if (output_partial_domain) then
                    write (dbroot) x_output_idx%end - x_output_idx%beg, 0, 0, dbvars
                else
                    write (dbroot) m_root, 0, 0, dbvars
                end if
            end if
        end if

    end subroutine s_open_formatted_database_file

    !> Open the interface data file for appending extracted interface coordinates.
    impure subroutine s_open_intf_data_file()

        character(LEN=path_len + 3*name_len) :: file_path

        write (file_path, '(A)') '/intf_data.dat'
        file_path = trim(case_dir) // trim(file_path)

        open (211, FILE=trim(file_path), form='formatted', POSITION='append', STATUS='unknown')

    end subroutine s_open_intf_data_file

    !> Open the energy data file for appending volume-integrated energy budget quantities.
    impure subroutine s_open_energy_data_file()

        character(LEN=path_len + 3*name_len) :: file_path

        write (file_path, '(A)') '/eng_data.dat'
        file_path = trim(case_dir) // trim(file_path)

        open (251, FILE=trim(file_path), form='formatted', POSITION='append', STATUS='unknown')

    end subroutine s_open_energy_data_file

    !> Write the computational grid (cell-boundary coordinates) to the formatted database slave and master files.
    impure subroutine s_write_grid_to_formatted_database_file(t_step)

        integer, intent(in) :: t_step

        ! NAG compiler requires these to be statically sized
        character(LEN=4*name_len), dimension(num_procs) :: meshnames
        integer, dimension(num_procs)                   :: meshtypes
        integer                                         :: i
        integer                                         :: ierr

        if (format == 1) then
            ! For multidimensional data sets, the spatial extents of all of the grid(s) handled by the local processor(s) are
            ! recorded so that they may be written, by root processor, to the formatted database master file.
            if (num_procs > 1) then
                call s_mpi_gather_spatial_extents(spatial_extents)
            else if (p > 0) then
                spatial_extents(:,0) = (/minval(x_cb), minval(y_cb), minval(z_cb), maxval(x_cb), maxval(y_cb), maxval(z_cb)/)
            else if (n > 0) then
                spatial_extents(:,0) = (/minval(x_cb), minval(y_cb), maxval(x_cb), maxval(y_cb)/)
            else
                spatial_extents(:,0) = (/minval(x_cb), maxval(x_cb)/)
            end if

            ! Next, the root processor proceeds to record all of the spatial extents in the formatted database master file. In
            ! addition, it also records a sub-domain connectivity map so that the entire grid may be reassembled by looking at the
            ! master file.
            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (meshnames(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:rectilinear_grid'
                end do

                meshtypes = DB_QUAD_RECT

                err = DBSET2DSTRLEN(len(meshnames(1)))
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIOPT(optlist, DBOPT_EXTENTS_SIZE, size(spatial_extents, 1))
                err = DBADDDOPT(optlist, DBOPT_EXTENTS, spatial_extents)
                err = DBPUTMMESH(dbroot, 'rectilinear_grid', 16, num_procs, meshnames, len_trim(meshnames), meshtypes, optlist, &
                                 & ierr)
                err = DBFREEOPTLIST(optlist)
            end if

            ! Finally, the local quadrilateral mesh, either 2D or 3D, along with its offsets that indicate the presence and size of
            ! ghost zone layer(s), are put in the formatted database slave file.

            if (p > 0) then
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIAOPT(optlist, DBOPT_LO_OFFSET, size(lo_offset), lo_offset)
                err = DBADDIAOPT(optlist, DBOPT_HI_OFFSET, size(hi_offset), hi_offset)
                err = DBPUTQM(dbfile, 'rectilinear_grid', 16, 'x', 1, 'y', 1, 'z', 1, x_cb, y_cb, z_cb, dims, 3, DB_DOUBLE, &
                              & DB_COLLINEAR, optlist, ierr)
                err = DBFREEOPTLIST(optlist)
            else if (n > 0) then
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIAOPT(optlist, DBOPT_LO_OFFSET, size(lo_offset), lo_offset)
                err = DBADDIAOPT(optlist, DBOPT_HI_OFFSET, size(hi_offset), hi_offset)
                err = DBPUTQM(dbfile, 'rectilinear_grid', 16, 'x', 1, 'y', 1, 'z', 1, x_cb, y_cb, DB_F77NULL, dims, 2, DB_DOUBLE, &
                              & DB_COLLINEAR, optlist, ierr)
                err = DBFREEOPTLIST(optlist)
            else
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIAOPT(optlist, DBOPT_LO_OFFSET, size(lo_offset), lo_offset)
                err = DBADDIAOPT(optlist, DBOPT_HI_OFFSET, size(hi_offset), hi_offset)
                err = DBPUTQM(dbfile, 'rectilinear_grid', 16, 'x', 1, 'y', 1, 'z', 1, x_cb, DB_F77NULL, DB_F77NULL, dims, 1, &
                              & DB_DOUBLE, DB_COLLINEAR, optlist, ierr)
                err = DBFREEOPTLIST(optlist)
            end if
        else if (format == 2) then
            ! Multidimensional local grid data is written to the formatted database slave file. Recall that no master file to
            ! maintained in multidimensions.
            if (p > 0) then
                if (precision == 1) then
                    write (dbfile) real(x_cb, sp), real(y_cb, sp), real(z_cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x_cb(x_output_idx%beg - 1:x_output_idx%end), y_cb(y_output_idx%beg - 1:y_output_idx%end), &
                               & z_cb(z_output_idx%beg - 1:z_output_idx%end)
                    else
                        write (dbfile) x_cb, y_cb, z_cb
                    end if
                end if
            else if (n > 0) then
                if (precision == 1) then
                    write (dbfile) real(x_cb, sp), real(y_cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x_cb(x_output_idx%beg - 1:x_output_idx%end), y_cb(y_output_idx%beg - 1:y_output_idx%end)
                    else
                        write (dbfile) x_cb, y_cb
                    end if
                end if

                ! One-dimensional local grid data is written to the formatted database slave file. In addition, the local grid data
                ! is put together by the root process and written to the master file.
            else
                if (precision == 1) then
                    write (dbfile) real(x_cb, sp)
                else
                    if (output_partial_domain) then
                        write (dbfile) x_cb(x_output_idx%beg - 1:x_output_idx%end)
                    else
                        write (dbfile) x_cb
                    end if
                end if

                if (num_procs > 1) then
                    call s_mpi_defragment_1d_grid_variable()
                else
                    x_root_cb(:) = x_cb(:)
                end if

                if (proc_rank == 0) then
                    if (precision == 1) then
                        write (dbroot) real(x_root_cb, wp)
                    else
                        if (output_partial_domain) then
                            write (dbroot) x_root_cb(x_output_idx%beg - 1:x_output_idx%end)
                        else
                            write (dbroot) x_root_cb
                        end if
                    end if
                end if
            end if
        end if

    end subroutine s_write_grid_to_formatted_database_file

    !> Write a single flow variable field to the formatted database slave and master files for a given time step.
    impure subroutine s_write_variable_to_formatted_database_file(varname, t_step)

        character(LEN=*), intent(in) :: varname
        integer, intent(in)          :: t_step

        ! NAG compiler requires these to be statically sized
        character(LEN=4*name_len), dimension(num_procs) :: varnames
        integer, dimension(num_procs)                   :: vartypes
        integer                                         :: i, j, k
        integer                                         :: ierr

        if (format == 1) then
            ! Determining the extents of the flow variable on each local process and gathering all this information on root process
            if (num_procs > 1) then
                call s_mpi_gather_data_extents(q_sf, data_extents)
            else
                data_extents(:,0) = (/minval(q_sf), maxval(q_sf)/)
            end if

            if (proc_rank == 0) then
                do i = 1, num_procs
                    write (varnames(i), '(A,I0,A,I0,A)') '../p', i - 1, '/', t_step, '.silo:' // trim(varname)
                end do

                vartypes = DB_QUADVAR

                err = DBSET2DSTRLEN(len(varnames(1)))
                err = DBMKOPTLIST(2, optlist)
                err = DBADDIOPT(optlist, DBOPT_EXTENTS_SIZE, 2)
                err = DBADDDOPT(optlist, DBOPT_EXTENTS, data_extents)
                err = DBPUTMVAR(dbroot, trim(varname), len_trim(varname), num_procs, varnames, len_trim(varnames), vartypes, &
                                & optlist, ierr)
                err = DBFREEOPTLIST(optlist)
            end if

            if (wp == dp) then
                if (precision == 1) then
                    do i = -offset_x%beg, m + offset_x%end
                        do j = -offset_y%beg, n + offset_y%end
                            do k = -offset_z%beg, p + offset_z%end
                                q_sf_s(i, j, k) = real(q_sf(i, j, k), sp)
                            end do
                        end do
                    end do
                end if
            else if (wp == sp) then
                do i = -offset_x%beg, m + offset_x%end
                    do j = -offset_y%beg, n + offset_y%end
                        do k = -offset_z%beg, p + offset_z%end
                            q_sf_s(i, j, k) = q_sf(i, j, k)
                        end do
                    end do
                end do
            end if

            #:for PRECISION, SFX, DBT in [(1,'_s','DB_FLOAT'),(2,'',"DB_DOUBLE")]
                if (precision == ${PRECISION}$) then
                    if (p > 0) then
                        err = DBPUTQV1(dbfile, trim(varname), len_trim(varname), 'rectilinear_grid', 16, q_sf${SFX}$, dims - 1, &
                                       & 3, DB_F77NULL, 0, ${DBT}$, DB_ZONECENT, DB_F77NULL, ierr)
                    else if (n > 0) then
                        err = DBPUTQV1(dbfile, trim(varname), len_trim(varname), 'rectilinear_grid', 16, q_sf${SFX}$, dims - 1, &
                                       & 2, DB_F77NULL, 0, ${DBT}$, DB_ZONECENT, DB_F77NULL, ierr)
                    else
                        err = DBPUTQV1(dbfile, trim(varname), len_trim(varname), 'rectilinear_grid', 16, q_sf${SFX}$, dims - 1, &
                                       & 1, DB_F77NULL, 0, ${DBT}$, DB_ZONECENT, DB_F77NULL, ierr)
                    end if
                end if
            #:endfor
        else
            ! Writing the name of the flow variable and its data, associated with the local processor, to the formatted database
            ! slave file
            if (precision == 1) then
                write (dbfile) varname, real(q_sf, wp)
            else
                write (dbfile) varname, q_sf
            end if

            ! In 1D, the root process also takes care of gathering the flow variable data from all of the local processor(s) and
            ! writes it to the formatted database master file.
            if (n == 0) then
                if (num_procs > 1) then
                    call s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)
                else
                    q_root_sf(:,:,:) = q_sf(:,:,:)
                end if

                if (proc_rank == 0) then
                    if (precision == 1) then
                        write (dbroot) varname, real(q_root_sf, wp)
                    else
                        write (dbroot) varname, q_root_sf
                    end if
                end if
            end if
        end if

    end subroutine s_write_variable_to_formatted_database_file

    !> Extract the volume-fraction interface contour from primitive fields and write the coordinates to the interface data file.
    impure subroutine s_write_intf_data_file(q_prim_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf
        integer :: i, j, k, l, cent
        integer :: counter, root  !< number of data points extracted to fit shape to SH perturbations
        real(wp), allocatable :: x_td(:), y_td(:), x_d1(:), y_d1(:), y_d(:), x_d(:)
        real(wp) :: axp, axm, ayp, aym, tgp, euc_d, thres, maxalph_loc, maxalph_glb

        allocate (x_d1(m*n))
        allocate (y_d1(m*n))
        counter = 0
        maxalph_loc = 0._wp
        do k = 0, p
            do j = 0, n
                do i = 0, m
                    if (q_prim_vf(E_idx + 2)%sf(i, j, k) > maxalph_loc) then
                        maxalph_loc = q_prim_vf(E_idx + 2)%sf(i, j, k)
                    end if
                end do
            end do
        end do

        call s_mpi_allreduce_max(maxalph_loc, maxalph_glb)
        if (p > 0) then
            do l = 0, p
                if (z_cc(l) < dz(l) .and. z_cc(l) > 0) then
                    cent = l
                end if
            end do
        else
            cent = 0
        end if

        thres = 0.9_wp*maxalph_glb
        do k = 0, n
            do j = 0, m
                axp = q_prim_vf(E_idx + 2)%sf(j + 1, k, cent)
                axm = q_prim_vf(E_idx + 2)%sf(j, k, cent)
                ayp = q_prim_vf(E_idx + 2)%sf(j, k + 1, cent)
                aym = q_prim_vf(E_idx + 2)%sf(j, k, cent)
                if ((axp > thres .and. axm < thres) .or. (axp < thres .and. axm > thres) .or. (ayp > thres .and. aym < thres) &
                    & .or. (ayp < thres .and. aym > thres)) then
                    if (counter == 0) then
                        counter = counter + 1
                        x_d1(counter) = x_cc(j)
                        y_d1(counter) = y_cc(k)
                    else
                        tgp = sqrt(dx(j)**2 + dy(k)**2)
                        do i = 1, counter
                            euc_d = sqrt((x_cc(j) - x_d1(i))**2 + (y_cc(k) - y_d1(i))**2)
                            if (euc_d < tgp) then
                                exit
                            else if (i == counter) then
                                counter = counter + 1
                                x_d1(counter) = x_cc(j)
                                y_d1(counter) = y_cc(k)
                            end if
                        end do
                    end if
                end if
            end do
        end do

        allocate (x_d(counter), y_d(counter))

        do i = 1, counter
            y_d(i) = y_d1(i)
            x_d(i) = x_d1(i)
        end do
        root = 0

        call s_mpi_gather_data(x_d, counter, x_td, root)
        call s_mpi_gather_data(y_d, counter, y_td, root)
        if (proc_rank == 0) then
            do i = 1, size(x_td)
                if (i == size(x_td)) then
                    write (211, '(F12.9,1X,F12.9,1X,I4)') x_td(i), y_td(i), size(x_td)
                else
                    write (211, '(F12.9,1X,F12.9,1X,F3.1)') x_td(i), y_td(i), 0._wp
                end if
            end do
        end if

    end subroutine s_write_intf_data_file

    !> Compute volume-integrated kinetic, potential, and internal energies and write the energy budget to the energy data file.
    impure subroutine s_write_energy_data_file(q_prim_vf, q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in) :: q_prim_vf, q_cons_vf
        real(wp) :: Elk, Egk, Elp, Egint, Vb, Vl, pres_av, Et
        real(wp) :: rho, pres, dV, tmp, gamma, pi_inf, MaxMa, MaxMa_glb, maxvel, c, Ma, H, qv
        real(wp), dimension(num_vels) :: vel
        real(wp), dimension(num_fluids) :: adv
        integer :: i, j, k, l, s  !< looping indices

        Egk = 0._wp
        Elp = 0._wp
        Egint = 0._wp
        Vb = 0._wp
        maxvel = 0._wp
        MaxMa = 0._wp
        Vl = 0._wp
        Elk = 0._wp
        Et = 0._wp
        Vb = 0._wp
        dV = 0._wp
        pres_av = 0._wp
        pres = 0._wp
        c = 0._wp

        do k = 0, p
            do j = 0, n
                do i = 0, m
                    pres = 0._wp
                    dV = dx(i)*dy(j)*dz(k)
                    rho = 0._wp
                    gamma = 0._wp
                    pi_inf = 0._wp
                    qv = 0._wp
                    pres = q_prim_vf(E_idx)%sf(i, j, k)
                    Egint = Egint + q_prim_vf(E_idx + 2)%sf(i, j, k)*(gammas(2)*pres)*dV
                    do s = 1, num_vels
                        vel(s) = q_prim_vf(num_fluids + s)%sf(i, j, k)
                        Egk = Egk + 0.5_wp*q_prim_vf(E_idx + 2)%sf(i, j, k)*q_prim_vf(2)%sf(i, j, k)*vel(s)*vel(s)*dV
                        Elk = Elk + 0.5_wp*q_prim_vf(E_idx + 1)%sf(i, j, k)*q_prim_vf(1)%sf(i, j, k)*vel(s)*vel(s)*dV
                        if (abs(vel(s)) > maxvel) then
                            maxvel = abs(vel(s))
                        end if
                    end do
                    do l = 1, adv_idx%end - E_idx
                        adv(l) = q_prim_vf(E_idx + l)%sf(i, j, k)
                        gamma = gamma + adv(l)*gammas(l)
                        pi_inf = pi_inf + adv(l)*pi_infs(l)
                        rho = rho + adv(l)*q_prim_vf(l)%sf(i, j, k)
                        qv = qv + adv(l)*q_prim_vf(l)%sf(i, j, k)*qvs(l)
                    end do

                    H = ((gamma + 1._wp)*pres + pi_inf + qv)/rho

                    call s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, 0._wp, 0._wp, c, qv)

                    Ma = maxvel/c
                    if (Ma > MaxMa .and. (adv(1) > (1.0_wp - 1.0e-10_wp))) then
                        MaxMa = Ma
                    end if
                    Vl = Vl + adv(1)*dV
                    Vb = Vb + adv(2)*dV
                    pres_av = pres_av + adv(1)*pres*dV
                    Et = Et + q_cons_vf(E_idx)%sf(i, j, k)*dV
                end do
            end do
        end do

        tmp = pres_av
        call s_mpi_allreduce_sum(tmp, pres_av)
        tmp = Vl
        call s_mpi_allreduce_sum(tmp, Vl)

        call s_mpi_allreduce_max(MaxMa, MaxMa_glb)
        tmp = Elk
        call s_mpi_allreduce_sum(tmp, Elk)
        tmp = Egint
        call s_mpi_allreduce_sum(tmp, Egint)
        tmp = Egk
        call s_mpi_allreduce_sum(tmp, Egk)
        tmp = Vb
        call s_mpi_allreduce_sum(tmp, Vb)
        tmp = Et
        call s_mpi_allreduce_sum(tmp, Et)

        Elp = pres_av/Vl*Vb
        if (proc_rank == 0) then
            write (251, '(10X, 8F24.8)') Elp, Egint, Elk, Egk, Et, Vb, Vl, MaxMa_glb
        end if

    end subroutine s_write_energy_data_file

    !> Close the formatted database slave file and, for the root process, the master file.
    impure subroutine s_close_formatted_database_file()

        integer :: ierr

        if (format == 1) then
            ierr = DBCLOSE(dbfile)
            if (proc_rank == 0) ierr = DBCLOSE(dbroot)
        else
            close (dbfile)
            if (n == 0 .and. proc_rank == 0) close (dbroot)
        end if

    end subroutine s_close_formatted_database_file

    !> Close the interface data file.
    impure subroutine s_close_intf_data_file()

        close (211)

    end subroutine s_close_intf_data_file

    !> Close the energy data file.
    impure subroutine s_close_energy_data_file()

        close (251)

    end subroutine s_close_energy_data_file

    !> Deallocate module arrays and release all data-output resources.
    impure subroutine s_finalize_data_output_module()

        deallocate (q_sf)
        if (n == 0) deallocate (q_root_sf)

        ! Deallocating spatial and data extents and also the variables for the offsets and the one bookkeeping the number of
        ! cell-boundaries in each active coordinate direction. Note that all these variables were only needed by Silo-HDF5 format
        ! for multidimensional data.
        if (format == 1) then
            deallocate (spatial_extents)
            deallocate (data_extents)
            deallocate (lo_offset)
            deallocate (hi_offset)
            deallocate (dims)
        end if

    end subroutine s_finalize_data_output_module

end module m_data_output
