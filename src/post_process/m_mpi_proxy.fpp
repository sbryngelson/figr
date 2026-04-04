
module m_mpi_proxy

    use mpi  !< Message passing interface (MPI) module

    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use ieee_arithmetic

    implicit none

    !! from all processes to the root process
    integer, allocatable, dimension(:) :: recvcounts
    integer, allocatable, dimension(:) :: displs

contains

    !> Computation of parameters, allocation procedures, and/or any other tasks needed to properly setup the module
    impure subroutine s_initialize_mpi_proxy_module

        integer :: i     !< Generic loop iterator
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        ! Allocating and configuring the receive counts and the displacement vector variables used in variable-gather communication
        ! procedures. Note that these are only needed for either multidimensional runs that utilize the Silo database file format or
        ! for 1D simulations.

        if ((format == 1 .and. n > 0) .or. n == 0) then
            allocate (recvcounts(0:num_procs - 1))
            allocate (displs(0:num_procs - 1))

            if (n == 0) then
                call MPI_GATHER(m + 1, 1, MPI_INTEGER, recvcounts(0), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            else if (proc_rank == 0) then
                recvcounts = 1
            end if

            if (proc_rank == 0) then
                displs(0) = 0

                do i = 1, num_procs - 1
                    displs(i) = displs(i - 1) + recvcounts(i - 1)
                end do
            end if
        end if

    end subroutine s_initialize_mpi_proxy_module

    !> Since only processor with rank 0 is in charge of reading and checking the consistency of the user provided inputs, these are
    !! not available to the remaining processors. This subroutine is then in charge of broadcasting the required information.
    impure subroutine s_mpi_bcast_user_inputs

        integer :: i     !< Generic loop iterator
        integer :: ierr  !< Generic flag used to identify and report MPI errors
        ! Logistics

        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in [ 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',               &
            & 't_step_start', 't_step_stop', 't_step_save',                    &
            & 'model_eqns', 'num_fluids', 'bc_x%beg', 'bc_x%end', 'bc_y%beg',  &
            & 'bc_y%end', 'bc_z%beg', 'bc_z%end', 'flux_lim', 'format',        &
            & 'precision', 'fd_order',                                          &
            & 'n_start' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'parallel_io',                                          &
            & 'rho_wrt', 'E_wrt', 'pres_wrt', 'gamma_wrt', 'sim_data',         &
            & 'heat_ratio_wrt', 'pi_inf_wrt', 'pres_inf_wrt', 'cons_vars_wrt', &
            & 'prim_vars_wrt', 'c_wrt', 'qm_wrt','schlieren_wrt',              &
            & 'file_per_process',                                               &
            & 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt',                          &
            & 'output_partial_domain', 'bc_io',                                 &
            & 'down_sample', 'double_mach']
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        call MPI_BCAST(flux_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(omega_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mom_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(vel_wrt(1), 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_rho_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_rho_e_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(alpha_wrt(1), num_fluids_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

        do i = 1, num_fluids_max
            call MPI_BCAST(fluid_pp(i)%gamma, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%pi_inf, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%cv, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%qv, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(fluid_pp(i)%qvp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        #:for VAR in [ 'pref', 'rhoref',                                        &
            & 't_save', 't_stop',                                               &
            & 'x_output%beg', 'x_output%end', 'y_output%beg', &
            & 'y_output%end', 'z_output%beg', 'z_output%end', 'dt']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        call MPI_BCAST(schlieren_alpha(1), num_fluids_max, mpi_p, 0, MPI_COMM_WORLD, ierr)

    end subroutine s_mpi_bcast_user_inputs

    !> Gather spatial extents from all ranks for Silo database metadata
    impure subroutine s_mpi_gather_spatial_extents(spatial_extents)

        real(wp), dimension(1:,0:), intent(inout) :: spatial_extents

        integer  :: ierr  !< Generic flag used to identify and report MPI errors
        real(wp) :: ext_temp(0:num_procs - 1)

        ! Simulation is 3D

        if (p > 0) then
            ! Minimum spatial extent in the x-direction
            call MPI_GATHERV(minval(x_cb), 1, mpi_p, spatial_extents(1, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Minimum spatial extent in the y-direction
            call MPI_GATHERV(minval(y_cb), 1, mpi_p, spatial_extents(2, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Minimum spatial extent in the z-direction
            call MPI_GATHERV(minval(z_cb), 1, mpi_p, spatial_extents(3, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum spatial extent in the x-direction
            call MPI_GATHERV(maxval(x_cb), 1, mpi_p, spatial_extents(4, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum spatial extent in the y-direction
            call MPI_GATHERV(maxval(y_cb), 1, mpi_p, spatial_extents(5, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum spatial extent in the z-direction
            call MPI_GATHERV(maxval(z_cb), 1, mpi_p, spatial_extents(6, 0), recvcounts, 6*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)
            ! Simulation is 2D
        else if (n > 0) then
            ! Minimum spatial extent in the x-direction
            call MPI_GATHERV(minval(x_cb), 1, mpi_p, spatial_extents(1, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Minimum spatial extent in the y-direction
            call MPI_GATHERV(minval(y_cb), 1, mpi_p, spatial_extents(2, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum spatial extent in the x-direction
            call MPI_GATHERV(maxval(x_cb), 1, mpi_p, spatial_extents(3, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum spatial extent in the y-direction
            call MPI_GATHERV(maxval(y_cb), 1, mpi_p, spatial_extents(4, 0), recvcounts, 4*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)
            ! Simulation is 1D
        else
            ! For 1D, recvcounts/displs are sized for grid defragmentation (m+1 per rank), not for scalar gathers. Use MPI_GATHER
            ! instead.

            ! Minimum spatial extent in the x-direction
            call MPI_GATHER(minval(x_cb), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) spatial_extents(1,:) = ext_temp

            ! Maximum spatial extent in the x-direction
            call MPI_GATHER(maxval(x_cb), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) spatial_extents(2,:) = ext_temp
        end if

    end subroutine s_mpi_gather_spatial_extents

    !> Collect the sub-domain cell-boundary or cell-center location data from all processors and put back together the grid of the
    !! entire computational domain on the rank 0 processor. This is only done for 1D simulations.
    impure subroutine s_mpi_defragment_1d_grid_variable

        integer :: ierr  !< Generic flag used to identify and report MPI errors
        ! Silo-HDF5 database format

        if (format == 1) then
            call MPI_GATHERV(x_cc(0), m + 1, mpi_p, x_root_cc(0), recvcounts, displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Binary database format
        else
            call MPI_GATHERV(x_cb(0), m + 1, mpi_p, x_root_cb(0), recvcounts, displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            if (proc_rank == 0) x_root_cb(-1) = x_cb(-1)
        end if

    end subroutine s_mpi_defragment_1d_grid_variable

    !> Gather the Silo database metadata for the flow variable's extents to boost performance of the multidimensional visualization.
    impure subroutine s_mpi_gather_data_extents(q_sf, data_extents)

        real(wp), dimension(:,:,:), intent(in)                  :: q_sf
        real(wp), dimension(1:2,0:num_procs - 1), intent(inout) :: data_extents

        integer  :: ierr  !< Generic flag used to identify and report MPI errors
        real(wp) :: ext_temp(0:num_procs - 1)

        if (n > 0) then
            ! Multi-D: recvcounts = 1, so strided MPI_GATHERV works correctly Minimum flow variable extent
            call MPI_GATHERV(minval(q_sf), 1, mpi_p, data_extents(1, 0), recvcounts, 2*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

            ! Maximum flow variable extent
            call MPI_GATHERV(maxval(q_sf), 1, mpi_p, data_extents(2, 0), recvcounts, 2*displs, mpi_p, 0, MPI_COMM_WORLD, ierr)
        else
            ! 1D: recvcounts/displs are sized for grid defragmentation (m+1 per rank), not for scalar gathers. Use MPI_GATHER
            ! instead.

            ! Minimum flow variable extent
            call MPI_GATHER(minval(q_sf), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) data_extents(1,:) = ext_temp

            ! Maximum flow variable extent
            call MPI_GATHER(maxval(q_sf), 1, mpi_p, ext_temp, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            if (proc_rank == 0) data_extents(2,:) = ext_temp
        end if

    end subroutine s_mpi_gather_data_extents

    !> Gather the sub-domain flow variable data from all processors and reassemble it for the entire computational domain on the
    !! rank 0 processor. This is only done for 1D simulations.
    impure subroutine s_mpi_defragment_1d_flow_variable(q_sf, q_root_sf)

        real(wp), dimension(0:m), intent(in)    :: q_sf
        real(wp), dimension(0:m), intent(inout) :: q_root_sf

        integer :: ierr  !< Generic flag used to identify and report MPI errors
        ! Gathering the sub-domain flow variable data from all the processes and putting it back together for the entire
        ! computational domain on the process with rank 0

        call MPI_GATHERV(q_sf(0), m + 1, mpi_p, q_root_sf(0), recvcounts, displs, mpi_p, 0, MPI_COMM_WORLD, ierr)

    end subroutine s_mpi_defragment_1d_flow_variable

    !> Deallocation procedures for the module
    impure subroutine s_finalize_mpi_proxy_module

        ! Deallocating the receive counts and the displacement vector variables used in variable-gather communication procedures
        if ((format == 1 .and. n > 0) .or. n == 0) then
            deallocate (recvcounts)
            deallocate (displs)
        end if

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
