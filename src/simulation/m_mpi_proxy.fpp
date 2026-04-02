!>
!! @file
!! @brief Contains module m_mpi_proxy

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief MPI halo exchange, domain decomposition, and buffer packing/unpacking for the simulation solver
module m_mpi_proxy

    use mpi  !< Message passing interface (MPI) module

    use m_helper_basic
    use m_helper
    use m_derived_types
    use m_global_parameters
    use m_mpi_common
    use m_nvtx
    use ieee_arithmetic

    implicit none

    integer, private, allocatable, dimension(:) :: ib_buff_send  !< IB marker send buffer for halo exchange
    integer, private, allocatable, dimension(:) :: ib_buff_recv  !< IB marker receive buffer for halo exchange
    integer                                     :: i_halo_size
    $:GPU_DECLARE(create='[i_halo_size]')

contains

    !> Initialize the MPI proxy module
    subroutine s_initialize_mpi_proxy_module()

    end subroutine s_initialize_mpi_proxy_module

    !> Since only the processor with rank 0 reads and verifies the consistency of user inputs, these are initially not available to
    !! the other processors. Then, the purpose of this subroutine is to distribute the user inputs to the remaining processors in
    !! the communicator.
    impure subroutine s_mpi_bcast_user_inputs()

        integer :: i, j  !< Generic loop iterator
        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['t_step_old', 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',  &
            & 't_step_start','t_step_stop','t_step_save','t_step_print',       &
            & 'model_eqns','time_stepper',                                     &
            & 'precision', 'bc_x%beg', 'bc_x%end',               &
            & 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end',     &
            & 'n_start',    &
            & 'num_bc_patches', 'num_igr_iters', 'num_igr_warm_start_iters' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'run_time_info',     &
            & 'mp_weno', 'rdma_mpi', 'bc_io', &
            & 'mixture_err', 'parallel_io', &
            & 'prim_vars_wrt', 'weno_avg', 'file_per_process',   &
            & 'cfl_adap_dt', 'cfl_const_dt', 'cfl_dt',       &
            & 'down_sample']
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'dt','weno_eps','teno_CT','pref','rhoref', &
            & 'bc_x%vb1','bc_x%vb2','bc_x%vb3','bc_x%ve1','bc_x%ve2','bc_x%ve3', &
            & 'bc_y%vb1','bc_y%vb2','bc_y%vb3','bc_y%ve1','bc_y%ve2','bc_y%ve3', &
            & 'bc_z%vb1','bc_z%vb2','bc_z%vb3','bc_z%ve1','bc_z%ve2','bc_z%ve3', &
            & 'x_domain%beg', 'x_domain%end', 'y_domain%beg', 'y_domain%end',    &
            & 'z_domain%beg', 'z_domain%end', 'x_a', 'x_b', 'y_a', 'y_b', 'z_a', &
            & 'z_b', 't_stop', 't_save', 'cfl_target', 'alf_factor' ]
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:if not MFC_CASE_OPTIMIZATION
            call MPI_BCAST(mapped_weno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(teno, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(weno_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(num_fluids, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(wenoz_q, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_pres_lim, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(igr_iter_solver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(viscous, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(recon_type, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(muscl_order, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(muscl_lim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endif

        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf','cv','qv','qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
            call MPI_BCAST(fluid_pp(i)%Re(1), 2, mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        do i = 1, num_fluids_max
            #:for VAR in ['bc_x%alpha_rho_in','bc_x%alpha_in','bc_y%alpha_rho_in','bc_y%alpha_in','bc_z%alpha_rho_in', &
                & 'bc_z%alpha_in']
                call MPI_BCAST(${VAR}$ (i), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        ! NVIDIA UVM variables
        call MPI_BCAST(nv_uvm_out_of_core, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nv_uvm_igr_temps_on_gpu, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nv_uvm_pref_gpu, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    end subroutine s_mpi_bcast_user_inputs

    !> Broadcast random phase numbers from rank 0 to all MPI processes
    impure subroutine s_mpi_send_random_number(phi_rn, num_freq)

        integer, intent(in)                            :: num_freq
        real(wp), intent(inout), dimension(1:num_freq) :: phi_rn

        integer :: ierr  !< Generic flag used to identify and report MPI errors

        call MPI_BCAST(phi_rn, num_freq, mpi_p, 0, MPI_COMM_WORLD, ierr)

    end subroutine s_mpi_send_random_number

    !> Finalize the MPI proxy module
    subroutine s_finalize_mpi_proxy_module()

    end subroutine s_finalize_mpi_proxy_module

end module m_mpi_proxy
