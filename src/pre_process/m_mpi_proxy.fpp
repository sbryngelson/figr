
module m_mpi_proxy

    use mpi

    use m_helper
    use m_derived_types
    use m_global_parameters
    use m_mpi_common

    implicit none

contains
    !> Since only processor with rank 0 is in charge of reading and checking the consistency of the user provided inputs, these are
    !! not available to the remaining processors. This subroutine is then in charge of broadcasting the required information.
    impure subroutine s_mpi_bcast_user_inputs

        integer :: i, j
        integer :: ierr

        call MPI_BCAST(case_dir, len(case_dir), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

        #:for VAR in ['t_step_old', 't_step_start', 'm', 'n', 'p', 'm_glb', 'n_glb', 'p_glb',  &
            & 'loops_x', 'loops_y', 'loops_z', 'model_eqns', 'num_fluids',     &
            & 'weno_order', 'precision', 'num_patches',                         &
            & 'n_start', 'num_bc_patches', 'recon_type',                        &
            & 'muscl_order', 'igr_order' ]
            call MPI_BCAST(${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        #:for VAR in [ 'old_grid','old_ic','stretch_x','stretch_y','stretch_z',&
            & 'parallel_io', 'file_per_process', 'cfl_adap_dt',               &
            & 'cfl_const_dt', 'cfl_dt', 'viscous',                            &
            & 'bc_io', 'down_sample', 'double_mach' ]
            call MPI_BCAST(${VAR}$, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
        #:endfor
        #:for VAR in [ 'x_domain%beg', 'x_domain%end', 'y_domain%beg',         &
            & 'y_domain%end', 'z_domain%beg', 'z_domain%end', 'a_x', 'a_y',    &
            & 'a_z', 'x_a', 'x_b', 'y_a', 'y_b', 'z_a', 'z_b', 'bc_x%beg',     &
            & 'bc_x%end', 'bc_y%beg', 'bc_y%end', 'bc_z%beg', 'bc_z%end',      &
            & 'pref', 'rhoref']
            call MPI_BCAST(${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
        #:endfor

        do i = 1, num_bc_patches_max
            #:for VAR in ['geometry', 'type', 'dir', 'loc']
                call MPI_BCAST(patch_bc(i)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_bc(i)%radius, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in ['centroid', 'length']
                call MPI_BCAST(patch_bc(i)%${VAR}$, size(patch_bc(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

        do i = 1, num_patches_max
            #:for VAR in [ 'geometry', 'smooth_patch_id']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%smoothen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%non_axis_sym, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
            call MPI_BCAST(patch_icpp(i)%alter_patch(0), num_patches_max, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

            #:for VAR in [ 'x_centroid', 'y_centroid', 'z_centroid',           &
                & 'length_x', 'length_y', 'length_z', 'radius', 'epsilon',     &
                & 'beta', 'smooth_coeff', 'rho',                                &
                & 'pres', 'gamma', 'pi_inf', 'hcid', 'cv', 'qv', 'qvp']
                call MPI_BCAST(patch_icpp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ '2', '3', '4', '5', '6', '7', '8', '9']
                call MPI_BCAST(patch_icpp(i)%a(${VAR}$), 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            #:for VAR in [ 'normal', 'radii', 'vel', 'alpha_rho', 'alpha' ]
                call MPI_BCAST(patch_icpp(i)%${VAR}$, size(patch_icpp(i)%${VAR}$), mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor

            call MPI_BCAST(patch_icpp(i)%sph_har_coeff, size(patch_icpp(i)%sph_har_coeff), mpi_p, 0, MPI_COMM_WORLD, ierr)
        end do

        ! Fluid physical parameters
        do i = 1, num_fluids_max
            #:for VAR in [ 'gamma','pi_inf', 'cv', 'qv', 'qvp' ]
                call MPI_BCAST(fluid_pp(i)%${VAR}$, 1, mpi_p, 0, MPI_COMM_WORLD, ierr)
            #:endfor
        end do

    end subroutine s_mpi_bcast_user_inputs

end module m_mpi_proxy
