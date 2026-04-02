#:include 'macros.fpp'

!>
!! @file
!! @brief  Contains module m_start_up

!> @brief Reads and validates user inputs, allocates variables, and configures MPI decomposition and I/O for post-processing

module m_start_up

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_mpi_common
    use m_boundary_common
    use m_variables_conversion
    use m_data_input
    use m_data_output
    use m_derived_variables
    use m_helper
    use m_compile_specific
    use m_checker_common
    use m_checker
    use m_finite_differences
    ! m_chemistry removed (IGR-only build)

#ifdef MFC_MPI
    use mpi
#endif

    implicit none

    integer                                 :: ierr

contains

    !> Reads the configuration file post_process.inp, in order to populate parameters in module m_global_parameters.f90 with the
    !! user provided inputs
    impure subroutine s_read_input_file

        character(LEN=name_len) :: file_loc
        logical                 :: file_check
        integer                 :: iostatus
        character(len=1000)     :: line

        namelist /user_inputs/ case_dir, m, n, p, t_step_start, t_step_stop, t_step_save, model_eqns, num_fluids, bc_x, bc_y, &
            & bc_z, fluid_pp, format, precision, output_partial_domain, x_output, y_output, z_output, alpha_rho_wrt, rho_wrt, &
            & mom_wrt, vel_wrt, E_wrt, pres_wrt, alpha_wrt, gamma_wrt, heat_ratio_wrt, pi_inf_wrt, pres_inf_wrt, cons_vars_wrt, &
            & prim_vars_wrt, c_wrt, omega_wrt, qm_wrt, schlieren_wrt, schlieren_alpha, fd_order, flux_lim, flux_wrt, &
            & parallel_io, rhoref, pref, file_per_process, cfl_adap_dt, cfl_const_dt, t_save, t_stop, n_start, cfl_target, &
            & sim_data, num_bc_patches, igr_order, down_sample, alpha_rho_e_wrt

        file_loc = 'post_process.inp'
        inquire (FILE=trim(file_loc), EXIST=file_check)

        if (file_check) then
            open (1, FILE=trim(file_loc), form='formatted', STATUS='old', ACTION='read')
            read (1, NML=user_inputs, iostat=iostatus)

            if (iostatus /= 0) then
                backspace (1)
                read (1, fmt='(A)') line
                print *, 'Invalid line in namelist: ' // trim(line)
                call s_mpi_abort('Invalid line in post_process.inp. It is ' // 'likely due to a datatype mismatch. Exiting.')
            end if

            close (1)

            call s_update_cell_bounds(cells_bounds, m, n, p)

            if (down_sample) then
                m = int((m + 1)/3) - 1
                n = int((n + 1)/3) - 1
                p = int((p + 1)/3) - 1
            end if

            m_glb = m
            n_glb = n
            p_glb = p

            nGlobal = int(m_glb + 1, kind=8)*int(n_glb + 1, kind=8)*int(p_glb + 1, kind=8)

            if (cfl_adap_dt .or. cfl_const_dt) cfl_dt = .true.

            if (any((/bc_x%beg, bc_x%end, bc_y%beg, bc_y%end, bc_z%beg, bc_z%end/) == -17) .or. num_bc_patches > 0) then
                bc_io = .true.
            end if
        else
            call s_mpi_abort('File post_process.inp is missing. Exiting.')
        end if

    end subroutine s_read_input_file

    !> Checking that the user inputs make sense, i.e. that the individual choices are compatible with the code's options and that
    !! the combination of these choices results into a valid configuration for the post-process
    impure subroutine s_check_input_file

        character(LEN=len_trim(case_dir)) :: file_loc
        logical                           :: dir_check

        case_dir = adjustl(case_dir)

        file_loc = trim(case_dir) // '/.'

        call my_inquire(file_loc, dir_check)

        if (dir_check .neqv. .true.) then
            call s_mpi_abort('Unsupported choice for the value of ' // 'case_dir. Exiting.')
        end if

        call s_check_inputs_common()
        call s_check_inputs()

    end subroutine s_check_input_file

    !> Load grid and conservative data for a time step, fill ghost-cell buffers, and convert to primitive variables.
    impure subroutine s_perform_time_step(t_step)

        integer, intent(inout) :: t_step

        if (proc_rank == 0) then
            if (cfl_dt) then
                print '(" [", I3, "%]  Saving ", I8, " of ", I0, " Time Avg = ", ES16.6,  " Time/step = ", ES12.6, "")', &
                    & int(ceiling(100._wp*(real(t_step - n_start)/(n_save)))), t_step, n_save, wall_time_avg, wall_time
            else
                print '(" [", I3, "%]  Saving ", I8, " of ", I0, " @ t_step = ", I8, " Time Avg = ", ES16.6,  " Time/step = ", ES12.6, "")', &
                    & int(ceiling(100._wp*(real(t_step - t_step_start)/(t_step_stop - t_step_start + 1)))), &
                    & (t_step - t_step_start)/t_step_save + 1, (t_step_stop - t_step_start)/t_step_save + 1, t_step, &
                    & wall_time_avg, wall_time
            end if
        end if

        call s_read_data_files(t_step)

        if (buff_size > 0) then
            call s_populate_grid_variables_buffers()
            call s_populate_variables_buffers(bc_type, q_cons_vf)
        end if

        call s_convert_conservative_to_primitive_variables(q_cons_vf, q_prim_vf, idwbuff)

    end subroutine s_perform_time_step

    !> Derive requested flow quantities from primitive variables and write them to the formatted database files.
    impure subroutine s_save_data(t_step, varname, pres, c, H)

        integer, intent(inout)                 :: t_step
        character(LEN=name_len), intent(inout) :: varname
        real(wp), intent(inout)                :: pres, c, H
        real(wp)                               :: theta1, theta2

        integer       :: i, j, k, l
        integer       :: x_beg, x_end, y_beg, y_end, z_beg, z_end

        if (output_partial_domain) then
            call s_define_output_region
            x_beg = -offset_x%beg + x_output_idx%beg
            x_end = offset_x%end + x_output_idx%end
            y_beg = -offset_y%beg + y_output_idx%beg
            y_end = offset_y%end + y_output_idx%end
            z_beg = -offset_z%beg + z_output_idx%beg
            z_end = offset_z%end + z_output_idx%end
        else
            x_beg = -offset_x%beg
            x_end = offset_x%end + m
            y_beg = -offset_y%beg
            y_end = offset_y%end + n
            z_beg = -offset_z%beg
            z_end = offset_z%end + p
        end if

        call s_open_formatted_database_file(t_step)

        if (sim_data .and. proc_rank == 0) then
            call s_open_intf_data_file()
            call s_open_energy_data_file()
        end if

        if (sim_data) then
            call s_write_intf_data_file(q_prim_vf)
            call s_write_energy_data_file(q_prim_vf, q_cons_vf)
        end if

        call s_write_grid_to_formatted_database_file(t_step)

        if (omega_wrt(2) .or. omega_wrt(3) .or. qm_wrt .or. schlieren_wrt) then
            call s_compute_finite_difference_coefficients(m, x_cc, fd_coeff_x, buff_size, fd_number, fd_order, offset_x)
        end if

        if (omega_wrt(1) .or. omega_wrt(3) .or. qm_wrt .or. (n > 0 .and. schlieren_wrt)) then
            call s_compute_finite_difference_coefficients(n, y_cc, fd_coeff_y, buff_size, fd_number, fd_order, offset_y)
        end if

        if (omega_wrt(1) .or. omega_wrt(2) .or. qm_wrt .or. (p > 0 .and. schlieren_wrt)) then
            call s_compute_finite_difference_coefficients(p, z_cc, fd_coeff_z, buff_size, fd_number, fd_order, offset_z)
        end if

        if (model_eqns == 2) then
            do i = 1, num_fluids
                if (alpha_rho_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    q_sf(:,:,:) = q_cons_vf(i)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I0)') 'alpha_rho', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '
                end if
            end do
        end if

        if (rho_wrt) then
            q_sf(:,:,:) = rho_sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'rho'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        do i = 1, E_idx - mom_idx%beg
            if (mom_wrt(i) .or. cons_vars_wrt) then
                q_sf(:,:,:) = q_cons_vf(i + cont_idx%end)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A,I0)') 'mom', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        do i = 1, E_idx - mom_idx%beg
            if (vel_wrt(i) .or. prim_vars_wrt) then
                q_sf(:,:,:) = q_prim_vf(i + cont_idx%end)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                write (varname, '(A,I0)') 'vel', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        do i = 1, E_idx - mom_idx%beg
            if (flux_wrt(i)) then
                call s_derive_flux_limiter(i, q_prim_vf, q_sf)

                write (varname, '(A,I0)') 'flux', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        if (E_wrt .or. cons_vars_wrt) then
            q_sf(:,:,:) = q_cons_vf(E_idx)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'E'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (pres_wrt .or. prim_vars_wrt) then
            q_sf(:,:,:) = q_prim_vf(E_idx)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'pres'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (model_eqns == 2) then
            do i = 1, num_fluids - 1
                if (alpha_wrt(i) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                    q_sf(:,:,:) = q_cons_vf(i + E_idx)%sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
                    write (varname, '(A,I0)') 'alpha', i
                    call s_write_variable_to_formatted_database_file(varname, t_step)

                    varname(:) = ' '
                end if
            end do

            if (alpha_wrt(num_fluids) .or. (cons_vars_wrt .or. prim_vars_wrt)) then
                do k = z_beg, z_end
                    do j = y_beg, y_end
                        do i = x_beg, x_end
                            q_sf(i, j, k) = 1._wp
                            do l = 1, num_fluids - 1
                                q_sf(i, j, k) = q_sf(i, j, k) - q_cons_vf(E_idx + l)%sf(i, j, k)
                            end do
                        end do
                    end do
                end do
                write (varname, '(A,I0)') 'alpha', num_fluids
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end if

        if (gamma_wrt) then
            q_sf(:,:,:) = gamma_sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'gamma'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (heat_ratio_wrt) then
            call s_derive_specific_heat_ratio(q_sf)

            write (varname, '(A)') 'heat_ratio'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (pi_inf_wrt) then
            q_sf(:,:,:) = pi_inf_sf(x_beg:x_end,y_beg:y_end,z_beg:z_end)
            write (varname, '(A)') 'pi_inf'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (pres_inf_wrt) then
            call s_derive_liquid_stiffness(q_sf)

            write (varname, '(A)') 'pres_inf'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (c_wrt) then
            do k = -offset_z%beg, p + offset_z%end
                do j = -offset_y%beg, n + offset_y%end
                    do i = -offset_x%beg, m + offset_x%end
                        do l = 1, adv_idx%end - E_idx
                            adv(l) = q_prim_vf(E_idx + l)%sf(i, j, k)
                        end do

                        pres = q_prim_vf(E_idx)%sf(i, j, k)

                        H = ((gamma_sf(i, j, k) + 1._wp)*pres + pi_inf_sf(i, j, k) + qv_sf(i, j, k))/rho_sf(i, j, k)

                        call s_compute_speed_of_sound(pres, rho_sf(i, j, k), gamma_sf(i, j, k), pi_inf_sf(i, j, k), H, adv, &
                                                      & 0._wp, 0._wp, c, qv_sf(i, j, k))

                        q_sf(i, j, k) = c
                    end do
                end do
            end do

            write (varname, '(A)') 'c'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        do i = 1, 3
            if (omega_wrt(i)) then
                call s_derive_vorticity_component(i, q_prim_vf, q_sf)

                write (varname, '(A,I0)') 'omega', i
                call s_write_variable_to_formatted_database_file(varname, t_step)

                varname(:) = ' '
            end if
        end do

        if (p > 0 .and. qm_wrt) then
            call s_derive_qm(q_prim_vf, q_sf)

            write (varname, '(A)') 'qm'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (schlieren_wrt) then
            call s_derive_numerical_schlieren_function(q_cons_vf, q_sf)

            write (varname, '(A)') 'schlieren'
            call s_write_variable_to_formatted_database_file(varname, t_step)

            varname(:) = ' '
        end if

        if (sim_data .and. proc_rank == 0) then
            call s_close_intf_data_file()
            call s_close_energy_data_file()
        end if

        call s_close_formatted_database_file()

    end subroutine s_save_data

    !> Initialize all post-process sub-modules and set up I/O pointers.
    impure subroutine s_initialize_modules

        call s_initialize_global_parameters_module()
        if (num_procs > 1) then
            call s_initialize_mpi_proxy_module()
            call s_initialize_mpi_common_module()
        end if
        call s_initialize_boundary_common_module()
        call s_initialize_variables_conversion_module()
        call s_initialize_data_input_module()
        call s_initialize_derived_variables_module()
        call s_initialize_data_output_module()

        if (parallel_io .neqv. .true.) then
            s_read_data_files => s_read_serial_data_files
        else
            s_read_data_files => s_read_parallel_data_files
        end if

    end subroutine s_initialize_modules

    !> Set up the MPI environment, read and broadcast user inputs, and decompose the computational domain.
    impure subroutine s_initialize_mpi_domain

        num_dims = 1 + min(1, n) + min(1, p)

        call s_mpi_initialize()

        if (proc_rank == 0) then
            call s_assign_default_values_to_user_inputs()
            call s_read_input_file()
            call s_check_input_file()

            print '(" Post-processing a ", I0, "x", I0, "x", I0, " case on ", I0, " rank(s)")', m, n, p, num_procs
        end if

        call s_mpi_bcast_user_inputs()
        call s_initialize_parallel_io()
        call s_mpi_decompose_computational_domain()

    end subroutine s_initialize_mpi_domain

    !> Finalize all post-process sub-modules.
    impure subroutine s_finalize_modules

        s_read_data_files => null()

        call s_finalize_data_output_module()
        call s_finalize_derived_variables_module()
        call s_finalize_data_input_module()
        call s_finalize_variables_conversion_module()
        if (num_procs > 1) then
            call s_finalize_mpi_proxy_module()
            call s_finalize_mpi_common_module()
        end if
        call s_finalize_global_parameters_module()

        call s_mpi_finalize()

    end subroutine s_finalize_modules

end module m_start_up
