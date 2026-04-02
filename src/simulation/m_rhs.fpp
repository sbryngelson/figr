!>
!! @file
!! @brief Contains module m_rhs

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Assembles the right-hand side of the governing equations using IGR Riemann solvers and physical source terms
module m_rhs

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_variables_conversion
    use m_nvtx
    use m_boundary_common
    use m_helper
    use m_igr

    implicit none

    private; public :: s_initialize_rhs_module, s_compute_rhs, s_finalize_rhs_module

contains

    !> Initialize the RHS module
    impure subroutine s_initialize_rhs_module

        $:GPU_ENTER_DATA(copyin='[idwbuff]')
        $:GPU_UPDATE(device='[idwbuff]')

    end subroutine s_initialize_rhs_module

    !> Compute the right-hand side of the semi-discrete governing equations for a single time stage
    impure subroutine s_compute_rhs(q_cons_vf, q_prim_vf, bc_type, rhs_vf, t_step, &

        & time_avg, stage)

        type(scalar_field), dimension(sys_size), intent(inout)     :: q_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout)     :: q_prim_vf
        type(integer_field), dimension(1:num_dims,1:2), intent(in) :: bc_type
        type(scalar_field), dimension(sys_size), intent(inout)     :: rhs_vf
        integer, intent(in) :: t_step
        real(wp), intent(inout) :: time_avg
        integer, intent(in) :: stage
        real(wp) :: t_start, t_finish
        integer :: id
        integer(kind=8) :: i, j, k, l, q  !< Generic loop iterators

        ! RHS: halo exchange -> IGR Riemann solve -> source terms

        call nvtxStartRange("COMPUTE-RHS")

        call cpu_time(t_start)

        call nvtxStartRange("RHS-COMMUNICATION")
        call s_populate_variables_buffers(bc_type, q_cons_vf)
        call nvtxEndRange

        if (cfl_dt) then
            if (mytime >= t_stop) return
        else
            if (t_step == t_step_stop) return
        end if

        ! Loop over coordinate directions for dimensional splitting
        do id = 1, num_dims
            if (id == 1) then
                $:GPU_PARALLEL_LOOP(private='[i, j, k, l]', collapse=4)
                do l = -1, p + 1
                    do k = -1, n + 1
                        do j = -1, m + 1
                            do i = 1, sys_size
                                rhs_vf(i)%sf(j, k, l) = 0._stp
                            end do
                        end do
                    end do
                end do
                $:END_GPU_PARALLEL_LOOP()
            end if

            call nvtxStartRange("IGR_RIEMANN")
            call s_igr_riemann_solver(q_cons_vf, rhs_vf, id)
            call nvtxEndRange

            if (id == 1) then
                call nvtxStartRange("IGR_Jacobi")
                call s_igr_iterative_solve(q_cons_vf, bc_type, t_step)
                call nvtxEndRange

                call nvtxStartRange("IGR_SIGMA")
                call s_igr_sigma_x(q_cons_vf, rhs_vf)
                call nvtxEndRange
            end if
        end do
        ! END: Dimensional Splitting Loop

        ! END: Additional physics and source terms

        call cpu_time(t_finish)

        if (t_step >= 2) then
            time_avg = (abs(t_finish - t_start) + (t_step - 2)*time_avg)/(t_step - 1)
        else
            time_avg = 0._wp
        end if

        call nvtxEndRange

    end subroutine s_compute_rhs

    !> Module deallocation and/or disassociation procedures
    impure subroutine s_finalize_rhs_module

    end subroutine s_finalize_rhs_module

end module m_rhs
