!>
!! @file
!! @brief Contains module m_assign_variables

#:include 'case.fpp'
#:include 'macros.fpp'

!> @brief Assigns initial primitive variables to computational cells based on patch geometry
module m_assign_variables

    use m_derived_types
    use m_global_parameters
    use m_variables_conversion
    use m_helper_basic

    implicit none

    type(scalar_field) :: alf_sum

    !> Pointer to mixture or species patch assignment routine
    procedure(s_assign_patch_xxxxx_primitive_variables), pointer :: s_assign_patch_primitive_variables => null()
    !> Abstract interface to the two subroutines that assign the patch primitive variables, either mixture or species, depending on
    !! the subroutine, to a particular cell in the computational domain
    abstract interface

        !> Skeleton of s_assign_patch_species_primitive_variables
        subroutine s_assign_patch_xxxxx_primitive_variables(patch_id, j, k, l, eta, q_prim_vf, patch_id_fp)

            import :: scalar_field, sys_size, n, m, p, wp

            integer, intent(in)                                      :: patch_id
            integer, intent(in)                                      :: j, k, l
            real(wp), intent(in)                                     :: eta
            type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

#ifdef MFC_MIXED_PRECISION
            integer(kind=1), dimension(0:m,0:n,0:p), intent(inout) :: patch_id_fp
#else
            integer, dimension(0:m,0:n,0:p), intent(inout) :: patch_id_fp
#endif

        end subroutine s_assign_patch_xxxxx_primitive_variables
    end interface

    private
    public :: s_initialize_assign_variables_module, s_assign_patch_primitive_variables, &
        & s_assign_patch_species_primitive_variables, s_finalize_assign_variables_module

contains

    !> Allocate volume fraction sum and set the patch primitive variable assignment procedure pointer.
    impure subroutine s_initialize_assign_variables_module

        ! Select procedure pointer based on multicomponent flow model

        ! Volume fraction model (model_eqns == 2 only in IGR-only build)
        s_assign_patch_primitive_variables => s_assign_patch_species_primitive_variables

    end subroutine s_initialize_assign_variables_module

    !> Assign the species primitive variables for the volume fraction model
    impure subroutine s_assign_patch_species_primitive_variables(patch_id, j, k, l, eta, q_prim_vf, patch_id_fp)

        $:GPU_ROUTINE(parallelism='[seq]')

        integer, intent(in)  :: patch_id
        integer, intent(in)  :: j, k, l
        real(wp), intent(in) :: eta
#ifdef MFC_MIXED_PRECISION
        integer(kind=1), dimension(0:m,0:n,0:p), intent(inout) :: patch_id_fp
#else
        integer, dimension(0:m,0:n,0:p), intent(inout) :: patch_id_fp
#endif
        type(scalar_field), dimension(1:sys_size), intent(inout) :: q_prim_vf

        ! Density, gamma, and liquid stiffness from current and smoothing patches
        real(wp)                       :: rho           !< density
        real(wp)                       :: gamma
        real(wp)                       :: pi_inf        !< stiffness from SEOS
        real(wp)                       :: qv            !< reference energy from SEOS
        real(wp)                       :: orig_rho
        real(wp)                       :: orig_gamma
        real(wp)                       :: orig_pi_inf
        real(wp)                       :: orig_qv
        real(wp)                       :: Ys(1:1)
        real(stp), dimension(sys_size) :: orig_prim_vf  !< Vector to hold original values of cell for smoothing purposes
        integer                        :: i
        integer                        :: smooth_patch_id

        smooth_patch_id = patch_icpp(patch_id)%smooth_patch_id

        do i = 1, sys_size
            orig_prim_vf(i) = q_prim_vf(i)%sf(j, k, l)
        end do

        call s_convert_to_mixture_variables(q_prim_vf, j, k, l, orig_rho, orig_gamma, orig_pi_inf, orig_qv)

        if (num_fluids > 1) then
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha(i - E_idx)
            end do
        end if

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(patch_id)%alpha_rho(i)
        end do

        call s_convert_to_mixture_variables(q_prim_vf, j, k, l, patch_icpp(patch_id)%rho, patch_icpp(patch_id)%gamma, &
                                            & patch_icpp(patch_id)%pi_inf, patch_icpp(patch_id)%qv)

        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha_rho(i)
        end do

        if (num_fluids > 1) then
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j, k, l) = patch_icpp(smooth_patch_id)%alpha(i - E_idx)
            end do
        end if

        call s_convert_to_mixture_variables(q_prim_vf, j, k, l, patch_icpp(smooth_patch_id)%rho, &
                                            & patch_icpp(smooth_patch_id)%gamma, patch_icpp(smooth_patch_id)%pi_inf, &
                                            & patch_icpp(smooth_patch_id)%qv)

        q_prim_vf(E_idx)%sf(j, k, l) = (eta*patch_icpp(patch_id)%pres + (1._wp - eta)*orig_prim_vf(E_idx))

        if (num_fluids > 1) then
            do i = adv_idx%beg, adv_idx%end
                q_prim_vf(i)%sf(j, k, l) = eta*patch_icpp(patch_id)%alpha(i - E_idx) + (1._wp - eta)*orig_prim_vf(i)
            end do
        end if

        ! mixture density is an input
        do i = 1, cont_idx%end
            q_prim_vf(i)%sf(j, k, l) = eta*patch_icpp(patch_id)%alpha_rho(i) + (1._wp - eta)*orig_prim_vf(i)
        end do

        call s_convert_to_mixture_variables(q_prim_vf, j, k, l, rho, gamma, pi_inf, qv)

        do i = 1, E_idx - mom_idx%beg
            q_prim_vf(i + cont_idx%end)%sf(j, k, &
                      & l) = (eta*patch_icpp(patch_id)%vel(i) + (1._wp - eta)*orig_prim_vf(i + cont_idx%end))
        end do

        ! Set streamwise velocity to hyperbolic tangent function of y
        if (mixlayer_vel_profile) then
            q_prim_vf(1 + cont_idx%end)%sf(j, k, &
                      & l) = (eta*patch_icpp(patch_id)%vel(1)*tanh(y_cc(k)*mixlayer_vel_coef) + (1._wp - eta)*orig_prim_vf(1 &
                      & + cont_idx%end))
        end if

        if (1._wp - eta < 1.e-16_wp) patch_id_fp(j, k, l) = patch_id

    end subroutine s_assign_patch_species_primitive_variables

    !> Nullify the patch primitive variable assignment procedure pointer.
    impure subroutine s_finalize_assign_variables_module

        s_assign_patch_primitive_variables => null()

    end subroutine s_finalize_assign_variables_module

end module m_assign_variables
