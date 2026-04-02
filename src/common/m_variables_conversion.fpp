!>
!! @file
!! @brief Contains module m_variables_conversion

#:include 'macros.fpp'
#:include 'case.fpp'

!> @brief Conservative-to-primitive variable conversion, mixture property evaluation, and pressure computation
module m_variables_conversion

    use m_derived_types
    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic
    use m_helper

    implicit none

    private
    public :: s_initialize_variables_conversion_module, &
              s_convert_to_mixture_variables, &
              s_convert_species_to_mixture_variables, &
              s_convert_species_to_mixture_variables_acc, &
              s_convert_conservative_to_primitive_variables, &
              s_convert_primitive_to_conservative_variables, &
              s_convert_primitive_to_flux_variables, &
              s_compute_pressure, &
              s_compute_species_fraction, &
#ifndef MFC_PRE_PROCESS
    s_compute_speed_of_sound, &
#endif
    s_finalize_variables_conversion_module

    ! In simulation, gammas, pi_infs, and qvs are already declared in m_global_variables
#ifndef MFC_SIMULATION
    real(wp), allocatable, public, dimension(:) :: gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps
    $:GPU_DECLARE(create='[gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps]')
#endif

    real(wp), allocatable, dimension(:,:) :: Res_vc
    $:GPU_DECLARE(create='[Res_vc]')

    integer :: is1b, is2b, is3b, is1e, is2e, is3e
    $:GPU_DECLARE(create='[is1b, is2b, is3b, is1e, is2e, is3e]')

    real(wp), allocatable, dimension(:,:,:), public :: rho_sf     !< Scalar density function
    real(wp), allocatable, dimension(:,:,:), public :: gamma_sf   !< Scalar sp. heat ratio function
    real(wp), allocatable, dimension(:,:,:), public :: pi_inf_sf  !< Scalar liquid stiffness function
    real(wp), allocatable, dimension(:,:,:), public :: qv_sf      !< Scalar liquid energy reference function

contains

    !> Dispatch to s_convert_species_to_mixture_variables.
    subroutine s_convert_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv, Re_K)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in)                                 :: i, j, k
        real(wp), intent(out), target                       :: rho, gamma, pi_inf, qv
        real(wp), optional, dimension(2), intent(out)       :: Re_K

        call s_convert_species_to_mixture_variables(q_vf, i, j, k, rho, gamma, pi_inf, qv, Re_K)

    end subroutine s_convert_to_mixture_variables

    !> Compute the pressure from the appropriate equation of state
    subroutine s_compute_pressure(energy, alf, dyn_p, pi_inf, gamma, rho, qv, pres, stress, mom, pres_mag)

        $:GPU_ROUTINE(function_name='s_compute_pressure',parallelism='[seq]', cray_noinline=True)

        real(stp), intent(in)           :: energy, alf
        real(wp), intent(in)            :: dyn_p
        real(wp), intent(in)            :: pi_inf, gamma, rho, qv
        real(wp), intent(out)           :: pres
        real(stp), intent(in), optional :: stress, mom
        real(wp), intent(in), optional  :: pres_mag

        pres = (energy - dyn_p - pi_inf - qv)/gamma

    end subroutine s_compute_pressure

    !> Convert species volume fractions and partial densities to mixture density, gamma, pi_inf, and qv. Given conservative or
    !! primitive variables, computes the density, the specific heat ratio function and the liquid stiffness function from q_vf and
    !! stores the results into rho, gamma and pi_inf.
    subroutine s_convert_species_to_mixture_variables(q_vf, k, l, r, rho, gamma, pi_inf, qv, Re_K)

        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in)                                 :: k, l, r
        real(wp), intent(out), target                       :: rho
        real(wp), intent(out), target                       :: gamma
        real(wp), intent(out), target                       :: pi_inf
        real(wp), intent(out), target                       :: qv
        real(wp), optional, dimension(2), intent(out)       :: Re_K
        real(wp), dimension(num_fluids)                     :: alpha_rho_K, alpha_K
        integer                                             :: i, j  !< Generic loop iterator
        ! Computing the density, the specific heat ratio function and the liquid stiffness function, respectively

        call s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)

        ! Calculating the density, the specific heat ratio function, the liquid stiffness function, and the energy reference
        ! function, respectively, from the species analogs
        rho = 0._wp; gamma = 0._wp; pi_inf = 0._wp; qv = 0._wp
        do i = 1, num_fluids
            rho = rho + alpha_rho_K(i)
            gamma = gamma + alpha_K(i)*gammas(i)
            pi_inf = pi_inf + alpha_K(i)*pi_infs(i)
            qv = qv + alpha_rho_K(i)*qvs(i)
        end do

#ifdef MFC_SIMULATION
        ! Computing the shear and bulk Reynolds numbers from species analogs
        if (viscous) then
            do i = 1, 2
                Re_K(i) = dflt_real; if (Re_size(i) > 0) Re_K(i) = 0._wp

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/fluid_pp(Re_idx(i, j))%Re(i) + Re_K(i)
                end do

                Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
            end do
        end if
#endif

        ! Post process requires rho_sf/gamma_sf/pi_inf_sf/qv_sf to also be updated
#ifdef MFC_POST_PROCESS
        rho_sf(k, l, r) = rho
        gamma_sf(k, l, r) = gamma
        pi_inf_sf(k, l, r) = pi_inf
        qv_sf(k, l, r) = qv
#endif

    end subroutine s_convert_species_to_mixture_variables

    !> GPU-accelerated conversion of species volume fractions and partial densities to mixture density, gamma, pi_inf, and qv.
    subroutine s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K)

        $:GPU_ROUTINE(function_name='s_convert_species_to_mixture_variables_acc', parallelism='[seq]', cray_noinline=True)

        real(wp), intent(out) :: rho_K, gamma_K, pi_inf_K, qv_K
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(inout) :: alpha_rho_K, alpha_K
        #:else
            real(wp), dimension(num_fluids), intent(inout) :: alpha_rho_K, alpha_K
        #:endif
        real(wp), dimension(2), intent(out) :: Re_K
        real(wp)                            :: alpha_K_sum
        integer                             :: i, j  !< Generic loop iterators
        rho_K = 0._wp; gamma_K = 0._wp; pi_inf_K = 0._wp; qv_K = 0._wp
        do i = 1, num_fluids
            rho_K = rho_K + alpha_rho_K(i)
            gamma_K = gamma_K + alpha_K(i)*gammas(i)
            pi_inf_K = pi_inf_K + alpha_K(i)*pi_infs(i)
            qv_K = qv_K + alpha_rho_K(i)*qvs(i)
        end do

        if (viscous) then
            do i = 1, 2
                Re_K(i) = dflt_real

                if (Re_size(i) > 0) Re_K(i) = 0._wp

                do j = 1, Re_size(i)
                    Re_K(i) = alpha_K(Re_idx(i, j))/Res_vc(i, j) + Re_K(i)
                end do

                Re_K(i) = 1._wp/max(Re_K(i), sgm_eps)
            end do
        end if

    end subroutine s_convert_species_to_mixture_variables_acc

    !> Initialize the variables conversion module.
    impure subroutine s_initialize_variables_conversion_module

        integer :: i, j

        $:GPU_ENTER_DATA(copyin='[is1b, is1e, is2b, is2e, is3b, is3e]')

        @:ALLOCATE(gammas (1:num_fluids))
        @:ALLOCATE(gs_min (1:num_fluids))
        @:ALLOCATE(pi_infs(1:num_fluids))
        @:ALLOCATE(ps_inf(1:num_fluids))
        @:ALLOCATE(cvs    (1:num_fluids))
        @:ALLOCATE(qvs    (1:num_fluids))
        @:ALLOCATE(qvps    (1:num_fluids))

        do i = 1, num_fluids
            gammas(i) = fluid_pp(i)%gamma
            gs_min(i) = 1.0_wp/gammas(i) + 1.0_wp
            pi_infs(i) = fluid_pp(i)%pi_inf
            ps_inf(i) = pi_infs(i)/(1.0_wp + gammas(i))
            cvs(i) = fluid_pp(i)%cv
            qvs(i) = fluid_pp(i)%qv
            qvps(i) = fluid_pp(i)%qvp
        end do
        $:GPU_UPDATE(device='[gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps]')

#ifdef MFC_SIMULATION
        if (viscous) then
            @:ALLOCATE(Res_vc(1:2, 1:Re_size_max))
            do i = 1, 2
                do j = 1, Re_size(i)
                    Res_vc(i, j) = fluid_pp(Re_idx(i, j))%Re(i)
                end do
            end do

            $:GPU_UPDATE(device='[Res_vc, Re_idx, Re_size]')
        end if
#endif

#ifdef MFC_POST_PROCESS
        ! Allocating the density, the specific heat ratio function and the liquid stiffness function, respectively

        ! Simulation is at least 2D
        if (n > 0) then
            ! Simulation is 3D
            if (p > 0) then
                allocate (rho_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                allocate (gamma_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                allocate (pi_inf_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))
                allocate (qv_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,-buff_size:p + buff_size))

                ! Simulation is 2D
            else
                allocate (rho_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (gamma_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (pi_inf_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
                allocate (qv_sf(-buff_size:m + buff_size,-buff_size:n + buff_size,0:0))
            end if

            ! Simulation is 1D
        else
            allocate (rho_sf(-buff_size:m + buff_size,0:0,0:0))
            allocate (gamma_sf(-buff_size:m + buff_size,0:0,0:0))
            allocate (pi_inf_sf(-buff_size:m + buff_size,0:0,0:0))
            allocate (qv_sf(-buff_size:m + buff_size,0:0,0:0))
        end if
#endif

    end subroutine s_initialize_variables_conversion_module

    !> Convert conserved variables (rho*alpha, rho*u, E, alpha) to primitives (rho, u, p, alpha). Conversion depends on model_eqns:
    !! each model has different variable sets and EOS.
    subroutine s_convert_conservative_to_primitive_variables(qK_cons_vf, qK_prim_vf, ibounds)

        type(scalar_field), dimension(sys_size), intent(in)    :: qK_cons_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: qK_prim_vf
        type(int_bounds_info), dimension(1:3), intent(in)      :: ibounds

        #:if USING_AMD and not MFC_CASE_OPTIMIZATION
            real(wp), dimension(3) :: alpha_K, alpha_rho_K
            real(wp)               :: rhoYks(1:10)
        #:else
            real(wp), dimension(num_fluids) :: alpha_K, alpha_rho_K
            real(wp)                        :: rhoYks(1:1)
        #:endif
        real(wp), dimension(2) :: Re_K
        real(wp)               :: rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K
        real(wp)               :: pres
        integer                :: i, j, k, l  !< Generic loop iterators
        real(wp)               :: T
        real(wp)               :: pres_mag

        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_K, alpha_rho_K, Re_K, rho_K, gamma_K, pi_inf_K, qv_K, dyn_pres_K, rhoYks, &
                            & pres, T, pres_mag]')
        do l = ibounds(3)%beg, ibounds(3)%end
            do k = ibounds(2)%beg, ibounds(2)%end
                do j = ibounds(1)%beg, ibounds(1)%end
                    dyn_pres_K = 0._wp

                    call s_compute_species_fraction(qK_cons_vf, j, k, l, alpha_rho_K, alpha_K)

#ifdef MFC_SIMULATION
                    call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K)
#else
                    call s_convert_to_mixture_variables(qK_cons_vf, j, k, l, rho_K, gamma_K, pi_inf_K, qv_K)
#endif

                    ! Non-reacting: partial densities are directly primitive (alpha_i * rho_i)
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, contxe
                        qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                    end do

#ifdef MFC_SIMULATION
                    rho_K = max(rho_K, sgm_eps)
#endif

                    ! Recover velocity from momentum: u = rho*u / rho, and accumulate dynamic pressure 0.5*rho*|u|^2
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = momxb, momxe
                        qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)/rho_K
                        dyn_pres_K = dyn_pres_K + 5.e-1_wp*qK_cons_vf(i)%sf(j, k, l)*qK_prim_vf(i)%sf(j, k, l)
                    end do

                    pres_mag = 0._wp

                    call s_compute_pressure(qK_cons_vf(E_idx)%sf(j, k, l), qK_cons_vf(alf_idx)%sf(j, k, l), dyn_pres_K, pi_inf_K, &
                                            & gamma_K, rho_K, qv_K, pres, pres_mag=pres_mag)

                    qK_prim_vf(E_idx)%sf(j, k, l) = pres

                    if (num_fluids > 1) then
                        $:GPU_LOOP(parallelism='[seq]')
                        do i = advxb, advxe
                            qK_prim_vf(i)%sf(j, k, l) = qK_cons_vf(i)%sf(j, k, l)
                        end do
                    end if
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()

    end subroutine s_convert_conservative_to_primitive_variables

    !> Convert primitives (rho, u, p, alpha) to conserved variables (rho*alpha, rho*u, E, alpha).
    impure subroutine s_convert_primitive_to_conservative_variables(q_prim_vf, q_cons_vf)

        type(scalar_field), dimension(sys_size), intent(in)    :: q_prim_vf
        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf

        ! Density, specific heat ratio function, liquid stiffness function and dynamic pressure, as defined in the incompressible
        ! flow sense, respectively
        real(wp)               :: rho
        real(wp)               :: gamma
        real(wp)               :: pi_inf
        real(wp)               :: qv
        real(wp)               :: dyn_pres
        real(wp), dimension(2) :: Re_K
        integer                :: i, j, k, l  !< Generic loop iterators
        real(wp), dimension(1) :: Ys
        real(wp)               :: e_mix, mix_mol_weight, T

#ifndef MFC_SIMULATION
        ! Converting the primitive variables to the conservative variables
        do l = 0, p
            do k = 0, n
                do j = 0, m
                    ! Obtaining the density, specific heat ratio function and the liquid stiffness function, respectively
                    call s_convert_to_mixture_variables(q_prim_vf, j, k, l, rho, gamma, pi_inf, qv, Re_K)

                    if (num_fluids > 1) then
                        ! Transferring the advection equation(s) variable(s)
                        do i = adv_idx%beg, adv_idx%end
                            q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                        end do
                    end if

                    ! Transferring the continuity equation(s) variable(s)
                    do i = 1, contxe
                        q_cons_vf(i)%sf(j, k, l) = q_prim_vf(i)%sf(j, k, l)
                    end do

                    ! Zeroing out the dynamic pressure since it is computed iteratively by cycling through the velocity equations
                    dyn_pres = 0._wp

                    ! Computing momenta and dynamic pressure from velocity
                    do i = momxb, momxe
                        q_cons_vf(i)%sf(j, k, l) = rho*q_prim_vf(i)%sf(j, k, l)
                        dyn_pres = dyn_pres + q_cons_vf(i)%sf(j, k, l)*q_prim_vf(i)%sf(j, k, l)/2._wp
                    end do

                    ! Five-equation model: E = Gamma*p + 0.5*rho*|u|^2 + pi_inf + qv
                    q_cons_vf(E_idx)%sf(j, k, l) = gamma*q_prim_vf(E_idx)%sf(j, k, l) + dyn_pres + pi_inf + qv
                end do
            end do
        end do
#else
        if (proc_rank == 0) then
            call s_mpi_abort('Conversion from primitive to ' // 'conservative variables not ' // 'implemented. Exiting.')
        end if
#endif

    end subroutine s_convert_primitive_to_conservative_variables

    !> Convert primitive variables to Eulerian flux variables.
    subroutine s_convert_primitive_to_flux_variables(qK_prim_vf, FK_vf, FK_src_vf, is1, is2, is3, s2b, s3b)

        integer, intent(in)                                                           :: s2b, s3b
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(in)        :: qK_prim_vf
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,1:), intent(inout)     :: FK_vf
        real(wp), dimension(0:,idwbuff(2)%beg:,idwbuff(3)%beg:,advxb:), intent(inout) :: FK_src_vf
        type(int_bounds_info), intent(in)                                             :: is1, is2, is3

        ! Partial densities, density, velocity, pressure, energy, advection variables, the specific heat ratio and liquid stiffness
        ! functions, the shear and volume Reynolds numbers and the Weber numbers

        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3)  :: alpha_rho_K
            real(wp), dimension(3)  :: alpha_K
            real(wp), dimension(3)  :: vel_K
            real(wp), dimension(10) :: Y_K
        #:else
            real(wp), dimension(num_fluids) :: alpha_rho_K
            real(wp), dimension(num_fluids) :: alpha_K
            real(wp), dimension(num_vels)   :: vel_K
            real(wp), dimension(1)          :: Y_K
        #:endif
        real(wp)               :: rho_K
        real(wp)               :: vel_K_sum
        real(wp)               :: pres_K
        real(wp)               :: E_K
        real(wp)               :: gamma_K
        real(wp)               :: pi_inf_K
        real(wp)               :: qv_K
        real(wp), dimension(2) :: Re_K
        real(wp)               :: T_K, mix_mol_weight, R_gas
        integer                :: i, j, k, l  !< Generic loop iterators

        is1b = is1%beg; is1e = is1%end
        is2b = is2%beg; is2e = is2%end
        is3b = is3%beg; is3e = is3%end

        $:GPU_UPDATE(device='[is1b, is2b, is3b, is1e, is2e, is3e]')

        ! Computing the flux variables from the primitive variables, without accounting for the contribution of either viscosity or
        ! capillarity
#ifdef MFC_SIMULATION
        $:GPU_PARALLEL_LOOP(collapse=3, private='[alpha_rho_K, vel_K, alpha_K, Re_K, Y_K, rho_K, vel_K_sum, pres_K, E_K, gamma_K, &
                            & pi_inf_K, qv_K, T_K, mix_mol_weight, R_gas]')
        do l = is3b, is3e
            do k = is2b, is2e
                do j = is1b, is1e
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, contxe
                        alpha_rho_K(i) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advxb, advxe
                        alpha_K(i - E_idx) = qK_prim_vf(j, k, l, i)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K(i) = qK_prim_vf(j, k, l, contxe + i)
                    end do

                    vel_K_sum = 0._wp
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        vel_K_sum = vel_K_sum + vel_K(i)**2._wp
                    end do

                    pres_K = qK_prim_vf(j, k, l, E_idx)
                    call s_convert_species_to_mixture_variables_acc(rho_K, gamma_K, pi_inf_K, qv_K, alpha_K, alpha_rho_K, Re_K)

                    ! Computing the energy from the pressure

                    ! Computing the energy from the pressure
                    E_K = gamma_K*pres_K + pi_inf_K + 5.e-1_wp*rho_K*vel_K_sum + qv_K

                    ! mass flux, this should be \alpha_i \rho_i u_i
                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, contxe
                        FK_vf(j, k, l, i) = alpha_rho_K(i)*vel_K(dir_idx(1))
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = 1, num_vels
                        FK_vf(j, k, l, contxe + dir_idx(i)) = rho_K*vel_K(dir_idx(1))*vel_K(dir_idx(i)) + pres_K*dir_flg(dir_idx(i))
                    end do

                    ! energy flux, u(E+p)
                    FK_vf(j, k, l, E_idx) = vel_K(dir_idx(1))*(E_K + pres_K)

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advxb, advxe
                        FK_vf(j, k, l, i) = vel_K(dir_idx(1))*alpha_K(i - E_idx)
                    end do

                    $:GPU_LOOP(parallelism='[seq]')
                    do i = advxb, advxe
                        FK_src_vf(j, k, l, i) = vel_K(dir_idx(1))
                    end do
                end do
            end do
        end do
        $:END_GPU_PARALLEL_LOOP()
#endif

    end subroutine s_convert_primitive_to_flux_variables

    !> Compute partial densities and volume fractions
    subroutine s_compute_species_fraction(q_vf, k, l, r, alpha_rho_K, alpha_K)

        $:GPU_ROUTINE(function_name='s_compute_species_fraction', parallelism='[seq]', cray_noinline=True)
        type(scalar_field), dimension(sys_size), intent(in) :: q_vf
        integer, intent(in)                                 :: k, l, r
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(out) :: alpha_rho_K, alpha_K
        #:else
            real(wp), dimension(num_fluids), intent(out) :: alpha_rho_K, alpha_K
        #:endif
        integer  :: i
        real(wp) :: alpha_K_sum

        if (num_fluids == 1) then
            alpha_rho_K(1) = q_vf(contxb)%sf(k, l, r)
            alpha_K(1) = 1._wp
        else
            do i = 1, num_fluids - 1
                alpha_rho_K(i) = q_vf(i)%sf(k, l, r)
                alpha_K(i) = q_vf(advxb + i - 1)%sf(k, l, r)
            end do
            alpha_rho_K(num_fluids) = q_vf(num_fluids)%sf(k, l, r)
            alpha_K(num_fluids) = 1._wp - sum(alpha_K(1:num_fluids - 1))
        end if

    end subroutine s_compute_species_fraction

    !> Deallocate fluid property arrays and post-processing fields allocated during module initialization.
    impure subroutine s_finalize_variables_conversion_module()

        ! Deallocating the density, the specific heat ratio function and the liquid stiffness function
#ifdef MFC_POST_PROCESS
        deallocate (rho_sf, gamma_sf, pi_inf_sf, qv_sf)
#endif

        @:DEALLOCATE(gammas, gs_min, pi_infs, ps_inf, cvs, qvs, qvps)

    end subroutine s_finalize_variables_conversion_module

#ifndef MFC_PRE_PROCESS
    !> Compute the speed of sound from thermodynamic state variables, supporting multiple equation-of-state models.
    subroutine s_compute_speed_of_sound(pres, rho, gamma, pi_inf, H, adv, vel_sum, c_c, c, qv)

        $:GPU_ROUTINE(parallelism='[seq]')

        real(wp), intent(in) :: pres
        real(wp), intent(in) :: rho, gamma, pi_inf, qv
        real(wp), intent(in) :: H
        #:if not MFC_CASE_OPTIMIZATION and USING_AMD
            real(wp), dimension(3), intent(in) :: adv
        #:else
            real(wp), dimension(num_fluids), intent(in) :: adv
        #:endif
        real(wp), intent(in)  :: vel_sum
        real(wp), intent(in)  :: c_c
        real(wp), intent(out) :: c

        c = (H - 5.e-1*vel_sum - qv/rho)/gamma
        c = sqrt(c)

    end subroutine s_compute_speed_of_sound
#endif
end module m_variables_conversion
