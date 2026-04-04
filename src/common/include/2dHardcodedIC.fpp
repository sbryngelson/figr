#:def Hardcoded2DVariables()
    ! Place any declaration of intermediate variables here
    real(wp) :: eps
    real(wp) :: r, rmax, gam, umax, p0
    real(wp) :: rhoH, rhoL, pRef, pInt, h, lam, wl, amp, intH, intL, alph

    ! # 207
    real(wp) :: sigma, gauss1, gauss2
    ! # 208
    real(wp) :: ei, d, fsm, alpha_air, alpha_sf6
    ! #284
    real(wp) :: ux_tl, ux_tr, ux_bl, ux_br, uy_tl, uy_tr, uy_bl, uy_br, rho_tl, rho_tr, rho_bl, rho_br, En_tl, En_tr, En_bl, &
         & En_br, p_tl, p_tr, p_bl, p_br
    ! #285
    real(wp) :: sb

    eps = 1.e-9_wp
#:enddef

#:def Hardcoded2D()
    select case (patch_icpp(patch_id)%hcid)  ! 2D_hardcoded_ic example case
    case (200)  ! Two-fluid cubic interface
        if (y_cc(j) <= (-x_cc(i)**3 + 1)**(1._wp/3._wp)) then
            ! Volume Fractions
            q_prim_vf(advxb)%sf(i, j, 0) = eps
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - eps
            q_prim_vf(contxb)%sf(i, j, 0) = eps*1000._wp
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - eps)*1._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = 1000._wp
        end if
    case (202)  ! Gresho vortex (Gouasmi et al 2022 JCP)
        r = ((x_cc(i) - 0.5_wp)**2 + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2._wp/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4*(1 - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*(-2 + 4*log(2._wp))
        end if
    case (203)  ! Gresho vortex (Gouasmi et al 2022 JCP) with density correction
        r = ((x_cc(i) - 0.5_wp)**2._wp + (y_cc(j) - 0.5_wp)**2)**0.5_wp
        rmax = 0.2_wp

        gam = 1._wp + 1._wp/fluid_pp(1)%gamma
        umax = 2*pi*rmax*patch_icpp(patch_id)%vel(2)
        p0 = umax**2*(1._wp/(gam*patch_icpp(patch_id)%vel(2)**2) - 0.5_wp)

        if (r < rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -(y_cc(j) - 0.5_wp)*umax/rmax
            q_prim_vf(momxe)%sf(i, j, 0) = (x_cc(i) - 0.5_wp)*umax/rmax
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2._wp/2._wp)
        else if (r < 2*rmax) then
            q_prim_vf(momxb)%sf(i, j, 0) = -((y_cc(j) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(momxe)%sf(i, j, 0) = ((x_cc(i) - 0.5_wp)/r)*umax*(2._wp - r/rmax)
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2*((r/rmax)**2/2._wp + 4._wp*(1._wp - (r/rmax) + log(r/rmax)))
        else
            q_prim_vf(momxb)%sf(i, j, 0) = 0._wp
            q_prim_vf(momxe)%sf(i, j, 0) = 0._wp
            q_prim_vf(E_idx)%sf(i, j, 0) = p0 + umax**2._wp*(-2._wp + 4*log(2._wp))
        end if

        q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(E_idx)%sf(i, j, 0)**(1._wp/gam)
    case (204)  ! Rayleigh-Taylor instability
        rhoH = 3._wp
        rhoL = 1._wp
        pRef = 1.e5_wp
        pInt = pRef
        h = 0.7_wp
        lam = 0.2_wp
        wl = 2._wp*pi/lam
        amp = 0.05_wp/wl

        intH = amp*sin(2._wp*pi*x_cc(i)/lam - pi/2._wp) + h

        alph = 0.5_wp*(1._wp + tanh((y_cc(j) - intH)/2.5e-3_wp))

        if (alph < eps) alph = eps
        if (alph > 1._wp - eps) alph = 1._wp - eps

        if (y_cc(j) > intH) then
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            q_prim_vf(E_idx)%sf(i, j, 0) = pref + rhoH*9.81_wp*(1.2_wp - y_cc(j))
        else
            q_prim_vf(advxb)%sf(i, j, 0) = alph
            q_prim_vf(advxe)%sf(i, j, 0) = 1._wp - alph
            q_prim_vf(contxb)%sf(i, j, 0) = alph*rhoH
            q_prim_vf(contxe)%sf(i, j, 0) = (1._wp - alph)*rhoL
            pInt = pref + rhoH*9.81_wp*(1.2_wp - intH)
            q_prim_vf(E_idx)%sf(i, j, 0) = pInt + rhoL*9.81_wp*(intH - y_cc(j))
        end if
    case (205)  ! 2D lung wave interaction problem
        h = 0.0_wp  ! non dim origin y
        lam = 1.0_wp  ! non dim lambda
        amp = patch_icpp(patch_id)%a(2)  ! to be changed later!       !non dim amplitude

        intH = amp*sin(2*pi*x_cc(i)/lam - pi/2) + h

        if (y_cc(j) > intH) then
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if
    case (206)  ! 2D lung wave interaction problem - horizontal domain
        h = 0.0_wp  ! non dim origin y
        lam = 1.0_wp  ! non dim lambda
        amp = patch_icpp(patch_id)%a(2)

        intL = amp*sin(2*pi*y_cc(j)/lam - pi/2) + h

        if (x_cc(i) > intL) then  ! this is the liquid
            q_prim_vf(contxb)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(1)
            q_prim_vf(contxe)%sf(i, j, 0) = patch_icpp(1)%alpha_rho(2)
            q_prim_vf(E_idx)%sf(i, j, 0) = patch_icpp(1)%pres
            q_prim_vf(advxb)%sf(i, j, 0) = patch_icpp(1)%alpha(1)
            q_prim_vf(advxe)%sf(i, j, 0) = patch_icpp(1)%alpha(2)
        end if
    case (207)  ! Kelvin Helmholtz Instability
        sigma = 0.05_wp/sqrt(2.0_wp)
        gauss1 = exp(-(y_cc(j) - 0.75_wp)**2/(2.0_wp*sigma**2))
        gauss2 = exp(-(y_cc(j) - 0.25_wp)**2/(2.0_wp*sigma**2))
        q_prim_vf(momxb + 1)%sf(i, j, 0) = 0.1_wp*sin(4.0_wp*pi*x_cc(i))*(gauss1 + gauss2)
    case (208)  ! Richtmeyer Meshkov Instability
        lam = 1.0_wp
        eps = 1.0e-6_wp
        ei = 5.0_wp
        ! Smoothening function to smooth out sharp discontinuity in the interface
        if (x_cc(i) <= 0.7_wp*lam) then
            d = x_cc(i) - lam*(0.4_wp - 0.1_wp*sin(2.0_wp*pi*(y_cc(j)/lam + 0.25_wp)))
            fsm = 0.5_wp*(1.0_wp + erf(d/(ei*sqrt(dx*dy))))
            alpha_air = eps + (1.0_wp - 2.0_wp*eps)*fsm
            alpha_sf6 = 1.0_wp - alpha_air
            q_prim_vf(contxb)%sf(i, j, 0) = alpha_sf6*5.04_wp
            q_prim_vf(contxe)%sf(i, j, 0) = alpha_air*1.0_wp
            q_prim_vf(advxb)%sf(i, j, 0) = alpha_sf6
            q_prim_vf(advxe)%sf(i, j, 0) = alpha_air
        end if
    case (270)  ! 2D extrusion of 1D profile from external data
        ! This hardcoded case extrudes a 1D profile to initialize a 2D simulation domain
        @: HardcodedReadValues()
    case (280)  ! Isentropic vortex
        ! This is patch is hard-coded for test suite optimization used in the 2D_isentropicvortex case: This analytic patch uses
        ! geometry 2
        if (patch_id == 1) then
            q_prim_vf(E_idx)%sf(i, j, &
                      & 0) = 1.0*(1.0 - (1.0/1.0)*(5.0/(2.0*pi))*(5.0/(8.0*1.0*(1.4 + 1.0)*pi))*exp(2.0*1.0*(1.0 - (x_cc(i) &
                      & - patch_icpp(1)%x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0)))**(1.4 + 1.0)
            q_prim_vf(contxb + 0)%sf(i, j, &
                      & 0) = 1.0*(1.0 - (1.0/1.0)*(5.0/(2.0*pi))*(5.0/(8.0*1.0*(1.4 + 1.0)*pi))*exp(2.0*1.0*(1.0 - (x_cc(i) &
                      & - patch_icpp(1)%x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0)))**1.4
            q_prim_vf(momxb + 0)%sf(i, j, &
                      & 0) = 0.0 + (y_cc(j) - patch_icpp(1)%y_centroid)*(5.0/(2.0*pi))*exp(1.0*(1.0 - (x_cc(i) - patch_icpp(1) &
                      & %x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0))
            q_prim_vf(momxb + 1)%sf(i, j, &
                      & 0) = 0.0 - (x_cc(i) - patch_icpp(1)%x_centroid)*(5.0/(2.0*pi))*exp(1.0*(1.0 - (x_cc(i) - patch_icpp(1) &
                      & %x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0))
        end if
    case (281)  ! Acoustic pulse
        ! This is patch is hard-coded for test suite optimization used in the 2D_acoustic_pulse case: This analytic patch uses
        ! geometry 2
        if (patch_id == 2) then
            q_prim_vf(E_idx)%sf(i, j, &
                      & 0) = 101325*(1 - 0.5*(1.4 - 1)*(0.4)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1.4/(1.4 - 1))
            q_prim_vf(contxb + 0)%sf(i, j, &
                      & 0) = 1*(1 - 0.5*(1.4 - 1)*(0.4)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1/(1.4 - 1))
        end if
    case (282)  ! Zero-circulation vortex
        ! This is patch is hard-coded for test suite optimization used in the 2D_zero_circ_vortex case: This analytic patch uses
        ! geometry 2
        if (patch_id == 2) then
            q_prim_vf(E_idx)%sf(i, j, &
                      & 0) = 101325*(1 - 0.5*(1.4 - 1)*(0.1/0.3)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1.4/(1.4 - 1))
            q_prim_vf(contxb + 0)%sf(i, j, &
                      & 0) = 1*(1 - 0.5*(1.4 - 1)*(0.1/0.3)**2*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2))))**(1/(1.4 - 1))
            q_prim_vf(momxb + 0)%sf(i, j, &
                      & 0) = 112.99092883944267*(1 - (0.1/0.3))*y_cc(j)*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2)))
            q_prim_vf(momxb + 1)%sf(i, j, 0) = 112.99092883944267*((0.1/0.3))*x_cc(i)*exp(0.5*(1 - sqrt(x_cc(i)**2 + y_cc(j)**2)))
        end if
    case (283)  ! IGR Isentropic vortex
        ! Modified initial condition so as to match the initial condition used in IGR User Guide
        if (patch_id == 1) then
            q_prim_vf(E_idx)%sf(i, j, &
                      & 0) = 1.0*(1.0 - (1.0/1.0)*(5.0**2*0.4/(8.0*(1.4)*pi**2))*exp(1.0*(1.0 - (x_cc(i) - patch_icpp(1) &
                      & %x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0)))**(1.4/0.4)
            q_prim_vf(contxb + 0)%sf(i, j, &
                      & 0) = 1.0*(1.0 - (1.0/1.0)*(5.0**2*0.4/(8.0*(1.4)*pi**2))*exp(1.0*(1.0 - (x_cc(i) - patch_icpp(1) &
                      & %x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0)))**(1/0.4)
            q_prim_vf(momxb + 0)%sf(i, j, &
                      & 0) = 0.1 + (y_cc(j) - patch_icpp(1)%y_centroid)*(5.0/(2.0*pi))*exp(0.5*(1.0 - (x_cc(i) - patch_icpp(1) &
                      & %x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0))
            q_prim_vf(momxb + 1)%sf(i, j, &
                      & 0) = 0.0 - (x_cc(i) - patch_icpp(1)%x_centroid)*(5.0/(2.0*pi))*exp(0.5*(1.0 - (x_cc(i) - patch_icpp(1) &
                      & %x_centroid)**2.0 - (y_cc(j) - patch_icpp(1)%y_centroid)**2.0))
        end if
    case(284)  ! 2D IGR riemann test
        ! Modified initial condition to include entropy wave with fluctuations in density and energy
        if (patch_id == 1) then
            gam = 1.4_wp
            ux_tl = 1.206_wp
            ux_tr = 0.0_wp
            ux_bl = 1.206_wp
            ux_br = 0.0_wp

            uy_tl = 0.0_wp
            uy_tr = 0.0_wp
            uy_bl = 1.206_wp
            uy_br = 1.206_wp

            rho_tl = 0.5323_wp
            rho_tr = 1.5_wp
            rho_bl = 0.138_wp
            rho_br = 0.5323_wp

            p_tl = 0.3_wp
            p_tr = 1.5_wp
            p_bl = 0.029_wp
            p_br = 0.3_wp

            En_tl = p_tl/((gam - 1.0)*rho_tl)
            En_tr = p_tr/((gam - 1.0)*rho_tr)
            En_bl = p_bl/((gam - 1.0)*rho_bl)
            En_br = p_br/((gam - 1.0)*rho_br)

            call symm_2d_blending(q_prim_vf(momxb)%sf(i, j, 0), ux_tl, ux_tr, ux_bl, ux_br, 0.75_wp, 0.75_wp, 0.0025_wp, i, j)
            call symm_2d_blending(q_prim_vf(momxb + 1)%sf(i, j, 0), uy_tl, uy_tr, uy_bl, uy_br, 0.75_wp, 0.75_wp, 0.0025_wp, i, j)

            call symm_2d_blending(q_prim_vf(contxb)%sf(i, j, 0), rho_tl, rho_tr, rho_bl, rho_br, 0.75_wp, 0.75_wp, 0.0025_wp, i, j)
            q_prim_vf(contxb)%sf(i, j, 0) = q_prim_vf(contxb)%sf(i, j, 0) + 0.025 + 0.025*sin(2*pi*x_cc(i)*25)*sin(2*pi*y_cc(j)*25)

            call symm_2d_blending(q_prim_vf(E_idx)%sf(i, j, 0), p_tl, p_tr, p_bl, p_br, 0.75_wp, 0.75_wp, 0.0025_wp, i, j)
        end if
    case (285)
        ! This patch is hard-coded for 2D_IGR_double_mach case with Mach number set to 10
        gam = 1._wp + 1._wp/fluid_pp(1)%gamma

        if (patch_id == 1) then
            sb = y_cc(j) - tan(theta_dm)*(x_cc(i) - xr_dm)
            pshock = (2*gam*Mach**2 - (gam - 1._wp))/(gam + 1._wp)
            rhoshock = ((gam + 1._wp)*Mach**2)*rho0_dm/((gam - 1._wp)*Mach**2 + 2)
            velshock = (Mach - (Mach*rho0_dm/rhoshock))/2

            q_prim_vf(contxb)%sf(i, j, 0) = 0.5_wp*(1._wp + tanh(cf*sb))*rhoshock + 0.5_wp*(1._wp - tanh(cf*sb))*rho0_dm
            q_prim_vf(E_idx)%sf(i, j, 0) = 0.5_wp*(1._wp + tanh(cf*sb))*pshock + 0.5_wp*(1._wp - tanh(cf*sb))*p0_dm
            q_prim_vf(momxb)%sf(i, j, &
                      & 0) = 0.5_wp*(1._wp + tanh(cf*sb))*(velshock*tan(theta_dm)) + 0.5_wp*(1._wp - tanh(cf*sb))*u0_dm
            q_prim_vf(momxb + 1)%sf(i, j, 0) = 0.5_wp*(1._wp + tanh(cf*sb))*(-velshock) + 0.5_wp*(1._wp - tanh(cf*sb))*v0_dm
        end if
    case default
        if (proc_rank == 0) then
            call s_int_to_str(patch_id, iStr)
            call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
        end if
    end select
#:enddef
