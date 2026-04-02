#:include 'case.fpp'
#:include 'macros.fpp'

!>
!! @file
!! @brief Contains module m_helper

!> @brief Utility routines for coordinate transforms, array sampling, and special functions
module m_helper

    use m_derived_types
    use m_global_parameters
    use ieee_arithmetic  !< For checking NaN

    implicit none

    private
    public :: s_int_to_str, f_xor, f_logical_to_int, associated_legendre, real_ylm, &
        & double_factorial, factorial, f_cut_on, f_cut_off, s_downsample_data, s_upsample_data

contains

    !> Convert an integer to its trimmed string representation.
    elemental subroutine s_int_to_str(i, res)

        integer, intent(in)             :: i
        character(len=*), intent(inout) :: res

        write (res, '(I0)') i
        res = trim(res)

    end subroutine s_int_to_str

    !> Perform XOR on lhs and rhs.
    elemental function f_xor(lhs, rhs) result(res)

        logical, intent(in) :: lhs, rhs
        logical             :: res

        res = (lhs .and. .not. rhs) .or. (.not. lhs .and. rhs)

    end function f_xor

    !> Convert a logical to 1 or 0.
    elemental function f_logical_to_int(predicate) result(int)

        logical, intent(in) :: predicate
        integer             :: int

        if (predicate) then
            int = 1
        else
            int = 0
        end if

    end function f_logical_to_int

    !> Real spherical harmonic Y_lm(theta, phi). theta = polar angle from +z (acos(z/r)), phi = atan2(y,x). Uses associated Legendre
    !! P_l^|m|(cos theta). Standard normalisation.
    function real_ylm(theta, phi, l, m) result(Y)

        integer, intent(in)  :: l, m
        real(wp), intent(in) :: theta, phi
        real(wp)             :: Y, x, prefac
        integer              :: m_abs

        m_abs = abs(m)
        if (m_abs > l) then
            Y = 0._wp
            return
        end if
        x = cos(theta)
        prefac = sqrt((2*l + 1)*real(factorial(l - m_abs), wp)/real(factorial(l + m_abs), wp)/(4._wp*pi))
        if (m == 0) then
            Y = prefac*associated_legendre(x, l, 0)
        else if (m > 0) then
            Y = prefac*sqrt(2._wp)*associated_legendre(x, l, m_abs)*cos(m*phi)
        else
            Y = prefac*sqrt(2._wp)*associated_legendre(x, l, m_abs)*sin(m_abs*phi)
        end if

    end function real_ylm

    !> Associated Legendre polynomial P_l^m(x) (Ferrers function, Condon-Shortley phase). Valid for integer l >= 0, 0 <= m <= l, and
    !! x in [-1,1]. Returns 0 for |m| > l or l < 0. Formulas: DLMF 14.10.3 (recurrence in degree), Wikipedia "Associated Legendre
    !! polynomials" (P_l^l and P_l^{l-1} identities). Recurrence: (l-m)P_l^m = (2l-1)x P_{l-1}^m - (l+m-1)P_{l-2}^m.
    !! @param x argument (typically cos(theta)), should be in [-1,1]
    !! @param l degree (>= 0)
    !! @param m_order order (0 <= m_order <= l)
    recursive function associated_legendre(x, l, m_order) result(result_P)

        integer, intent(in)  :: l, m_order
        real(wp), intent(in) :: x
        real(wp)             :: result_P
        real(wp)             :: one_minus_x2

        ! Out-of-domain: P_l^m = 0 for |m| > l or l < 0 (standard convention)

        if (l < 0 .or. m_order < 0 .or. m_order > l) then
            result_P = 0._wp
            return
        end if

        if (m_order <= 0 .and. l <= 0) then
            result_P = 1._wp
        else if (l == 1 .and. m_order <= 0) then
            result_P = x
        else if (l == 1 .and. m_order == 1) then
            one_minus_x2 = max(0._wp, 1._wp - x**2)
            result_P = -sqrt(one_minus_x2)
        else if (m_order == l) then
            ! P_l^l(x) = (-1)^l (2l-1)!! (1-x^2)^(l/2). Use real exponent for odd l
            one_minus_x2 = max(0._wp, 1._wp - x**2)
            result_P = (-1)**l*real(double_factorial(2*l - 1), wp)*one_minus_x2**(0.5_wp*real(l, wp))
        else if (m_order == l - 1) then
            result_P = x*(2*l - 1)*associated_legendre(x, l - 1, l - 1)
        else
            result_P = ((2*l - 1)*x*associated_legendre(x, l - 1, m_order) - (l + m_order - 1)*associated_legendre(x, l - 2, &
                        & m_order))/(l - m_order)
        end if

    end function associated_legendre

    !> Calculate the double factorial of an integer
    elemental function double_factorial(n_in) result(R_result)

        integer, intent(in)      :: n_in
        integer, parameter       :: int64_kind = selected_int_kind(18)  !< 18 bytes for 64-bit integer
        integer(kind=int64_kind) :: R_result
        integer                  :: i

        R_result = product((/(i, i=n_in, 1, -2)/))

    end function double_factorial

    !> Calculate the factorial of an integer
    elemental function factorial(n_in) result(R_result)

        integer, intent(in)      :: n_in
        integer, parameter       :: int64_kind = selected_int_kind(18)  !< 18 bytes for 64-bit integer
        integer(kind=int64_kind) :: R_result
        integer                  :: i

        R_result = product((/(i, i=n_in, 1, -1)/))

    end function factorial

    !> Calculate a smooth cut-on function that is zero for x values smaller than zero and goes to one, for generating smooth initial
    !! conditions
    function f_cut_on(x, eps) result(fx)

        real(wp), intent(in) :: x, eps
        real(wp)             :: fx

        fx = 1 - f_gx(x/eps)/(f_gx(x/eps) + f_gx(1 - x/eps))

    end function f_cut_on

    !> Calculate a smooth cut-off function that is one for x values smaller than zero and goes to zero, for generating smooth
    !! initial conditions
    function f_cut_off(x, eps) result(fx)

        real(wp), intent(in) :: x, eps
        real(wp)             :: fx

        fx = f_gx(x/eps)/(f_gx(x/eps) + f_gx(1 - x/eps))

    end function f_cut_off

    !> Helper function for f_cut_on and f_cut_off
    function f_gx(x) result(gx)

        real(wp), intent(in) :: x
        real(wp)             :: gx

        if (x > 0) then
            gx = exp(-1._wp/x)
        else
            gx = 0._wp
        end if

    end function f_gx

    !> Downsample conservative variable fields by a factor of 3 in each direction using volume averaging.
    subroutine s_downsample_data(q_cons_vf, q_cons_temp, m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds)

        type(scalar_field), dimension(sys_size), intent(inout) :: q_cons_vf, q_cons_temp

        ! Down sampling variables
        integer                :: i, j, k, l
        integer                :: ix, iy, iz, x_id, y_id, z_id
        integer, intent(inout) :: m_ds, n_ds, p_ds, m_glb_ds, n_glb_ds, p_glb_ds

        m_ds = int((m + 1)/3) - 1
        n_ds = int((n + 1)/3) - 1
        p_ds = int((p + 1)/3) - 1

        m_glb_ds = int((m_glb + 1)/3) - 1
        n_glb_ds = int((n_glb + 1)/3) - 1
        p_glb_ds = int((p_glb + 1)/3) - 1

        do l = -1, p_ds + 1
            do k = -1, n_ds + 1
                do j = -1, m_ds + 1
                    x_id = 3*j + 1
                    y_id = 3*k + 1
                    z_id = 3*l + 1
                    do i = 1, sys_size
                        q_cons_temp(i)%sf(j, k, l) = 0

                        do iz = -1, 1
                            do iy = -1, 1
                                do ix = -1, 1
                                    q_cons_temp(i)%sf(j, k, l) = q_cons_temp(i)%sf(j, k, &
                                                & l) + (1._wp/27._wp)*q_cons_vf(i)%sf(x_id + ix, y_id + iy, z_id + iz)
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do

    end subroutine s_downsample_data

    !> Upsample conservative variable fields from a coarsened grid back to the original resolution using interpolation.
    subroutine s_upsample_data(q_cons_vf, q_cons_temp)

        type(scalar_field), intent(inout), dimension(sys_size) :: q_cons_vf, q_cons_temp
        integer                                                :: i, j, k, l
        integer                                                :: ix, iy, iz
        integer                                                :: x_id, y_id, z_id
        real(wp), dimension(4)                                 :: temp

        do l = 0, p
            do k = 0, n
                do j = 0, m
                    do i = 1, sys_size
                        ix = int(j/3._wp)
                        iy = int(k/3._wp)
                        iz = int(l/3._wp)

                        x_id = j - int(3*ix) - 1
                        y_id = k - int(3*iy) - 1
                        z_id = l - int(3*iz) - 1

                        temp(1) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy, iz) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, iy, iz)
                        temp(2) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy + y_id, iz) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, &
                             & iy + y_id, iz)
                        temp(3) = (2._wp/3._wp)*temp(1) + (1._wp/3._wp)*temp(2)

                        temp(1) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy, iz + z_id) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, &
                             & iy, iz + z_id)
                        temp(2) = (2._wp/3._wp)*q_cons_temp(i)%sf(ix, iy + y_id, &
                             & iz + z_id) + (1._wp/3._wp)*q_cons_temp(i)%sf(ix + x_id, iy + y_id, iz + z_id)
                        temp(4) = (2._wp/3._wp)*temp(1) + (1._wp/3._wp)*temp(2)

                        q_cons_vf(i)%sf(j, k, l) = (2._wp/3._wp)*temp(3) + (1._wp/3._wp)*temp(4)
                    end do
                end do
            end do
        end do

    end subroutine s_upsample_data

end module m_helper
