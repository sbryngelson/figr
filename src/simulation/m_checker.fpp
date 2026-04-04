
#:include 'macros.fpp'

module m_checker

    use m_global_parameters
    use m_mpi_common
    use m_helper_basic

    implicit none

    private; public :: s_check_inputs

contains

    impure subroutine s_check_inputs

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs

end module m_checker
