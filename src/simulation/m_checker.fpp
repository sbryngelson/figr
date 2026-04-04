!!@file
!!@brief Contains module m_checker

#:include 'macros.fpp'
#:include 'case.fpp'

module m_checker

    use m_global_parameters
    use m_mpi_proxy
    use m_helper
    use m_helper_basic

    implicit none

    private; public :: s_check_inputs

contains

    !> Checks compatibility of parameters in the input file. Used by the simulation stage
    impure subroutine s_check_inputs

        call s_check_inputs_compilers
        call s_check_inputs_time_stepping

    end subroutine s_check_inputs

    !> Checks constraints on compiler options
    impure subroutine s_check_inputs_compilers

    end subroutine s_check_inputs_compilers

    !> Checks constraints on time stepping parameters
    impure subroutine s_check_inputs_time_stepping

        if (.not. cfl_dt) then
            @:PROHIBIT(dt <= 0)
        end if

    end subroutine s_check_inputs_time_stepping

end module m_checker
