
#:include 'macros.fpp'

module m_checker_common

    use m_global_parameters
    use m_mpi_proxy
    use m_helper_basic

    implicit none

    private; public :: s_check_inputs_common, wp

contains

    impure subroutine s_check_inputs_common

#ifndef FIGR_SIMULATION
        call s_check_total_cells
#endif

    end subroutine s_check_inputs_common

#ifndef FIGR_SIMULATION
    impure subroutine s_check_total_cells

        character(len=18) :: numStr
        integer(kind=8)   :: min_cells

        min_cells = int(2, kind=8)**int(min(1, m) + min(1, n) + min(1, p), kind=8)*int(num_procs, kind=8)
        call s_int_to_str(2**(min(1, m) + min(1, n) + min(1, p))*num_procs, numStr)

        @:PROHIBIT(nGlobal < min_cells, &
                   & "Total number of cells must be at least (2^[number of dimensions])*num_procs, " // "which is currently " &
                   & // trim(numStr))

    end subroutine s_check_total_cells
#endif

end module m_checker_common
