#:def Hardcoded1DVariables()
    ! Place any declaration of intermediate variables here
#:enddef

#:def Hardcoded1D()
    select case (patch_icpp(patch_id)%hcid)
    case (170)  ! 1D profile from external data
        @: HardcodedReadValues()
    case (180)  ! Shu-Osher problem
        if (patch_id == 2) then
            q_prim_vf(contxb + 0)%sf(i, 0, 0) = 1 + 0.2*sin(5*x_cc(i))
        end if
    case (181)  ! Titarev-Torro problem
        q_prim_vf(contxb + 0)%sf(i, 0, 0) = 1 + 0.1*sin(20*x_cc(i)*pi)
    case default
        call s_int_to_str(patch_id, iStr)
        call s_mpi_abort("Invalid hcid specified for patch " // trim(iStr))
    end select
#:enddef
