
module m_compile_specific

    implicit none

contains

    impure subroutine s_create_directory(dir_name)

        character(LEN=*), intent(in) :: dir_name

#ifdef _WIN32
        call system('mkdir "' // dir_name // '" 2> NUL')
#else
        call system('mkdir -p "' // dir_name // '"')
#endif

    end subroutine s_create_directory

    impure subroutine s_delete_directory(dir_name)

        character(LEN=*), intent(in) :: dir_name

#ifdef _WIN32
        call system('rmdir "' // dir_name // '" /s /q')
#else
        call system('rm -r "' // dir_name // '"')
#endif

    end subroutine s_delete_directory

    impure subroutine my_inquire(fileloc, dircheck)

        character(LEN=*), intent(in) :: fileloc
        logical, intent(inout)       :: dircheck

#ifdef __INTEL_COMPILER
        inquire (DIRECTORY=trim(fileloc), EXIST=dircheck)
#else
        inquire (FILE=trim(fileloc), EXIST=dircheck)
#endif

    end subroutine my_inquire

end module m_compile_specific
