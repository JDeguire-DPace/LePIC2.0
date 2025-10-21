module mod_functionsText
    use iso_fortran_env
    implicit none
    PRIVATE

    public parse_path

    contains
        subroutine parse_path(line, filename)
        character(len=*),  intent(in)  :: line
        character(len=*),  intent(out) :: filename
        integer :: pos
        pos = index(line,'=')
        if (pos > 0) then
        filename = adjustl(trim(line(pos+1:)))
        else
        filename = ''
        end if
    end subroutine parse_path
end module mod_functionsText
