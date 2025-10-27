module mod_functionsText
    use iso_fortran_env
    implicit none
    PRIVATE

    public parse_path, parse_value_after_equal, split_ws

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

    subroutine parse_value_after_equal(line, outvalue)
        character(len=*), intent(in)  :: line
        character(len=*), intent(out) :: outvalue
        integer :: peq, i1, i2
        character(len=:), allocatable :: raw
        outvalue = ''
        peq = index(line, '=')
        if (peq == 0) return
        raw = adjustl( line(peq+1:) )
        i1 = 1
        do while (i1 <= len(raw) .and. raw(i1:i1) == ' ') ; i1=i1+1 ; end do
        i2 = len_trim(raw); if (i1>i2) return
        if (raw(i1:i1) == '"' .and. raw(i2:i2) == '"') then
            outvalue = raw(i1+1:i2-1)
        else
            outvalue = raw(i1:i2)
        end if
    end subroutine parse_value_after_equal


    subroutine split_ws(s, tokens, ntokens)
        implicit none
        character(len=*), intent(in)  :: s
        character(len=*), intent(out) :: tokens(:)
        integer,          intent(out) :: ntokens
        integer :: i, n, L, start, finish
        character(len=:), allocatable :: ss

        ss = adjustl(s)
        L = len_trim(ss)
        n = 0
        i = 1

        do while (i <= L)
            ! skip spaces
            do while (i <= L .and. ss(i:i) == ' ')
                i = i + 1
            end do
            if (i > L) exit

            start = i
            ! advance to next space
            do while (i <= L .and. ss(i:i) /= ' ')
                i = i + 1
            end do
            finish = i - 1

            n = n + 1
            if (n <= size(tokens)) then
                tokens(n) = ss(start:finish)
            end if
        end do

        ntokens = n
        ! blank out unused tokens (nice for debugging)
        do i = n+1, size(tokens)
            tokens(i) = ''
        end do
    end subroutine split_ws


end module mod_functionsText
