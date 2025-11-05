module mod_functionsText
    use iso_fortran_env
    implicit none
    PRIVATE

    public parse_path, parse_value_after_equal, split_ws, to_lower
    public next_data_line, is_single_integer,strip_inline_comment

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

    subroutine next_data_line(iu, ios)
        implicit none
        integer, intent(in)  :: iu
        integer, intent(out) :: ios

        character(len=1024) :: line
        ios = 0

        do
            read(iu,'(A)',iostat=ios) line
            if (ios /= 0) return    ! EOF or error -> exit

            ! Trim whitespace
            line = trim(adjustl(line))

            ! Skip blank or comment-only lines
            if (line == '' .or. line(1:1) == '!') cycle

            ! We found a useful data line.
            ! Rewind the read pointer one line so caller can read it properly.
            backspace(iu)
            exit
        end do
    end subroutine next_data_line

    pure function to_lower(str) result(out)
        implicit none
        character(*), intent(in) :: str
        character(len(str))      :: out
        integer :: i, c

        do i = 1, len(str)
            c = iachar(str(i:i))
            if (c >= iachar('A') .and. c <= iachar('Z')) then
            out(i:i) = achar(c + 32)   ! convert ASCII upper → lower
            else
            out(i:i) = str(i:i)
            end if
        end do
    end function to_lower

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

    function is_single_integer(s)
        logical :: is_single_integer
        character(len=*), intent(in) :: s
        integer :: a, b, ios1, ios2
        character(len=:), allocatable :: t
        t = trim(adjustl(s))
        if (len_trim(t) == 0) then
        is_single_integer = .false.; return
        end if
        ! Try to read TWO ints: if that works, not a single integer line
        ios2 = 0; read(t,*,iostat=ios2) a, b
        if (ios2 == 0) then
        is_single_integer = .false.; return
        end if
        ! Try to read ONE int: if that works, treat as sentinel
        ios1 = 0; read(t,*,iostat=ios1) a
        is_single_integer = (ios1 == 0)
    end function is_single_integer

    subroutine strip_inline_comment(line)
        character(len=*), intent(inout) :: line
        integer :: b
        line = adjustl(line)
        b = index(line,'!')
        if (b == 1) then
        line = ''
        else if (b > 1) then
        line = trim(line(:b-1))
        end if
    end subroutine strip_inline_comment


end module mod_functionsText
