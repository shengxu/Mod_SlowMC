module BiSearch

contains

    subroutine binarySearch (a, value, ind, stat)
    implicit none

        real(8), intent(in), dimension(:), target :: a
        real(8), intent(in)         :: value
        integer, intent(out)        :: ind
        logical, intent(out)        :: stat  !stat: true if value falls on a grid point
        real(8), pointer            :: p(:)
        integer                  :: n, mid, offset

!        print *, a
        stat = .false.
        p => a
        n = size(p)
!        write(11, *) n
        ind = 0
        if (value < p(1))  then
            ind = 0
            return
        else if (value > p(n))  then
            ind = n
            return
        else if (value == p(1)) then
            ind = 1
            stat = .true.
            return
        else if (value == p(n)) then
            ind = n
            stat = .true.
            return
        else
            offset = 0
            mid = 0
            do while (size(p) > 0)
                if (size(p) == 2) then
                    ind = offset + 1
                    return
                else
                    mid = size(p)/2 + 1
                    if (p(mid) > value) then
                        p => p(:mid)
                    else if (p(mid) < value) then
                        offset = offset + mid - 1
                        p => p(mid:)
                    else
                        ind = offset + mid   ! SUCCESS!!
                        stat = .true.
                        return
                    end if
                end if
            end do
        end if
    end subroutine binarySearch

end module BiSearch
