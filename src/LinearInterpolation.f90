module LinearInterpolation
    use BiSearch, only: binarySearch
    implicit none

contains

    subroutine LinInterp (x, y, x_int, y_int)

        real(8), intent(in)         :: x(:)
        real(8), intent(in)         :: y(:)
        real(8), intent(in)         :: x_int
        real(8), intent(out)        :: y_int
        integer                     :: n
        integer                     :: ind
        logical                     :: stat  !stat: true if value falls on a grid point

        n = size(x)
        if (x_int < x(1)) then
          y_int = y(1)+(y(2)-y(1))/(x(2)-x(1))*(x_int-x(1))
        else if (x_int == x(1)) then
          y_int = y(1)
        else if (x_int > x(n)) then
          y_int = y(n)+(y(n)-y(n-1))/(x(n)-x(n-1))*(x_int-x(n))
        else if (x_int == x(n)) then
          y_int = y(n)    
        else
          call binarySearch(x, x_int, ind, stat)
          if (stat) then  ! x_int falls on a grid point
            y_int = y(ind)
          else            ! x_int falls in the interval
            y_int = y(ind)+(y(ind+1)-y(ind))/(x(ind+1)-x(ind))*(x_int-x(ind))
          end if
        end if      

    end subroutine LinInterp

end module LinearInterpolation
