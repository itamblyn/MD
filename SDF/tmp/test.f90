program test
  use vectorfun ! this makes all declarations and definitions of module available
  implicit none
  integer :: i,n
  real    :: r 
  real, allocatable :: x(:)
  print *,'n ='
  read *, n
  allocate(x(n))
  do i=1,n
    print *,'x(',i,') ='
    read *, x(i)
  end do
  r=length(x)
  print *, 'length =',r
end program
