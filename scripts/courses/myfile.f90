subroutine myfunc(vari, varo, nx)

  real(8),intent(in)  :: vari(nx)
  integer, intent(in) :: nx
  real(8),intent(out) :: varo(nx)
  varo = vari**2 + 3

end subroutine myfunc

