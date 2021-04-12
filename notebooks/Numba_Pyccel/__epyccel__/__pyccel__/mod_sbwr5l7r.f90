module mod_sbwr5l7r

use, intrinsic :: ISO_C_BINDING

implicit none

contains

!........................................
function solve_1d_burger_pyccel(u, un, nt, nx, dt, dx, nu) result( &
      Out_0001)

  implicit none

  integer(C_INT32_T) :: Out_0001
  real(C_DOUBLE), intent(inout) :: u(0:)
  real(C_DOUBLE), intent(inout) :: un(0:)
  integer(C_INT32_T), value :: nt
  integer(C_INT32_T), value :: nx
  real(C_DOUBLE), value :: dt
  real(C_DOUBLE), value :: dx
  real(C_DOUBLE), value :: nu
  integer(C_INT32_T) :: n
  integer(C_INT32_T) :: i

  !fill the update of u
  do n = 0_C_INT32_T, nt-1_C_INT32_T, 1_C_INT32_T
    do i = 0_C_INT32_T, nx-1_C_INT32_T, 1_C_INT32_T
      un(i) = u(i)
    end do
    do i = 1_C_INT32_T, nx - 1_C_INT32_T-1_C_INT32_T, 1_C_INT32_T
      u(i) = un(i) - un(i) * dt / dx * (un(i) - un(i - 1_C_INT32_T)) + &
      nu * dt / dx ** 2_C_INT32_T * (un(i + 1_C_INT32_T) - 2_C_INT32_T &
      * un(i) + un(i - 1_C_INT32_T))
    end do
  end do
  Out_0001 = 0_C_INT32_T
  return

end function solve_1d_burger_pyccel
!........................................

end module mod_sbwr5l7r
