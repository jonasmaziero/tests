!-----------------------------------------------------------------------------------------------------------------------------------
program nonlocality_tests

call min_test() ;   stop ! Tests for the measurement induced nonlocality: ok

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine min_test()  ! Tests for the MEASUREMENT INDUCED NONLOCALITY
implicit none
integer :: nqb  ! For the number of qubits
real(8) :: nn  ! Number of qubits (real version)
integer :: d  ! For the dimensions
real(8) :: theta, dtheta  ! For the weight of W and GHZ, and its variation
complex(8), allocatable :: GHZ(:), W(:), psi(:), rho(:,:), rho_ab(:,:)
real(8) :: x  ! Auxiliary variable for cos(theta)
real(8) :: min_hs  ! For the ``numerical'' calculation of MIN
real(8) :: c1, c2, c3  ! Auxiliary variables for the enalytical formula
real(8) :: MINa, MINn  ! For the analytical and 'numerical' formulas of MIN
open(unit=13, file='min_hs_n10.dat', status='unknown')

! Tests for the superposition of n-qubit W and GHZ states
nqb = 10 ;   d = 2**nqb ;   allocate( GHZ(1:d), W(1:d), psi(1:d), rho(1:d,1:d), rho_ab(1:4,1:4) )
nn = dble(nqb) ; c1 = -nn**3.d0+6.d0*nn**2.d0-32.d0*nn+16.d0 ; c2 = 3.d0*nn**3.d0+16.d0*nn ; c3 = 3.d0*nn**3.d0+2.d0*nn**2.d0
dtheta = 1.d-2 ;   theta = 0.d0 - dtheta 
do
  theta = theta + dtheta ;   x = dcos(theta)
  MINa = ( c1*x**6.d0 + c2*x**4.d0 - c3*x**2.d0 + nn**3.d0 )/( 2.d0*(nn**2.d0)*( (nn**2.d0 - 6.d0*nn + 4.d0)*x**2.d0 + 2.d0*nn ) )
  call state_GHZ(nqb, GHZ) ;   call state_W(nqb, W) ;   psi = x*W + dsqrt(1.d0-x**2.d0)*GHZ ;   call projector(psi, d, rho)
  call partial_trace_b_he(rho, 2**2, 2**(nqb-2), rho_ab) ; MINn = min_hs('a', 2, 2, rho_ab) ; write(13,*) x**2.d0, MINn, MINa
  if(theta > 1.570796326794897d0) exit
enddo  
deallocate( GHZ, W, psi, rho, rho_ab )
         
end
!-----------------------------------------------------------------------------------------------------------------------------------