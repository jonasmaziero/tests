!-----------------------------------------------------------------------------------------------------------------------------------
program entanglement_tests

!call test_eentropy()
call test_schcoeff()

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine test_eentropy()  ! Tests for the entanglement entropy function
implicit none
real(8) :: w, dw  ! Weight for controling the degree of entanglement of a pure state, and the step for its variation
complex(8) :: psi(1:4)  ! The state vector
real(8) :: eentropy  ! For the entanglement entropy function
open(unit=13, file='ee.dat', status='unknown')

dw = 0.01d0 ;   w = 0 - dw
do
  w = w + dw ;   psi(1) = w ;   psi(2) = 0.d0 ;   psi(3) = 0.d0 ;   psi(4) = 1.d0 - w
  write(13,*) w, eentropy(psi, 4, 2, 2, 'a')
  if ( w > 1.d0 ) exit
enddo


end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine test_schcoeff()  ! Tests for the Schmidt decomposition subroutine
complex(8) :: psi(1:4)  ! The bipartite state vector
real(8) :: schcoeff(1:2)  ! For the Schmidt coefficients
complex(8) :: eigvec_a(1:2,1:2), eigvec_b(1:2,1:2)  ! For the eigenvectors of the reduced density matrices
real(8) :: w  ! Weight for controling the degree of entanglement of a pure state

w = 0.5d0
psi(1) = sqrt(w) ;   psi(2) = 0.d0 ;   psi(3) = 0.d0 ;   psi(4) = sqrt(1.d0 - w)
call schmidt_coefficients(psi, 4, 2, 2, schcoeff, eigvec_a, eigvec_b)
write(*,*) 'schcoeff = ', schcoeff
write(*,*) 'eigvec_a = ' ;   call array_display_cr(2, 2, eigvec_a)
write(*,*) 'eigvec_b = ' ;   call array_display_cr(2, 2, eigvec_b)


end
!-----------------------------------------------------------------------------------------------------------------------------------