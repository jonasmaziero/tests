!###################################################################################################################################
!                                           Basic tests for the PARTIAL TRACE functions
!###################################################################################################################################
program thamiltonian

call tIsing1D()

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine tIsing1D()
implicit none
integer :: d  ! Dimension of the density matrix and hamiltonian
complex(8), allocatable :: rho(:,:)  ! density matrix
complex(8), allocatable :: ham(:,:)  ! Hamiltonian
complex(8), allocatable :: A(:,:) ! matrix we send to Lapack
real(8), allocatable :: W(:)  ! eigenvalues of rho
real(8) :: h, dh  ! transverse magnetic field and its variation
integer :: ns  ! number of spins
integer :: j  ! for counters
open(unit=13, file='energies.dat', status='unknown')

  ns = 4
  d = 2**ns ;   allocate(ham(1:d,1:d), A(1:d,1:d), W(1:d)) 
  dh = 1.d-2 ;   h = 0.d0 - dh + 1.d-15
  doh: do
    h = h + dh ;  if(h > 2.d0) exit doh
    call Ising_1D(ns, h, ham) ;   A = ham !;   call array_display_cr(4, 4, ham) ;   stop
    call lapack_zheevd('V', d, A, W) ;   write(13,*) h, (W(j), j=1,d)
  enddo doh
  deallocate(ham, A, W)

end
!###################################################################################################################################
