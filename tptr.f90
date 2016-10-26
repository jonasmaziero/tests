!###################################################################################################################################
!                                           Basic tests for the PARTIAL TRACE functions
!###################################################################################################################################
program tpartial_trace

call tptr()

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine tptr()
implicit none
integer :: d, da, db, dc  ! dimensions
complex(8), allocatable :: rho(:,:), rho_a(:,:), rho_c(:,:), rho_ac(:,:)  ! density matrix
complex(8), allocatable :: ham(:,:)  ! Hamiltonian
complex(8), allocatable :: A(:,:) ! matrix we send to Lapack
real(8), allocatable :: W(:)  ! eigenvalues of rho
real(8) :: h, dh  ! transverse magnetic field and its variation
integer :: ns  ! number of spins
integer :: j  ! for counters
real(8) :: SvN, neumann, coh, coh_l1n
real(8) :: Ca, Cc, Cac, Cnl
real(8) :: t1, t2  ! times
open(unit=13, file='Cnl10.dat', status='unknown')
!open(unit=14, file='dcoh2.dat', status='unknown')

call cpu_time(t1)

!dons: do ns = 3, 10
  ns = 10 ;  da = 2 ; dc = 2 ; db = 2**(ns-2);  d = 2**ns  
  allocate(ham(1:d,1:d), A(1:d,1:d), W(1:d), rho(1:d,1:d), rho_a(1:da,1:da), rho_c(1:dc,1:dc), rho_ac(1:da*dc,1:da*dc)) 
  dh = 5.d-3 ;   h = 0.2d0 - dh + 1.d-5
  doh: do
    !if (h <= 0.45d0 .or. h >= 0.55d0) then ;  dh = 1.d-2 ; else ; dh = 1.d-3 ; endif
       h = h + dh ;  if(h > 0.7d0) exit doh
       call Ising_1D(ns, h, ham) 
       !;   A = ham !;   call array_display_cr(4, 4, ham) ;   stop
       !call lapack_zheevd('V', d, A, W) ;   write(13,*) h, (W(j)-W(1), j=1,d)
       call state_thermal(1.d-15, ham, d, rho)
!call partial_trace_3(rho,da,db,dc,rho_ac); call partial_trace_a(rho_ac,da,dc,rho_c); call partial_trace_b(rho_ac,da,dc,rho_a) 
call ptr_3(rho, da, db, dc, rho_ac) ;    call ptr_a(rho_ac, da, dc, rho_c) ;   call ptr_b(rho_ac, da, dc, rho_a)
       Ca = coh_l1n(da, rho_a) ;   Cc = coh_l1n(dc, rho_c) ;   Cac = coh_l1n(da*dc, rho_ac) ;   Cnl = Cac - Ca - Cc
       write(13,*) 1.d0/h, Cnl, Ca, Cc, Cac  
  enddo doh
  deallocate(ham, A, W, rho, rho_a, rho_c, rho_ac)
!enddo dons

call cpu_time(t2) ; write(*,*) 'ns=',ns,'  t=', t2-t1

end
!###################################################################################################################################