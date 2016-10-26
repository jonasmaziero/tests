!###################################################################################################################################
!                                           Basic tests for the PARTIAL TRANSPOSE functions
!###################################################################################################################################
program tptransp

call tptranspose()
call tptranspose_gen()

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine tptranspose()
implicit none
complex(8) :: psi_p(1:4), psi_m(1:4), phi_p(1:4), phi_m(1:4), A(1:4,1:4)
complex(8) :: rho(1:4,1:4), rho_pt(1:4,1:4), proj(1:4,1:4), identity(1:4,1:4)
real(8) :: w, dw, negativity, log_negativity, Ehs, En, Eln, arg
open(unit=13, file='werner2.dat', status='unknown')

call bell_basis(psi_p, psi_m, phi_p, phi_m) ;   call identity_c(4, identity) ;   call projector(psi_m, 4, proj)
dw = 1.d0/100.d0 ;   w = 0.d0 - dw  ! Sets the width of the step in w
do
  w = w + dw ;   if( w > 1.d0 ) exit ;   rho = (1.d0-w)*(identity/4.d0) + w*proj
  call partial_transpose_b(2, 2, rho, rho_pt) !;   Eln = log_negativity(4, rho_pt)
  A = rho_pt ;   call entanglement_hs(4, A, Ehs, 'n')
  write(13,*) w, negativity(4, rho_pt), Ehs
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine tptranspose_gen()   ! Tests for a 3-qubit Werner state
implicit none
complex(8) :: GHZ(1:8), rho(1:8,1:8), rho_pt(1:8,1:8), proj(1:8,1:8), identity(1:8,1:8), A(1:8,1:8)
real(8) :: w, dw, negativity, log_negativity, Ehs, En, Eln
integer :: di(1:3), ssys(1:3)
open(unit=13, file='werner3.dat', status='unknown')

di(1) = 2 ; di(2) = 2 ; di(3) = 2  ! Dimension of the sub-systems
ssys(1) = 1 ; ssys(2) = 1 ;  ssys(3) = 0  ! determines over which sub-systems to take the PT
call state_GHZ(3, GHZ) ;   call identity_c(8, identity) ;   call projector(GHZ, 8, proj)
dw = 1.d0/100.d0 ;   w = 0.d0 - dw  ! Sets the width of the step in w
do
  w = w + dw ;   if( w > 1.d0 ) exit ;   rho = (1.d0-w)*(identity/8.d0) + w*proj
  call partial_transpose(8, rho, rho_pt, 3, di, ssys) !;   En = negativity(8, rho_pt) !;   Eln = log_negativity(8, rho_pt)
  A = rho_pt ;   call entanglement_hs(8, A, Ehs, 'y') ;   write(13,*) w, negativity(8, rho_pt), Ehs
enddo
call partial_transpose(8, A, rho_pt, 3, di, ssys)

end
!###################################################################################################################################