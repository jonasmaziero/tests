!-----------------------------------------------------------------------------------------------------------------------------------
program dynamics_tests

!call corr_ad() ! Tests for the dynamics of correlations under local amplitude damping: 
call bds_lad() ! Test for the dynamics of Bell-diagonal states under local aplitude damping channels

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine corr_ad()
implicit none
real(8) :: c1, c2, c3, p, dp, discord_tr_xs, discord_hs, discord_hsa, concurrence_2qb
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)
open(unit=13, file='dynad.dat', status='unknown')
 
 !c1 = 0.91d0; c2 = 0.5d0 ; c3 = -0.43d0
 c1 = 0.6d0 ; c2 = 0.6d0 ; c3 = -0.4d0 
 call rho_bds(c1, c2, c3, rho)

dp = 0.005d0 ;   p = 1.d-15 - dp
do
  p = p + dp ;   if( p > 1.d0 ) exit
  call rho_ad(p, rho, rhop)
  write(13,*) p, concurrence_2qb(rhop), discord_tr_xs('a', rhop), discord_hs('a', 2, 2, rhop), discord_hsa('a', 2, 2, rhop)
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine bds_lad()  ! Test for the dynamics of Bell-diagonal states under local aplitude damping channels
implicit none
real(8) :: r11i, r33i, r30i, r11, r30, r33, r33t, s11s, s30s, s33ts
real(8) :: Ec, concurrence_2qb, E2_cone, E2, En, negativity, D2, D1, discord_tr_xs, Dhs, discord_hs
complex(8) :: rho(1:4,1:4), rhop(1:4,1:4)
real(8) :: stokes_parameter
integer :: ord_pm(1:2)
real(8) :: p, dp
open(unit=13, file='bds_lad.dat', status='unknown')
open(unit=14, file='bds_lad_.dat', status='unknown')

!call rdm_bds(rho)  ! gets a random Bell-diagonal state
!ord_pm(1) = 1 ;   ord_pm(2) = 1 ;   r11 = stokes_parameter(rho, ord_pm, 2)
!ord_pm(1) = 3 ;   ord_pm(2) = 0 ;   r30 = stokes_parameter(rho, ord_pm, 2)
!ord_pm(1) = 3 ;   ord_pm(2) = 3 ;   r33 = stokes_parameter(rho, ord_pm, 2)
!write(*,*) 'r11=', r11, '  r30=', r30, '  r33=', r33

call rho_bds(0.5d0, 0.5d0, -0.6d0, rho) ;   call rho_tests(4, rho)
!call rho_bds(-0.5d0, -0.5d0, -0.6d0, rho) ;   call rho_tests(4, rho)
!call rho_bds(0.4d0, 0.4d0, -0.3d0, rho) ;   call rho_tests(4, rho)
!call rho_bds(-0.8d0, -0.8d0, -0.7d0, rho) ;   call rho_tests(4, rho)

dp = 0.005d0 ;   p = 1.d-15 - dp
do
  p = p + dp ;   if( p > 1.d0 ) exit
  call rho_ad(p, rho, rhop)
  Ec = concurrence_2qb(rhop) ;   En = negativity(2, 2, 'a', rhop) ;   D1 = discord_tr_xs('a', rhop)
  Dhs = discord_hs('a', 2, 2, rhop)
  ord_pm(1) = 1 ;   ord_pm(2) = 1 ;   r11 = stokes_parameter(rhop, ord_pm, 2)
  ord_pm(1) = 3 ;   ord_pm(2) = 0 ;   r30 = stokes_parameter(rhop, ord_pm, 2)
  ord_pm(1) = 3 ;   ord_pm(2) = 3 ;   r33 = stokes_parameter(rhop, ord_pm, 2) ;   r33t = (1.d0+r33)/2.d0
  !write(*,*) r11, r30, r33
  ! with min d_2(r,s)
  s11s = (r11/3.d0)*(1.d0+2.d0*r33t/sqrt(r11**2.d0+r30**2.d0))
  s30s = (r30/3.d0)*(1.d0+2.d0*r33t/sqrt(r11**2.d0+r30**2.d0))
  s33ts = r33t*(1.d0+2.d0*r33t/sqrt(r11**2.d0+r30**2.d0))/(-1.d0+4.d0*r33t/sqrt(r11**2.d0+r30**2.d0))
  if ( (1.d0+r33)**2.d0 < 4.d0*(r11**2.d0+r30**2.d0) ) then
    E2 = (1.d0/sqrt(2.d0))*sqrt((r11-s11s)**2.d0 + (r30-s30s)**2.d0 + 2.d0*(r33t-s33ts)**2.d0)
  else
    E2 = 0.d0
  endif
  !write(14,*) p, r11, r30, r33t, s11s, s30s, s33ts
  ! with min d(cone)
  s11s = r11*(1.d0/2.d0)*(1.d0+r33t/sqrt(r11**2.d0+r30**2.d0))
  s30s = r30*(1.d0/2.d0)*(1.d0+r33t/sqrt(r11**2.d0+r30**2.d0))
  s33ts = (1.d0/2.d0)*(sqrt(r11**2.d0+r30**2.d0)+r33t)
  if ( (1.d0+r33)**2.d0 < 4.d0*(r11**2.d0+r30**2.d0) ) then
    E2_cone = (1.d0/sqrt(2.d0))*sqrt((r11-s11s)**2.d0 + (r30-s30s)**2.d0 + 2.d0*(r33t-s33ts)**2.d0)
  else
    E2_cone = 0.d0
  endif
  !D2 = (1.d0/2.d0)*min(sqrt(2.d0)*abs(r11) , sqrt(r33**2.d0+2.d0*r30**2.d0))
  D2 = (1.d0/2.d0)*sqrt(2.d0)*abs(r11)
  write(13,*) p, E2, E2_cone, Ec, 2.d0*En, D2, D1, Dhs
enddo


end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine rho_tests(d, rho)
implicit none
integer :: d  ! Dimension of the density matrix
complex(8) :: rho(1:d,1:d)  ! The random density matrix
complex(8), allocatable :: rhoa(:,:)  ! For the adjoint of rho
real(8), allocatable :: W(:)  ! For the eigenvalues of rho
real(8) :: trace_he  ! For the trace function
real(8) :: norm_hs  ! For the Hilbert-Schmidt norm

allocate( rhoa(1:d,1:d), W(1:d) )

write(*,*) 'Verifying if the trace is equal to one'
write(*,*) 'Tr(rho) = ', trace_he(d, rho)  

write(*,*) 'Verifying Hermiticity'
call adjoint(d, d, rho, rhoa) ;   write(*,*) '||rho-rhoa||_2 = ', norm_hs(d, d, rho-rhoa) 

write(*,*) 'Verifying positivity'
call lapack_zheevd('N', d, rho, W) ;   write(*,*) 'Eigvals = ', W

deallocate( rhoa, W )

end
