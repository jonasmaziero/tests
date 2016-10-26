!-----------------------------------------------------------------------------------------------------------------------------------
program discord_tests

!call gellmann_test() ;   stop  ! Tests for the generalized Gell Mann matrices: ok
!call bloch_gellmann_test() ;   stop  ! Tests for the Bloch vector with GGMM: ok
!call corrmat_gellmann_test() ;   stop ! Tests for the correlation matrix with GGMM: ok
call discord_hsa_test() ;   stop ! Tests for the ameliorated 2-norm discord: ok
!call discord_mid_test() ;   stop ! Tests for the discord MID: ok
!call discord_easy_test() ;   stop ! Tests for the 'easy' discord: ?

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine gellmann_test()  ! ok
implicit none
integer :: d  ! Dimension of SU(d)
complex(8), allocatable :: ggmm(:,:)

d = 4 ;   allocate( ggmm(1:d,1:d) )
call gellmann(d, 1, 3, 4, ggmm) !gellmann(d, group, k, l, ggmm) 
call array_display(d, d, ggmm)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine bloch_gellmann_test()  ! ok
implicit none
integer :: d  ! Dimension of SU(d)
character(10), dimension(5) :: optg  ! Options for the generators
complex(8), allocatable :: rho(:,:)  ! For the random density matrix to be used in the tests
real(8), allocatable :: A(:,:)  ! Auxiliary matrix for visual comparison of the Bloch vectors
real(8), allocatable :: bv(:)  ! For the Bloch vectors
integer :: j ! Ayxiliary variable for counters
integer :: dd  ! Auxiliary variable for the dimensions

d = 4 ;   dd = d**2-1 ;   allocate( rho(1:d,1:d), bv(1:dd), A(1:2,1:dd) )
optg = 'std' ;   call rng_init(optg) ;   call rdmg(optg, d, rho)  ! Gets a random density matrix
!write(*,*) "Density matrix" ;   call array_display(d, d, rho)

call bloch_vector_gellmann(d, rho, bv) ;   do j = 1, dd ;   A(1,j) = bv(j) ;   enddo  ! Optmized calculation
call bloch_vector_gellmann_unopt(d, rho, bv) ;   do j = 1, dd ;   A(2,j) = bv(j) ;   enddo  ! Un-optimized calculation

write(*,*) "Bloch's vectors" ;   call array_display_rr(2, dd, A)

deallocate(rho, bv, A)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine corrmat_gellmann_test()  ! ok
implicit none
integer :: d, da, db  ! Dimensions
character(10), dimension(5) :: optg  ! Options for the generators
complex(8), allocatable :: rho(:,:)  ! For the random density matrix to be used in the tests
real(8), allocatable :: A(:,:)  ! Auxiliary matrix for visual comparison of the Bloch vectors
real(8), allocatable :: cm1(:,:), cm2(:,:)  ! For the Bloch vectors
integer :: j ! Ayxiliary variable for counters
integer :: dda, ddb  ! Auxiliary variable for the dimensions

da = 4 ;   dda = da**2-1 ;   db = da ; ddb = db**2-1 ;   d = da*db ;  allocate( rho(1:d,1:d), cm1(1:dda,1:ddb), cm2(1:dda,1:ddb) )
optg = 'std' ;   call rng_init(optg) ;   call rdmg(optg, d, rho)  ! Gets a random density matrix
!write(*,*) "Density matrix" ;   call array_display(d, d, rho)

call corrmat_gellmann(da, db, rho, cm1)  ! Optmized calculation
call corrmat_gellmann_unopt(da, db, rho, cm2)  ! Un-optimized calculation

write(*,*) "cm_opt=" ;   call array_display_rr(dda, ddb, cm1)
write(*,*) "cm_unopt=" ;   call array_display_rr(dda, ddb, cm2)
write(*,*) "cm_opt - cm_unopt =" ;   call array_display_rr(dda, ddb, cm1-cm2)

deallocate(rho, cm1, cm2)

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine discord_hsa_test()  ! ok
implicit none
real(8) :: x, dx  ! Parameter for the 'mixedness' of the Werner state, and for its variation
integer :: d, da, db  ! For the dimension of the whole Hilbert space, and of the subsystems
complex(8), allocatable :: rho(:,:)  ! For the composite density matrix
real(8) :: adhsa  ! For the analytical ameliorated Hilbert-Schmidt discord
real(8) :: ndhsa, discord_hsa, discord_hsa_unopt, ndhs, discord_hs, discord_hs_2qb  ! For the numerical (ameliorated) Hilbert-Schmidt discord
real(8) :: t1, t2  ! For the cpu times
open(unit=13, file='dhsa_werner_d4_uo.dat', status='unknown')

call cpu_time(t1)

d = 4 ;   allocate(rho(1:d,1:d)) ;   da = anint(sqrt(dble(d))) ;   db = da ;   write(*,*) 'd, da:', d, da
dx = 0.05d0 ;   x = -1.d0 - dx
do
  x = x + dx
  call state_werner(da, x, rho)
  adhsa = ((x*dble(da) - 1.d0)**2.d0)/(dble((da-1)*(da+1)**2))  ! For Werner states
  !call state_isotropic(da, x, rho)
  !adhsa = ((x*dble(da*da) - 1.d0)**2.d0)/(dble(da*(da-1)*(da+1)**2))  ! For Isotropic states
  !write(*,*) "Eigenvalues if Xi:", ((x*dble(da)-1.d0)**2.d0)/(dble(da**2)*(dble(da**2)-1.d0)**2.d0)  ! ok
  !ndhsa = discord_hsa('a', da, db, rho) !;   ndhs = discord_hs('a', da, db, rho)
  ndhsa = discord_hsa_unopt('a', da, db, rho)
  write(13,*) x, adhsa, ndhsa
  if(x > 1.d0) exit
enddo

deallocate( rho )

call cpu_time(t2) ;   write(*,*) 't2-t1=',t2-t1

end
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
real(8) function discord_hsa_unopt(ssys, da, db, rho)  ! Returns the AMENDED HILBERT-SCHMIDT discord (unopt calc)
! Ref: S. J. Akhtarshenas et al., Computable measure of quantum correlation, QIP 14, 247 (2015).
implicit none
character(1), intent(in) :: ssys  ! Tells if sub-system a or b is classical one, in the minimization
integer, intent(in) :: da, db  ! Dimension of the subsystems
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! The bipartite density matrix, represented in the global computational basis
complex(8), allocatable :: rho_a(:,:), rho_b(:,:)  ! The reduced states of the subsystems
real(8), allocatable :: bv(:)  ! For the Bloch vector, of rho_a or of rho_b
real(8), allocatable :: proj_bv(:,:)  ! For the projector on the Bloch vector: x*x^t
real(8), allocatable :: corrmat(:,:)  ! For the correlation matrix
real(8), allocatable :: Xi(:,:)  ! For the S-correlation matrix, with S = A, B
real(8), allocatable :: W(:)  ! For the eigenvalues of e.g. Xi = (2/da*db)*( a*a^t + (2/db)*C*C^t)
real(8) :: purity  ! For the purity function
integer :: dda, ddb  ! Auxiliary variable for the dimensions

dda = da**2 - 1 ;   ddb = db**2 - 1

if (ssys == 'a' ) then ! CQ states are 'classical'
  allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a)
  allocate( bv(1:dda) ) ;   call bloch_vector_gellmann_unopt(da, rho_a, bv) ;   deallocate( rho_a )  ! Computes the Bloch vector
  allocate( proj_bv(1:dda,1:dda) ) ;   call projector_re(bv, dda, proj_bv) ;   deallocate( bv )  
  allocate( corrmat(1:dda,1:ddb) ) ;   call corrmat_gellmann_unopt(da, db, rho, corrmat)  ! Computes the correlation matrix
  ! A-correlation matrix
  allocate( Xi(1:dda,1:dda) ) ;   Xi = (2.d0/dble(da*da*db))*( proj_bv + (2.d0/dble(db))*matmul(corrmat,transpose(corrmat)) )
  deallocate( proj_bv, corrmat )
  allocate( W(1:dda) ) ;   call lapack_dsyevd('N', dda, Xi, W) ;   deallocate( Xi )
  allocate( rho_b(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b)
  discord_hsa_unopt = (1.d0/purity(db, rho_b))*sum(W(1:(da*(da-1)))) ;   deallocate( W, rho_b )  ! for the ameliorated 2-norm discord
  
else if ( ssys == 'b' ) then  ! QC states are 'classical'
  allocate( rho_b(1:db,1:db) ) ;   call partial_trace_a_he(rho, da, db, rho_b)
  allocate( bv(1:ddb) ) ;   call bloch_vector_gellmann_unopt(db, rho_b, bv) ;   deallocate( rho_b )  ! Computes the Bloch vector
  allocate( proj_bv(1:ddb,1:ddb) ) ;   call projector_re(bv, ddb, proj_bv) ;   deallocate( bv ) 
  allocate( corrmat(1:dda,1:ddb) ) ;   call corrmat_gellmann_unopt(da, db, rho, corrmat)  ! Computes the correlation matrix
  ! B-correlation matrix
  allocate( Xi(1:ddb,1:ddb) ) ;   Xi = (2.d0/dble(da*db*db))*( proj_bv + (2.d0/dble(da))*matmul(transpose(corrmat),corrmat) ) 
  deallocate( proj_bv, corrmat )
  allocate( W(1:ddb) ) ;   call lapack_dsyevd('N', ddb, Xi, W) ;   deallocate( Xi )
  allocate( rho_a(1:da,1:da) ) ;   call partial_trace_b_he(rho, da, db, rho_a)
  discord_hsa_unopt = (1.d0/purity(da, rho_a))*sum(W(1:(db*(db-1)))) ;   deallocate( W, rho_a )  ! for the ameliorated 2-norm discord
endif

end
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine discord_mid_test()  ! ok
implicit none
real(8) :: x, dx  ! Parameter for the 'mixedness' of the Werner state, and for its variation
integer :: d, da, db  ! For the dimension of the whole Hilbert space, and of the subsystems
complex(8), allocatable :: rho(:,:)  ! For the composite density matrix
real(8) :: admid  ! For the analytical discord mid
real(8) :: ndmid, discord_mid  ! For the numerical discord MID
real(8) :: log2  ! For the log base two function
open(unit=13, file='dmid.dat', status='unknown')

d = 36 ;   allocate(rho(1:d,1:d)) ;   da = anint(sqrt(dble(d))) ;   db = da ;   write(*,*) 'd, da:', d, da
dx = 0.005d0 ;   x = -1.d0 + 1.d-15 - dx
do
  x = x + dx
  call state_werner(da, x, rho)  ! For Werner states
  admid = ((1.d0-x)/2.d0)*log2((1.d0-x)/(dble(da*(da-1)))) & 
          + (((1.d0+x)*(dble(da-1)))/(dble(2*(da+1))))*log2((1.d0+x)/(dble(da*(da+1)))) &
          - ((dble(da)-x)/(dble(da+1)))*log2((dble(da)-x)/(dble(da*(da**2-1))))
  !call state_isotropic(da, x, rho)  ! For Isotropic states
  !admid = ((1.d0-x)/(dble(da)+1.d0))*log2((1.d0-x)/(dble(da*da)-1.d0)) + x*log2(x) & 
  !        - ((1.d0+x*dble(da))/(dble(da)+1.d0))*log2((1.d0+x*dble(da))/(dble(da*(da+1))))
  ndmid = discord_mid(da, db, rho)
  write(13,*) x, admid, ndmid+0.01d0
  if(x > 1.d0) exit
enddo

deallocate( rho )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine discord_easy_test()  ! ?
implicit none
real(8) :: p, dp  ! Parameter for the 'mixedness' of state
integer :: d, da, db  ! For the dimension of the whole Hilbert space, and of the subsystems
complex(8), allocatable :: rho(:,:)  ! For the composite density matrix
real(8) :: adb  ! For the analytical block discord
real(8) :: ndb, discord_easy  ! For the numerical block discord
open(unit=13, file='deasy.dat', status='unknown')

d = 6 ;   allocate(rho(1:d,1:d)) ;   da = 2 ;   db = 3
dp = 0.005d0 ;   p = 1.d-15 - dp
do
  p = p + dp
  adb = sqrt(2.d0 - 2.d0*sqrt(1.d0 - (40.d0*p**4.d0)/81.d0))
  rho = 0.d0
  rho(1,1) = (1.d0+p)/6.d0 ;    rho(1,5) = p/3.d0 ;    rho(1,6) = p/3.d0 ;    rho(2,2) = (1.d0-p)/6.d0
  rho(3,3) = (1.d0-p)/6.d0 ;   rho(4,4) = (1.d0-p)/6.d0 ;   rho(5,1) = p/3.d0 ;   rho(5,5) = (1.d0+p)/6.d0 ;   rho(5,6) = p/3.d0
  rho(6,1) = p/3.d0 ;   rho(6,5) = p/3.d0 ;   rho(6,6) = (1.d0+p)/6.d0
  ndb = discord_easy('b', da, db, rho)
  write(13,*) p, adb, ndb
  if( p > 1.d0 ) exit
enddo
deallocate( rho )

end
!-----------------------------------------------------------------------------------------------------------------------------------