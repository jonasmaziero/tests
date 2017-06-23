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
! (UNOPTIMIZED) Implementation of the partial trace directly from its definition
! Ref. J. Maziero, Computing partial traces and reduced density matrices, Int. J. Mod. Phys. C 28, 1750005 (2017)
!###################################################################################################################################
subroutine ptr(rho, d, di, nss, ssys, dr, rhor)  ! Returns the partial trace for general multipartite systems
implicit none
integer :: nss  ! Number of subsystems
integer :: di(1:nss)  ! Vector specifying the dimensions of the sub-systems
integer :: ssys(1:nss)  ! Vector (with components equal to 0 or 1) specifying the subsystems to be traced out. 
                        ! If ssys(j) = 0 the j-th subsystem is traced out. If ssys(j) = 1 it is NOT traced out.
integer :: d ! Total dimension (is the product of the sub-systems dimensions)
complex(8) :: rho(1:d,1:d)  ! Total matrix (given as input)
integer :: dr  ! Dimension of the reduced matrix (is the product of the dimensions of the sub-systems we shall not trace out)
complex(8) :: rhor(1:dr,1:dr)  ! Reduced matrix
complex(8), allocatable :: mat1(:,:), mat2(:,:), rhored(:,:)  ! Auxiliary matrices
integer :: j, k, l  ! Auxiliary variables for counters
integer :: da, db, dc  ! Auxiliary variables for the dimensions

! For bipartite systems
if ( nss == 2 ) then
  if ( ssys(1) == 0 ) then ; call ptr_a(rho, di(1), di(2), rhor)
  else if ( ssys(2) == 0 ) then ; call ptr_b(rho, di(1), di(2), rhor) 
  endif
  rhored = rhor
! For multipartite systems
else if ( nss >= 3 ) then
  ! Left partial traces
    l = 0 ;   do ;   if ( ssys(l+1) == 1 ) exit ;   l = l + 1 ;   enddo  ! l defines up to which position we shall trace out
    if ( l == 0 ) then
      mat1 = rho  ! This matrix shall be used below if l = 0
    else if ( l > 0 ) then
      if ( l == 1 ) then ;   da = di(1) ;   else ;   da = product(di(1:l)) ;  endif ;   db = d/da
      allocate( mat1(1:db,1:db) ) ;   call ptr_a(rho, da, db, mat1)
      d = db  ! After this operation, the matrix left over is mat1, whose dimension is d = db  
      rhored = mat1    
    endif
  ! Right partial traces
    k = nss+1 ;   do ;   if ( ssys(k-1) == 1 ) exit ;   k = k - 1 ;   enddo  ! k defines up to which position we shall trace out
    if ( k == (nss+1) ) then
      mat2 = mat1  ! This matrix shall be used below if k = nss+1
      rhored = mat2
    else if ( k < (nss+1) ) then
      if ( k == nss ) then ;   db = di(nss) ;   else ;   db = product(di(k:nss)) ;  endif ;   da = d/db
      allocate( mat2(1:da,1:da) ) ;   call ptr_b(mat1, da, db, mat2) ;   deallocate( mat1 )
      d = da  ! After this operation, the matrix left over is mat2, whose dimension is d = da
      rhored = mat2
    endif
    if ( l == 0 .and. k == (nss+1) ) deallocate( mat1 )  ! To avoid allocating mat1 two times in these cases
  ! Inner partial traces
    if ( (k-l) > 3 ) then  ! If (k-l) <= 3 there is no need to take inner partial traces
    do j = (l+2), (k-2)
      if ( ssys(j) == 0 ) then
        db = di(j) ;   if ( j == (k-2) ) then ;   dc = di(k-1) ;   else ;   dc = product(di(j+1:k-1)) ;  endif ;   da = d/(db*dc)
        allocate( mat1(1:da*dc,1:da*dc) ) ;   call ptr_3(mat2, da, db, dc, mat1) ;   rhored = mat1
        d = da*dc ;   deallocate( mat2 ) ;   allocate( mat2(1:d,1:d) ) ;   mat2 = mat1 ;   deallocate( mat1 )
      endif
    enddo
    endif
endif

rhor = rhored

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ptr_a(rho, da, db, rho_b)  ! Returns the left partial trace for a bi-partite matrix
implicit none
integer, intent(in) :: da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
complex(8), intent(out) :: rho_b(1:db,1:db)  !  Reduced matrix
integer :: j  ! Auxiliary variable for counters
complex(8), allocatable :: ida(:,:), idb(:,:), kp1(:,:), kp2(:,:)  ! For the indentities and Kronecker products

allocate( ida(1:da,1:da), idb(1:db,1:db), kp1(1:db,1:da*db), kp2(1:da*db,1:db) )
ida = 0.d0 ;   forall (j=1:da) ida(j,j) = 1.d0
if(db == da)then ; idb = ida ; else ; idb = 0.d0 ; forall (j=1:db) idb(j,j) = 1.d0 ;  endif
rho_b = 0.d0
do j = 1, da  ! we use the indentity to get the computational basis 
  call kronecker_product_c(ida(j,:), 1, da, idb, db, db, kp1) ;   call kronecker_product_c(ida(:,j), da, 1, idb, db, db, kp2)
  rho_b = rho_b + matmul(kp1, matmul(rho,kp2))
enddo
deallocate( ida, idb, kp1, kp2 )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ptr_b(rho, da, db, rho_a)  ! Returns the right partial trace for a bi-partite matrix
implicit none
integer, intent(in) :: da, db ! Dimensions of the subsystems (the dimension of the whole system is d = da*db)
complex(8), intent(in) :: rho(1:da*db,1:da*db)  ! Bipartite matrix (computational basis representation of the ragarded operator)
complex(8), intent(out) :: rho_a(1:da,1:da)  !  Reduced matrix
integer :: j  ! Auxiliary variable for counters
complex(8), allocatable :: ida(:,:), idb(:,:), kp1(:,:), kp2(:,:)  ! For the indentities and Kronecker products

allocate( ida(1:da,1:da), idb(1:db,1:db), kp1(1:da,1:da*db), kp2(1:da*db,1:da) )
ida = 0.d0 ;   forall (j=1:da) ida(j,j) = 1.d0
if(db == da)then ; idb = ida ; else ; idb = 0.d0 ; forall (j=1:db) idb(j,j) = 1.d0 ;  endif
rho_a = 0.d0
do j = 1, db  ! we use the indentity to get the computational basis 
  call kronecker_product_c(ida, da, da, idb(j,:), 1, db, kp1) ;   call kronecker_product_c(ida, da, da, idb(:,j), db, 1, kp2)
  rho_a = rho_a + matmul(kp1, matmul(rho,kp2))
enddo
deallocate( ida, idb, kp1, kp2 )

end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine ptr_3(rho, da, db, dc, rho_ac)  ! Returns the inner partial trace for a three-partite matrix
implicit none
integer, intent(in) :: da, db, dc  ! Dimension of the three sub-systems (the dimension of the whole system is d = da*db*dc)
complex(8), intent(in) :: rho(1:da*db*dc,1:da*db*dc)  ! Three-partite matrix (computational basis representation of the ragarded operator)
complex(8), intent(out) :: rho_ac(1:da*dc,1:da*dc)  ! Bipartite reduced matrix
integer :: j  ! Auxiliary variable for counters
complex(8), allocatable :: ida(:,:), idb(:,:), idc(:,:), kp1(:,:), kp2(:,:), kp3(:,:), kp4(:,:)  ! For the indentities and Kronecker products

allocate(ida(1:da,1:da),idb(1:db,1:db),idc(1:dc,1:dc))
allocate(kp1(1:da,1:da*db),kp2(1:da*dc,1:da*db*dc),kp3(1:da*db,1:da),kp4(1:da*db*dc,1:da*dc))
ida = 0.d0 ;   forall (j=1:da) ida(j,j) = 1.d0 
if(db == da)then ; idb = ida ; else ; idb = 0.d0 ; forall (j=1:db) idb(j,j) = 1.d0 ; endif
if(dc == da)then ; idc = ida ; else if(dc == db)then ; idc = idb ; else ; idc = 0.d0 ; forall (j=1:dc) idc(j,j) = 1.d0 ; endif
rho_ac = 0.d0
do j = 1, db  ! we use the indentity to get the computational basis 
  call kronecker_product_c(ida, da, da, idb(j,:), 1, db, kp1) ;   call kronecker_product_c(kp1, da, da*db, idc, dc, dc, kp2)
  call kronecker_product_c(ida, da, da, idb(:,j), db, 1, kp3) ;   call kronecker_product_c(kp3, da*db, da, idc, dc, dc, kp4)
  rho_ac = rho_ac + matmul(kp2, matmul(rho,kp4))
enddo
deallocate( ida, idb, idc, kp1, kp2, kp3, kp4 )

end
!###################################################################################################################################