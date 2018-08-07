!-----------------------------------------------------------------------------------------------------------------------------------
program tests
 implicit none
 !call test_mat_sqrt()
 call test_werner()
end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine test_mat_sqrt()
  implicit none
  complex(8) :: sigma_0(2,2), sigma_1(2,2), sigma_2(2,2), sigma_3(2,2)
  complex(8) :: Asr(2,2)
  call pauli_group(sigma_0, sigma_1, sigma_2, sigma_3)
  call mat_sqrt(2, sigma_1, Asr)
  call array_display(2, 2, Asr)
end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine test_werner()
  implicit none
  complex(8) :: rho(4,4)
  real(8) :: w, Dhe
  w = 1.d0
  call rho_bds(-w, -w, -w, rho)!; call array_display(4,4,rho)
  call discord_he_test(2,2,rho,Dhe);  write(*,*) Dhe
end
!-----------------------------------------------------------------------------------------------------------------------------------
subroutine discord_he_test(da,db,rho,Dhe)
  ! Hellinger discord of qubit-qudit systems
  ! Ref: J. Phys. A: Math. Theor. 49 (2016) 235301
  integer :: da, db, d
  complex(8) :: rho(da*db,da*db)
  real(8) :: Dhe
  !complex(8) :: rhosr(da*db,da*db), rhosrA(da,da), rhosrB(db,db)
  !real(8) :: W(da**2-1), bvB(db**2-1), bvA(da**2-1), corrmat(da**2-1,db**2-1), op(da**2-1,da**2-1)
  d = da*db
  call mat_sqrt(d,rho,rhosr); call array_display(4,4,rhosr)
  !call partial_trace_a_he(rhosr, da, db, rhosrB)
  !call bloch_vector_gellmann(db, rhosrB, bvB)
  !call partial_trace_b_he(rhosr, da, db, rhosrA)
  !call bloch_vector_gellmann(da, rhosrA, bvA)
  !call corrmat_gellmann(da, db, rhosr, corrmat)
  !call outer_product_re(da**2-1, bvA, bvA, op)
  !call lapack_zheevd('N', da**2-1, op+matmul(corrmat,transpose(corrmat)), W)
  !discord_he = 2.d0-2.d0*sqrt((trace_he(d,rhosr))**2.d0 + (norm_r(da**2-1,bvB))**2.d0 + maxval(W))
end
