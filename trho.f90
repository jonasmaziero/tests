!###################################################################################################################################
!                                                   Basic tests for density matrices
!###################################################################################################################################
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
!###################################################################################################################################