!###################################################################################################################################
!                                               Tests - LAPACK eigensolvers
!###################################################################################################################################
program lapack_tests  ! Some simple tests on LAPACK eigensolvers
implicit none
integer, parameter :: nqb = 8  ! No. of qubits
integer :: ord_pm(1:nqb)  ! Vector for the order of the Pauli matrices
integer, parameter :: d = 2**nqb  ! Dimension of the matrices and vectors
complex(8) :: kp(1:d,1:d), A(1:d,1:d)  ! Matrix for the Kronecker product and the one to be diagonalized
real(8) :: W(1:d)  ! Vector of real eigenvalues
complex(8) :: Wc(1:d)  ! Vector of complex eigenvalues
real(8) :: rn(1:nqb)  ! For the random numbers used to choose ord_pm
integer :: j  ! Auxiliary variable for counters
real(8) :: t1, t2  ! Times
character(10), dimension(5) :: optg  ! Options for the generators

call cpu_time(t1)

optg = 'std'  ! Sets the methods to be used by the generators

write(*,*) "## Performing some tests for the LAPACK eigensolvers"

call rng(optg, d, rn)  ! Gets the random numbers
do j = 1, nqb  ! Choose the indexes for the order of the Pauli matrices
     if ( (rn(j) >= 0.d0) .and. (rn(j) < 0.25d0) ) then; ord_pm(j) = 0
else if ( (rn(j) >= 0.25d0) .and. (rn(j) < 0.5d0) ) then; ord_pm(j) = 1
else if ( (rn(j) >= 0.5d0) .and. (rn(j) < 0.75d0) ) then; ord_pm(j) = 2
else if ( (rn(j) >= 0.75d0) .and. (rn(j) <= 1.d0) ) then; ord_pm(j) = 3
     endif
enddo
write(*,*) "ord_pm", ord_pm

call kron_prod_pauli_mat(ord_pm, nqb, kp)  ! Gets the tensor product of the choosed Pauli matrices
A = kp ;   call lapack_zheevd('N', d, A, W)  ! Calls LAPACK zheevd, which computes the eigenvalues of Hermitian matrices
write(*,*) "zheevd", sum(W), ". The sum should be zero."

A = kp ;   call lapack_zgeev('N', d, A, Wc)  ! Calls LAPACK zgeev, which computes the eigenvalues of general matrices
W = real(Wc) ;   call qsort(W,d)  ! We define A again because it could be changed by zheevd
write(*,*) "zgeev", sum(W), ". The sum should be zero."

call cpu_time(t2) ;   write(*,*) 'time=',t2-t1,'seconds'  ! Writes the time taken by this program

end
!###################################################################################################################################