!-----------------------------------------------------------------------------------------------------------------------------------
program test
integer :: j, k, l
integer :: vec_bin(1:4)
integer :: nd = 4
integer :: dec

do j = 0, 1 ;   do k = 0, 1 ;   do l = 0, 1 ;  do m = 0, 1
  vec_bin(1) = j ;   vec_bin(2) = k ;   vec_bin(3) = l ;   vec_bin(4) = m
  call bin2dec(vec_bin, nd, dec) ;   write(*,*) j, k, l, m, dec
enddo ;   enddo ;   enddo ;   enddo

write(*,*) "---------------------------------------"

do j = 0, 2**nd-1
  dec = j ;   call dec2bin(dec, nd, vec_bin) ;   write(*,*) j, vec_bin(1), vec_bin(2), vec_bin(3), vec_bin(4)
enddo

end
!-----------------------------------------------------------------------------------------------------------------------------------