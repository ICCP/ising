   subroutine calcsize(r,N)
   implicit none
   integer, intent(in) :: N
   real*8, intent(in) :: r(N,3)
   real*8 xcm,ycm,zcm,rcm,rg,ri
   integer j
   
   ! calculate the radius of gyration
   ! R_g^2 = (1/N) \sum <(r_i - R_cm)^2>
   
   ! calculate the center of the distribution, R_cm
   xcm = 0
   ycm = 0
   zcm = 0
   do j=1,N
      xcm = xcm + r(j,1)
      ycm = ycm + r(j,2)
      zcm = zcm + r(j,3)
   enddo
   
   rcm = sqrt(xcm**2 + ycm**2 + zcm**2)/N
   
   ! R_g calculation
   rg = 0
   do j=1,N
      ri = sqrt(r(j,1)**2 + r(j,2)**2 + r(j,3)**2)
      rg = rg + (ri - rcm)**2
   enddo
   rg = sqrt(rg/N)
   open(unit=20,file="RvsNT2010000.txt",access="append")
   write(20,*) log(rg), log(real(N))
   close(20)
      
   end subroutine calcsize
