subroutine subwavefieldresample(n1,n2,npml,nt,dt,ir,nr,sis,wavefield)
  use kwsinc
  implicit none
  integer:: n1,n2,npml,nt,i1,i2,it,nr
  integer:: ir(nr,5)
  real:: wavefield(n1,n2,nt)
  real:: sis(nt,nr)
  integer:: irec
  real:: dt

  !hicks parameter
  real                              :: coeff
  integer                           :: irad1min,irad1max,irad2min,irad2max
  integer                           :: ialpha1,ialpha2,nalpha_frees,ix1,ix2

  sis=0

  do it=1, nt
  do irec=1, nr
     ix1=ir(irec,1)-npml
     ix2=ir(irec,2)-npml

     ialpha1=ir(irec,3)
     ialpha2=ir(irec,4)
     nalpha_frees=ir(irec,5)

     irad1min=MAX(1-ix1,   -irad(1))
     irad1max=MIN(n1-ix1, irad(1))
     irad2min=MAX(1-ix2,   -irad(2))
     irad2max=MIN(n2-ix2, irad(2))

     do i1=irad1min,irad1max
     do i2=irad2min,irad2max
        coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
        sis(it,irec)=sis(it,irec) + coeff*sngl(wavefield(ix1 + i1,ix2 + i2, it))
     end do
     end do
  end do
  enddo

end subroutine
