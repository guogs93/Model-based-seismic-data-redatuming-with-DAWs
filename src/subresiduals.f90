! ==============================================================
! SUBROUTINE SUBMISFIT
! Compute the data residuals (res) and the l2 misfit function (cost)
! ==============================================================

subroutine subresiduals(is,dataobs,datacal,res,nt,nr)
  implicit none
  integer                           :: nt,nr,ir,it,is,ii
  real                              :: dataobs(nt,nr),datacal(nt,nr),res(nt,nr)
  real                              :: cost

  do ir=1,nr
     do it=1,nt
        res(nt-it+1,ir)=dataobs(it,ir)-datacal(it,ir)
     end do
  end do

end subroutine subresiduals
