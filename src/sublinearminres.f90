!***********************************************************************************
! linear minization residuals (MINRES) to solve Ax=b
! solve the normal system H_d delta d^e = delta d^*
! gaoshan guo at CNRS geoazur
!***********************************************************************************

subroutine sublinearminres(istot,grid,inv,ir)
  use parameters
  use kwsinc
  use minresModule
  implicit none
  type (grid_type)                  :: grid
  type (inv_type)                   :: inv
  integer                           :: niter,istot,nn,i,unit
  real, allocatable                 :: deltade01d(:),deltade1d(:),deltadr1d(:)
  character(len=160)                :: namei
  integer                           :: ir(grid%nrec,5)
  real                              :: deltadr_norm, deltade_norm, fac

  ! Local variables
  logical             :: normal, precon, ifLS
  integer             :: n, nout
  real                :: shift, pertM

  logical   :: checkA 
  integer   :: itnlim, istop, itn, nprint
  real      :: Anorm, Acond, Arnorm, rnorm, rtol, r1norm, ynorm
  real      :: enorm, etol, wnorm, xnorm

  nn=grid%nt*grid%nrec
 
  allocate(deltade1d(nn),deltade01d(nn),deltadr1d(nn))
 
  ! transform the initial data residuals into 1D
  !call sub2d21d(deltade0,grid%nrec,grid%nt,nn,deltade01d)

  deltade01d=0

  ! transform the observed data residuals into 1D
  call sub2d21d(inv%deltadr,grid%nrec,grid%nt,nn,deltadr1d) 

  normal=.false.
  precon=.true.

  shift=0

  checkA=.false.
  itnlim = grid%maxiter
  rtol   = grid%error

  if(grid%oio.eq.1)then
    nout = 111
  else
    nout = -111
  endif

  call subname(istot,namei)

  if(nout.gt.0)then
    open(nout,file='MINRES_'//namei(1:LEN_TRIM(namei))//'.txt',status='unknown')
  endif

  ! adjoint simulation
  call submodeling_adjoint_esource_data(istot,grid%c0,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,&
              inv%deltadr,grid%type_rec,grid%nrec,ir,inv%deltade)

  call datanorm(inv%deltadr,grid%nt,grid%nrec,deltadr_norm)
  call datanorm(inv%deltade,grid%nt,grid%nrec,deltade_norm)

  fac=deltadr_norm/deltade_norm*grid%mu

  shift=-fac

  call MINRES(nn, deltadr1d, shift, checkA, normal, &
             deltade01d, deltade1d, itnlim, nout, rtol,istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm,&
             istot,grid,inv,ir)
 
  if(nout.gt.0)then
    close(nout)
  endif

  call sub1d22d(deltade1d,grid%nrec,grid%nt,nn,inv%deltade)

 deallocate(deltade1d,deltade01d,deltadr1d)

end subroutine

!subroutine MINRES( n, b, shift, checkA, precon, x0, x, itnlim, nout, rtol, &
!                     istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm, &
!                     istot,grid,inv,ir)
!    use parameters
!    integer    :: n, itnlim, nout
!    logical    :: checkA, precon
!    real       :: b(n), x0(n)
!    real       :: shift, rtol
!    real       :: x(n)
!    integer    :: istop, itn
!    real       :: Anorm, Acond, rnorm, Arnorm, ynorm

    ! wavefield modeing
!    type (grid_type)                  :: grid
!    type (inv_type)                   :: inv
!    integer                           :: !istot
!    integer                           :: ir(grid%nrec,5)

!    x0=1

!    b=1

!    x=1

!    inv%deltade=0

!end subroutine


subroutine datanorm(datao,nt,nr,data_norm)
 implicit none
 integer:: nt,nr,i,j
 real:: datao(nt,nr)
 real:: data_norm

 data_norm=0

 do i=1, nt
 do j=1, nr
    data_norm=data_norm+datao(i,j)**2
 enddo
 enddo

 data_norm=sqrt(data_norm)

end subroutine

!***************************************
! The transform form 2d to 1d
!***************************************
subroutine sub2d21d(deltade,nr,nt,nn,deltade1d)
  implicit none
  integer :: nr,nn,nt,it,ir,k
  real :: deltade1d(nn),deltade(nt,nr)

  k=0
  do ir=1, nr
     do it=1, nt
        k=k+1
        deltade1d(k)=deltade(it,ir)
     end do
  end do

end subroutine sub2d21d

!******************************************
! The transform form 1d to 2d
!******************************************
subroutine sub1d22d(deltade1d,nr,nt,nn,deltade)
  implicit none
  integer :: nr,nn,nt,it,ir,k
  real :: deltade1d(nn),deltade(nt,nr)

  k=0
  do ir=1, nr
     do it=1, nt
        k=k+1
        deltade(it,ir)=deltade1d(k)
     end do
  end do

end subroutine sub1d22d

