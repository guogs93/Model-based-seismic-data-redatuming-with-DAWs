
module minresModule

   use parameters
   implicit none
   public   :: MINRES, blurring, precodition

contains

  subroutine MINRES( n, b, shift, checkA, precon, x0, x, itnlim, nout, rtol, &
                     istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm, &
                     istot,grid,inv,ir)

    integer :: n, itnlim, nout
    logical :: checkA, precon
    real :: b(n), x0(n)
    real :: shift, rtol
    real :: x(n)
    integer :: istop, itn
    real :: Anorm, Acond, rnorm, Arnorm, ynorm

    ! wavefield modeing
    type (grid_type)                  :: grid
    type (inv_type)                   :: inv
    integer                           :: istot
    integer                           :: ir(grid%nrec,5)
    !
    !     Local arrays and variables
      real  :: r1(n), r2(n), v(n), w(n), w1(n), w2(n), y(n)
      real  :: alfa  , beta  , beta1 , cs    ,          &
                   dbar  , delta , denom , diag  ,          &
                   eps   , epsa  , epsln , epsr  , epsx  ,  &
                   gamma , gbar  , gmax  , gmin  ,          &
                   oldb  , oldeps, qrnorm, phi   , phibar,  &
                   rhs1  , rhs2  , rnorml, rootl ,          &
                   Arnorml,        relArnorml,              &
                   s     , sn    , t     , tnorm2, ynorm2, z
      integer   :: i
      logical   :: debug, prnt

    ! Local constants
      real,         parameter :: zero =  0.0,  one = 1.0
      real,         parameter :: ten  = 10.0
      character(len=*), parameter :: enter = ' Enter MINRES.  '
      character(len=*), parameter :: exitt = ' Exit  MINRES.  '
      character(len=*), parameter :: msg(-1:8) =                  &
        (/ 'beta2 = 0.  If M = I, b and x are eigenvectors of A', & ! -1
           'beta1 = 0.  The exact solution is  x = 0           ', & !  0
           'Requested accuracy achieved, as determined by rtol ', & !  1
           'Reasonable accuracy achieved, given eps            ', & !  2
           'x has converged to an eigenvector                  ', & !  3
           'Acond has exceeded 0.1/eps                         ', & !  4
           'The iteration limit was reached                    ', & !  5
           'Aprod  does not define a symmetric matrix          ', & !  6
           'Msolve does not define a symmetric matrix          ', & !  7
           'Msolve does not define a pos-def preconditioner    ' /) !  8
    !-------------------------------------------------------------------

    intrinsic       :: abs, dot_product, epsilon, min, max, sqrt

    ! Print heading and initialize.
    
    debug = .false.
    eps   = epsilon(eps)

    if (nout > 0) then
       write(nout, 1000) enter, n, checkA, precon, itnlim, rtol, shift
    end if

    istop  = 0
    itn    = 0
    Anorm  = zero
    Acond  = zero
    rnorm  = zero
    ynorm  = zero
    x(1:n) = 0 ! x0(1:n)

    !-------------------------------------------------------------------
    ! Set up y and v for the first Lanczos vector v1.
    ! y = beta1 P' v1, where P = C**(-1).
    ! v is really P' v1.
    !-------------------------------------------------------------------
    y      = b
    r1     = b
    if ( precon ) call precodition(istot,b,y,n,grid)

    beta1  = dot_product(b,y)

    if (beta1 < zero) then     ! M must be indefinite.
       istop = 8
       go to 900
    end if

    if (beta1 == zero) then    ! b = 0 exactly.  Stop with x = 0.
       istop = 0
       go to 900
    end if

    beta1  = sqrt( beta1 )     ! Normalize y to get v1 later.

    !-------------------------------------------------------------------
    ! See if Msolve is symmetric.
    !-------------------------------------------------------------------
    if (checkA  .and.  precon) then
       call precodition(istot,y,r2,n,grid)
       s      = dot_product(y ,y )
       t      = dot_product(r1,r2)
       z      = abs(s - t)
       epsa   = (s + eps) * eps**0.33333
       if (z > epsa) then
          istop = 7
          go to 900
       end if
    end if

    !-------------------------------------------------------------------
    ! See if Aprod  is symmetric.  Initialize Arnorm.
    !-------------------------------------------------------------------
    if (checkA) then
       call blurring (istot, grid, inv, ir, n, y, w )
       call blurring (istot, grid, inv, ir, n, w, r2 )
       s      = dot_product(w,w )
       t      = dot_product(y,r2)
       z      = abs(s - t)
       epsa   = (s + eps) * eps**0.33333
       if (z > epsa) then
          istop = 6
          go to 900
       end if
       Arnorml = sqrt(s);
    else
       call blurring (istot,grid,inv,ir, n, y, w )
       Arnorml = sqrt( dot_product(w,w) )
    end if

    !-------------------------------------------------------------------
    ! Initialize other quantities.
    !-------------------------------------------------------------------
    oldb   = zero
    beta   = beta1
    dbar   = zero
    epsln  = zero
    qrnorm = beta1
    phibar = beta1
    rhs1   = beta1
    rhs2   = zero
    tnorm2 = zero
    ynorm2 = zero
    cs     = - one
    sn     = zero
    w(1:n) = zero
    w2(1:n)= zero
    r2(1:n)= r1

    if (debug) then
       write(*,*) ' '
       write(*,*) 'b    ', b
       write(*,*) 'beta ', beta
       write(*,*) ' '
    end if

    !===================================================================
    ! Main iteration loop.
    !===================================================================
    do
       itn = itn + 1               ! k = itn = 1 first time through

       !----------------------------------------------------------------
       ! Obtain quantities for the next Lanczos vector vk+1, k = 1, 2,...
       ! The general iteration is similar to the case k = 1 with v0 = 0:
       !
       !   p1      = Operator * v1  -  beta1 * v0,
       !   alpha1  = v1'p1,
       !   q2      = p2  -  alpha1 * v1,
       !   beta2^2 = q2'q2,
       !   v2      = (1/beta2) q2.
       !
       ! Again, y = betak P vk,  where  P = C**(-1).
       ! .... more description needed.
       !----------------------------------------------------------------
       s      = one / beta            ! Normalize previous vector (in y).
       v      = s*y(1:n)              ! v = vk if P = I

       call blurring (istot, grid, inv, ir, n, v, y )
       y      = y - shift*v           ! call daxpy ( n, (- shift), v, 1, y, 1 )
       if (itn >= 2) then
          y   = y - (beta/oldb)*r1    ! call daxpy ( n, (- beta/oldb), r1, 1, y,1 )
       end if

       alfa   = dot_product(v,y)      ! alphak
       y      = y - (alfa/beta)*r2    ! call daxpy ( n, (- alfa/beta), r2, 1, y, 1 )
       r1     = r2
       r2     = y
       if ( precon ) call precodition(istot,r2,y,n,grid)

       oldb   = beta                  ! oldb = betak
       beta   = dot_product(r2,y)     ! beta = betak+1^2
       if (beta < zero) then
          istop = 6
          go to 900
       end if

       beta   = sqrt( beta )          ! beta = betak+1
       tnorm2 = tnorm2 + alfa**2 + oldb**2 + beta**2

       if (itn == 1) then                   ! Initialize a few things.
          if (beta/beta1 <= ten*eps) then   ! beta2 = 0 or ~ 0.
             istop = -1                     ! Terminate later.
          end if
         !tnorm2 = alfa**2
          gmax   = abs( alfa )              ! alpha1
          gmin   = gmax                     ! alpha1
       end if

       ! Apply previous rotation Qk-1 to get
       !   [deltak epslnk+1] = [cs  sn][dbark    0   ]
       !   [gbar k dbar k+1]   [sn -cs][alfak betak+1].

       oldeps = epsln
       delta  = cs * dbar  +  sn * alfa ! delta1 = 0         deltak
       gbar   = sn * dbar  -  cs * alfa ! gbar 1 = alfa1     gbar k
       epsln  =               sn * beta ! epsln2 = 0         epslnk+1
       dbar   =            -  cs * beta ! dbar 2 = beta2     dbar k+1

       ! Compute the next plane rotation Qk

       gamma  = sqrt( gbar**2 + beta**2 )   ! gammak
       cs     = gbar / gamma                ! ck
       sn     = beta / gamma                ! sk
       phi    = cs * phibar                 ! phik
       phibar = sn * phibar                 ! phibark+1

       if (debug) then
          write(*,*) ' '
          write(*,*) 'v    ', v
          write(*,*) 'alfa ', alfa
          write(*,*) 'beta ', beta
          write(*,*) 'gamma', gamma
          write(*,*) 'delta', delta
          write(*,*) 'gbar ', gbar
          write(*,*) 'epsln', epsln
          write(*,*) 'dbar ', dbar
          write(*,*) 'phi  ', phi
          write(*,*) 'phiba', phibar
          write(*,*) ' '
       end if

       ! Update  x.

       denom = one/gamma

       do i = 1, n
          w1(i) = w2(i)
          w2(i) = w(i)
          w(i)  = ( v(i) - oldeps*w1(i) - delta*w2(i) ) * denom
          x(i)  =   x(i) +   phi * w(i)
       end do

       ! Go round again.

       gmax   = max( gmax, gamma )
       gmin   = min( gmin, gamma )
       z      = rhs1 / gamma
       ynorm2 = z**2  +  ynorm2
       rhs1   = rhs2  -  delta * z
       rhs2   =       -  epsln * z

       ! Estimate various norms and test for convergence.

       Anorm  = sqrt( tnorm2 )
       ynorm  = sqrt( ynorm2 )
       epsa   = Anorm * eps
       epsx   = Anorm * ynorm * eps
       epsr   = Anorm * ynorm * rtol
       diag   = gbar
       if (diag == zero) diag = epsa

       qrnorm = phibar
       rnorml = rnorm
       rnorm  = qrnorm
       rootl       = sqrt( gbar**2 +dbar**2  )  ! norm([gbar; dbar]);
       Arnorml     = rnorml*rootl               ! ||A r_{k-1} ||
       relArnorml  = rootl  /  Anorm;           ! ||Ar|| / (||A|| ||r||)     
       !relArnorml = Arnorml / Anorm;           ! ||Ar|| / ||A|| 

       ! Estimate  cond(A).
       ! In this version we look at the diagonals of  R  in the
       ! factorization of the lower Hessenberg matrix,  Q * H = R,
       ! where H is the tridiagonal matrix from Lanczos with one
       ! extra row, beta(k+1) e_k^T.

       Acond  = gmax / gmin

       ! See if any of the stopping criteria are satisfied.
       ! In rare cases, istop is already -1 from above (Abar = const*I).

       if (istop == 0) then
          if (itn    >= itnlim    ) istop = 5
          !if (Acond  >= 0.1d+0/eps) istop = 4
          !if (epsx   >= beta1     ) istop = 3
          !if (qrnorm <= epsx  .or.  relArnorml <= epsx) istop = 2
          if (qrnorm/beta1 <= rtol) istop = 1
       end if


       ! See if it is time to print something.

       if (nout > 0) then
          prnt   = .false.
          if (n      <= 40         ) prnt = .true.
          if (itn    <= 10         ) prnt = .true.
          if (itn    >= itnlim - 10) prnt = .true.
          if (mod(itn,10)  ==     0) prnt = .true.
          if (qrnorm <=  ten * epsx) prnt = .true.
          if (qrnorm <=  ten * epsr) prnt = .true.
          if (relArnorml<= ten*epsx) prnt = .true.
          if (relArnorml<= ten*epsr) prnt = .true.
          if (Acond  >= 1.0d-2/eps ) prnt = .true.
          if (istop  /=  0         ) prnt = .true.

          if ( prnt ) then
             if (    itn     == 1) write(nout, 1200)
             write(nout, 1300) itn, x(1), qrnorm/beta1, Anorm, Acond
             if (mod(itn,10) == 0) write(nout, 1500)
          end if
       end if
       if (istop /= 0) exit

    end do
    !===================================================================
    ! End of iteration loop.
    !===================================================================

    ! Display final status.

900 Arnorm = Arnorml
    if (nout  > 0) then
       write(nout, 2000) exitt, istop, itn,   &
                         exitt, Anorm, Acond, &
                         exitt, rnorm, ynorm, Arnorm
       write(nout, 3000) exitt, msg(istop)
    end if

    return

 1000 format(// 1p,    a, 5x, 'Solution of symmetric   Ax = b'    &
              / ' n      =', i7, 5x, 'checkA =', l4, 12x,         &
                 'precon =', l4                                   &
              / ' itnlim =', i7, 5x, 'rtol   =', e11.2, 5x,       &
                 'shift  =', e23.14)
 1200 format(// 5x, 'itn', 8x, 'x(1)', 10x,                       &
                'norm(r)', 3x, 'norm(A)', 3X, 'cond(A)')
 1300 format(1p, i8, e19.10, 3e10.2)
 1500 format(1x)
 2000 format(/ 1p, a, 5x, 'istop =', i3,   14x, 'itn   =', i8     &
             /     a, 5x, 'Anorm =', e12.4, 5x, 'Acond =', e12.4  &
             /     a, 5x, 'rnorm =', e12.4, 5x, 'ynorm =', e12.4, 5x, 'Arnorml=', e12.4)
 3000 format(      a, 5x, a )

end subroutine MINRES 

!*******************************
subroutine blurring(istot,grid,inv,ir,nn,deltade1d,Hdeltade1d)
  use parameters
  use kwsinc
  implicit none
  type (grid_type)                  :: grid
  type (inv_type)                   :: inv
  integer                           :: nn,i,j
  real                              :: deltade1d(nn),Hdeltade1d(nn)
  integer                           :: istot,irec0
  real                              :: alpha,cost,deltade2_norm,deltade_norm,fac
  real, allocatable                 :: dataesc(:,:),deltade(:,:),deltade2(:,:),Hdeltade(:,:),deltadew(:,:)
  integer*8                         :: bwd,fwd
  integer                           :: ir(grid%nrec,5)
  integer                           :: unit

  allocate(deltade(grid%nt,grid%nrec),deltadew(grid%nt,grid%nrec))
  allocate(Hdeltade(grid%nt,grid%nrec),deltade2(grid%nt,grid%nrec))

  call sub1d22d(deltade1d,grid%nrec,grid%nt,nn,deltade)

  !*******************************************************************
  ! calculate gradient g=F*{ (S(m)S(m)^T+mu I) deltade - mu deltadr }
  !*******************************************************************

  ! adjoint simulation
  call submodeling_adjoint_esource_data(istot,grid%c0,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,&
              deltade,grid%type_rec,grid%nrec,ir,deltade2)

  call datanorm(deltade2,grid%nt,grid%nrec,deltade2_norm)
  call datanorm(deltade,grid%nt,grid%nrec,deltade_norm)

  do j=1, grid%nt
     Hdeltade(j,:)=deltade2(grid%nt-j+1,:) !+ fac*deltade(j,:)
  enddo

  ! transform the data residuals into 1D
  call sub2d21d(Hdeltade,grid%nrec,grid%nt,nn,Hdeltade1d)

  deallocate(deltade,Hdeltade,deltade2,deltadew)

end subroutine

subroutine precodition(is,res1d,res_preco1d,nn,grid)
  use parameters
  implicit none
  type (grid_type)            :: grid
  integer                     :: nn,it,ir,unit,is
  real                        :: res1d(nn), res_preco1d(nn)
  real, allocatable           :: res2d(:,:), res_preco2d(:,:), res_precom2d(:,:)

  allocate(res2d(grid%nt,grid%nrec), res_preco2d(grid%nt,grid%nrec), res_precom2d(grid%nt,grid%nrec))

  call sub1d22d(res1d,grid%nrec,grid%nt,nn,res2d)

  !call mute(is,res2d,grid%nrec,grid%nt,grid%dt,res_preco2d)

  call filmartin2d(res2d,grid%nrec,grid%nt,grid%dt,grid%flow,grid%fhig,grid%llow,grid%lhig,res_precom2d)

  !open(NEWUNIT=unit,file='res.bin',access='stream',form='unformatted')
  !do ir=1, grid%nrec
  !do it=1, grid%nt
  !   write(unit) res2d(grid%nt-it+1,ir)
  !enddo
  !enddo
  !close(unit)

  !open(NEWUNIT=unit,file='res_preco.bin',access='stream',form='unformatted')
  !do ir=1, grid%nrec
  !do it=1, grid%nt
  !   write(unit) res_preco2d(grid%nt-it+1,ir)
  !enddo
  !enddo
  !close(unit)  

  !open(NEWUNIT=unit,file='res_precom.bin',access='stream',form='unformatted')
  !do ir=1, grid%nrec
  !do it=1, grid%nt
  !   write(unit) res_precom2d(grid%nt-it+1,ir)
  !enddo
  !enddo
  !close(unit)

  call sub2d21d(res_precom2d,grid%nrec,grid%nt,nn,res_preco1d)

  deallocate(res2d,res_preco2d,res_precom2d)

end subroutine

subroutine mute(is,res_preco2d,nr,nt,dt,res_precom2d)
 implicit none
 integer:: nr, nt, is, ir, it
 real:: dt
 real:: res_preco2d(nt,nr), res_precom2d(nt,nr) 
 real, allocatable:: t0(:)
 integer, allocatable:: it0(:)
 character(len=160):: namei

 call subname(is,namei)

 allocate(t0(nr),it0(nr))

 open(unit=11,file='ftraveltime_'//namei(1:LEN_TRIM(namei))//'.txt')
 do ir=1, nr
    read(11,*) t0(ir)
 enddo
 close(11)

 it0(:)=nint(t0(:)/dt)+1

 do ir=1, nr
    do it=1, nt
       if(it.ge.nt-it0(ir))then
         res_precom2d(it,ir)=0
       else
         res_precom2d(it,ir)=res_preco2d(it,ir)
       endif
    enddo
 enddo

 deallocate(t0,it0)

end subroutine

subroutine filmartin2d(deltadeo,nr,nt,dt,flow,fhig,llow,lhig,deltade)
  implicit none
  integer:: nr,nt,nn,i,j,k,nto,i2,nf,nord
  real:: deltadeo(nt,nr),deltade(nt,nr)
  real:: flow, fhig,llow,lhig
  real:: dt, fe, fmax
  integer :: conv
  integer :: iflow, ifhig, illow, ilhig
  real, allocatable               :: trace(:)
  real, allocatable               :: han1(:),han2(:)
  real, allocatable               :: fmartin(:)
  real*8, allocatable             :: dmartin(:)
  real*8, allocatable             :: preal(:),pimag(:)
  real*8, allocatable             :: ybuf(:)

  integer, parameter              :: idnum=14
  integer inum(idnum)
  data inum /8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536/

do i=1,idnum
     if (nt.le.inum(i)) then
         nto=inum(i)
         go to 100
     end if
  end do

100  write(*,*) "nto = ",nto

  ! exponent of the Hanning taper
  nord=2

  !max frequency

  nf=nto
  fe=1./dt
  fmax=fe/2

  if (fhig.gt.fmax) then
     write(*,*) 'WARNING: fhig > fmax; fhig forced to fmax'
     fhig=fmax
     lhig=0.
  end if

  conv=float(nto)/fe
  iflow=int(flow*conv)+1
  ifhig=int(fhig*conv)+1
  illow=int(llow*conv)+1
  ilhig=int(lhig*conv)+1

  if (iflow-illow.lt.1.and.iflow.gt.0) then
     illow=iflow-1
  end if

  if (ifhig+ilhig.gt.nf) then
     ilhig=nf-ifhig
  end if

  write(*,*) 'iflow ifhig illow ilhig = ',iflow,ifhig,illow,ilhig

  allocate(trace(nto))
  allocate(han1(illow+1))
  allocate(han2(ilhig+1))
  allocate(fmartin(nf+1))
  allocate(dmartin(nf+1))
  allocate(preal(nto))
  allocate(pimag(nto))
  allocate(ybuf(nto))

  trace(:)=0.
  han1(:)=0.
  han2(:)=0.

! Building filter

  fmartin(:)=0.

! Low-cut branch

  if (illow.gt.0) then
     call sphanning(han1,illow,nord)
     do j=1,illow
        fmartin(iflow-j+1)=han1(j)
     end do
  end if

! High-cut branch

  han2(:)=0.

  if (ilhig.gt.0) then
     call sphanning(han2,ilhig,nord)
     do j=1,ilhig
        fmartin(ifhig+j-1)=han2(j)
     end do
  end if

  do j=iflow,ifhig
     fmartin(j)=1.
  end do

  dmartin(:)=dble(fmartin(:)+0.0001)

  do i2=1, nr
     pimag(:)=0.
     do j=1, nt
        trace(j)=deltadeo(j,i2)
     enddo
     preal(:)=dble(trace(:))
     call realft2s(nto,2,preal,pimag,ybuf)
     do j=1, nto
        preal(j)=preal(j)*dmartin(j)
        pimag(j)=pimag(j)*dmartin(j)
     end do
     call realft2s(nto,1,preal,pimag,ybuf)
     do j=1, nt
        deltade(j,i2)=sngl(preal(j))
     enddo
  end do

  deallocate(trace)
  deallocate(han1)
  deallocate(han2)
  deallocate(fmartin)
  deallocate(dmartin)
  deallocate(preal)
  deallocate(pimag)
  deallocate(ybuf)

end subroutine

!**************************************
! hanning window
!**************************************
subroutine sphanning(han,nf,nord)
  implicit none
  integer:: j, nf, nord
  real:: han(nf+1), han1, han2, han3
  real,parameter:: pi=3.14159265

  do j=1, nf+1
     han1=2.*pi*float(j-1)/float(2*nf)
     han2=(cos(han1)+1.)/2.+0.000001
     han3=log(han2)*float(nord)
     han(j)=exp(han3)
  end do

end subroutine sphanning

end module minresModule
