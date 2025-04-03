module minresModule

   use parameters
   implicit none
   public   :: MINRES, blurring, precodition

contains

  subroutine MINRES( n, b, shift, checkA, precon, x0, x, itnlim, nout, rtol, &
                     istop, itn, Anorm, Acond, rnorm, Arnorm, ynorm, &
                     istot,grid) !,inv,ir)

    integer,  intent(in)    :: n, itnlim, nout
    logical,  intent(in)    :: checkA, precon
    real, intent(in)    :: b(n), x0(n)
    real, intent(in)    :: shift, rtol
    real, intent(out)   :: x(n)
    integer,  intent(out)   :: istop, itn
    real, intent(out)   :: Anorm, Acond, rnorm, Arnorm, ynorm

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
    
   write(*,*)'error2'

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
    if ( precon ) call precodition(b,y,n,grid)

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
       call precodition(y,r2,n,grid)
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
       if ( precon ) call precodition(r2,y,n,grid)

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
          if (qrnorm/beta1 <= 0.00001) istop = 1
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
  include 'fftw3.f'
  type (grid_type)                  :: grid
  type (inv_type)                   :: inv
  integer                           :: nn,i,j
  real                              :: deltade1d(nn),Hdeltade1d(nn)
  integer                           :: istot,irec0
  real                              :: alpha,cost,deltade2_norm,deltade_norm,fac
  real, allocatable                 :: dataesc(:,:),deltade(:,:),deltade2(:,:),Hdeltade(:,:),deltadew(:,:)
  integer*8                         :: bwd,fwd
  integer                           :: ir(grid%nrec,5)

  allocate(deltade(grid%nt,grid%nrec),deltadew(grid%nt,grid%nrec))
  allocate(Hdeltade(grid%nt,grid%nrec),deltade2(grid%nt,grid%nrec))

  call sub1d22d(deltade1d,grid%nrec,grid%nt,nn,deltade)

  !*******************************************************************
  ! calculate gradient g=F*{ (S(m)S(m)^T+mu I) deltade - mu deltadr }
  !*******************************************************************

  write(*,*)'Apply Hessian'

  ! adjoint simulation
  call submodeling_adjoint_esource_data(istot,grid%c0,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,&
              deltade,grid%type_rec,grid%nrec,ir,deltade2)

  call datanorm(deltade2,grid%nt,grid%nrec,deltade2_norm)
  call datanorm(deltade,grid%nt,grid%nrec,deltade_norm)

  fac=deltade2_norm/deltade_norm*grid%mu

  do j=1, grid%nt
     Hdeltade(j,:)=deltade2(grid%nt-j+1,:) + fac*deltadew(j,:)
  enddo

  ! transform the data residuals into 1D
  call sub2d21d(Hdeltade,grid%nrec,grid%nt,nn,Hdeltade1d)

  deallocate(deltade,Hdeltade,deltade2,deltadew)

end subroutine


subroutine precodition(res1d,res_preco1d,nn,grid)
  use parameters
  implicit none
  type (grid_type)            :: grid
  integer                     :: nn
  real                        :: res1d(nn), res_preco1d(nn)

  res_preco1d=res1d

end subroutine

end module minresModule
