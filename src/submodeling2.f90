! ================================================================================
! ================================================================================
! ================================================================================
! SUBROUTINE SUBMODELING: explicit time-stepping scheme with finite differences.
! Velocity-stress formulation. O(Delta t**2,Delta x**4).
! Forward simulation.
! ================================================================================
! ================================================================================
! ================================================================================

subroutine submodeling(is0,vp,n1,n2,npml,h,ofs,nt,dt,s,type_src,type_rec,is,ir,nr,sis,wavefield)
  use kwsinc
  implicit none

  integer,intent(in)                :: nt,n1,n2,npml
  real,intent(in)                   :: h,dt
  real                              :: s(nt)
  integer                           :: nr,is0,is(5),ir(nr,5),is1,is2

  ! wavefield arrays
  real, allocatable                 :: p(:,:),px(:,:),pz(:,:),vx(:,:),vz(:,:)
  real                              :: der1,der2
  real                              :: wavefield(n1,n2,nt)
  real                              :: sis(nt,nr)

  ! subsurface properties
  real                              :: vp(n1,n2)
  real, allocatable                 :: rho(:,:),vpe(:,:),rhoe(:,:),bue(:,:),bu1(:,:),bu2(:,:),buw1(:,:),buw2(:,:)
  real, allocatable                 :: kappa(:,:),kappaw1(:,:),kappaw2(:,:)

  !PML arrays
  real, allocatable                 :: spg(:),spg1(:)
  real, allocatable                 :: damp1a(:),damp1a1(:),damp2a(:),damp2a1(:),damp1b(:),damp1b1(:),damp2b(:),damp2b1(:)

  real, parameter                   :: xxfac=0.12

  integer                           :: n1e,n2e,i1,i2,i1min,i1max,i2min,i2max,it,imax,idir,ishift,irec
  integer                           :: type_src,type_rec
  integer                           :: ofs
  real                              :: vmax,xmax,tt
  real                              :: dtdx
  
  !hicks parameter
  real                              :: coeff
  integer                           :: irad1min,irad1max,irad2min,irad2max
  integer                           :: ialpha1,ialpha2,nalpha_frees,ix1,ix2

  real, parameter                   :: a0=1.125
  real, parameter                   :: a1=-1./24.

  logical, parameter                :: kws=.true.


  ! ############ Grid parameters #############################################
  n2e=n2+2*npml
  n1e=n1+2*npml
  dtdx=dt/h

  ! ##########################################################
  ! PML sponges (split PML implementation)
  allocate (damp2a(n2e))
  allocate (damp2a1(n2e))
  allocate (damp2b(n2e))
  allocate (damp2b1(n2e))
  allocate (damp1a(n1e))
  allocate (damp1a1(n1e))
  allocate (damp1b(n1e))
  allocate (damp1b1(n1e))

  allocate (spg(npml+1))
  allocate (spg1(npml+1))
  call cspongea(spg,spg1,npml,xxfac,damp2a,damp2a1,n2)
  call cspongeb(spg,spg1,npml,xxfac,damp2b,damp2b1,n2)
  call cspongea(spg,spg1,npml,xxfac,damp1a,damp1a1,n1)
  call cspongeb(spg,spg1,npml,xxfac,damp1b,damp1b1,n1)
  deallocate (spg)
  deallocate (spg1)

  ! ############### Staggered grids ###########################################
  allocate (kappa(n1e,n2e))
  allocate (kappaw2(n1e,n2e))
  allocate (kappaw1(n1e,n2e))
  allocate (buw2(n1e,n2e))
  allocate (buw1(n1e,n2e))

  allocate (bu1(n1e,n2e))
  allocate (bu2(n1e,n2e))
  allocate (rho(n1,n2))
  allocate (vpe(n1e,n2e))
  allocate (rhoe(n1e,n2e))
  allocate (bue(n1e,n2e))

  rho(:,:)=1.
 
  call submodext(vp,n1,n2,npml,vpe)
  call submodext(rho,n1,n2,npml,rhoe)
  call substagger(vpe,rhoe,bue,kappa,bu2,bu1,n1e,n2e)

  deallocate (vpe)
  deallocate (rho)
  deallocate (rhoe)

  do i2=1,n2e
     do i1=1,n1e
        kappaw2(i1,i2)=dtdx*kappa(i1,i2)*damp2a1(i2)
        kappaw1(i1,i2)=dtdx*kappa(i1,i2)*damp1a1(i1)
        buw2(i1,i2)=dtdx*bu2(i1,i2)*damp2b1(i2)
        buw1(i1,i2)=dtdx*bu1(i1,i2)*damp1b1(i1)
     end do
  end do

  deallocate (bu2)
  deallocate (bu1)

  ! #################################################################
  ! Memory allocation for wavefiels (p=px+pz due to split pml)
  allocate (p(n1e,n2e))
  allocate (px(n1e,n2e))
  allocate (pz(n1e,n2e))
  allocate (vx(n1e,n2e))
  allocate (vz(n1e,n2e))

  ! ==========================================================
  ! LOOP OVER SHOTS
  ! ==========================================================

  wavefield(:,:,:)=0.
  sis(:,:)=0.

  is1=is(1)
  is2=is(2)

  p(:,:)=0.	 ! Initial conditions
  px(:,:)=0.
  pz(:,:)=0.
  vx(:,:)=0.
  vz(:,:)=0.
  
  ! ==========================================================
  ! LOOP OVER TIME
  ! ==========================================================

  tt=0.
  do it=1,nt
     tt=tt+dt

     ! ############ Pressure wavefield ##############################################
     do i2=3,n2e-2
        do i1=3,n1e-2
           der2=a0*(vx(i1,i2)-vx(i1,i2-1))+a1*(vx(i1,i2+1)-vx(i1,i2-2))
           der1=a0*(vz(i1,i2)-vz(i1-1,i2))+a1*(vz(i1+1,i2)-vz(i1-2,i2))
           px(i1,i2)=damp2a(i2)*px(i1,i2)+kappaw2(i1,i2)*der2
           pz(i1,i2)=damp1a(i1)*pz(i1,i2)+kappaw1(i1,i2)*der1
        end do
     end do

     if ( type_src.eq.0 ) then
        ix1=is(1)
        ix2=is(2)

        if(kws)then

          ialpha1=is(3)
          ialpha2=is(4)
          nalpha_frees=is(5)

          irad1min=MAX(1-ix1,   -irad(1))
          irad1max=MIN(n1e-ix1, irad(1))
          irad2min=MAX(1-ix2,   -irad(2))
          irad2max=MIN(n2e-ix2, irad(2))

          do i1=irad1min,irad1max
             do i2=irad2min,irad2max
                coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                px(ix1+i1,ix2+i2) = px(ix1+i1,ix2+i2) + coeff*dt*s(it)*kappa(ix1+i1,ix2+i2)
             end do
          end do
        else
          px(ix1,ix2)=dt*s(it)*kappa(ix1,ix2)
        endif
     endif
     
     p(:,:)=px(:,:)+pz(:,:)

     if (ofs.eq.1) then						!free surface boundary condition if used
        do i1=1,npml+1
           p(i1,:)=0.
        end do
     end if

     ! ######### Store pressure seismograms and wavefields ##################################
     if (type_rec.eq.0) then
        do i2=1,n2
           do i1=1,n1
              wavefield(i1,i2,it)=p(i1+npml,i2+npml)
           end do
        end do
        
        do irec=1, nr
          ix1=ir(irec,1)
          ix2=ir(irec,2)
 
          if(kws)then

            ialpha1=ir(irec,3)
            ialpha2=ir(irec,4)
            nalpha_frees=ir(irec,5)

            irad1min=MAX(1-ix1,   -irad(1))
            irad1max=MIN(n1e-ix1, irad(1))
            irad2min=MAX(1-ix2,   -irad(2))
            irad2max=MIN(n2e-ix2, irad(2))

            do i1=irad1min,irad1max
               do i2=irad2min,irad2max
                  coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                  sis(it,irec)=sis(it,irec) + coeff*sngl(p(ix1 + i1,ix2 + i2))
               end do
            end do
          else
            sis(it,irec)=sis(it,irec)+sngl(p(ix1,ix2))
          endif

        end do
        
     end if

     ! ########  Particle velocity wavefields ######################################
     do i2=3,n2e-2
        do i1=3,n1e-2
           der2=a0*(p(i1,i2+1)-p(i1,i2))+a1*(p(i1,i2+2)-p(i1,i2-1))
           der1=a0*(p(i1+1,i2)-p(i1,i2))+a1*(p(i1+2,i2)-p(i1-1,i2))
           vx(i1,i2)=damp2b(i2)*vx(i1,i2)+buw2(i1,i2)*der2
           vz(i1,i2)=damp1b(i1)*vz(i1,i2)+buw1(i1,i2)*der1
        end do
     end do

     if (type_src.eq.1) then
        ix1=is(1)
        ix2=is(2)
        if(kws)then
          ialpha1=is(3)
          ialpha2=is(4)
          nalpha_frees=is(5)
          irad1min=MAX(1-ix1,   -irad(1))
          irad1max=MIN(n1e-ix1, irad(1))
          irad2min=MAX(1-ix2,   -irad(2))
          irad2max=MIN(n2e-ix2, irad(2))
          do i1=irad1min,irad1max
             do i2=irad2min,irad2max
                coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                vz(ix1+i1,ix2+i2) = vz(ix1+i1,ix2+i2) + coeff*dt*s(it)*bue(ix1+i1,ix2+i2)
             end do
          end do
        else
          vz(ix1,ix2) = dt*s(it)*bue(ix1,ix2)
        endif
     end if

     if (ofs.eq.1) then				!free surface boundary condition if used
        do i2=1,n2e
           vz(npml,i2)=vz(npml+1,i2)
        end do
     end if

  end do						!end of loop over time

  deallocate (vx)
  deallocate (vz)
  deallocate (p)
  deallocate (px)
  deallocate (pz)
  deallocate (damp1a)
  deallocate (damp1a1)
  deallocate (damp1b)
  deallocate (damp1b1)
  deallocate (damp2a)
  deallocate (damp2a1)
  deallocate (damp2b)
  deallocate (damp2b1)
  deallocate (buw1)
  deallocate (buw2)
  deallocate (kappa)
  deallocate (kappaw2)
  deallocate (kappaw1)
  deallocate (bue)
 
end subroutine submodeling

! ================================================================================
! ================================================================================
! ================================================================================
! SUBROUTINE submodeling_adjoint: explicit time-stepping scheme with finite differences.
! Velocity-stress formulation. O(Delta t**2,Delta x**4).
! Adjoint simulation.
! ================================================================================
! ================================================================================
! ================================================================================

subroutine submodeling_adjoint(is,vp,n1,n2,npml,h,ofs,nt,dt,res,nr,ir,wavefieldp)
  USE kwsinc
  implicit none

  integer                           :: nt,n1,n2,npml
  real                              :: h,dt

  integer                           :: nr,ir(nr,5)

  ! wavefield arrays
  real, allocatable                 :: p(:,:),px(:,:),pz(:,:),vx(:,:),vz(:,:)
  real                              :: der1,der2
  real                              :: wavefieldp(n1,n2,nt),res(nt,nr)

  ! subsurface properties
  real                              :: vp(n1,n2)
  real, allocatable                 :: rho(:,:),vpe(:,:),rhoe(:,:),bue(:,:),bu1(:,:),bu2(:,:),buw1(:,:),buw2(:,:)
  real, allocatable                 :: kappa(:,:),kappaw1(:,:),kappaw2(:,:)

  !PML arrays
  real, allocatable                 :: spg(:),spg1(:)
  real, allocatable                 :: damp1a(:),damp1a1(:),damp2a(:),damp2a1(:),damp1b(:),damp1b1(:),damp2b(:),damp2b1(:)

  real, parameter                   :: xxfac=0.12

  integer                           :: n1e,n2e,is,i1,i2,i1min,i1max,i2min,i2max,it,imax,idir,ishift,irec
  integer                           :: ofs
  real                              :: tt
  real                              :: dtdx

  !hicks parameter
  real                              :: coeff
  integer                           :: irad1min,irad1max,irad2min,irad2max
  integer                           :: ialpha1,ialpha2,nalpha_frees,ix1,ix2

  real, parameter                   :: a0=1.125
  real, parameter                   :: a1=-1./24.

  logical, parameter                :: kws=.true.


  ! ############ Grid parameters #############################################
  n2e=n2+2*npml
  n1e=n1+2*npml
  dtdx=dt/h

  ! ##########################################################
  ! PML sponges (split PML implementation)
  allocate (damp2a(n2+2*npml))
  allocate (damp2a1(n2+2*npml))
  allocate (damp2b(n2+2*npml))
  allocate (damp2b1(n2+2*npml))
  allocate (damp1a(n1+2*npml))
  allocate (damp1a1(n1+2*npml))
  allocate (damp1b(n1+2*npml))
  allocate (damp1b1(n1+2*npml))

  allocate (spg(npml+1))
  allocate (spg1(npml+1))
  call cspongea(spg,spg1,npml,xxfac,damp2a,damp2a1,n2)
  call cspongeb(spg,spg1,npml,xxfac,damp2b,damp2b1,n2)
  call cspongea(spg,spg1,npml,xxfac,damp1a,damp1a1,n1)
  call cspongeb(spg,spg1,npml,xxfac,damp1b,damp1b1,n1)
  deallocate (spg)
  deallocate (spg1)

  ! ############### Staggered grids ###########################################
  allocate (kappa(n1e,n2e))
  allocate (kappaw2(n1e,n2e))
  allocate (kappaw1(n1e,n2e))
  allocate (buw2(n1e,n2e))
  allocate (buw1(n1e,n2e))

  allocate (bu1(n1e,n2e))
  allocate (bu2(n1e,n2e))
  allocate (rho(n1,n2))
  allocate (vpe(n1e,n2e))
  allocate (rhoe(n1e,n2e))
  allocate (bue(n1e,n2e))

  rho(:,:)=1.
  ! write(*,*) 'VMAX = ',vmax
  call submodext(vp,n1,n2,npml,vpe)
  call submodext(rho,n1,n2,npml,rhoe)
  call substagger(vpe,rhoe,bue,kappa,bu2,bu1,n1e,n2e)
  deallocate (vpe)
  deallocate (rho)
  deallocate (bue)
  deallocate (rhoe)
  do i2=1,n2e
     do i1=1,n1e
        kappaw2(i1,i2)=dtdx*kappa(i1,i2)*damp2a1(i2)
        kappaw1(i1,i2)=dtdx*kappa(i1,i2)*damp1a1(i1)
        buw2(i1,i2)=dtdx*bu2(i1,i2)*damp2b1(i2)
        buw1(i1,i2)=dtdx*bu1(i1,i2)*damp1b1(i1)
     end do
  end do
  deallocate (bu2)
  deallocate (bu1)


  ! #################################################################
  ! Memory allocation for wavefiels (p=px+pz due to split pml)
  allocate (p(n1e,n2e))
  allocate (px(n1e,n2e))
  allocate (pz(n1e,n2e))
  allocate (vx(n1e,n2e))
  allocate (vz(n1e,n2e))

  ! ==========================================================
  ! LOOP OVER SHOTS
  ! ==========================================================

  wavefieldp(:,:,:)=0.
  
  p(:,:)=0.	 ! Initial conditions
  px(:,:)=0.
  pz(:,:)=0.
  vx(:,:)=0.
  vz(:,:)=0.
    
  ! ==========================================================
  ! LOOP OVER TIME
  ! ==========================================================

  tt=0.
  do it=1,nt
     tt=tt+dt

     ! ############ Pressure wavefield ##############################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(vx(i1,i2)-vx(i1,i2-1))+a1*(vx(i1,i2+1)-vx(i1,i2-2))
           der1=a0*(vz(i1,i2)-vz(i1-1,i2))+a1*(vz(i1+1,i2)-vz(i1-2,i2))
           px(i1,i2)=damp2a(i2)*px(i1,i2)-kappaw2(i1,i2)*der2
           pz(i1,i2)=damp1a(i1)*pz(i1,i2)-kappaw1(i1,i2)*der1
        end do
     end do
     
     ! pressure adjoint source is applied

     do irec=1, nr
        ix1=ir(irec,1)
        ix2=ir(irec,2)

        if(kws)then
          ialpha1=ir(irec,3)
          ialpha2=ir(irec,4)
          nalpha_frees=ir(irec,5)

          irad1min=MAX(1-ix1,   -irad(1))
          irad1max=MIN(n1e-ix1, irad(1))
          irad2min=MAX(1-ix2,   -irad(2))
          irad2max=MIN(n2e-ix2, irad(2))

          do i1=irad1min,irad1max
             do i2=irad2min,irad2max
                coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                px(ix1+i1,ix2+i2) = px(ix1+i1,ix2+i2) - coeff*dt*res(it,irec)*kappa(ix1+i1,ix2+i2)
             end do
          end do
        else
          px(ix1,ix2)=px(ix1,ix2)-dt*res(it,irec)*kappa(ix1,ix2)
        endif
     end do
     
     ! unsplit pml
     p(:,:)=px(:,:)+pz(:,:)

     if (ofs.eq.1) then						!free surface boundary condition if used
        do i1=1,npml+1
           p(i1,:)=0.
        end do
     end if

     ! ######### Store pressure data and wavefields ##################################
     do i2=1,n2
        do i1=1,n1
           wavefieldp(i1,i2,it)=p(i1+npml,i2+npml)
        end do
     end do

     ! ########  Particle velocity wavefields ######################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(p(i1,i2+1)-p(i1,i2))+a1*(p(i1,i2+2)-p(i1,i2-1))
           der1=a0*(p(i1+1,i2)-p(i1,i2))+a1*(p(i1+2,i2)-p(i1-1,i2))
           vx(i1,i2)=damp2b(i2)*vx(i1,i2)-buw2(i1,i2)*der2
           vz(i1,i2)=damp1b(i1)*vz(i1,i2)-buw1(i1,i2)*der1
        end do
     end do

     if (ofs.eq.1) then				!free surface boundary condition if used
        do i2=1,n2e
           vz(npml,i2)=vz(npml+1,i2)
        end do
     end if

  end do						!end of loop over time

  deallocate (vx)
  deallocate (vz)
  deallocate (p)
  deallocate (px)
  deallocate (pz)
  deallocate (damp1a)
  deallocate (damp1a1)
  deallocate (damp1b)
  deallocate (damp1b1)
  deallocate (damp2a)
  deallocate (damp2a1)
  deallocate (damp2b)
  deallocate (damp2b1)
  deallocate (buw1)
  deallocate (buw2)
  deallocate (kappaw2)
  deallocate (kappaw1)
  deallocate (kappa)
 
end subroutine submodeling_adjoint



subroutine submodeling_adjoint_esource(is,vp,n1,n2,npml,h,ofs,nt,dt,res,type_rec,nr,ir,sis,wavefieldp)
  USE kwsinc
  implicit none

  integer                           :: nt,n1,n2,npml,icount
  real                              :: h,dt,ss

  integer                           :: nr,ir(nr,5)

  ! wavefield arrays
  real, allocatable                 :: p(:,:),px(:,:),pz(:,:),vx(:,:),vz(:,:)
  real                              :: der1,der2
  real                              :: wavefieldp(n1,n2,nt),res(nt,nr),sis(nt,nr)
  real, allocatable                 :: sp(:,:,:),svx(:,:,:),svz(:,:,:)

  ! subsurface properties
  real                              :: vp(n1,n2)
  real, allocatable                 :: rho(:,:),vpe(:,:),rhoe(:,:),bue(:,:),bu1(:,:),bu2(:,:),buw1(:,:),buw2(:,:)
  real, allocatable                 :: kappa(:,:),kappaw1(:,:),kappaw2(:,:)

  !PML arrays
  real, allocatable                 :: spg(:),spg1(:)
  real, allocatable                 :: damp1a(:),damp1a1(:),damp2a(:),damp2a1(:),damp1b(:),damp1b1(:),damp2b(:),damp2b1(:)

  real, parameter                   :: xxfac=0.12

  integer                           :: n1e,n2e,is,i1,i2,it,imax,idir,ishift,irec
  integer                           :: ofs,type_rec
  real                              :: tt,cost
  real                              :: dtdx

  !hicks parameter
  real                              :: coeff
  integer                           :: irad1min,irad1max,irad2min,irad2max
  integer                           :: ialpha1,ialpha2,nalpha_frees,ix1,ix2

  real, parameter                   :: a0=1.125
  real, parameter                   :: a1=-1./24.

  logical, parameter                :: kws=.true.


  ! ############ Grid parameters #############################################
  n2e=n2+2*npml
  n1e=n1+2*npml
  dtdx=dt/h

  ! ##########################################################
  ! PML sponges (split PML implementation)
  allocate (damp2a(n2+2*npml))
  allocate (damp2a1(n2+2*npml))
  allocate (damp2b(n2+2*npml))
  allocate (damp2b1(n2+2*npml))
  allocate (damp1a(n1+2*npml))
  allocate (damp1a1(n1+2*npml))
  allocate (damp1b(n1+2*npml))
  allocate (damp1b1(n1+2*npml))

  allocate (spg(npml+1))
  allocate (spg1(npml+1))
  call cspongea(spg,spg1,npml,xxfac,damp2a,damp2a1,n2)
  call cspongeb(spg,spg1,npml,xxfac,damp2b,damp2b1,n2)
  call cspongea(spg,spg1,npml,xxfac,damp1a,damp1a1,n1)
  call cspongeb(spg,spg1,npml,xxfac,damp1b,damp1b1,n1)
  deallocate (spg)
  deallocate (spg1)

  ! ############### Staggered grids ###########################################
  allocate (kappa(n1e,n2e))
  allocate (kappaw2(n1e,n2e))
  allocate (kappaw1(n1e,n2e))
  allocate (buw2(n1e,n2e))
  allocate (buw1(n1e,n2e))

  allocate (bu1(n1e,n2e))
  allocate (bu2(n1e,n2e))
  allocate (rho(n1,n2))
  allocate (vpe(n1e,n2e))
  allocate (rhoe(n1e,n2e))
  allocate (bue(n1e,n2e))

  rho(:,:)=1.
  call submodext(vp,n1,n2,npml,vpe)
  call submodext(rho,n1,n2,npml,rhoe)
  call substagger(vpe,rhoe,bue,kappa,bu2,bu1,n1e,n2e)

  deallocate (vpe)
  deallocate (rho)
  deallocate (rhoe)

  do i2=1,n2e
     do i1=1,n1e
        kappaw2(i1,i2)=dtdx*kappa(i1,i2)*damp2a1(i2)
        kappaw1(i1,i2)=dtdx*kappa(i1,i2)*damp1a1(i1)
        buw2(i1,i2)=dtdx*bu2(i1,i2)*damp2b1(i2)
        buw1(i1,i2)=dtdx*bu1(i1,i2)*damp1b1(i1)
     end do
  end do
  deallocate (bu2)
  deallocate (bu1)

  ! #################################################################
  ! Memory allocation for wavefiels (p=px+pz due to split pml)
  allocate (p(n1e,n2e))
  allocate (px(n1e,n2e))
  allocate (pz(n1e,n2e))
  allocate (vx(n1e,n2e))
  allocate (vz(n1e,n2e))

  allocate (sp(n1,n2,nt),svx(n1,n2,nt),svz(n1,n2,nt))

  ! ==========================================================
  ! LOOP OVER SHOTS
  ! ==========================================================

  sp(:,:,:)=0.
  svx(:,:,:)=0.
  svz(:,:,:)=0.
  
  p(:,:)=0.      ! Initial conditions
  px(:,:)=0.
  pz(:,:)=0.
  vx(:,:)=0.
  vz(:,:)=0.

  ! ==========================================================
  ! LOOP OVER TIME
  ! ==========================================================

  tt=0.
  do it=1,nt
     tt=tt+dt

     ! ############ Pressure wavefield  ##############################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(vx(i1,i2)-vx(i1,i2-1))+a1*(vx(i1,i2+1)-vx(i1,i2-2))
           der1=a0*(vz(i1,i2)-vz(i1-1,i2))+a1*(vz(i1+1,i2)-vz(i1-2,i2))
           px(i1,i2)=damp2a(i2)*px(i1,i2)-kappaw2(i1,i2)*der2
           pz(i1,i2)=damp1a(i1)*pz(i1,i2)-kappaw1(i1,i2)*der1
        end do
     end do

     ! pressure adjoint source is applied

     do irec=1, nr
        ix1=ir(irec,1)
        ix2=ir(irec,2)

        if(kws)then
          ialpha1=ir(irec,3)
          ialpha2=ir(irec,4)
          nalpha_frees=ir(irec,5)

          irad1min=MAX(1-ix1,   -irad(1))
          irad1max=MIN(n1e-ix1, irad(1))
          irad2min=MAX(1-ix2,   -irad(2))
          irad2max=MIN(n2e-ix2, irad(2))

          do i1=irad1min,irad1max
             do i2=irad2min,irad2max
                coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                px(ix1+i1,ix2+i2) = px(ix1+i1,ix2+i2) - coeff*dt*res(it,irec)*kappa(ix1+i1,ix2+i2)
             end do
          end do
        else
          px(ix1,ix2)=px(ix1,ix2)-dt*res(it,irec)*kappa(ix1,ix2)
        endif
     end do

     ! unsplit pml
     p(:,:)=px(:,:)+pz(:,:)

     if (ofs.eq.1) then                                         !free surface boundary condition if used
        do i1=1,npml+1
           p(i1,:)=0.
        end do
     end if

     ! ######### Store pressure data and wavefields ##################################
     do i2=1,n2
        do i1=1,n1
           sp(i1,i2,it)=p(i1+npml,i2+npml)
        end do
     end do

     ! ########  Particle velocity wavefields ######################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(p(i1,i2+1)-p(i1,i2))+a1*(p(i1,i2+2)-p(i1,i2-1))
           der1=a0*(p(i1+1,i2)-p(i1,i2))+a1*(p(i1+2,i2)-p(i1-1,i2))
           vx(i1,i2)=damp2b(i2)*vx(i1,i2)-buw2(i1,i2)*der2
           vz(i1,i2)=damp1b(i1)*vz(i1,i2)-buw1(i1,i2)*der1
        end do
     end do

     if (ofs.eq.1) then                         !free surface boundary condition if used
        do i2=1,n2e
           vz(npml,i2)=vz(npml+1,i2)
        end do
     end if

     do i2=1,n2
        do i1=1,n1
           svx(i1,i2,it)=vx(i1+npml,i2+npml)
           svz(i1,i2,it)=vz(i1+npml,i2+npml)
        end do
     end do

  end do                                        !end of loop over time

  cost=0

  do it=1, nt
  do i2=1, n2
  do i1=1, n1
     cost=cost+0.5*sp(i1,i2,it)**2
  enddo
  enddo
  enddo

  write(*,*)'cost=',cost

  wavefieldp(:,:,:)=0.
  sis(:,:)=0.

  p(:,:)=0.      ! Initial conditions
  px(:,:)=0.
  pz(:,:)=0.
  vx(:,:)=0.
  vz(:,:)=0.

  ! ==========================================================
  ! LOOP OVER TIME
  ! ==========================================================

  tt=0.
  do it=1,nt
     tt=tt+dt

     ! ############ Pressure wavefield ##############################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(vx(i1,i2)-vx(i1,i2-1))+a1*(vx(i1,i2+1)-vx(i1,i2-2))
           der1=a0*(vz(i1,i2)-vz(i1-1,i2))+a1*(vz(i1+1,i2)-vz(i1-2,i2))
           px(i1,i2)=damp2a(i2)*px(i1,i2)+kappaw2(i1,i2)*der2
           pz(i1,i2)=damp1a(i1)*pz(i1,i2)+kappaw1(i1,i2)*der1
        end do
     end do

     icount=nt-it+1
     ss=-1.

    ! Apply pressure extended source
     do i2=2,n2-1
        do i1=2,n1-1
           px(i1+npml,i2+npml)=px(i1+npml,i2+npml)+ss*dt*sp(i1,i2,icount)*kappa(i1+npml,i2+npml)
        end do
     end do

     p(:,:)=px(:,:)+pz(:,:)

     if (ofs.eq.1) then                                         !free surface boundary condition if used
        do i1=1,npml+1
           p(i1,:)=0.
        end do
     end if

     ! ######### Store seismograms and wavefields
     if (type_rec.eq.0) then
        do i2=1,n2
           do i1=1,n1
              wavefieldp(i1,i2,it)=p(i1+npml,i2+npml)
           end do
        end do

        do irec=1, nr
           ix1=ir(irec,1)
           ix2=ir(irec,2)
           if(kws)then
             ialpha1=ir(irec,3)
             ialpha2=ir(irec,4)
             nalpha_frees=ir(irec,5)

             irad1min=MAX(1-ix1,   -irad(1))
             irad1max=MIN(n1e-ix1, irad(1))
             irad2min=MAX(1-ix2,   -irad(2))
             irad2max=MIN(n2e-ix2, irad(2))

             do i1=irad1min,irad1max
                do i2=irad2min,irad2max
                   coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                   sis(it,irec) = sis(it,irec) + coeff*sngl(p(ix1 + i1,ix2 + i2))
                end do
             end do
           else
             sis(it,irec)=sis(it,irec) + p(ix1,ix2)
           endif
        end do

     end if
 ! ########  Particle velocity wavefields ######################################

     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(p(i1,i2+1)-p(i1,i2))+a1*(p(i1,i2+2)-p(i1,i2-1))
           der1=a0*(p(i1+1,i2)-p(i1,i2))+a1*(p(i1+2,i2)-p(i1-1,i2))
           vx(i1,i2)=damp2b(i2)*vx(i1,i2)+buw2(i1,i2)*der2
           vz(i1,i2)=damp1b(i1)*vz(i1,i2)+buw1(i1,i2)*der1
        end do
     end do

     ! Apply vx extended source
     do i2=1,n2
        do i1=1,n1
           vx(i1+npml,i2+npml)=vx(i1+npml,i2+npml)+ss*dt*svx(i1,i2,icount)*bue(i1+npml,i2+npml)
        end do
     end do

     ! Apply vz extended source
     do i2=1,n2
        do i1=1,n1
           vz(i1+npml,i2+npml)=vz(i1+npml,i2+npml)+ss*dt*svz(i1,i2,icount)*bue(i1+npml,i2+npml)
        end do
     end do

     if (ofs.eq.1) then                         !free surface boundary condition if used
        do i2=1,n2e
           vz(npml,i2)=vz(npml+1,i2)
        end do
     end if
 
  end do                                                !end of loop over time
  
  deallocate (vx)
  deallocate (vz)
  deallocate (p)
  deallocate (px)
  deallocate (pz)
  deallocate (damp1a)
  deallocate (damp1a1)
  deallocate (damp1b)
  deallocate (damp1b1)
  deallocate (damp2a)
  deallocate (damp2a1)
  deallocate (damp2b)
  deallocate (damp2b1)
  deallocate (buw1)
  deallocate (buw2)
  deallocate (kappaw2)
  deallocate (kappaw1)
  deallocate (kappa)
  deallocate (sp,svx,svz)
  deallocate (bue)

end subroutine submodeling_adjoint_esource

subroutine submodeling_adjoint_esource_data(is,vp,n1,n2,npml,h,ofs,nt,dt,res,type_rec,nr,ir,sis)
  USE kwsinc
  implicit none

  integer                           :: nt,n1,n2,npml,icount
  real                              :: h,dt,ss

  integer                           :: nr,ir(nr,5)

  ! wavefield arrays
  real, allocatable                 :: p(:,:),px(:,:),pz(:,:),vx(:,:),vz(:,:)
  real                              :: der1,der2
  real                              :: res(nt,nr),sis(nt,nr)
  real, allocatable                 :: sp(:,:,:),svx(:,:,:),svz(:,:,:)

  ! subsurface properties
  real                              :: vp(n1,n2)
  real, allocatable                 :: rho(:,:),vpe(:,:),rhoe(:,:),bue(:,:),bu1(:,:),bu2(:,:),buw1(:,:),buw2(:,:)
  real, allocatable                 :: kappa(:,:),kappaw1(:,:),kappaw2(:,:)

  !PML arrays
  real, allocatable                 :: spg(:),spg1(:)
  real, allocatable                 :: damp1a(:),damp1a1(:),damp2a(:),damp2a1(:),damp1b(:),damp1b1(:),damp2b(:),damp2b1(:)

  real, parameter                   :: xxfac=0.12

  integer                           :: n1e,n2e,is,i1,i2,it,imax,idir,ishift,irec
  integer                           :: ofs,type_rec
  real                              :: tt,cost
  real                              :: dtdx

  !hicks parameter
  real                              :: coeff
  integer                           :: irad1min,irad1max,irad2min,irad2max
  integer                           :: ialpha1,ialpha2,nalpha_frees,ix1,ix2

  real, parameter                   :: a0=1.125
  real, parameter                   :: a1=-1./24.

  logical, parameter                :: kws=.true.


  ! ############ Grid parameters #############################################
  n2e=n2+2*npml
  n1e=n1+2*npml
  dtdx=dt/h

  ! ##########################################################
  ! PML sponges (split PML implementation)
  allocate (damp2a(n2+2*npml))
  allocate (damp2a1(n2+2*npml))
  allocate (damp2b(n2+2*npml))
  allocate (damp2b1(n2+2*npml))
  allocate (damp1a(n1+2*npml))
  allocate (damp1a1(n1+2*npml))
  allocate (damp1b(n1+2*npml))
  allocate (damp1b1(n1+2*npml))

  allocate (spg(npml+1))
  allocate (spg1(npml+1))
  call cspongea(spg,spg1,npml,xxfac,damp2a,damp2a1,n2)
  call cspongeb(spg,spg1,npml,xxfac,damp2b,damp2b1,n2)
  call cspongea(spg,spg1,npml,xxfac,damp1a,damp1a1,n1)
  call cspongeb(spg,spg1,npml,xxfac,damp1b,damp1b1,n1)
  deallocate (spg)
  deallocate (spg1)

  ! ############### Staggered grids ###########################################
  allocate (kappa(n1e,n2e))
  allocate (kappaw2(n1e,n2e))
  allocate (kappaw1(n1e,n2e))
  allocate (buw2(n1e,n2e))
  allocate (buw1(n1e,n2e))

  allocate (bu1(n1e,n2e))
  allocate (bu2(n1e,n2e))
  allocate (rho(n1,n2))
  allocate (vpe(n1e,n2e))
  allocate (rhoe(n1e,n2e))
  allocate (bue(n1e,n2e))

  rho(:,:)=1.
  call submodext(vp,n1,n2,npml,vpe)
  call submodext(rho,n1,n2,npml,rhoe)
  call substagger(vpe,rhoe,bue,kappa,bu2,bu1,n1e,n2e)

  deallocate (vpe)
  deallocate (rho)
  deallocate (rhoe)

  do i2=1,n2e
     do i1=1,n1e
        kappaw2(i1,i2)=dtdx*kappa(i1,i2)*damp2a1(i2)
        kappaw1(i1,i2)=dtdx*kappa(i1,i2)*damp1a1(i1)
        buw2(i1,i2)=dtdx*bu2(i1,i2)*damp2b1(i2)
        buw1(i1,i2)=dtdx*bu1(i1,i2)*damp1b1(i1)
     end do
  end do
  deallocate (bu2)
  deallocate (bu1)

  ! #################################################################
  ! Memory allocation for wavefiels (p=px+pz due to split pml)
  allocate (p(n1e,n2e))
  allocate (px(n1e,n2e))
  allocate (pz(n1e,n2e))
  allocate (vx(n1e,n2e))
  allocate (vz(n1e,n2e))

  allocate (sp(n1,n2,nt),svx(n1,n2,nt),svz(n1,n2,nt))

  ! ==========================================================
  ! LOOP OVER SHOTS
  ! ==========================================================

  sp(:,:,:)=0.
  svx(:,:,:)=0.
  svz(:,:,:)=0.
  
  p(:,:)=0.      ! Initial conditions
  px(:,:)=0.
  pz(:,:)=0.
  vx(:,:)=0.
  vz(:,:)=0.

  ! ==========================================================
  ! LOOP OVER TIME
  ! ==========================================================

  tt=0.
  do it=1,nt
     tt=tt+dt

     ! ############ Pressure wavefield  ##############################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(vx(i1,i2)-vx(i1,i2-1))+a1*(vx(i1,i2+1)-vx(i1,i2-2))
           der1=a0*(vz(i1,i2)-vz(i1-1,i2))+a1*(vz(i1+1,i2)-vz(i1-2,i2))
           px(i1,i2)=damp2a(i2)*px(i1,i2)-kappaw2(i1,i2)*der2
           pz(i1,i2)=damp1a(i1)*pz(i1,i2)-kappaw1(i1,i2)*der1
        end do
     end do

     ! pressure adjoint source is applied

     do irec=1, nr
        ix1=ir(irec,1)
        ix2=ir(irec,2)

        if(kws)then
          ialpha1=ir(irec,3)
          ialpha2=ir(irec,4)
          nalpha_frees=ir(irec,5)

          irad1min=MAX(1-ix1,   -irad(1))
          irad1max=MIN(n1e-ix1, irad(1))
          irad2min=MAX(1-ix2,   -irad(2))
          irad2max=MIN(n2e-ix2, irad(2))

          do i1=irad1min,irad1max
             do i2=irad2min,irad2max
                coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                px(ix1+i1,ix2+i2) = px(ix1+i1,ix2+i2) - coeff*dt*res(it,irec)*kappa(ix1+i1,ix2+i2)
             end do
          end do
        else
          px(ix1,ix2)=px(ix1,ix2)-dt*res(it,irec)*kappa(ix1,ix2)
        endif
     end do

     ! unsplit pml
     p(:,:)=px(:,:)+pz(:,:)

     if (ofs.eq.1) then                                         !free surface boundary condition if used
        do i1=1,npml+1
           p(i1,:)=0.
        end do
     end if

     ! ######### Store pressure data and wavefields ##################################
     do i2=1,n2
        do i1=1,n1
           sp(i1,i2,it)=p(i1+npml,i2+npml)
        end do
     end do

! ########  Particle velocity wavefields ######################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(p(i1,i2+1)-p(i1,i2))+a1*(p(i1,i2+2)-p(i1,i2-1))
           der1=a0*(p(i1+1,i2)-p(i1,i2))+a1*(p(i1+2,i2)-p(i1-1,i2))
           vx(i1,i2)=damp2b(i2)*vx(i1,i2)-buw2(i1,i2)*der2
           vz(i1,i2)=damp1b(i1)*vz(i1,i2)-buw1(i1,i2)*der1
        end do
     end do

     if (ofs.eq.1) then                         !free surface boundary condition if used
        do i2=1,n2e
           vz(npml,i2)=vz(npml+1,i2)
        end do
     end if

     do i2=1,n2
        do i1=1,n1
           svx(i1,i2,it)=vx(i1+npml,i2+npml)
           svz(i1,i2,it)=vz(i1+npml,i2+npml)
        end do
     end do

  end do                                        !end of loop over time

  sis(:,:)=0.

  p(:,:)=0.      ! Initial conditions
  px(:,:)=0.
  pz(:,:)=0.
  vx(:,:)=0.
  vz(:,:)=0.

  ! ==========================================================
  ! LOOP OVER TIME
  ! ==========================================================

  tt=0.
  do it=1,nt
     tt=tt+dt

     ! ############ Pressure wavefield ##############################################
     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(vx(i1,i2)-vx(i1,i2-1))+a1*(vx(i1,i2+1)-vx(i1,i2-2))
           der1=a0*(vz(i1,i2)-vz(i1-1,i2))+a1*(vz(i1+1,i2)-vz(i1-2,i2))
           px(i1,i2)=damp2a(i2)*px(i1,i2)+kappaw2(i1,i2)*der2
           pz(i1,i2)=damp1a(i1)*pz(i1,i2)+kappaw1(i1,i2)*der1
        end do
     end do

     icount=nt-it+1
     ss=-1.

    ! Apply pressure extended source
     do i2=1,n2
        do i1=1,n1
           px(i1+npml,i2+npml)=px(i1+npml,i2+npml)+ss*dt*sp(i1,i2,icount)*kappa(i1+npml,i2+npml)
        end do
     end do

     p(:,:)=px(:,:)+pz(:,:)

     if (ofs.eq.1) then                                         !free surface boundary condition if used
        do i1=1,npml+1
           p(i1,:)=0.
        end do
     end if

     ! ######### Store seismograms and wavefields
     if (type_rec.eq.0) then

        do irec=1, nr
           ix1=ir(irec,1)
           ix2=ir(irec,2)
           if(kws)then
             ialpha1=ir(irec,3)
             ialpha2=ir(irec,4)
             nalpha_frees=ir(irec,5)

             irad1min=MAX(1-ix1,   -irad(1))
             irad1max=MIN(n1e-ix1, irad(1))
             irad2min=MAX(1-ix2,   -irad(2))
             irad2max=MIN(n2e-ix2, irad(2))

             do i1=irad1min,irad1max
                do i2=irad2min,irad2max
                   coeff = coef_hicks2D_monopole(i1,i2,ialpha1,ialpha2,nalpha_frees)
                   sis(it,irec) = sis(it,irec) + coeff*sngl(p(ix1 + i1,ix2 + i2))
                end do
             end do
           else
             sis(it,irec)=sis(it,irec) + p(ix1,ix2)
           endif
        end do

     end if
 ! ########  Particle velocity wavefields ######################################

     do i2=3,n2e-3
        do i1=3,n1e-3
           der2=a0*(p(i1,i2+1)-p(i1,i2))+a1*(p(i1,i2+2)-p(i1,i2-1))
           der1=a0*(p(i1+1,i2)-p(i1,i2))+a1*(p(i1+2,i2)-p(i1-1,i2))
           vx(i1,i2)=damp2b(i2)*vx(i1,i2)+buw2(i1,i2)*der2
           vz(i1,i2)=damp1b(i1)*vz(i1,i2)+buw1(i1,i2)*der1
        end do
     end do

     ! Apply vx extended source
     do i2=1,n2
        do i1=1,n1
           vx(i1+npml,i2+npml)=vx(i1+npml,i2+npml)+ss*dt*svx(i1,i2,icount)*bue(i1+npml,i2+npml)
        end do
     end do

     ! Apply vz extended source
     do i2=1,n2
        do i1=1,n1
           vz(i1+npml,i2+npml)=vz(i1+npml,i2+npml)+ss*dt*svz(i1,i2,icount)*bue(i1+npml,i2+npml)
        end do
     end do

     if (ofs.eq.1) then                         !free surface boundary condition if used
        do i2=1,n2e
           vz(npml,i2)=vz(npml+1,i2)
        end do
     end if
 
  end do                                                !end of loop over time
  
  deallocate (vx)
  deallocate (vz)
  deallocate (p)
  deallocate (px)
  deallocate (pz)
  deallocate (damp1a)
  deallocate (damp1a1)
  deallocate (damp1b)
  deallocate (damp1b1)
  deallocate (damp2a)
  deallocate (damp2a1)
  deallocate (damp2b)
  deallocate (damp2b1)
  deallocate (buw1)
  deallocate (buw2)
  deallocate (kappaw2)
  deallocate (kappaw1)
  deallocate (kappa)
  deallocate (sp,svx,svz)
  deallocate (bue)

end subroutine submodeling_adjoint_esource_data


! ================================================================================
! ================================================================================
! SUBROUTINE substagger: Build staggered grids
! ================================================================================
! ================================================================================
! ================================================================================

subroutine substagger(vp,rho,bu,kappa,bu2,bu1,n1,n2)
  real vp(n1,n2),rho(n1,n2),bu(n1,n2),bu1(n1,n2),bu2(n1,n2)
  real(kind=4) kappa(n1,n2)

  ! #########################################################################

  do i2=1,n2
     do i1=1,n1
        kappa(i1,i2)=dble(rho(i1,i2))*dble(vp(i1,i2))**2
        if (rho(i1,i2).gt.0.) then
           bu(i1,i2)=1./rho(i1,i2)
        else
           bu(i1,i2)=1e10
        end if
     end do
  end do

  ! #########################################################################

  ! bu2 (grid coincident with Tauzy)

  do i2=1,n2
     do i1=1,n1-1
	bu1(i1,i2)=0.5*(bu(i1+1,i2)+bu(i1,i2))
     end do
     bu1(n1,i2)=bu2(n1-1,i2)
  end do

  ! #########################################################################

  ! bu1 (grid coincident with Tauxy)

  do i2=1,n2-1
     do i1=1,n1
	bu2(i1,i2)=0.5*(bu(i1,i2+1)+bu(i1,i2))
     end do
  end do
  do i1=1,n1
     bu2(i1,n2)=bu1(i1,n2-1)
  end do

  ! #########################################################################

end subroutine substagger

! ================================================================================
! ================================================================================
! ================================================================================
! SUBROUTINE subaugment: Augment grids with PMLs
! Extend a n1 x n2 grid with npml-wide layers
! ================================================================================
! ================================================================================
! ================================================================================

subroutine submodext(v,n1,n2,npml,ve)
  real v(n1,n2),ve(n1+2*npml,n2+2*npml)

  ! kernel

  do i2=1,n2
     do i1=1,n1
	ve(i1+npml,i2+npml)=v(i1,i2)
     end do
  end do

  ! right- and left- sides

  do i2=1,npml
     do i1=1,n1
	ve(i1+npml,i2)=v(i1,1)
	ve(i1+npml,i2+n2+npml)=v(i1,n2)
     end do
  end do

  ! upper- and lower- sides

  do i2=1,n2
     do i1=1,npml
	ve(i1,i2+npml)=v(1,i2)
	ve(npml+n1+i1,i2+npml)=v(n1,i2)
     end do
  end do

  ! corners

  do i2=1,npml
     do i1=1,npml
	ve(i1,i2)=v(1,1)
	ve(i1,npml+n2+i2)=v(1,n2)
	ve(i1+npml+n1,i2)=v(n1,1)
	ve(i1+npml+n1,i2+npml+n2)=v(n1,n2)
     end do
  end do

end subroutine submodext

! ================================================================================
! ================================================================================
! ================================================================================
! SUBROUTINE cspongea: PML functions a
! ================================================================================
! ================================================================================
! ================================================================================

subroutine cspongea(spg,spg1,npml,fac,damp,damp1,n)
  real spg(npml+1),spg1(npml+1),fac,damp(n+2*npml),damp1(n+2*npml)

  do i=1,npml
     x=float(npml-i)
     spg(i)=exp(-(fac*x)**2)
     spg1(i)=exp(-0.5*(fac*x)**2)
  end do


  do i=1,npml
     damp(i)=spg(i)
     damp1(i)=spg1(i)
     damp(npml+n+i)=spg(npml-i+1)
     damp1(npml+n+i)=spg1(npml-i+1)
  end do

  do i=1,n
     damp(i+npml)=1.
     damp1(i+npml)=1.
  end do

end subroutine cspongea

! ================================================================================
! ================================================================================
! ================================================================================
! SUBROUTINE cspongea: PML functions b
! ================================================================================
! ================================================================================
! ================================================================================

subroutine cspongeb(spg,spg1,npml,fac,damp,damp1,n)
  real spg(npml+1),spg1(npml+1),fac,damp(n+2*npml),damp1(n+2*npml)

  do i=1,npml-1
     x=float(npml-1-i)+0.5
     spg(i)=exp(-(fac*x)**2)
     spg1(i)=exp(-0.5*(fac*x)**2)
  end do

  do i=1,npml-1
     damp(i)=spg(i)
     damp1(i)=spg1(i)
  end do

  do i=1,n+1
     damp(i+npml-1)=1.
     damp1(i+npml-1)=1.
  end do

  do i=1,npml
     x=float(i-1)+0.5
     spg(i)=exp(-(fac*x)**2)
     spg1(i)=exp(-0.5*(fac*x)**2)
  end do

  do i=1,npml
     damp(npml+n+i)=spg(i)
     damp1(npml+n+i)=spg1(i)
  end do

end subroutine cspongeb

