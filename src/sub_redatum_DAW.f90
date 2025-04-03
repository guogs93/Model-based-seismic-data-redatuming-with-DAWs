subroutine sub_redatum_DAW(grid,inv,src_tab,rec_tab)
  use mpi
  use parameters
  use kwsinc
  use minresModule
  implicit none
  include 'common.h'

  TYPE (grid_type)                  :: grid
  TYPE (inv_type)                   :: inv
  integer                           :: i,j,k
  real, allocatable                 :: wavefieldr(:,:,:),wavefieldesc(:,:,:)
  integer                           :: src_tab(grid%nsrc,5),rec_tab(grid%nrec,5)
  integer                           :: unit
  character(len=160)                :: namei
  integer                           :: isrc
  integer                           :: ista,iend,ierr
 
  allocate(wavefieldr(grid%n1,grid%n2,grid%nt))
  allocate(wavefieldesc(grid%n1,grid%n2,grid%nt))

  ! data samping for Nyquist theory
  
  grid%nto=nint(float(grid%nt-1)*grid%dt/grid%dto)+1

  if(mype==0) write(*,*)'Niquist samping number', grid%nto

  if(mype==0) write(*,*)'start computation'

  call para_range(1,grid%nsrc,nproc,mype,ista,iend)
 
  write(*,*)ista,iend

  do isrc=ista,iend

        allocate(inv%dataobs(grid%nt,grid%nrec))
        allocate(inv%datar(grid%nt,grid%nrec))
        allocate(inv%deltadr(grid%nt,grid%nrec))
        allocate(inv%deltade(grid%nt,grid%nrec))
        allocate(inv%datae(grid%nt,grid%nrec))
        allocate(inv%dataesc(grid%nt,grid%nrec))
        allocate(inv%dataescf(grid%nt,grid%nrec))
        allocate(inv%dataobsf(grid%nt,grid%nrec))

        inv%dataobs(:,:)=0.
        inv%datar(:,:)=0.
        inv%deltadr(:,:)=0
        inv%deltade(:,:)=0
        inv%datae(:,:)=0
        inv%dataesc(:,:)=0
      
        ! =============================================================================
        ! MODE 0: forward simulation
        ! =============================================================================
        if (grid%mode==0) then

           if (mype==0) write(*,*) 'forward simulation'

           call submodeling(isrc,grid%ctrue,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,grid%s, &
                grid%type_src,grid%type_rec,src_tab(isrc,:),rec_tab,grid%nrec,inv%dataobs,wavefieldr)

           call subname(isrc,namei)
           open(NEWUNIT=unit,file='dataobs_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           write(unit) inv%dataobs
           close(unit)
          
           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='wavefield_obs'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) wavefieldr
              close(unit)
           end if
            
        elseif (grid%mode==1) then

           if (mype==0) write(*,*) 'compute data-assimilated wavefields'

           !**********************************************************
           ! step 1: calculate reduced data residuals (delta d^r)
           !**********************************************************

           ! step 1.1: Forward -> reduced wavefield and data (u^r and d^r) 
           call submodeling(isrc,grid%c0,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt,grid%s, &
                grid%type_src,grid%type_rec,src_tab(isrc,:),rec_tab,grid%nrec,inv%datar,wavefieldr)
                
           call subname(isrc,namei)
        
           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='datar_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) inv%datar
              close(unit)

              open(NEWUNIT=unit,file='wavefield_ini'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) wavefieldr
              close(unit)
           end if

           open(NEWUNIT=unit,file='dataobs_org_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           read(unit) inv%dataobs
           close(unit)

           ! step 1.2: calculate reduced data residuals (time is reversed at here) (delta d^r=d^obs-d^r)

           call subresiduals(isrc,inv%dataobs,inv%datar,inv%deltadr,grid%nt,grid%nrec)

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='deltadr_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do i=1, grid%nrec
              do j=1, grid%nt
                  write(unit) inv%deltadr(grid%nt-j+1,i)
              enddo
              enddo
              close(unit)
           end if

           call sublinearminres(isrc,grid,inv,rec_tab)

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='deltade_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do i=1, grid%nrec
              do j=1, grid%nt
                  write(unit) inv%deltade(grid%nt-j+1,i)
              enddo
              enddo
              close(unit)
           end if

           !**********************************************************************************
           ! step 3: adjoint -> calculate source extension  (deltab^e=1\mu*G^T delta d^e)
           !**********************************************************************************

           call submodeling_adjoint_esource(isrc,grid%c0,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt, &
                inv%deltade,grid%type_rec,grid%nrec,rec_tab,inv%dataesc,wavefieldesc)

           call filmartin2d(inv%dataesc,grid%nrec,grid%nt,grid%dt,grid%flow,grid%fhig,grid%llow,grid%lhig,inv%dataescf)

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='dataesc_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) inv%dataesc(:,:)
              close(unit)
           end if

           ! summmation of reduced wave and scatter wave
           wavefieldr(:,:,:)=wavefieldr(:,:,:) + wavefieldesc(:,:,:)
         
           inv%datae(:,:)=inv%datar(:,:) + inv%dataesc(:,:)

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='datae_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) inv%datae
              close(unit)

              open(NEWUNIT=unit,file='wavefield_exd'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) wavefieldr
              close(unit)
           end if

        elseif (grid%mode==3) then

           if (mype==0) write(*,*) 'compute data-assimilated wavefields with Symes method'

           !**********************************************************
           ! step 1: calculate reduced data residuals (delta d^r)
           !**********************************************************

           call subname(isrc,namei)

           open(NEWUNIT=unit,file='dataobs_org_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           read(unit) inv%dataobs
           close(unit)

           ! step 1.2: calculate reduced data residuals (time is reversed at here) (delta d^r=d^obs-d^r)

           !call subresiduals(isrc,inv%dataobs,inv%datar,inv%deltadr,grid%nt,grid%nrec)

           do i=1, grid%nrec
              do j=1, grid%nt
                 inv%deltadr(grid%nt-j+1,i)=inv%dataobs(j,i)
              enddo
           enddo

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='deltadr_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do i=1, grid%nrec
              do j=1, grid%nt
                  write(unit) inv%deltadr(grid%nt-j+1,i)
              enddo
              enddo
              close(unit)
           end if

           call sublinearminres(isrc,grid,inv,rec_tab)

            if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='deltade_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              do i=1, grid%nrec
              do j=1, grid%nt
                  write(unit) inv%deltade(grid%nt-j+1,i)
              enddo
              enddo
              close(unit)
           end if

           !**********************************************************************************
           ! step 3: adjoint -> calculate source extension  (deltab^e=1\mu*G^T delta d^e)
           !**********************************************************************************

           call submodeling_adjoint_esource(isrc,grid%c0,grid%n1,grid%n2,grid%npml,grid%h,grid%ofs,grid%nt,grid%dt, &
                inv%deltade,grid%type_rec,grid%nrec,rec_tab,inv%dataesc,wavefieldesc)

           call filmartin2d(inv%dataesc,grid%nrec,grid%nt,grid%dt,grid%flow,grid%fhig,grid%llow,grid%lhig,inv%dataescf)

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='dataesc_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) inv%dataesc(:,:)
              close(unit)
           end if

           ! summmation of reduced wave and scatter wave
           wavefieldr(:,:,:) = wavefieldesc(:,:,:)

           !inv%datae(:,:) = inv%dataesc(:,:)

           if (grid%oio.eq.1) then
              open(NEWUNIT=unit,file='wavefield_exd'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
              write(unit) wavefieldr
              close(unit)
           end if


        elseif (grid%mode==2)then
 
           call subname(isrc,namei)

           open(NEWUNIT=unit,file='wavefield_exd'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           read(unit) wavefieldr
           close(unit)

           call subwavefieldresample(grid%n1,grid%n2,grid%npml,grid%nt,grid%dt,rec_tab,grid%nrec,inv%dataobs,wavefieldr)

           call filmartin2d(inv%dataobs,grid%nrec,grid%nt,grid%dt,grid%flow,grid%fhig,grid%llow,grid%lhig,inv%dataobsf)

           open(NEWUNIT=unit,file='dataobs_redatum_'//namei(1:LEN_TRIM(namei))//'.bin',access='stream',form='unformatted')
           write(unit) inv%dataobsf
           close(unit)

        end if
 
        deallocate(inv%dataobs)
        deallocate(inv%datar)
        deallocate(inv%deltadr)
        deallocate(inv%deltade)
        deallocate(inv%dataesc)
        deallocate(inv%datae)
        deallocate(inv%dataobsf)
        deallocate(inv%dataescf) 

  end do !end do is

  ! ============================================================================================
  !    Deallocation
  ! ============================================================================================

  deallocate(wavefieldr)
  deallocate(wavefieldesc)

end subroutine sub_redatum_DAW
