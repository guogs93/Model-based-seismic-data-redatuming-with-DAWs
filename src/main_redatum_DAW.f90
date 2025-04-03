!****************************************
! global parameter
!****************************************

module parameters
! ==============================================================================
! PARAMETERS RELATED TO THE FORWARD PROBLEM (grid, acquisition)
! ==============================================================================

    TYPE grid_type
    INTEGER                     :: mode
    INTEGER                     :: n1,n2,npml,nt,nto
    REAL                        :: h,dt,t0,dto,fhig,flow,llow,lhig
    REAL                        :: peakfreq,maxfreq,mu,maxiter,error
    REAL,ALLOCATABLE            :: ctrue(:,:),c0(:,:),s(:)
    INTEGER                     :: type_src,type_rec
    INTEGER                     :: ofs
    INTEGER                     :: oio

    INTEGER                     :: nsrc,nrec
    REAL, ALLOCATABLE           :: xs(:,:),xr(:,:)
    INTEGER, ALLOCATABLE        :: is(:,:),ir(:,:)

    END TYPE grid_type

! ==============================================================================
! PARAMETERS RELATED TO THE INVERSE PROBLEM
! ==============================================================================

    TYPE inv_type
    REAL, ALLOCATABLE           :: dataobs(:,:),datar(:,:),deltadr(:,:)
    REAL, ALLOCATABLE           :: deltade(:,:),datae(:,:),dataesc(:,:)
    REAL, ALLOCATABLE           :: dataescf(:,:),dataobsf(:,:)

    END TYPE inv_type

end module parameters

!**************************************************************************************************
! Seismic data redatuming with data-assimilated wavefield reconstruction (DAW)
! author: Gaoshan guo at geoazur
! Reference: 
! G. Guo, S. Operto, H. Aghamiry. Model-based seismic redatuming with time-domain data-assimilated
! wavefield reconstruction, In Proceedings of the 85th Annual EAGE Meeting (Oslo).
!**************************************************************************************************

PROGRAM main_redatum_DAW
  USE mpi
  USE parameters
  USE kwsinc
  implicit none
  INCLUDE 'common.h'

  TYPE (grid_type)                  :: grid
  TYPE (inv_type)                   :: inv
  integer                           :: infompi,unit
  character(len=160)                :: name_true,name_init,name_acqui,name_src
  real*8                            :: tbegin,tend
  
  ! parameters for hicks
  integer, allocatable              :: src_tab(:,:),rec_tab(:,:)

  call MPI_INIT(infompi)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,infompi)
  call MPI_COMM_RANK(MPI_COMM_WORLD,mype,infompi)

  tbegin=MPI_WTIME()

  ! ###############################################################################
  ! Read input parameters
  ! ###############################################################################

  call subreadpar(name_true,name_init,name_acqui,name_src,grid)

  if (mype.eq.0) then
     write(*,*) ''
     write(*,*) '=============================================================================='
     write(*,*) 'Run mode: ',grid%mode
     write(*,*) 'Initial velocity model file name (name_init)',name_init(1:LEN_TRIM(name_init))
     write(*,*) 'Source wavelets', name_src(1:LEN_TRIM(name_src))
     write(*,*) 'Acquisition file name (name_init)',name_acqui(1:LEN_TRIM(name_acqui))
     write(*,*) 'Space & time grid dimension and interval (n1 n2 h nt dt)',grid%n1,grid%n2,grid%h,grid%nt,grid%dt
     write(*,*) 'Peak frequency of Ricker wavelet (peakfreq)',grid%peakfreq
     write(*,*) 'frequency band', grid%flow,grid%fhig,grid%llow,grid%lhig
     write(*,*) 'Penalty term', grid%mu
  end if

  ! ############################################################################
  ! SET CONSTANTS. CALL SUBHICKS TO COMPUTE coef_hicks
  ! ############################################################################

  call allocatekws2d
  call kwsinc2d(grid%ofs)

  ! ###############################################################################
  ! dynamic memory allocation
  ! ###############################################################################
  
  allocate(grid%c0(grid%n1,grid%n2))
  allocate(grid%ctrue(grid%n1,grid%n2))
  allocate(grid%s(grid%nt))

  ! #####################################  
  ! read in the wavelet
  ! #####################################

  if(mype.eq.0)then
    open(NEWUNIT=unit,file=name_src,access='direct',recl=grid%nt*4)
    read(unit,rec=1) grid%s
    close(unit)
    grid%s(:)=grid%s(:)/grid%h**2
  endif
  call MPI_BCAST(grid%s,grid%nt,MPI_REAL,0,MPI_COMM_WORLD,infompi)

  ! ###############################################################################
  ! open and read input models
  ! ###############################################################################

  if (mype.eq.0) then
     open(NEWUNIT=unit,file=name_true,access='direct',recl=grid%n1*grid%n2*4)
     read(unit,rec=1) grid%ctrue
     close(unit)
     open(NEWUNIT=unit,file=name_init,access='direct',recl=grid%n1*grid%n2*4)
     read(unit,rec=1) grid%c0
     close(unit)
  end if

  call MPI_BCAST(grid%ctrue,grid%n1*grid%n2,MPI_REAL,0,MPI_COMM_WORLD,infompi)
  call MPI_BCAST(grid%c0,grid%n1*grid%n2,MPI_REAL,0,MPI_COMM_WORLD,infompi)

  ! ###############################################################################
  ! read acquisition
  ! ###############################################################################
  
  OPEN(unit=20,file=name_acqui)
  read(20,*)grid%nsrc,grid%nrec
  close(20)

  write(*,*)'Number of sources and receivers',grid%nsrc,grid%nrec

  allocate(src_tab(grid%nsrc,5),rec_tab(grid%nrec,5))

  call subacqui(name_acqui,grid,src_tab,rec_tab)

  ! ###############################################################################
  ! seismic data redatuming with data-assimilated wavefield
  ! ###############################################################################
 
  call sub_redatum_DAW(grid,inv,src_tab,rec_tab)
 
  deallocate(grid%ctrue)
  deallocate(grid%c0)
  deallocate(grid%s)
  deallocate(src_tab,rec_tab)

  ! ###########################################################
  ! Deallocate memory for kws coefficients
  ! ###########################################################

  call deallocatekws2d 

  tend=MPI_WTIME()

  if(mype.eq.0) write(*,*) "ELAPSED TIME = ",tend-tbegin

  call MPI_FINALIZE(infompi)

end program main_redatum_DAW
