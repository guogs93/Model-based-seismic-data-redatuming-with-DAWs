subroutine subreadpar(name_true,name_init,name_acqui,name_src,grid)
  USE parameters
  IMPLICIT NONE
          
  TYPE (grid_type)   :: grid
  CHARACTER(LEN=160) :: name_true,name_init,name_acqui,name_src
  INTEGER, parameter :: unit=10

  open(unit,file='redatum.par')
  read(unit,*) grid%mode
  read(unit,*) name_true,name_init
  read(unit,*) name_acqui
  read(unit,*) name_src
  read(unit,*) grid%flow,grid%fhig,grid%llow,grid%lhig
  read(unit,*) grid%n1,grid%n2,grid%npml,grid%h,grid%nt,grid%dt,grid%dto
  read(unit,*) grid%type_src,grid%type_rec
  read(unit,*) grid%ofs
  read(unit,*) grid%mu
  read(unit,*) grid%oio
  read(unit,*) grid%error, grid%maxiter
  close(unit)

end subroutine subreadpar

