!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!        /\        CYBO-Unstuctured Euler Solver  !!!
!!!       /  \                                      !!!
!!!      /    \      Chris M. Yu & Britton J. Olson !!!
!!!     /______\     Department on Aero/Astro       !!!
!!!    /\      /\    Stanford University            !!!
!!!   /  \    /  \   Spring 2010                    !!!
!!!  /    \  /    \                                 !!! 
!!! /______\/______\                                !!!
!!!                                                 !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!             

PROGRAM cybo
USE inputs
IMPLICIT NONE
NAMELIST /INPUT/ mesh_name
character(len=90) :: inputFile
integer :: funit

! Check to make sure an input file was given
if (iargc() .eq. 0) stop 'Usage: ./cybo inputfile'
call getarg(1,inputFile)

! Read in namelist file and use to read in grid and setup
! arrays for the solver
funit=11
open(unit=funit,file=TRIM(inputFile),form='FORMATTED',status='OLD')
  read(unit=funit,nml=INPUT)
close(funit)

call read_mesh
call write_tec

END PROGRAM


SUBROUTINE allocate_mesh
USE globals

allocate(x(numpts))
allocate(y(numpts))
allocate(u(numpts))
allocate(tri(3,numtri))

END SUBROUTINE
