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
USE mesh
IMPLICIT NONE
NAMELIST /INPUT/ mesh_name,out_file,tsmax
CHARACTER(len=90) :: inputFile
INTEGER :: funit

! Check to make sure an input file was given
IF (iargc() .EQ. 0) STOP 'Usage: ./cybo inputfile'
CALL GETARG(1,inputFile)

! Read in namelist file and use to read in grid and setup
! arrays for the solver
funit=11
open(unit=funit,file=TRIM(inputFile),form='FORMATTED',status='OLD')
  read(unit=funit,nml=INPUT)
close(funit)

call read_mesh

u = x**2 - y**2
call solver

call write_tec

END PROGRAM

