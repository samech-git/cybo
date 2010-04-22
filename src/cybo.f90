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
OPEN(UNIT=funit,FILE=TRIM(inputFile),FORM='FORMATTED',STATUS='OLD')
  READ(UNIT=funit,NML=INPUT)
CLOSE(funit)

CALL read_mesh

u = x**2 - y**2
CALL init_field
CALL solver

CALL write_tec

END PROGRAM


SUBROUTINE init_field
USE euler
IMPLICIT NONE
DOUBLE PRECISION :: u_in,v_in,P_in,rho_in

CALL allocate_euler

rho_in = 1.0d0
P_in = 1.0d0
u_in = .001d0
v_in = 0.0d0

inlet(1) = rho_in
inlet(2) = rho_in*u_in
inlet(3) = rho_in*v_in
inlet(4) = p_in

rho = rho_in
rhou = rho_in*u_in
rhov = rho_in*v_in
p = P_in
rhoE = P_in / gm1 + rho_in*u_in**2*v_in**2 / 2.0d0


CALL get_pressure



END SUBROUTINE
