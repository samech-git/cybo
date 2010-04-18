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

! Write out the unstructured mesh and variables
! in ASCII tecplot format
SUBROUTINE write_tec
USE globals
IMPLICIT NONE
integer :: funit,count,i
character(len=90) :: outfile,jobname

funit = 5
count = 0
jobname = 'out'
outfile = adjustr(trim(jobname)) // '_' // &
     adjustr(trim('000')) // '.tec'

open(funit,file=outfile,status='unknown')
write(funit,*) 'TITLE = "CYBO output" '
write(funit,*) 'VARIABLES="X","Y","u" '
write(funit,*) 'ZONE F=FEPOINT,ET=TRIANGLE'
write(funit,*) 'N=',numpts,',E=',numtri

u = x**2 - y**2
do i=1,numpts
   write(funit,*) x(i),y(i),u(i)  
end do
do i=1,numtri
   write(funit,*) tri(1,i),tri(2,i),tri(3,i)
end do

close(funit)


END SUBROUTINE



!! Read in a mesh file form Matlab routines
SUBROUTINE read_mesh
USE inputs, only: mesh_name
USE globals
IMPLICIT NONE
integer :: funit,i
character(len=30) :: comments

funit=2
open(unit=funit,file=mesh_name,status='old')
read(funit,*) comments
read(funit,*) numpts, numtri!, numbnd

call allocate_mesh

do i=1,numpts
   read(funit,*) x(i),y(i)
end do

do i=1,numtri
   read(funit,*) tri(1,i),tri(2,i),tri(3,i)
end do

close(funit)

END SUBROUTINE
