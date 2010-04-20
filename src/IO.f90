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
USE mesh
USE inputs
IMPLICIT NONE
integer :: funit,count,i
character(len=90) :: tecout, file_num
funit = 5

count = 2        ! Increment this counter to get a database of files
write(file_num,'(I4.4)') count
tecout = adjustr(trim(out_file)) // '_' // &
     adjustr(trim(file_num)) // '.tec'

open(funit,file=tecout,status='unknown')
write(funit,*) 'TITLE = "CYBO output" '
write(funit,*) 'VARIABLES="X","Y","u" '
write(funit,*) 'ZONE F=FEPOINT,ET=TRIANGLE'
write(funit,*) 'N=',numpts,',E=',numtri


do i=1,numpts
   !write(funit,*) real(x(i)),real(y(i)),real(u(i)) ! Single precision write
   write(funit,*) x(i),y(i),u(i)                   ! Double precision write 
end do
do i=1,numtri
   write(funit,*) tri(1,i),tri(2,i),tri(3,i)
end do

close(funit)


END SUBROUTINE



!! Read in a mesh file from Matlab routines
SUBROUTINE read_mesh
USE inputs, only: mesh_name
USE mesh
IMPLICIT NONE
integer :: funit,i
character(len=30) :: comments

funit=2
open(unit=funit,file=mesh_name,status='old')
read(funit,*) comments
read(funit,*) numpts, numtri, numedg

call allocate_mesh

! Read in the point locations
do i=1,numpts
   read(funit,*) x(i),y(i)
end do
! Read in the triangle node list (3 per volume)
do i=1,numtri
   read(funit,*) tri(1,i),tri(2,i),tri(3,i)
end do
! Read in the edge nodes (2 per edge + 1 for bc type)
do i=1,numedg
   read(funit,*) edg(1,i),edg(2,i),edg(3,i)
end do


close(funit)

END SUBROUTINE
