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
USE inputs, ONLY: mesh_name
USE mesh
IMPLICIT NONE
INTEGER :: funit,i,count
CHARACTER(LEN=30) :: comments

funit=2
OPEN(UNIT=funit,FILE=mesh_name,STATUS='old')
READ(funit,*) comments
READ(funit,*) numpts, numtri, numedg

CALL allocate_mesh

! Read in the point locations
DO i=1,numpts
   READ(funit,*) x(i),y(i)
END DO
! Read in the triangle node list (3 per volume)
DO i=1,numtri
   READ(funit,*) tri(1,i),tri(2,i),tri(3,i)
END DO
! Read in the edge nodes (2 per edge + 2 per tri + 1 for bc type)
DO i=1,numedg
   READ(funit,*) edg(1,i),edg(2,i),edg(3,i),edg(4,i),edg(5,i)
END DO
CLOSE(funit)

! Count the Boundary edges
count = 0
DO i=1,numedg
   IF (edg(5,i) .NE. unassigned) THEN
      count = count + 1
   END IF
END DO

ALLOCATE(bound(count))
ALLOCATE(inter(numedg-count))

! Keep track of Boundary edges
count = 0
DO i=1,numedg
   IF (edg(5,i) .NE. unassigned) THEN
      bound(count+1) = i
      count = count + 1
   END IF
END DO

! Keep track of Interior edges
count = 0
DO i=1,numedg
   IF (edg(5,i) .NE. unassigned) THEN
   ELSE
      inter(count+1) = i
      count = count + 1
   END IF
END DO

WRITE(*,*),'Read in:',mesh_name,':'
WRITE(*,*),'Nodes:',numpts
WRITE(*,*),'Triangles:',numtri
WRITE(*,*),'Edges:',numedg,'(',size(bound),'boundary +',size(inter),'interior)'

END SUBROUTINE
