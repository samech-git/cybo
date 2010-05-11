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
USE euler
USE mesh
USE inputs
IMPLICIT NONE
!DOUBLE PRECISION, DIMENSION(numpts) :: rhop,rhoup,rhovp
INTEGER, DIMENSION(numpts) :: sump
INTEGER :: funit,i
INTEGER :: n1,n2,n3
CHARACTER(LEN=90) :: tecout, file_num
funit = 5

count = count + 1        ! Increment this counter to get a database of files
WRITE(*,*) 'Writing tec file',count
write(file_num,'(I4.4)') count
tecout = adjustr(trim(out_file)) // '_' // &
     adjustr(trim(file_num)) // '.tec'

OPEN(funit,file=tecout,status='unknown')
WRITE(funit,*) 'TITLE = "CYBO output" '
WRITE(funit,*) 'VARIABLES="X","Y","rho","u","v","p" '
WRITE(funit,*) 'ZONE F=FEPOINT,ET=TRIANGLE'
WRITE(funit,*) 'N=',numpts,',E=',numtri

DO i=1,numpts
!   WRITE(funit,*) x(i),y(i),rho(i),rhou(i)/rho(i),rhov(i)/rho(i),p(i)      ! Double precision write 
   WRITE(funit,*) real(x(i)),real(y(i)),real(rho(i)),real(rhou(i)/rho(i)),real(rhov(i)/rho(i)),real(p(i))      ! Single precision write 
END DO
DO i=1,numtri
   WRITE(funit,*) tri(1,i),tri(2,i),tri(3,i)
END DO

CLOSE(funit)


END SUBROUTINE



!! Read in a mesh file from Matlab routines
SUBROUTINE read_mesh
USE inputs, ONLY: mesh_name
USE mesh
IMPLICIT NONE
INTEGER :: funit,i,count,n1,n2,n3
DOUBLE PRECISION :: tmp
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

! Calculate the summed area around each node (note: area spans numpts)
area = 0.d0
DO i=1,numtri
   n1 = tri(1,i)
   n2 = tri(2,i)
   n3 = tri(3,i)

   tmp = .5d0*( x(n1)*(y(n2)-y(n3)) & 
            & +x(n2)*(y(n3)-y(n1))+x(n3)*(y(n1)-y(n2)))

   area(n1) = area(n1) + tmp 
   area(n2) = area(n2) + tmp
   area(n3) = area(n3) + tmp
END DO

END SUBROUTINE
