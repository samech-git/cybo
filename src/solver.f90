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
SUBROUTINE solver
USE mesh
USE inputs
USE euler
IMPLICIT NONE
INTEGER :: k
DOUBLE PRECISION :: dt

CALL allocate_euler


DO k=1,tsmax
   CALL get_dt(dt)
   CALL rk4step(dt)
END DO


END SUBROUTINE 



SUBROUTINE get_dt(dt)
USE euler
USE mesh
IMPLICIT NONE
DOUBLE PRECISION :: dt

! Some routine that set the dt_max for the grid
dt = .001d0 !!!!!
!!!!!!!!!!!
!!!!!!!!!!!

END SUBROUTINE


SUBROUTINE rk4step(dt)
USE rk4
USE mesh, ONLY: numtri
USE euler
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: dt
DOUBLE PRECISION, DIMENSION(4,numtri) :: flux,phi,tmp1,tmp2
INTEGER :: i
!	Initialize some intermediate arrays
tmp1 = 0.0d0
tmp2 = 0.0d0
phi = 0.0d0

DO i=1,5
  	
   CALL get_flux(flux)
    
   tmp1 =  Ark(i)*phi
   PHI = -dt*flux + tmp1
   
   tmp2 =  Brk(i)*phi
   w(1:4,:) =  w(1:4,:) + tmp2 
   
   CALL get_pressure
   
END DO

END SUBROUTINE



SUBROUTINE get_pressure
USE euler, ONLY: w

w(5,:) = w(4,:)/w(1,:) !.... get pressure FIX IT
 
END SUBROUTINE

SUBROUTINE get_flux(flux)
USE euler
USE mesh
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(4,numtri), INTENT(OUT) :: flux

!! Routine to get the flux (residual) for each cell volume
flux = 0.0d0

! Some loop over the interior edges to get the flux balance of each
! triangle
!DO i=1,size(inter)
   !edg(3,inter(i))

!END DO


! Some loop over the boundary edges
!DO i=1,size(bound)

! bc_type = edge(5,bound(i))
! CALL bound_edge( 

!END DO



END SUBROUTINE

