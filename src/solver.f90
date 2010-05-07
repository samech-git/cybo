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


CALL allocate_rk4

CALL write_tec

DO k=1,tsmax
   CALL get_dt(dt)
   CALL rk4step(dt)
   WRITE(*,*) 'Step number:',k
   IF (MOD(k,10)==0) THEN
      WRITE(*,*) 'Writing tec file',count
      CALL write_tec
   END IF
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
USE mesh, ONLY: numpts
USE euler
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: dt
DOUBLE PRECISION, DIMENSION(4,numpts) :: flux,phi,tmp1,tmp2
INTEGER :: i

!	Initialize some intermediate arrays
tmp1 = 0.0d0
tmp2 = 0.0d0
phi = 0.0d0


! 5 Stage RK4 step
!!$DO i=1,5
!!$   flux = 0.0d0
!!$   CALL get_flux(flux)
!!$    
!!$   tmp1 =  Ark(i)*phi
!!$   PHI = -dt*flux + tmp1
!!$   
!!$   tmp2 =  Brk(i)*phi
!!$   w(1:4,:) =  w(1:4,:) + tmp2 
!!$   
!!$   CALL get_pressure
!!$   
!!$END DO

! Forward Euler Time stepping
flux = 0.0d0
CALL get_flux(flux)
w(1:4,:) = w(1:4,:) - dt*flux
CALL set_bc_points
CALL get_pressure

END SUBROUTINE



SUBROUTINE get_pressure
USE euler
IMPLICIT NONE

p = gm1*(rhoE/rho - (rhou**2 + rhov**2)/(2.0d0*rho) )
 
END SUBROUTINE

SUBROUTINE get_flux(flux)
USE euler
USE mesh
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(4,numtri), INTENT(INOUT) :: flux
DOUBLE PRECISION, DIMENSION(4) :: fs
INTEGER :: i,n1,n2,t1,t2
INTEGER :: bc,e
DOUBLE PRECISION :: dx,dy,qs1,qs2 
!! Routine to get the flux (residual) for each cell volume

!!         T1        Diagram of how the edge flux is used
!!         /\        to update the total flux of T1 and T2.
!!        /  \
!! Tri 1 /    \        
!!      /      \ 
!!   N1/__Edge__\N2
!!     \        /
!!      \      /  
!! Tri 2 \    /  
!!        \  /   
!!         \/
!!         T2

! Loop over the interior edges to get the flux balance of
! corresponding nodes

DO i=1,size(inter)
   t1 = edg(1,inter(i)) ! Node 1 of tri 1
   n1 = edg(2,inter(i)) ! Node 2 of tri 1/ node 1 of edge
   t2 = edg(3,inter(i)) ! Node 1 of tri 2
   n2 = edg(4,inter(i)) ! Node 2 of tri 2/ node 2 of edge
   
   dx = x(n2) - x(n1)   ! Get dx for edge
   dy = y(n2) - y(n1)   ! Get dy for edge
   
   qs1 = (rhou(n1)*dy - rhov(n1)*dx)/rho(n1)  ! Reused in flux
   qs2 = (rhou(n2)*dy - rhov(n2)*dx)/rho(n2)  ! Reused in flux

   fs(1) = .5d0*(qs1*rho(n1)  + qs2*rho(n2))
   fs(2) = .5d0*(qs1*rhou(n1) + qs2*rhou(n2)) + .5d0*(p(n1) + p(n2))*dy
   fs(3) = .5d0*(qs1*rhov(n1) + qs2*rhov(n2)) - .5d0*(p(n1) + p(n2))*dx
   fs(4) = .5d0*(qs1*(rhoE(n1)+p(n1)) + qs2*(rhoE(n2)+p(n2)))

   ! Add edge fluxes up for each triangle
   flux(:,t1) = flux(:,t1) - fs /area(1)  !** Need areas here ??
   flux(:,t2) = flux(:,t2) + fs /area(1)

END DO
!!$
! Some loop over the boundary edges
DO i=1,size(bound)
   
   e = bound(i)  ! Get the boundary edge index
   t1 = edg(1,e) ! Node 1 of tri 1
   t2 = edg(3,e) ! Node 1 of tri 2
   bc = edg(5,e) ! BC flag to set the flux here

   CALL bound_edge(e,bc,fs) 
   
   ! Add edge fluxes up for each triangle
   IF(t1 .NE. 0) flux(:,t1) = flux(:,t1) - fs /area(1)
   IF(t2 .NE. 0) flux(:,t2) = flux(:,t2) - fs /area(1)
    
END DO

END SUBROUTINE


SUBROUTINE bound_edge(e,bc,fs)
USE inputs
USE mesh, ONLY: x,y,edg,area
USE euler, ONLY: inlet,gm1
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(4), INTENT(OUT) :: fs
INTEGER,INTENT(IN) :: e,bc
INTEGER :: n1,n2,t1,t2,i
DOUBLE PRECISION :: rhoT,rhouT,rhovT,pT,rhoET,qs1,dx,dy

IF (bc == 1) THEN      ! Free stream
   
   t1 = edg(1,e) ! Node 1 of tri 1
   n1 = edg(2,e) ! Node 2 of tri 1/ node 1 of edge
   t2 = edg(3,e) ! Node 1 of tri 2
   n2 = edg(4,e) ! Node 2 of tri 2/ node 2 of edge
   
   dx = x(n2) - x(n1)   ! Get dx for edge
   dy = y(n2) - y(n1)   ! Get dy for edge
   
   !PRINT*,real(dx),real(dy)
   
   rhoT = inlet(1)
   rhouT = inlet(2)
   rhovT = inlet(3)
   pT = inlet(4)
   rhoET = pT/gm1 + (rhouT**2+rhovT**2)/(2.0d0*rhoT)
   
   qs1 = (rhouT*dy - rhovT*dx)/rhoT  ! Reused in flux
      
   fs(1) = qs1*rhoT
   fs(2) = qs1*rhouT + pT*dy
   fs(3) = qs1*rhovT - pT*dx
   fs(4) = qs1*(rhoET+pT)

ELSEIF (bc == 2) THEN  ! Slip wall

   fs = 0.0d0

ELSEIF (bc == 3) THEN  ! Outflow

   fs = 0.0d0

END IF


END SUBROUTINE


SUBROUTINE set_bc_points
USE euler
USE mesh

IMPLICIT NONE
DOUBLE PRECISION :: pT
INTEGER :: i ,n1,n2,e

DO i=1,size(bound)

   e = bound(i)
   n1 = edg(2,e)
   n2 = edg(4,e)

   rho(n1) = inlet(1)
   rhou(n1) = inlet(2)
   rhov(n1) = inlet(3)
   pT = inlet(4)
   rhoE(n1) = pT/gm1 + (rhou(n1)**2+rhov(n1)**2)/(2.0d0*rho(n1))
   
   rho(n2) = inlet(1)
   rhou(n2) = inlet(2)
   rhov(n2) = inlet(3)
   pT = inlet(4)
   rhoE(n2) = pT/gm1 + (rhou(n2)**2+rhov(n2)**2)/(2.0d0*rho(n2))
   
END DO



END SUBROUTINE
