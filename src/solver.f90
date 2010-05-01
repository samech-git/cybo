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
dt = .00001d0 !!!!!
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

!CALL get_flux(flux)
!w(1:4,:) = w(1:4,:) - dt*flux

CALL get_pressure


END SUBROUTINE



SUBROUTINE get_pressure
USE euler
IMPLICIT NONE
p = gm1*(rhoE - (rhou**2 + rhov**2)/(2.0d0*rho))
 
END SUBROUTINE

SUBROUTINE get_flux(flux)
USE euler
USE mesh
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(4,numtri), INTENT(OUT) :: flux
DOUBLE PRECISION, DIMENSION(4) :: fs
INTEGER :: i,n1,n2,t1,t2
INTEGER :: bc,e
DOUBLE PRECISION :: dx,dy,qs1,qs2 
!! Routine to get the flux (residual) for each cell volume
flux = 0.0d0

! Some loop over the interior edges to get the flux balance of each
! triangle
DO i=1,size(inter)
   n1 = edg(1,inter(i)) ! Node 1 of the edge
   n2 = edg(2,inter(i)) ! Node 2 of the edge
   t1 = edg(3,inter(i)) ! Tri 1 of the edge
   t2 = edg(4,inter(i)) ! Tri 2 of the edge

   dx = x(n2) - x(n1)   ! Get dx for edge
   dy = y(n2) - y(n1)   ! Get dy for edge

   qs1 = (rhou(t1)*dy - rhov(t1)*dx)/rho(t1)  ! Reused in flux
   qs2 = (rhou(t2)*dy - rhov(t2)*dx)/rho(t2)  ! Reused in flux

   
   fs(1) = .5d0*(qs1*rho(t1)  + qs2*rho(t2))
   fs(2) = .5d0*(qs1*rhou(t1) + qs2*rhou(t2)) + .5d0*(p(t1) + p(t2))*dy
   fs(3) = .5d0*(qs1*rhov(t1) + qs2*rhov(t2)) - .5d0*(p(t1) + p(t2))*dx
   fs(4) = .5d0*(qs1*(rhoE(t1)+p(t1)) + qs2*(rhoE(t2)+p(t2)))

   ! Add edge fluxes up for each triangle
   flux(:,t1) = flux(:,t1) + fs/area(t1)
   flux(:,t2) = flux(:,t2) + fs/area(t2)
END DO

! Some loop over the boundary edges
DO i=1,size(bound)
   
   e = bound(i)
   t1 = edg(3,e) ! Tri 1 of the edge
   t2 = edg(4,e) ! Tri 2 of the edge
   bc = edg(5,e)
   CALL bound_edge(e,bc,fs) 
   
   ! Add edge fluxes up for each triangle
    flux(:,t1) = flux(:,t1) + fs/area(t1)
       
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
   
   n1 = edg(1,e) ! Node 1 of the edge
   n2 = edg(2,e) ! Node 2 of the edge
   t1 = edg(3,e) ! Tri 1 of the edge
   t2 = edg(4,e) ! Tri 2 of the edge

   dx = x(n2) - x(n1)   ! Get dx for edge
   dy = y(n2) - y(n1)   ! Get dy for edge
   
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

ELSEIF (bc == 3) THEN  ! Outflow



END IF


END SUBROUTINE
