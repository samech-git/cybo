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
   IF (MOD(k,out_freq)==0) THEN
      CALL write_tec
   END IF
END DO


END SUBROUTINE 



SUBROUTINE get_dt(dt)
USE euler
USE mesh
USE inputs, ONLY: dt_fix
IMPLICIT NONE
DOUBLE PRECISION :: dt

! Some routine that set the dt_max for the grid
dt = dt_fix !!!!!
!!!!!!!!!!!
!!!!!!!!!!!

END SUBROUTINE


SUBROUTINE rk4step(dt)
USE rk4
USE mesh, ONLY: numpts
USE mesh, ONLY: bound
USE euler
IMPLICIT NONE
DOUBLE PRECISION, INTENT(IN) :: dt
DOUBLE PRECISION, DIMENSION(4,numpts) :: flux,phi,tmp1,tmp2
INTEGER :: i

!	Initialize some intermediate arrays
tmp1 = 0.0d0
tmp2 = 0.0d0
phi = 0.0d0

!!$! 5 Stage RK4 step
!!$DO i=1,5
!!$   flux = 0.0d0
!!$   CALL get_flux(flux)
!!$    
!!$   tmp1 =  Ark(i)*phi
!!$   PHI = -dt*flux + tmp1
!!$   
!!$   tmp2 =  Brk(i)*phi
!!$   w(1:4,:) =  w(1:4,:) + tmp2 
!!$   CALL set_bc_points
!!$   CALL get_pressure
!!$END DO
!!$
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

p = gm1*(rhoE - (rhou**2 + rhov**2)/(2.0d0*rho) )

END SUBROUTINE

SUBROUTINE get_flux(flux)
USE euler
USE mesh
USE inputs, ONLY: gamma
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(4,numpts), INTENT(INOUT) :: flux
DOUBLE PRECISION, DIMENSION(numpts) :: div
DOUBLE PRECISION, DIMENSION(4) :: fs,dfs
INTEGER :: i,n1,n2,t1,t2
INTEGER :: bc,e
DOUBLE PRECISION :: dx,dy,qs1,qs2,alpha,c1,c2,len
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
div = 0.0
CALL get_div(div)
DO i=1,size(inter)
   t1 = edg(1,inter(i)) ! Node 1 of tri 1
   n1 = edg(2,inter(i)) ! Node 2 of tri 1/ node 1 of edge
   t2 = edg(3,inter(i)) ! Node 1 of tri 2
   n2 = edg(4,inter(i)) ! Node 2 of tri 2/ node 2 of edge
   
   dx = x(n1) - x(n2)   ! Get dx for edge
   dy = y(n1) - y(n2)   ! Get dy for edge
   
   qs1 = (rhou(n1)*dy - rhov(n1)*dx)/rho(n1)  ! Reused in flux
   qs2 = (rhou(n2)*dy - rhov(n2)*dx)/rho(n2)  ! Reused in flux

   fs(1) = .5d0*(qs1*rho(n1)  + qs2*rho(n2))
   fs(2) = .5d0*(qs1*rhou(n1) + qs2*rhou(n2)) + .5d0*(p(n1) + p(n2))*dy
   fs(3) = .5d0*(qs1*rhov(n1) + qs2*rhov(n2)) - .5d0*(p(n1) + p(n2))*dx
   fs(4) = .5d0*(qs1*(rhoE(n1)+p(n1)) + qs2*(rhoE(n2)+p(n2)))

   ! Add scalar diffusion
   c1 = sqrt( p(n1)*gamma/rho(n1))
   c2 = sqrt( p(n2)*gamma/rho(n2))
   dx = x(t2) - x(t1)   ! Get dx for edge
   dy = y(t2) - y(t1)   ! Get dy for edge
   len = sqrt(dx**2 + dy**2)
   !alpha = ( abs(qs1 + qs2)/2.0d0 + (c1 + c2)/2.0d0 ) * len
   alpha = 50.0* abs( div(n1) + div(n2))/2.0d0 * len**2
   dfs = - alpha/2.0d0*(w(1:4,n1)-w(1:4,n2))
   
   ! Add edge fluxes up for each T point
   flux(:,t1) = flux(:,t1) - fs / area(t1) 
   flux(:,t2) = flux(:,t2) + fs / area(t2)

   ! Add diffusive fluxes for each N point
   flux(:,n1) = flux(:,n1) - dfs / area(n1)
   flux(:,n2) = flux(:,n2) + dfs / area(n2)
   

END DO
!!$
! Some loop over the boundary edges
DO i=1,size(bound)
   
   e = bound(i)  ! Get the boundary edge index
   t1 = edg(1,e) ! Node 1 of tri 1
   n1 = edg(2,e) ! Node 2 of tri 1/ node 1 of edge
   t2 = edg(3,e) ! Node 1 of tri 2
   n2 = edg(4,e) ! Node 2 of tri 2/ node 2 of edge
   bc = edg(5,e) ! BC flag to set the flux here

   CALL bound_edge(e,bc,fs) 
   
   ! Add edge fluxes up for each triangle
 
   IF(bc == 1 .or. bc == 3) THEN ! Free stream or outflow
      IF(t1 .NE. 0) flux(:,t1) = flux(:,t1) - fs/area(t1)
      IF(t2 .NE. 0) flux(:,t2) = flux(:,t2) + fs/area(t2)
   ELSE ! Wall
      IF(t1 .NE. 0) THEN
         flux(:,t1) = flux(:,t1) - fs/area(t1)   ! Contribution to interior node
         flux(:,n1) = flux(:,n1) - fs/area(n1)   ! Contribution to edge node 1
         flux(:,n2) = flux(:,n2) - fs/area(n2)   ! ...edge node 2
      END IF
      IF(t2 .NE. 0) THEN
         flux(:,t2) = flux(:,t2) + fs/area(t2)   ! Contribution to interior node
         flux(:,n1) = flux(:,n1) + fs/area(n1)   ! Contribution to edge node 1
         flux(:,n2) = flux(:,n2) + fs/area(n2)   ! ...edge node 2
      END IF
   END IF
   
    
END DO

END SUBROUTINE


SUBROUTINE bound_edge(e,bc,fs)
USE inputs
USE mesh, ONLY: x,y,edg,area
USE euler
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(4), INTENT(OUT) :: fs
INTEGER,INTENT(IN) :: e,bc
INTEGER :: n1,n2,t1,t2,i
DOUBLE PRECISION :: rhoT,rhouT,rhovT,pT,rhoET,qs1,qs2,dx,dy

IF (bc == 1 .or. bc == 2) THEN  ! Free stream or Slip wall
   t1 = edg(1,e)
   n1 = edg(2,e)
   t2 = edg(3,e)
   n2 = edg(4,e)

   dx = x(n1) - x(n2)   ! Get dx for edge
   dy = y(n1) - y(n2)   ! Get dy for edge
   
   !PRINT*,real(dx),real(dy)
   
   qs1 = (rhou(n1)*dy - rhov(n1)*dx)/rho(n1)  ! Reused in flux
   qs2 = (rhou(n2)*dy - rhov(n2)*dx)/rho(n2)  ! Reused in flux

   if (bc==2) then
      qs1 = 0.0
      qs2 = 0.0
   end if

   fs(1) = .5d0*(qs1*rho(n1)  + qs2*rho(n2))
   fs(2) = .5d0*(qs1*rhou(n1) + qs2*rhou(n2)) + .5d0*(p(n1) + p(n2))*dy
   fs(3) = .5d0*(qs1*rhov(n1) + qs2*rhov(n2)) - .5d0*(p(n1) + p(n2))*dx
   fs(4) = .5d0*(qs1*(rhoE(n1)+p(n1)) + qs2*(rhoE(n2)+p(n2)))

ELSEIF (bc == 3) THEN  ! Outflow

   fs = 0.0d0

END IF


END SUBROUTINE


SUBROUTINE set_bc_points
USE euler
USE mesh

IMPLICIT NONE
DOUBLE PRECISION :: pT
INTEGER :: i ,n1,n2,e,bc

DO i=1,size(bound)

   e = bound(i)
   n1 = edg(2,e)
   n2 = edg(4,e)
   bc = edg(5,e)

   IF (bc == 1) THEN
      
      rho(n1) = inlet(1)
      rhou(n1) = inlet(2)
      rhov(n1) = inlet(3)
      pT = inlet(4)
      p(n1) = pT
      rhoE(n1) = pT/gm1 + (rhou(n1)**2+rhov(n1)**2)/(2.0d0*rho(n1))
   
      rho(n2) = inlet(1)
      rhou(n2) = inlet(2)
      rhov(n2) = inlet(3)
      pT = inlet(4)
      p(n2) = pT
      rhoE(n2) = pT/gm1 + (rhou(n2)**2+rhov(n2)**2)/(2.0d0*rho(n2))
   END IF
   
END DO



END SUBROUTINE


SUBROUTINE get_div(div)
USE euler
USE mesh
IMPLICIT NONE
DOUBLE PRECISION, DIMENSION(numpts), INTENT(INOUT) :: div
DOUBLE PRECISION, DIMENSION(numpts) :: divT
DOUBLE PRECISION :: ds
INTEGER :: i,n1,n2,t1,t2
INTEGER :: bc,e
DOUBLE PRECISION :: dx,dy,qs1,qs2,alpha,c1,c2,len


DO i=1,size(inter)
   t1 = edg(1,inter(i)) ! Node 1 of tri 1
   n1 = edg(2,inter(i)) ! Node 2 of tri 1/ node 1 of edge
   t2 = edg(3,inter(i)) ! Node 1 of tri 2
   n2 = edg(4,inter(i)) ! Node 2 of tri 2/ node 2 of edge
   
   dx = x(n1) - x(n2)   ! Get dx for edge
   dy = y(n1) - y(n2)   ! Get dy for edge
   
   qs1 = (rhou(n1)*dy - rhov(n1)*dx)/rho(n1)  ! Reused in flux
   qs2 = (rhou(n2)*dy - rhov(n2)*dx)/rho(n2)  ! Reused in flux

   ds = .5d0*(qs1*rho(n1)  + qs2*rho(n2))
   
   ! Add edge fluxes up for each T point
   div(t1) = div(t1) - ds / area(t1) 
   div(t2) = div(t2) + ds / area(t2)   

END DO
!!$
! Some loop over the boundary edges
DO i=1,size(bound)
   
   e = bound(i)  ! Get the boundary edge index
   t1 = edg(1,e) ! Node 1 of tri 1
   n1 = edg(2,e) ! Node 2 of tri 1/ node 1 of edge
   t2 = edg(3,e) ! Node 1 of tri 2
   n2 = edg(4,e) ! Node 2 of tri 2/ node 2 of edge

   dx = x(n1) - x(n2)   ! Get dx for edge
   dy = y(n1) - y(n2)   ! Get dy for edge
   
   qs1 = (rhou(n1)*dy - rhov(n1)*dx)/rho(n1)  ! Reused in flux
   qs2 = (rhou(n2)*dy - rhov(n2)*dx)/rho(n2)  ! Reused in flux
   
   ds = .5d0*(qs1*rho(n1)  + qs2*rho(n2))
   
   ! Add edge dives up for each triangle
 
   !IF(bc == 1 .or. bc == 3) THEN ! Free stream or outflow
   !   IF(t1 .NE. 0) div(t1) = div(t1) - ds/area(t1)
   !   IF(t2 .NE. 0) div(t2) = div(t2) + ds/area(t2)
   !ELSE ! Wall
      IF(t1 .NE. 0) THEN
         div(t1) = div(t1) - ds/area(t1)   ! Contribution to interior node
         div(n1) = div(n1) - ds/area(n1)   ! Contribution to edge node 1
         div(n2) = div(n2) - ds/area(n2)   ! ...edge node 2
      END IF
      IF(t2 .NE. 0) THEN
         div(t2) = div(t2) + ds/area(t2)   ! Contribution to interior node
         div(n1) = div(n1) + ds/area(n1)   ! Contribution to edge node 1
         div(n2) = div(n2) + ds/area(n2)   ! ...edge node 2
      END IF
   !END IF

END DO

!!$!! Smooth div field
!!$divT = div
!!$div = 0.0d0
!!$DO i=1,size(inter)
!!$   e = inter(i)  ! Get the boundary edge index
!!$   t1 = edg(1,e) ! Node 1 of tri 1
!!$   n1 = edg(2,e) ! Node 2 of tri 1/ node 1 of edge
!!$   t2 = edg(3,e) ! Node 1 of tri 2
!!$   n2 = edg(4,e) ! Node 2 of tri 2/ node 2 of edge
!!$
!!$   div(t1) = div(t1) + ( divT(n1)*area(n1) + divT(n2)*area(n2) )/area(t1)
!!$   div(t2) = div(t2) + ( divT(n1)*area(n1) + divT(n2)*area(n2) )/area(t2)
!!$
!!$END DO



END SUBROUTINE get_div
