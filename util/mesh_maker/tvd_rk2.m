function tvd_rk2(alldata)

% This is the Navier-Stokes solver called by the GUI "Navier2d.m"
%
% The unsteady 2D incompressible Navier-Stokes equations are integrated
% using a 2nd order finite-volume (FV) method for the spatial discretisation 
% and a 1st/2nd order projection method for the time-stepping.
%
% The FV method is based on a non-staggered, vertex centred arrangement
% using the median dual mesh as the control volumes. 2nd order
% discretisations are developed for the divergence/gradient, Laplacian and
% non-linear operators:
%
%   - Divergence/gradient:  Green's theorem
%   - Laplacian          :  Green's theorem applied to face averaged
%                           gradients (unweighted least squares)
%   - Non-linear         :  Upwind piece-wise linear extrapolation with
%                           the HQUICK slope limiter
%
% The non-staggered arrangement was shown to support the common mesh scale
% pressure oscillation and hence the standard 2nd order incremental
% pressure correction method was modified to incorporate Rhie-Chow
% stabilisation. This reduced the order of the time integration to 1st/2nd.
%
% Unlike many other Navier-Stokes solvers this function makes use of a
% complete LU factorisation for the Poisson pressure equation. Using the
% fill reducing UMFPACK algorithm this approach was shown to be far more
% efficient than standard iterative approaches, although this could well be
% limited to MATLAB implementations and non-moving meshes.
%
% The momentum equations are solved in a semi-implicit manner. The
% non-linear terms are evaluated explicitly via a 2nd order TVD
% Runge-Kutta method. The BiCGSTAB method is used to solve the linear
% systems that arise due to the implicit treatment of the viscous terms.
%
%
%   AUTHOR : Darren Engwirda (2005-2006)
%   VERSION: 2.3 (08/05/2006)
%
% If you would like more information please contact me at:
%
%   d_engwirda@hotmail.com
%
% With many useful contributions from Steve Armfield at The University of
% Sydney, Australia.
%
% Naver2d is Copyright (C) 2005-2006 Darren Engwirda. See "copyright.m" for
% full details.
%
%
%   UPDATE 19/04/2006: A HEAP OF CHANGES!! 
%
% (1): The order of the stages in the TVD RK update was switched to put the
% implicit step last (with the pressure correction in the middle). This
% should improve the quality of Neumann BC's for the velocities (outflow,
% symmetry etc).
%
% (2): Inverse distance weighted least squares are now used in the Laplacian.
%
% (3): An advection/diffusion equation for a scalar tracer has been added.
% There is also the option to solve coupled thermal problems where the
% tracer is taken as the temperature. Thermal coupling is done via the
% Oberbeck-Boussinesq quasi-incompressible method.
%
% (4): A module has been added to calculate the [x,y] forces on boundaries
% within the flow, based on the pressure and skin friction effects.
%
% (5): A new pressure outflow boundary condition has been added. The values
% at boundary pressures nodes are extrapolated from the neighbours (setting 
% dP/dn=0). This hopefully should lead to a much less reflective boundary
% condition, although I have noticed some problems when inconsistent IC's
% are used (see readme).
%
% (6): General code structure changes. It's a bit faster and hopefully
% easier to look at??
%
%
%   UPDATE 22/04/2006: BETTER MEMORY OVERHEAD
%
% Based on some failed attempts to use really big meshes (200,000 tri), the
% memory usage has been improved. The most significant change is the use of
% only one Laplacian operator (as opposed to 4 previously). This operator
% is built with Neumann BC's. Dirichlet BC's are imposed approximately
% by adding large terms in the diagonal mass matricies. I have also tried
% to remove a heap of intermediate variables. Anyway, the memory usage is
% now about 1/2 that of previous versions. The maximum mesh size possible is
% computer dependent, although current versions of windows/MATLAB impose a 
% 1.2GB PF limit, independent of the amount of RAM you have...
%
% Big problems often can't fit in physical RAM and make use of virtual
% memory. This is bad. Generally I've seen that the CPU loading decreases
% dramatically (to like 5%) as soon as the memory use spills over into
% virtual. I'm currently trying to think about ways to fix this, but a 
% couple of obvious ones:
%   - get more RAM
%   - don't run any other programs in the background
%   - clear the MATLAB workspace before running navier2d
%
% One possible fix is to switch over to a fully iterative solver for the
% Poisson problem (like an SSOR preconditioned BiCGSTAB). For problems 
% that can use the RAM, I have seen that this can be about 3 times slower 
% than the UMFPACK LU method, BUT, because the LU information would no 
% longer need storage it has the potential to allow bigger meshes to fit 
% in the RAM, giving 100% CPU loads. I could switch over to the fully
% iterative solver once a certain size has been reached...
%
% At the moment, I've decided not to do this, because most computers these
% days have >=1GB of RAM, so MATLAB should run out of PF first. This is not
% always going to be the case (espec for 64bit systems) so if you think you
% are having memory problems, PLEASE CONTACT ME!
%
%
%   UPDATE 02/05/2006: NEW LAPLACIAN (FASTER!!)
%
% (1): A new discretisation of the Laplacian has been implemented. The new method 
% is very similar to a Galerkin FEM, although I (derived) and implemented
% it based on the FVM. The new method uses a much more compact stencil,
% only referencing the nearest-neighbours (as opposed to the next-to-nearest
% like the old method) so the resulting operator is much more sparse. This
% is good because:
%   - Less memory is used to store both the matrix and the LU information.
%   - SMVP's are completed faster.
%   - The LU forward/back solve is completed faster.
%   - The LU decomposition via UMFPACK seems to be stable for lower pivot
%     tolerances (gives sparser LU).
%
% The numerical results seem to be very similar compared to results obtained 
% using the old Laplacian, and based on a mesh refinement study of the steady 
% driven cavity flow the method shows very similar approx 2nd order convergence.
% Overall, the moral of this story is that (basically) identical results
% should now be obtained in about 75% of the CPU time taken previously and
% bigger meshes can be used before running out of RAM/memory.
%
% (2): Better boundary conditions. Extrapolated outflow BC's now use a piecewise
% linear extrapolation for all flow variables. This makes things much smoother,
% i.e the vorticity is now smoothly continuous through the boundary...


% Call the MATLAB memory garbage collection
pack


% Extract data
solver    = alldata.solver;         % Solver type: UVP, tracer, thermal
animation = alldata.animation;      % Animation settings
flag      = alldata.flag;           % Calculate lift/drag forces for edges with flag=true

% Mesh based data
p    = alldata.mesh.p;              % Nodes
t    = alldata.mesh.t;              % Triangulation
bnd  = alldata.mesh.bnd;            % True for bnd nodes
A    = alldata.mesh.A;              % Cell area
e    = alldata.mesh.e;              % Edges
pc   = alldata.mesh.pc;             % Centroids 
e2t  = alldata.mesh.e2t;            % Edge to triangle connectivity
n2n  = alldata.mesh.n2n;            % Node to node connectivity
hnx  = alldata.mesh.hnx;
hny  = alldata.mesh.hny;            
numn = size(p,1);                   % Number of nodes
nume = size(e,1);                   % Number of edges


% Sparse finite volume operators
[L,Gx,Gy,Lp,Up,pp,cp,bc] = build_operators(alldata.mesh,alldata.bc);


% ===================================================================
%                      INTEGRATION SETTINGS
% ===================================================================

nu    = alldata.settings(4);        % Viscosity
mu    = alldata.settings(6);        % Viscosity (tracer)
mmax  = alldata.settings(1);        % Max steps
tmax  = alldata.settings(2);        % Max time
CFL   = alldata.settings(3);        % CFL number
freq  = alldata.settings(5);        % Output frequency
cgtol = 1e-5;                       % Tolerance for BiCGSTAB solver
maxit = 25;                         % Max iters for BiCGSTAB solver  

if solver==3
    alpha = alldata.settings(7);    % Temperature/velocity coupling
else
    alpha = 0;                      % Tracer only - no coupling
end

% Flag BC's 
dbcu     = find(bc(:,1)==1);
dbcv     = find(bc(:,3)==1);
dbcp     = find(bc(:,5) >0);  
dbcs     = find(bc(:,7)==1);  
extrap_p = find(bc(:,5)==2);
extrap_u = extrap_p;                % At this stage only have "extrap"
extrap_v = extrap_p;                % BC's at outflows... I may extend
extrap_s = extrap_p;                % this in the future.

% Initial conditions
U = alldata.init.U;  U(dbcu) = bc(dbcu,2);
V = alldata.init.V;  V(dbcv) = bc(dbcv,4);
S = alldata.init.S;  S(dbcs) = bc(dbcs,8);
P = alldata.init.P;

% Viscous coefficients
diagL = full(diag(L));      % Diagonal terms

% Mass matrcies (use diagonal scaling to get Dirichlet BC's)
M_u = repmat(1,numn,1);  M_u = 2*M_u.*A/nu;  M_u(dbcu) = 1e3 * diagL(dbcu);  M_u(bc(:,1)==0) = 0; 
M_v = repmat(1,numn,1);  M_v = 2*M_v.*A/nu;  M_v(dbcv) = 1e3 * diagL(dbcv);  M_v(bc(:,3)==0) = 0;
M_s = repmat(1,numn,1);  M_s = 2*M_s.*A/mu;  M_s(dbcs) = 1e3 * diagL(dbcs);  M_s(bc(:,7)==0) = 0;


% ===================================================================
%                      NON-LINEAR COEFFICIENTS
% =================================================================== 

% Internal triangles
in = e2t(:,2)>0;

% Edge midpoints
xm = 0.5*(p(e(:,1),1)+p(e(:,2),1));
ym = 0.5*(p(e(:,1),2)+p(e(:,2),2));

% Median edge midpoints
xm1 = 0*xm;  xm1(in) = 0.5*(pc(e2t(in,1),1)+xm(in));  
xm2 = 0*xm;  xm2(in) = 0.5*(pc(e2t(in,2),1)+xm(in));  
ym1 = 0*ym;  ym1(in) = 0.5*(pc(e2t(in,1),2)+ym(in));
ym2 = 0*ym;  ym2(in) = 0.5*(pc(e2t(in,2),2)+ym(in));

% xy from nodes to median edge midpoints
dx = [xm1-p(e(:,1),1), xm1-p(e(:,2),1), xm2-p(e(:,1),1), xm2-p(e(:,2),1)];  
dy = [ym1-p(e(:,1),2), ym1-p(e(:,2),2), ym2-p(e(:,1),2), ym2-p(e(:,2),2)];

% 3rd node in associated triangles
% for each edge
N3 = repmat(0,nume,2);
for k = 1:nume
    if e2t(k,2)>0
        for q = 1:3
            if (t(e2t(k,1),q)~=e(k,1)) && (t(e2t(k,1),q)~=e(k,2))
                N3(k,1) = t(e2t(k,1),q);
            end
            if (t(e2t(k,2),q)~=e(k,1)) && (t(e2t(k,2),q)~=e(k,2))
                N3(k,2) = t(e2t(k,2),q);
            end
        end
    end
end

% COEFFCIENTS FOR LIFT/DRAG CALCS
if any(flag)
    he  = sqrt(hnx(flag,2).^2+hny(flag,2).^2);
    nxe = hnx(flag,2)./he;
    nye = hny(flag,2)./he;
end


% ===================================================================
%                   SETUP ANIMATION FIGURE WINDOWS
% ===================================================================

set(figure(2),'Name','Residuals','Units','Normalized','Position',[0.05 ,0.05 ,0.45,0.35]);
set(figure(3),'Name','UV Quiver','Units','Normalized','Position',[0.525,0.525,0.45,0.35]);
if animation(1)>0
    switch animation(1)
        case 1, set(figure(4),'Name','U Velocity'    ,'Units','Normalized','Position',[0.05,0.525,0.45,0.35]);
        case 2, set(figure(4),'Name','V Velocity'    ,'Units','Normalized','Position',[0.05,0.525,0.45,0.35]);
        case 3, set(figure(4),'Name','Total Velocity','Units','Normalized','Position',[0.05,0.525,0.45,0.35]);
        case 4, set(figure(4),'Name','Pressure'      ,'Units','Normalized','Position',[0.05,0.525,0.45,0.35]);
        case 5, set(figure(4),'Name','abs(Vorticity)','Units','Normalized','Position',[0.05,0.525,0.45,0.35]);
        case 6, set(figure(4),'Name','Tracer'        ,'Units','Normalized','Position',[0.05,0.525,0.45,0.35]);
    end
end
if animation(2)>0
    switch animation(2)
        case 1, set(figure(5),'Name','U Velocity'    ,'Units','Normalized','Position',[0.525,0.05,0.45,0.35]);
        case 2, set(figure(5),'Name','V Velocity'    ,'Units','Normalized','Position',[0.525,0.05,0.45,0.35]);
        case 3, set(figure(5),'Name','Total Velocity','Units','Normalized','Position',[0.525,0.05,0.45,0.35]);
        case 4, set(figure(5),'Name','Pressure'      ,'Units','Normalized','Position',[0.525,0.05,0.45,0.35]);
        case 5, set(figure(5),'Name','abs(Vorticity)','Units','Normalized','Position',[0.525,0.05,0.45,0.35]);
        case 6, set(figure(5),'Name','Tracer'        ,'Units','Normalized','Position',[0.525,0.05,0.45,0.35]);
    end
end
drawnow


% ===================================================================
%                          MAIN LOOP
% ===================================================================

% Memory garbage collection
clear i j pc xm ym xm1 ym1 xm2 ym2 he, pack

% Loop counters
tN     = 0; 
tV     = 0; 
tQ     = 0; 
tstart = cputime;
m      = 0;

% Initial grad(P)
Px = Gx*P;
Py = Gy*P;

% Residuals
resx     = repmat(0,mmax,1);
resy     = repmat(0,mmax,1);
time     = repmat(0,mmax,1);
residual = ~any(flag);

% Viscous time-step limit
if solver==1
    dtvisc = 2*min(A)/nu;
else
    dtvisc = 2*min(A)/max(nu,mu);    
end

% Estimate initial dt from BC's
sqrtA  = sqrt(A);
idt    = max(sqrt(U.^2+V.^2)./sqrtA)/CFL;
dt     = min(1/max(idt,eps), dtvisc);


Snew = S;
while (time(m+1)<tmax)&&(m<mmax)
    
    
    %     A NOTE ON THE HQUICK SCHEME: The non-linear eval is based 
    %     on a piecewise-linear upwind reconstruction at each 
    %     median edge midpoint. A "limiter" (HQUICK) is applied to the 
    %     reconstruction to try to make the scheme "positive" 
    %     (non-oscillatory) in regions of sharp gradient while 
    %     non-diffusive in smooth regions. The current implementation 
    %     is a short cut, and does not apply the limiter in a proper 
    %     multidimensional fashion. This means that the scheme is not 
    %     truly TVD and can admit oscillations. The scheme should be 
    %     bounded (TVB??).
    
    
    % TVD RK3 METHOD - 1ST RK STAGE:
 
    
    % ===================================================================
    %                       TRANSPORT EQUATIONS
    % ===================================================================
    tic
    
    if solver>1     %   Solve advection-diffusion equations 
                    %   for the [U,V] velocity components as 
                    %   well as the scalar tracer S.
        
                
        % NON-LINEAR EVAL (Limited piecewise-linear upwind)
        
        % smvp's
        ux = Gx*U;  uy = Gy*U;
        vx = Gx*V;  vy = Gy*V;
        sx = Gx*S;  sy = Gy*S;
        
        % Assemble non-linear vectors
        [Nu,Nv,Ns] = nonlinear_uvs(U,ux,uy,V,vx,vy,S,sx,sy,e,N3,in,hnx,hny,dx,dy);
        
        % Explicit Euler stage 1
        Ss      = S + dt*( (mu*(L*S) - Ns)./A                );  
        Us      = U + dt*( (nu*(L*U) - Nu)./A - Px           );  
        Vs      = V + dt*( (nu*(L*V) - Nv)./A - Py + alpha*S );  
        Ss(bnd) = S(bnd);
        Us(bnd) = U(bnd);
        Vs(bnd) = V(bnd);
        
        % Extrapolate boundary values where necessary
        if any(extrap_u), Us = bnd_extrap(Us,ux,uy,extrap_u,p,n2n); end
        if any(extrap_v), Vs = bnd_extrap(Vs,vx,vy,extrap_v,p,n2n); end
        if any(extrap_s), Ss = bnd_extrap(Ss,sx,sy,extrap_s,p,n2n); end
        
        
    else    %   Only solve advection-diffusion equations 
            %   for the [U,V] velocity components.

        
        % NON-LINEAR EVAL (Limited piecewise-linear upwind)
        
        % smvp's
        ux = Gx*U;  uy = Gy*U;
        vx = Gx*V;  vy = Gy*V;
        
        % Assemble non-linear vectors
        [Nu,Nv] = nonlinear_uv(U,ux,uy,V,vx,vy,e,N3,in,hnx,hny,dx,dy);
        
        % Explicit Euler stage 1
        Us      = U + dt*( (nu*(L*U) - Nu)./A - Px );  
        Vs      = V + dt*( (nu*(L*V) - Nv)./A - Py ); 
        Us(bnd) = U(bnd);
        Vs(bnd) = V(bnd);
        
        % Extrapolate boundary values where necessary
        if any(extrap_u), Us = bnd_extrap(Us,ux,uy,extrap_u,p,n2n); end
        if any(extrap_v), Vs = bnd_extrap(Vs,vx,vy,extrap_v,p,n2n); end
        
        
    end
    tN = tN + toc;
    
    
    % 2ND RK STAGE
    
    
    % ===================================================================
    %                       TRANSPORT EQUATIONS
    % ===================================================================
    
    if solver>1     %   Solve advection-diffusion equations 
                    %   for the [U,V] velocity components as 
                    %   well as the scalar tracer S.
        
                
        % NON-LINEAR EVAL (Limited piecewise-linear upwind)
        tic
        
        % smvp's
        ux = Gx*Us;  uy = Gy*Us;
        vx = Gx*Vs;  vy = Gy*Vs;
        sx = Gx*Ss;  sy = Gy*Ss;
        
        % Assemble non-linear vectors
        [Nu,Nv,Ns] = nonlinear_uvs(Us,ux,uy,Vs,vx,vy,Ss,sx,sy,e,N3,in,hnx,hny,dx,dy);
        
        tN = tN + toc;
        
        
        % Implicit viscous step (TVD RK2 for non-linear terms)
        tic 
        
        % Mass matrices
        M_u_dt = M_u/dt;
        M_v_dt = M_v/dt;
        M_s_dt = M_s/dt;
        
        % RHS vectors (put BC's in rhs)
        rhss      = ( A.*( S+Ss)/dt                  - Ns )/mu;  
        rhsu      = ( A.*((U+Us)/dt - Px)            - Nu )/nu;   
        rhsv      = ( A.*((V+Vs)/dt - Py + alpha*Ss) - Nv )/nu;  
        rhss(bnd) = M_s_dt(bnd).*bc(bnd,8);
        rhsu(bnd) = M_u_dt(bnd).*bc(bnd,2);
        rhsv(bnd) = M_v_dt(bnd).*bc(bnd,4);
        
        % Solve linear systems
        Us         = myBiCGSTAB(L,M_u_dt,rhsu,Us,cgtol,maxit,diagL);  
        Vs         = myBiCGSTAB(L,M_v_dt,rhsv,Vs,cgtol,maxit,diagL);  
        Snew       = myBiCGSTAB(L,M_s_dt,rhss,Ss,cgtol,maxit,diagL);  
        Us(dbcu)   = U(dbcu);
        Vs(dbcv)   = V(dbcv);
        Snew(dbcs) = S(dbcs);
        
        % Extrapolate boundary values where necessary
        if any(extrap_u), Us   = bnd_extrap(Us  ,ux,uy,extrap_u,p,n2n); end
        if any(extrap_v), Vs   = bnd_extrap(Vs  ,vx,vy,extrap_v,p,n2n); end
        if any(extrap_s), Snew = bnd_extrap(Snew,sx,sy,extrap_s,p,n2n); end
        
        tV = tV + toc;
        
        
    else    %   Only solve advection-diffusion equations 
            %   for the [U,V] velocity components.

        
        % NON-LINEAR EVAL (Limited piecewise-linear upwind)
        tic
        
        % smvp's
        ux = Gx*Us;  uy = Gy*Us;
        vx = Gx*Vs;  vy = Gy*Vs;
        
        % Assemble non-linear vectors
        [Nu,Nv] = nonlinear_uv(Us,ux,uy,Vs,vx,vy,e,N3,in,hnx,hny,dx,dy);
        
        tN = tN + toc;
        
        
        % Implicit viscous step (TVD RK2 for non-linear terms)
        tic 
        
        % Mass matrices
        M_u_dt = M_u/dt;
        M_v_dt = M_v/dt;
        
        % RHS vectors (put BC's in rhs)
        rhsu      = ( A.*((U+Us)/dt - Px) - Nu )/nu;   
        rhsv      = ( A.*((V+Vs)/dt - Py) - Nv )/nu;  
        rhsu(bnd) = M_u_dt(bnd).*bc(bnd,2);
        rhsv(bnd) = M_v_dt(bnd).*bc(bnd,4); 
        
        % Solve linear systems
        Us       = myBiCGSTAB(L,M_u_dt,rhsu,Us,cgtol,maxit,diagL);  
        Vs       = myBiCGSTAB(L,M_v_dt,rhsv,Vs,cgtol,maxit,diagL);  
        Us(dbcu) = U(dbcu);
        Vs(dbcv) = V(dbcv);
        
        % Extrapolate boundary values where necessary
        if any(extrap_u), Us = bnd_extrap(Us,ux,uy,extrap_u,p,n2n); end
        if any(extrap_v), Vs = bnd_extrap(Vs,vx,vy,extrap_v,p,n2n); end
        
        tV = tV + toc;
        
    end  


    % ===================================================================
    %                     PRESSURE POISSON EQUATION
    % ===================================================================
    tic
    
    % Pressure "free" velocity field
    Uh       = Us + dt*Px;  
    Vh       = Vs + dt*Py;  
    Uh(dbcu) = Us(dbcu);
    Vh(dbcv) = Vs(dbcv);
    
    % Divergence
    D       = A.*(Gx*Uh + Gy*Vh)/dt - L*P;  
    D(dbcp) = 0;
    
    % Solve Poisson for pressure correction
    Q  = cp*(Up\(Lp\(pp*D))); 
    Q  = Q-sum(Q)/numn;     % Set mean=0
    Qx = Gx*Q;
    Qy = Gy*Q;
 
    % ===================================================================
    %                        MOMENTUM CORRECTOR
    % ===================================================================
    
    Unew       = Us - dt*Qx;  
    Vnew       = Vs - dt*Qy;  
    Unew(dbcu) = Us(dbcu);
    Vnew(dbcv) = Vs(dbcv);
    
    % Pressure update
    P = P+Q;  Px = Px+Qx;  Py = Py+Qy;
    
    % Extrapolate boundary values where necessary
    if any(extrap_p)
        P  = bnd_extrap(P,Px,Py,extrap_p,p,n2n);
        Px = Gx*P;  
        Py = Gy*P;
    end
    
    
    % Extrapolate boundary values where necessary
    if any(extrap_u), Unew = bnd_extrap(Unew,ux,uy,extrap_u,p,n2n); end
    if any(extrap_v), Vnew = bnd_extrap(Vnew,vx,vy,extrap_v,p,n2n); end
    
    
    tQ = tQ + toc;
    
       
    % ===================================================================
    %                          FULL UPDATE
    % ===================================================================
    
    % Residuals
    if residual     % Velocity residuals
        
        resx(m+1) = sum(abs(Unew-U))/(dt*numn*max(norm(Unew,inf),eps));
        resy(m+1) = sum(abs(Vnew-V))/(dt*numn*max(norm(Vnew,inf),eps));
        
    else            % Lift/drag forces
        
        % Average velocity gradients on bnd edges
        uxe = 0.5*(ux(e(flag,1))+ux(e(flag,2))); uye = 0.5*(uy(e(flag,1))+uy(e(flag,2)));
        vxe = 0.5*(vx(e(flag,1))+vx(e(flag,2))); vye = 0.5*(vy(e(flag,1))+vy(e(flag,2)));
        
        % Form dvt/dn (normal gradient of tangential velocity) on bnd edges
        dvt_dn = nxe.*nye.*(uxe-vye) - (nxe.^2).*vxe + (nye.^2).*uye;
        
        % Average pressure on bnd edges
        Pe = 0.5*(P(e(flag,1))+P(e(flag,2)));
        
        % Drag/lift forces
        resx(m+1) = 2*sum( nu*hny(flag,2).*dvt_dn - hnx(flag,2).*Pe );          % Multiply by 2, remember that [hnx,hny]
        resy(m+1) = 2*sum( nu*hnx(flag,2).*dvt_dn + hny(flag,2).*Pe );          % on the boundary is *0.5...   
        
    end
    
    % Updates
    U = Unew;
    V = Vnew;
    S = Snew;
    
    % CFL based step-size
    idt = max(sqrt(U(~bnd).^2+V(~bnd).^2)./sqrtA(~bnd))/CFL;
    dt  = min(1/max(idt,eps), dtvisc);
    
    
    % ===================================================================
    %                           ANIMATION
    % ===================================================================

    if ~mod(m+1,freq)
        
        figure(2)       % Always plot residuals

        semilogy(time,abs(resx),'b',time,abs(resy),'r'), grid on 
        
        legend('x-residual','y-residual'); xlabel(['Time (sec) at step: ',num2str(m+1)])
        
        
        figure(3)       % Always plot quiver
        
        triquiver(U,V,alldata.mesh);
        
        
        if animation(1)>0, figure(4)    % Variable 1
            switch animation(1)
                case 1, H = U;
                case 2, H = V;
                case 3, H = sqrt(U.^2+V.^2);
                case 4, H = P;
                case 5, H = abs(Gx*V-Gy*U);
                case 6, H = S;
            end
            switch animation(3)
                case 1
                    trimesh(t,p(:,1),p(:,2),H);
                case 2
                    trisurf(t,p(:,1),p(:,2),H); shading interp
                case 3
                    tricontour(H,40,alldata.mesh);
            end
            if animation(4)==1
                view(2), axis equal
            else
                view(3)
            end
        end
        if animation(2)>0, figure(5)    % Variable 2
            switch animation(2)
                case 1, H = U;
                case 2, H = V;
                case 3, H = sqrt(U.^2+V.^2);
                case 4, H = P;
                case 5, H = abs(Gx*V-Gy*U);
                case 6, H = S;
            end
            switch animation(3)
                case 1
                    trimesh(t,p(:,1),p(:,2),H);
                case 2
                    trisurf(t,p(:,1),p(:,2),H); shading interp
                case 3
                    tricontour(H,40,alldata.mesh);
            end
            if animation(4)==1
                view(2), axis equal
            else
                view(3)
            end
        end

        drawnow
        
    end
       
    % Advance counters
    m         = m+1;
    time(m+1) = time(m)+dt;
    
end

% Save back to GUI data
alldata.init.U    = U;
alldata.init.V    = V;
alldata.init.S    = S;
alldata.init.P    = P;
alldata.init.W    = Gx*V-Gy*U;
alldata.init.time = time(1:mmax);
alldata.init.resx = resx(1:mmax);
alldata.init.resy = resy(1:mmax);

set(figure(1),'UserData',alldata);


% Print cost stats 
cost = {
       ['Total time: ',num2str(cputime-tstart),' seconds']
       ''
       ['Total steps: ',num2str(m)]
       ''
       ['Non-linear evals: ',num2str(tN),' seconds']
       ''
       ['Viscous evals: ',num2str(tV),' seconds']
       ''
       ['Pressure evals: ',num2str(tQ),' seconds']
       };
       
msgbox(cost,'Done','none')  

return



% SUB-FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nu,Nv] = nonlinear_uv(u,ux,uy,v,vx,vy,e,N3,in,hnx,hny,dx,dy);

% Assemble the non-linear vectors for the velocities

% Loop over edges
Nu = 0*u; Nv = Nu;
for k = 1:size(e,1)
    if in(k)
        
        % Nodes
        n1 = e(k,1);  n2 = e(k,2); 
        
        % Nodal velocities
        u1 = u(n1);  u2 = u(n2);  um = 5*(u1+u2)/12;
        v1 = v(n1);  v2 = v(n2);  vm = 5*(v1+v2)/12;
        
        % Velocity coeffs (median midpoint)
        um1 = um + u(N3(k,1))/6;  um2 = um + u(N3(k,2))/6;  
        vm1 = vm + v(N3(k,1))/6;  vm2 = vm + v(N3(k,2))/6;  
        
        % Integral face normal velocities
        un1 = hnx(k,1)*um1 + hny(k,1)*vm1;
        un2 = hnx(k,2)*um2 + hny(k,2)*vm2;
        
        % EDGE 1
        if un1<0
            % Extrapolation
            rGu = dx(k,2)*ux(n2)+dy(k,2)*uy(n2);
            rGv = dx(k,2)*vx(n2)+dy(k,2)*vy(n2);
            % Gradient ratios
            if rGu~=0,  Ru = (um1-u2)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm1-v2)/rGv;  else  Rv = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            % Upwind flux
            cu = un1*(u2+psiu);
            cv = un1*(v2+psiv);
        else
            % Extrapolation
            rGu = dx(k,1)*ux(n1)+dy(k,1)*uy(n1);
            rGv = dx(k,1)*vx(n1)+dy(k,1)*vy(n1);
            % Gradient ratios
            if rGu~=0,  Ru = (um1-u1)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm1-v1)/rGv;  else  Rv = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            % Upwind flux
            cu = un1*(u1+psiu);
            cv = un1*(v1+psiv);
        end
        
        % EDGE 2
        if un2<0
            % Extrapolation
            rGu = dx(k,4)*ux(n2)+dy(k,4)*uy(n2);
            rGv = dx(k,4)*vx(n2)+dy(k,4)*vy(n2);
            % Gradient ratios
            if rGu~=0,  Ru = (um2-u2)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm2-v2)/rGv;  else  Rv = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            % Upwind flux
            cu = cu + un2*(u2+psiu);
            cv = cv + un2*(v2+psiv);
        else
            % Extrapolation
            rGu = dx(k,3)*ux(n1)+dy(k,3)*uy(n1);
            rGv = dx(k,3)*vx(n1)+dy(k,3)*vy(n1);
            % Gradient ratios
            if rGu~=0,  Ru = (um2-u1)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm2-v1)/rGv;  else  Rv = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            % Upwind flux
            cu = cu + un2*(u1+psiu);
            cv = cv + un2*(v1+psiv);
        end
        
        % NODAL VECTORS
        Nu(n1) = Nu(n1)+cu;  Nu(n2) = Nu(n2)-cu;
        Nv(n1) = Nv(n1)+cv;  Nv(n2) = Nv(n2)-cv;
        
    end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Nu,Nv,Ns] = nonlinear_uvs(u,ux,uy,v,vx,vy,s,sx,sy,e,N3,in,hnx,hny,dx,dy);

% Assemble the non-linear vectors for the velocities and the scalar field

% Loop over edges
Nu = 0*u; Nv = Nu; Ns = Nu;
for k = 1:size(e,1)
    if in(k)
        
        % Nodes
        n1 = e(k,1);  n2 = e(k,2); 
        
        % Nodal velocities
        u1 = u(n1);  u2 = u(n2);  um = 5*(u1+u2)/12;
        v1 = v(n1);  v2 = v(n2);  vm = 5*(v1+v2)/12;
        s1 = s(n1);  s2 = s(n2);  sm = 5*(s1+s2)/12;
        
        % Velocity coeffs (median midpoint)
        um1 = um + u(N3(k,1))/6;  um2 = um + u(N3(k,2))/6;  
        vm1 = vm + v(N3(k,1))/6;  vm2 = vm + v(N3(k,2))/6;
        sm1 = sm + s(N3(k,1))/6;  sm2 = sm + s(N3(k,2))/6;
        
        % Integral face normal velocities
        un1 = hnx(k,1)*um1 + hny(k,1)*vm1;
        un2 = hnx(k,2)*um2 + hny(k,2)*vm2;
        
        % EDGE 1
        if un1<0
            % Extrapolation
            rGu = dx(k,2)*ux(n2)+dy(k,2)*uy(n2);
            rGv = dx(k,2)*vx(n2)+dy(k,2)*vy(n2);
            rGs = dx(k,2)*sx(n2)+dy(k,2)*sy(n2);
            % Gradient ratios
            if rGu~=0,  Ru = (um1-u2)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm1-v2)/rGv;  else  Rv = 0;  end
            if rGs~=0,  Rs = (sm1-s2)/rGs;  else  Rs = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            if Rs<=0,  psis = 0;  else  psis = 4*Rs*rGs/(3+Rs);  end
            % Upwind flux
            cu = un1*(u2+psiu);
            cv = un1*(v2+psiv);
            cs = un1*(s2+psis);
        else
            % Extrapolation
            rGu = dx(k,1)*ux(n1)+dy(k,1)*uy(n1);
            rGv = dx(k,1)*vx(n1)+dy(k,1)*vy(n1);
            rGs = dx(k,1)*sx(n1)+dy(k,1)*sy(n1);
            % Gradient ratios
            if rGu~=0,  Ru = (um1-u1)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm1-v1)/rGv;  else  Rv = 0;  end
            if rGs~=0,  Rs = (sm1-s1)/rGs;  else  Rs = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            if Rs<=0,  psis = 0;  else  psis = 4*Rs*rGs/(3+Rs);  end
            % Upwind flux
            cu = un1*(u1+psiu);
            cv = un1*(v1+psiv);
            cs = un1*(s1+psis);
        end
        
        % EDGE 2
        if un2<0
            % Extrapolation
            rGu = dx(k,4)*ux(n2)+dy(k,4)*uy(n2);
            rGv = dx(k,4)*vx(n2)+dy(k,4)*vy(n2);
            rGs = dx(k,4)*sx(n2)+dy(k,4)*sy(n2);
            % Gradient ratios
            if rGu~=0,  Ru = (um2-u2)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm2-v2)/rGv;  else  Rv = 0;  end
            if rGs~=0,  Rs = (sm2-s2)/rGs;  else  Rs = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            if Rs<=0,  psis = 0;  else  psis = 4*Rs*rGs/(3+Rs);  end
            % Upwind flux
            cu = cu + un2*(u2+psiu);
            cv = cv + un2*(v2+psiv);
            cs = cs + un2*(s2+psis);
        else
            % Extrapolation
            rGu = dx(k,3)*ux(n1)+dy(k,3)*uy(n1);
            rGv = dx(k,3)*vx(n1)+dy(k,3)*vy(n1);
            rGs = dx(k,3)*sx(n1)+dy(k,3)*sy(n1);
            % Gradient ratios
            if rGu~=0,  Ru = (um2-u1)/rGu;  else  Ru = 0;  end
            if rGv~=0,  Rv = (vm2-v1)/rGv;  else  Rv = 0;  end
            if rGs~=0,  Rs = (sm2-s1)/rGs;  else  Rs = 0;  end
            % HQUICK
            if Ru<=0,  psiu = 0;  else  psiu = 4*Ru*rGu/(3+Ru);  end
            if Rv<=0,  psiv = 0;  else  psiv = 4*Rv*rGv/(3+Rv);  end
            if Rs<=0,  psis = 0;  else  psis = 4*Rs*rGs/(3+Rs);  end
            % Upwind flux
            cu = cu + un2*(u1+psiu);
            cv = cv + un2*(v1+psiv);
            cs = cs + un2*(s1+psis);
        end
        
        % NODAL VECTORS
        Nu(n1) = Nu(n1)+cu;  Nu(n2) = Nu(n2)-cu;
        Nv(n1) = Nv(n1)+cv;  Nv(n2) = Nv(n2)-cv;
        Ns(n1) = Ns(n1)+cs;  Ns(n2) = Ns(n2)-cs;
        
    end
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = bnd_extrap(f,fx,fy,vec,p,n2n) 

% Set nodal values in VEC by average linear extrapolation 
% from neighbours

for iter = 1:4          % Coupled, so do some iters to converge
    for k = 1:length(vec)
        % Take average of neighbours
        z = 1; veck = vec(k); sumf = f(veck);
        while n2n(veck,z)>0
            % Neighbour
            nn = n2n(veck,z);
            % Extrapolation
            sumf = sumf + f(nn) + 0.5*( (p(veck,1)-p(nn,1))*(fx(veck)+fx(nn)) ...
                                      + (p(veck,2)-p(nn,2))*(fy(veck)+fy(nn)) ); 
            % Counter
            z = z+1;
        end
        f(veck) = sumf/z;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,iter] = myBiCGSTAB(V,M_dt,b,x,tol,maxit,diagV)

% BiCGSTAB iterative method for the viscous terms in the transport equations. 
% Jacobi pre-conditioned.

% Scale tolerance wrt RHS
tol = tol*norm(b,inf);

% MAIN LOOP
if tol==0                   % All zero rhs has all zero solution
    
    x = 0*x;  iter = 0;
    
else                        % BiCGSTAB method
    
    % Initial residual
    r     = x;              % Pre-alloc (important for speed...)
    r     = b - M_dt.*x + V*x;
    rdash = r'- sqrt(eps);  % Fudge by sqrt(eps) - gives better convergence??
    
    % Initialise
    a    = 0; p = r;
    v    = r; w = 1;
    rhol = 1;

    % Jacobi preconditioner
    M = 1./(M_dt-diagV);
    
    for iter = 1:maxit
    
        rho = rdash*r;

        rholw = rhol*w;
        if rholw==0 
            error('BiCGSTAB failure')
        end

        % Search vector
        p = r + (rho*a)*(p - w*v)/rholw;

        % HALF ITERATE STEP        
        % Jacobi Preconditioner
        ph = M.*p;     
        v  = M_dt.*ph-V*ph;        
        a  = rho/(rdash*v);
    
        % Residual for half iterate
        s = r - a*v;
    
        % FULL ITERATE STEP
        % Jacobi Preconditioner
        sh = M.*s;   
        t  = M_dt.*sh-V*sh;
        tt = t';
        w  = (tt*s)/(tt*t);
    
        % Solution & residual at full iterate
        x = x + a*ph + w*sh;
        r = s - w*t;
    
        % Break iteration if converged
        if norm(r,inf)<=tol
            return
        end
    
        % For the next step
        rhol = rho;
    
    end
    % iter==maxit to get here
    disp('WARNING: BiCGSTAB did not converge. Try decreasing the CFL number.')
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tricontour(Hn,num,data)

% This is the guts of my TRICONTOUR function from the FEX.

p    = data.p;      % Nodes
e    = data.e;      % Edges
t    = data.t;      % Triangulation
pt   = data.pc;     % Centroids
pe   = data.pe;     % Edge midpoints
e2t  = data.e2t;    % Edge to triangle connectivity    
eINt = data.eINt;   % Triangles as edges
numt = size(t,1);   % Number of triangles

% Endpoint nodes in edges
e1 = e(:,1); e2 = e(:,2);

% Interpolate to centroids and edge midpoints
Ht = (Hn(t(:,1))+Hn(t(:,2))+Hn(t(:,3)))/3;   
He = (Hn(e1)+Hn(e2))/2;  

% Data increment
dH = (max(Ht)-min(Ht))/(num+1);

% MAIN LOOP
x   = []; y = []; z = [];
in  = false(numt,1);
lev = max(Ht);
vec = 1:numt;
old = in;
for v = 1:num       % Loop over contouring levels

    % Find centroid values >= current level
    lev   = lev-dH;
    i     = vec(Ht>=lev);
    i     = i(~old(i));     % Don't need to check triangles from higher levels
    in(i) = true;

    % Locate boundary edges in group
    bnd  = [i; i; i];
    next = 1;
    for k = 1:length(i)
        ct    = i(k);
        count = 0;
        for q = 1:3     % Loop through edges in ct
            ce = eINt(ct,q);
            if ~in(e2t(ce,1)) || ((e2t(ce,2)>0)&&~in(e2t(ce,2)))    
                bnd(next) = ce;     % Found bnd edge
                next      = next+1;
            else
                count = count+1;    % Count number of non-bnd edges in ct
            end
        end
        if count==3                 % If 3 non-bnd edges ct must be in middle of group
            old(ct) = true;         % & doesn't need to be checked for the next level
        end
    end
    bnd = bnd(1:next-1);
    
    % Place nodes approximately on contours by interpolating across bnd edges
    numb  = next-1;
    penew = pe;
    cc    = repmat(0,2*numb,2);
    for k = 1:numb
        
        ce = bnd(k);        % Current edge
        t1 = e2t(ce,1);     % Associated
        t2 = e2t(ce,2);     % triangles
        
        if t2>0     % Move node based on neighbouring centroids
            dx = pt(t2,1)-pt(t1,1);
            dy = pt(t2,2)-pt(t1,2);
            r  = (lev-Ht(t1))/(Ht(t2)-Ht(t1));
        else        % Move node based on internal centroid & bnd midpoint
            dx = pe(ce,1)-pt(t1,1);
            dy = pe(ce,2)-pt(t1,2);
            r  = (lev-Ht(t1))/(He(ce)-Ht(t1));
        end
        
        % New node position
        penew(ce,1) = pt(t1,1)+r*dx;
        penew(ce,2) = pt(t1,2)+r*dy;
        
        % Do a temp connection between adjusted node & endpoint nodes in
        % ce so that the connectivity between neighbouring adjusted nodes
        % can be determined
        m         = 2*k-1;
        cc(m,1)   = e1(ce);
        cc(m,2)   = ce;
        cc(m+1,1) = e2(ce);
        cc(m+1,2) = ce;
        
    end
    
    % Sort connectivity to place connected edges in sucessive rows
    [j,i] = sort(cc(:,1));
    cc    = cc(i,1:2);
    
    % Connect adjacent adjusted nodes
    k    = 1;
    next = 1;
    while k<(2*numb)
        if cc(k,1)==cc(k+1,1)
            cc(next,1) = cc(k,2);
            cc(next,2) = cc(k+1,2);
            next       = next+1;
            k          = k+2;       % Skip over connected edge
        else
            k = k+1;                % Node has only 1 connection - will be picked up above
        end
    end
    cc = cc(1:next-1,1:2);
    
    % xzy data
    xc = [penew(cc(:,1),1),penew(cc(:,2),1)];
    yc = [penew(cc(:,1),2),penew(cc(:,2),2)];
    zc = lev + 0*xc;
    
    % Add contours to list
    x = [x; xc];
    y = [y; yc];
    z = [z; zc];

end

% Draw contours
newplot, patch('Xdata',x','Ydata',y','Zdata',z','Cdata',z','facecolor','none','edgecolor','flat');

axis([min(p(:,1)),max(p(:,1)),min(p(:,2)),max(p(:,2))])

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function triquiver(u,v,data)

% Draw a quiver plot of the velocity field, but do NOT scale the length of
% the arrows with magnitude. Instead colour the arrows with magnitude. Also
% plot a few contours of total velocity using tricontour.

x = data.p(:,1);
y = data.p(:,2);
A = data.A;

uv    = sqrt(u.^2+v.^2);                % Total velocity  
scale = 0.75*sqrt(A)./max(uv,eps);      % Scale arrow length with cell area

xend   = x+scale.*u;
yend   = y+scale.*v;
xarrow = 0.3*x+0.7*xend;
yarrow = 0.3*y+0.7*yend;

newplot, patch('xdata'    ,[x, xend; xarrow+0.2*scale.*v, xend; xarrow-0.2*scale.*v, xend]', ...
               'ydata'    ,[y, yend; yarrow-0.2*scale.*u, yend; yarrow+0.2*scale.*u, yend]', ...
               'cdata'    ,[uv,uv; uv,uv; uv,uv]'                                          , ... 
               'facecolor','none'                                                          , ...
               'edgecolor','flat');

axis equal, axis off, hold on

tricontour(uv,10,data); hold off

return
