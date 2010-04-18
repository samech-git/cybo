function [L,Gx,Gy,Lp,Up,pp,cp,bc] = build_operators(data,edgebc)

% Assemble all of the sparse FV operators and nodal boundary condition arrays
% used in the Navier-Stokes solver tvd_rk2.m.
%
% Darren Engwirda - 2006
%
% Naver2d is Copyright (C) 2005-2006 Darren Engwirda. See "copyright.m" for
% full details.
%
%
%   22/04/2006:
%
% The memory management has been substantially improved. Really big
% problems were starting to either exceed the 1.2GB PF limit imposed by
% MATLAB or wander into the virtual memory allocation (bringing the CPU 
% load down to 5%). The build process has been partitioned to try to mitigate
% this and some variables are written to disk, cleared and then re-loaded
% later, hopefully freeing the physical RAM and PF.
%
%
%   27/04/2006:
%
% A new Laplacian - see tvd_rk2 for more details.


tol = 0.02;     % Pivot tolerance for UMFPACK. Increase if
                % stability problems occur.


w = waitbar(0,'Building coefficient matrices and LU factors...');


% Setup nodal boundary condition arrays
[bc,corner] = boundary(edgebc,data);


% COEFFICIENT MATRICES

% Form divergence operators
[Gx,Gy] = trigrad(data,corner); waitbar(0.25,w); 

% Laplacian for viscous terms (build with Neumann BC's)
numn = size(data.p,1);
L    = del(data,repmat(0,numn,1)); waitbar(0.5,w);

% Save what we have
save temp.mat bc Gx Gy L

% Clean workspace
clear edgebc corner Gx Gy


% Modify boundary elements (if necessary) 
% for Poisson problem
if any(bc(:,5)>0)
    
    % Extract non-zeros and destroy old L
    [i,j,s] = find(L); clear L
    
    dbcp = bc(:,5)>0;   % True for Dirichlet nodes
    num  = find(dbcp);  % Node numbers
    
    % New Laplacian
    ok = ~dbcp(i);
    i  = [i(ok); num]; 
    j  = [j(ok); num];
    s  = [s(ok); repmat(1,length(num),1)];
    L  = sparse(i,j,s,numn,numn);
    
end

% LU factors for Poisson problem
[Lp,Up,pp,cp] = lu(L,tol); waitbar(0.75,w);  

% Clean workspace
clear L dbcp num ok i j s


% Re-load others
load temp.mat; waitbar(1,w);

% Remove temp file
delete temp.mat


close(w);  drawnow



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bc,corner] = boundary(edgebc,data)

% Setup nodal BC arrays. Flag "corner" nodes that have discontinuous BC's.

e    = data.e;      % Edges
be   = data.be;     % Boundary edges
numn = max(e(:));   % Number of nodes        

% Allocate BC and corner arrays
%
% bc = [Utype,Uval,Vtype,Vval,Ptype,Pval,Stype,Sval]
%
% where type==1 indicates a Dirichlet BC
%       type==0 indicates a Neumann BC
%
% corner = true for nodes with discontinous BC's

bc     = repmat(-1,numn,8);
corner = false(numn,1);

for k = 1:length(be)
    
    % Current edge
    ce = be(k);
    % End nodes
    n1 = e(ce,1); n2 = e(ce,2);
    
    % Node 1
    if bc(n1,1)==-1     % Unassigned
        % Take edge BC's
        bc(n1,:) = edgebc(ce,:);
    else
        for m = 1:4     % Loop through U,V,P,S BC's
            t = 2*m-1;  % Type
            v = 2*m;    % Value    
            if bc(n1,t)==0 && edgebc(ce,t)>0
                % Take Dirichlet BC from edge
                bc(n1,t) = edgebc(ce,t);
                bc(n1,v) = edgebc(ce,v);
            elseif bc(n1,t)==1 && edgebc(ce,t)>0
                if bc(n1,v)~=edgebc(ce,v)
                    % Take average Dirichlet
                    bc(n1,v)   = 0.5*(bc(n1,v)+edgebc(ce,v));
                    corner(n1) = true;
                end
            end
        end
    end
    
    % Node 2
    if bc(n2,1)==-1     % Unassigned
        % Take edge BC's
        bc(n2,:) = edgebc(ce,:);
    else
        for m = 1:4     % Loop through U,V,P,S BC's
            t = 2*m-1;  % Type
            v = 2*m;    % Value    
            if bc(n2,t)==0 && edgebc(ce,t)>0
                % Take Dirichlet BC from edge
                bc(n2,t) = edgebc(ce,t);
                bc(n2,v) = edgebc(ce,v);
            elseif bc(n2,t)==1 && edgebc(ce,t)>0
                if bc(n2,v)~=edgebc(ce,v)
                    % Take average Dirichlet
                    bc(n2,v)   = 0.5*(bc(n2,v)+edgebc(ce,v));
                    corner(n2) = true;
                end
            end
        end
    end
    
end

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Gx,Gy] = trigrad(data,corner)

% Build divergence/gradient operators via a 2nd order FV discretisation of
% Green's theorem on the median mesh. Take special care at "corner" 
% nodes to ensure that the correct velocity flux is taken along the boundary.

t    = data.t;          % Triangulation
e    = data.e;          % Edges
n2e  = data.n2e;        % Node to edge connectivity    
n2n  = data.n2n;        % Node to node connectivity
e2t  = data.e2t;        % Edge to triangle connectivity
hnx  = data.hnx;        % Median edge
hny  = data.hny;        % normals
A    = data.A;          % Cell area
numn = size(n2n,1);     % Number of nodes

% Max number of entries in [i,j,sx,sy]
num = 10*sum(sum(n2n>0,2));

% Pre-alloc
i = repmat(uint32(0),num,1); j = i; sx = double(i); sy = sx; 

next = 1;
for k = 1:numn
    m = 1;
    while n2n(k,m)>0    % Loop around neighbours of node k
        
        cn = n2n(k,m);  % Current node
        ce = n2e(k,m);  % Current edge
        Ac = A(k);      % Current area
        t1 = e2t(ce,1);
        t2 = e2t(ce,2);
        
        % Get median edge normals (ensure correct sign on boundary)
        if k==e(ce,1)
            hnx1 = hnx(ce,1); hny1 = hny(ce,1);
            hnx2 = hnx(ce,2); hny2 = hny(ce,2);
        else
            hnx1 = -hnx(ce,1); hny1 = -hny(ce,1);
            if e2t(ce,2)>0
                hnx2 = -hnx(ce,2); hny2 = -hny(ce,2); 
            else
                hnx2 = hnx(ce,2); hny2 = hny(ce,2);
            end
        end

        % EDGE 1
        % Self term
        i(next)  = k; 
        j(next)  = k;  
        sx(next) = 0.25*hnx1/Ac; 
        sy(next) = 0.25*hny1/Ac; 
        next     = next+1;
        % Neighbour term
        i(next)  = k; 
        j(next)  = cn; 
        sx(next) = 0.25*hnx1/Ac; 
        sy(next) = 0.25*hny1/Ac; 
        next     = next+1;
        % Centroid term
        for q = 1:3
            i(next)  = k; 
            j(next)  = t(t1,q); 
            sx(next) = hnx1/(6*Ac); 
            sy(next) = hny1/(6*Ac); 
            next     = next+1; 
        end
        
        % EDGE 2
        if t2>0         % ce not boundary
            % Self term
            i(next)  = k; 
            j(next)  = k;  
            sx(next) = 0.25*hnx2/Ac; 
            sy(next) = 0.25*hny2/Ac; 
            next     = next+1;
            % Neighbour term
            i(next)  = k; 
            j(next)  = cn; 
            sx(next) = 0.25*hnx2/Ac; 
            sy(next) = 0.25*hny2/Ac; 
            next     = next+1;
            % Centroid term
            for q = 1:3
                i(next)  = k; 
                j(next)  = t(t2,q); 
                sx(next) = hnx2/(6*Ac); 
                sy(next) = hny2/(6*Ac); 
                next     = next+1; 
            end
        else        % ce is boundary
            if corner(k)
                i(next)  = k; 
                j(next)  = cn; 
                sx(next) = hnx2/Ac; 
                sy(next) = hny2/Ac; 
                next     = next+1;
            elseif corner(cn)
                i(next)  = k; 
                j(next)  = k;  
                sx(next) = hnx2/Ac; 
                sy(next) = hny2/Ac; 
                next     = next+1;
            else
                % Self term
                i(next)  = k; 
                j(next)  = k;  
                sx(next) = 0.75*hnx2/Ac; 
                sy(next) = 0.75*hny2/Ac; 
                next     = next+1;
                % Neighbour term
                i(next)  = k; 
                j(next)  = cn; 
                sx(next) = 0.25*hnx2/Ac; 
                sy(next) = 0.25*hny2/Ac; 
                next     = next+1;
            end
        end

        m = m+1;

    end
end
% Can be over allocated
i = double(i(1:next-1)); j = double(j(1:next-1)); sx = sx(1:next-1); sy = sy(1:next-1);

% Sparse operators
Gy = sparse(i,j,sy,numn,numn);
Gx = sparse(i,j,sx,numn,numn);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = del(data,bc)

% Build the Laplacian operator using a hybrid FE/FV method. The gradients
% are obtained at the cell centroids and then used to construct the face
% fluxes on the median cell edges.
%
% Only uses nearest neighbours and seems to be fully 2nd order accuracte.

p    = data.p;          % Nodes
e    = data.e;          % Edges
t    = data.t;          % Triangles
n2n  = data.n2n;        % Node to node connectivity
n2e  = data.n2e;        % Node to edge connectivity
e2t  = data.e2t;        % Edge to triangle connectivity
bnd  = data.bnd;        % True for boundary nodes
hnx  = data.hnx;        % Median cell
hny  = data.hny;        % normals
numn = size(p,1);       % Number of nodes


% FORM FEM GRADIENTS
% Evaluate centroidal gradients (piecewise-linear interpolants)
x23 = p(t(:,2),1)-p(t(:,3),1);  y23 = p(t(:,2),2)-p(t(:,3),2);
x21 = p(t(:,2),1)-p(t(:,1),1);  y21 = p(t(:,2),2)-p(t(:,1),2);

% Denominators
detx = x23.*y21-x21.*y23;
dety = y23.*x21-y21.*x23;

% xy gradient coefficients at vertices [n1,n2,n3]
sx = [y23./detx, (y21-y23)./detx, -y21./detx];
sy = [x23./dety, (x21-x23)./dety, -x21./dety];


% FORM LAPLACIAN
% Max number of entries in [i,j,s]
num = 6*sum(sum(n2n>0,2));

% Pre-alloc
i = repmat(uint32(0),num,1); j = i; s = double(i); 

next = 1;
for k = 1:numn
    if ~bnd(k)||(bc(k)==0)      % Node k not Dirichlet
        
        % Loop around neighbours of node k
        m = 1;
        while n2e(k,m)>0
            
            cn = n2n(k,m);      % Neighbour
            ce = n2e(k,m);      % Current edge
            t1 = e2t(ce,1);
            t2 = e2t(ce,2);
            
            % ce has 2 median edges associated, deal with separately
            if k==e(ce,1)
                hnx1 = hnx(ce,1); hny1 = hny(ce,1);
                hnx2 = hnx(ce,2); hny2 = hny(ce,2);
            else
                hnx1 = -hnx(ce,1); hny1 = -hny(ce,1);
                hnx2 = -hnx(ce,2); hny2 = -hny(ce,2);
            end

            % Loop over associated median edges and assemble flux
            for z = 1:3
                % Edge 1
                i(next) = k;
                j(next) = t(t1,z);
                s(next) = hnx1*sx(t1,z) + hny1*sy(t1,z);
                next    = next+1;
                % Edge 2
                if t2>0
                    i(next) = k;
                    j(next) = t(t2,z);
                    s(next) = hnx2*sx(t2,z) + hny2*sy(t2,z);
                    next    = next+1;  
                end
            end
            
            % Counter
            m = m+1;
            
        end
        
    else    % Node k is Dirichlet
        i(next) = k;
        j(next) = k;
        s(next) = 1;
        next    = next+1;
    end
end

% Sparse operator
L = sparse(double(i(1:next-1)),double(j(1:next-1)),s(1:next-1),numn,numn);

return


% Old FV Laplacian
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function L = del(data,bc,A,B,C)
% 
% % Build Laplacian operators via a 2nd order FV discretisation on the median
% % mesh. Incorporate BC's directly into the coefficient matrix.
% 
% p    = data.p;          % Nodes
% e    = data.e;          % Edges
% n2n  = data.n2n;        % Node to node connectivity
% n2e  = data.n2e;        % Node to edge connectivity
% e2t  = data.e2t;        % Edge to triangle connectivity
% bnd  = data.bnd;        % True for boundary nodes
% d    = data.d;          % Edge length
% hnx  = data.hnx;        % Median cell
% hny  = data.hny;        % normals
% numn = size(p,1);       % Number of nodes
% 
% % Max number of entries in [i,j,s]
% numN = sum(n2n>0,2);
% num  = sum(numN.*(2*numN+2));
% 
% % Pre-alloc
% i = repmat(uint32(0),num,1); j = i; s = double(i); 
% 
% next = 1;
% for k = 1:numn
%     if ~bnd(k)||(bc(k)==0)      % Node k not Dirichlet
%         
%         % Loop around neighbours of node k
%         m = 1;
%         while n2e(k,m)>0
%             
%             cn = n2n(k,m);      % Neighbour
%             ce = n2e(k,m);      % Current edge
%             t2 = e2t(ce,2);
%             
%             % xy nodes
%             xk = p(k,1);  yk = p(k,2);
%             xc = p(cn,1); yc = p(cn,2);
%             
%             % Normal vector for ce (different from median normals)
%             nxw = (xc-xk)/d(ce);
%             nyw = (yc-yk)/d(ce);
%             
%             % ce has 2 median edges associated, deal with separately
%             if k==e(ce,1)
%                 hnx1 = hnx(ce,1); hny1 = hny(ce,1);
%                 hnx2 = hnx(ce,2); hny2 = hny(ce,2);
%             else
%                 hnx1 = -hnx(ce,1); hny1 = -hny(ce,1);
%                 hnx2 = -hnx(ce,2); hny2 = -hny(ce,2);
%             end
%             % Coefficients
%             ax1 = hnx1*(1-nxw^2); bx1 = hnx1*nxw*nyw;
%             ax2 = hnx2*(1-nxw^2); bx2 = hnx2*nxw*nyw;
%             ay1 = hny1*(1-nyw^2); by1 = hny1*nxw*nyw;
%             ay2 = hny2*(1-nyw^2); by2 = hny2*nxw*nyw;
%             
%             % Centre coefficients
%             sum1 = -(hnx1*nxw+hny1*nyw)/d(ce);
%             if t2>0
%                 sum1 = sum1 - (hnx2*nxw+hny2*nyw)/d(ce);    
%             end
%             sum2 = -sum1;
%             
%             % Gradient around k
%             Ak = 0.5*A(k); Bk = 0.5*B(k); Ck = 0.5*C(k);
%             z = 1;
%             while n2n(k,z)>0 
%                 % Neighbour
%                 nn = n2n(k,z);
%                 % xy increments
%                 x = xk-p(nn,1); 
%                 y = yk-p(nn,2);
%                 w = 1/(x^2+y^2);
%                 % Gradient coefficients
%                 gy = Bk*w*x-Ak*w*y;
%                 gx = Ck*w*x-Bk*w*y;
%                 % Centre term
%                 temp = (ax1*-gx)-(bx1*gy)+(ay1*gy)-(by1*-gx);
%                 sum1 = sum1 - temp;
%                 if t2>0
%                     cont = (ax2*gx)-(bx2*-gy)+(ay2*-gy)-(by2*gx);
%                     sum1 = sum1 + cont; 
%                     temp = temp - cont;
%                 end
%                 % Neighbour term
%                 i(next) = k;
%                 j(next) = nn;
%                 s(next) = temp;
%                 next    = next+1;
%                 % Counter
%                 z = z+1;
%             end
%             
%             % Gradient around cn
%             Ac = 0.5*A(cn); Bc = 0.5*B(cn); Cc = 0.5*C(cn);
%             z = 1;
%             while n2n(cn,z)>0
%                 % Neighbour
%                 nn = n2n(cn,z);
%                 % xy increment
%                 x = xc-p(nn,1); 
%                 y = yc-p(nn,2);
%                 w = 1/(x^2+y^2);
%                 % Gradient coefficients
%                 gy = Bc*w*x-Ac*w*y;
%                 gx = Cc*w*x-Bc*w*y;
%                 % Centre term
%                 temp = (ax1*-gx)-(bx1*gy)+(ay1*gy)-(by1*-gx);
%                 sum2 = sum2 - temp;
%                 if t2>0
%                     cont = (ax2*gx)-(bx2*-gy)+(ay2*-gy)-(by2*gx);
%                     sum2 = sum2 + cont; 
%                     temp = temp - cont;
%                 end
%                 % Neighbour term
%                 i(next) = k;
%                 j(next) = nn;
%                 s(next) = temp;
%                 next    = next+1;
%                 % Counter
%                 z = z+1;
%             end
%             
%             % Centre term
%             i(next) = k;
%             j(next) = k;
%             s(next) = sum1;
%             next    = next+1;
%             % Neighbour term
%             i(next) = k;
%             j(next) = cn;
%             s(next) = sum2;
%             next    = next+1;
%             
%             % Counter
%             m = m+1;
%             
%         end
%         
%     else    % Node k is Dirichlet
%         i(next) = k;
%         j(next) = k;
%         s(next) = 1;
%         next    = next+1;
%     end
% end
% 
% % Sparse operator
% L = sparse(double(i(1:next-1)),double(j(1:next-1)),s(1:next-1),numn,numn);
% 
% return
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [A,B,C] = least_squares(data)
% 
% % Form least squares gradient coefficients.
% 
% p    = data.p;          % Nodes
% n2n  = data.n2n;        % Node to node connectivity
% numn = size(p,1);       % Number of nodes
% 
% % LEAST SQUARES GRADIENT COEFFICIENTS (inverse distance weighted)
% A = repmat(0,numn,1); B = A; C = A;
% for k = 1:numn  
%     % Loop around neighbours of node k
%     a = 0; b = 0; c = 0; m = 1;
%     while n2n(k,m)>0
%         % xy increments
%         x = p(k,1)-p(n2n(k,m),1); 
%         y = p(k,2)-p(n2n(k,m),2);
%         w = 1/(x^2+y^2);
%         % Coefficients
%         a = a+w*x^2; b = b+w*x*y; c = c+w*y^2;
%         % Counter
%         m = m+1;
%     end
%     % Save
%     A(k) = a; B(k) = b; C(k) = c;
% end
% 
% % Divide by det
% detM = A.*C-B.^2; 
% A    = A./detM; 
% B    = B./detM; 
% C    = C./detM;
% 
% return