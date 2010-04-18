function data = data_structure(p,t)

% Build data structures used in the solver.
%
% This function builds the mesh parameters and connectivity that will be
% used in tvd_rk2.m. Building this initially saves HEAPS of CPU time later 
% in tvd_rk2.m, but of course, my meshes are fixed for the integration.
%
% Darren Engwirda - 2005-2006.
%
% Naver2d is Copyright (C) 2005-2006 Darren Engwirda. See "copyright.m" for
% full details.

w = waitbar(0,'Building data structures');

numn = size(p,1);
numt = size(t,1);
vect = 1:numt;

% DETERMINE UNIQUE EDGES IN MESH
 
e       = [t(:,[1,2]); t(:,[2,3]); t(:,[3,1])];             % Edges - not unique
vec     = (1:size(e,1))';                                   % List of edge numbers
[e,i,j] = unique(sort(e,2),'rows');                         % Unique edges
vec     = vec(j);                                           % Unique edge numbers
eINt    = [vec(vect), vec(vect+numt), vec(vect+2*numt)];    % Unique edges in each triangle

waitbar(0.2,w);

% DETERMINE EDGE TO TRIANGLE CONNECTIVITY

% Each row has two entries corresponding to the triangle numbers
% associated with each edge. Boundary edges have one entry = 0.
nume = size(e,1);
e2t  = repmat(0,nume,2);
ndx  = repmat(1,nume,1);
for k = 1:numt
    % Edge in kth triangle
    e1 = eINt(k,1); e2 = eINt(k,2); e3 = eINt(k,3);
    % Edge 1
    e2t(e1,ndx(e1)) = k; ndx(e1) = ndx(e1)+1;
    % Edge 2
    e2t(e2,ndx(e2)) = k; ndx(e2) = ndx(e2)+1;
    % Edge 3
    e2t(e3,ndx(e3)) = k; ndx(e3) = ndx(e3)+1;
end

waitbar(0.4,w);

% DETERMINE NODE TO EDGE CONNECTIVITY

% Determine maximum neighbours
j = e(:);
v = repmat(0,max(j),1);
for k = 1:length(j)
    jk = j(k); v(jk) = v(jk)+1;
end
maxN = max(v);

n2e = repmat(0,numn,maxN+1);
ndx = repmat(1,numn,1);
for k = 1:nume
    % End nodes
    n1 = e(k,1); n2 = e(k,2);
    % Connectivity
    n2e(n1,ndx(n1)) = k; ndx(n1) = ndx(n1)+1;
    n2e(n2,ndx(n2)) = k; ndx(n2) = ndx(n2)+1;
end

waitbar(0.6,w);

% DETERMINE NODE TO NODE CONNECTIVITY

n2n = repmat(0,numn,maxN+1);
for k = 1:numn
    next = 1; m = 1;
    while n2e(k,m)>0
        if e(n2e(k,m),1)==k
            n2n(k,next) = e(n2e(k,m),2); next = next+1;
        else
            n2n(k,next) = e(n2e(k,m),1); next = next+1; 
        end
        m = m+1;
    end
end

waitbar(0.8,w);

% FLAG BOUNDARY ELEMENTS

vec     = (1:nume)';            % Edge list
be      = vec(~all(e2t,2));     % Boundary edges
bn      = e(be,:);           
bn      = unique(bn(:));        % Boundary nodes
bnd     = false(numn,1);        
bnd(bn) = true;                 % True for boundary nodes

% IMPORTANT MESH DATA

d  = sqrt(sum((p(e(:,2),:)-p(e(:,1),:)).^2,2));     % Edge lengths
pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;       % Centroids
pe = 0.5*(p(e(:,1),:)+p(e(:,2),:));                 % Edge midpoints

% MEDIAN CELL PARAMETERS

hnx = repmat(0,nume,2); 
hny = hnx; 
A   = repmat(0,numn,1);
for k = 1:nume
    % Nodes and triangles
    n1 = e(k,1);   n2 = e(k,2);
    t1 = e2t(k,1); t2 = e2t(k,2);
    % xy nodes
    x1 = p(n1,1); y1 = p(n1,2);
    x2 = p(n2,1); y2 = p(n2,2);
    % Edge midpoint
    mx = 0.5*(x1+x2);
    my = 0.5*(y1+y2);

    % MEDIAN EDGE 1
    % Cell centroid
    cx = pc(t1,1);
    cy = pc(t1,2);
    % Area contribution
    ac    = 0.5*((mx-x1)*(cy-y1)-(my-y1)*(cx-x1));
    A(n1) = A(n1)+abs(ac);
    A(n2) = A(n2)+abs(0.5*((mx-x2)*(cy-y2)-(my-y2)*(cx-x2)));
    % Edge normal
    dx        = cx-mx; 
    dy        = cy-my;
    hnx(k,1) =  sign(ac)*dy;
    hny(k,1) = -sign(ac)*dx;
    
    % MEDIAN EDGE 2
    if t2>0
        % Cell centroid
        cx = pc(t2,1);
        cy = pc(t2,2);
        % Area contribution
        ac    = 0.5*((mx-x1)*(cy-y1)-(my-y1)*(cx-x1));
        A(n1) = A(n1)+abs(ac);
        A(n2) = A(n2)+abs(0.5*((mx-x2)*(cy-y2)-(my-y2)*(cx-x2)));
        % Edge normal
        dx        = cx-mx; 
        dy        = cy-my;
        hnx(k,2) =  sign(ac)*dy;
        hny(k,2) = -sign(ac)*dx;
    else        
        for q = 1:3
            if eINt(t1,q)==k
                if q==3, z = 1; else z = q+1; end
                % Boundary edge normals
                hnx(k,2) =  0.5*(p(t(t1,z),2)-p(t(t1,q),2));
                hny(k,2) = -0.5*(p(t(t1,z),1)-p(t(t1,q),1));
            end
        end        
    end
    
end

waitbar(1,w);

% CELL DATA STRUCTURE

eINt = uint32(eINt);
e2t  = uint32(e2t);
n2e  = uint32(n2e);
n2n  = uint32(n2n);

data = struct('p'     ,p      ,...      % Nodes
              'e'     ,e      ,...      % Edges
              't'     ,t      ,...      % Triangles
              'pc'    ,pc     ,...      % Centroids
              'pe'    ,pe     ,...      % Edge midpoints    
              'd'     ,d      ,...      % Edge length
              'hnx'   ,hnx    ,...      % Median cell face normals
              'hny'   ,hny    ,...      % (pre-multiplied by face length)
              'A'     ,A      ,...      % Median cell area
              'bnd'   ,bnd    ,...      % True for boundary nodes
              'be'    ,be     ,...      % Boundary edges
              'eINt'  ,eINt   ,...      % Triangles as edges
              'e2t'   ,e2t    ,...      % Edge to triangle connectivity
              'n2e'   ,n2e    ,...      % Node to edge connectivity
              'n2n'   ,n2n );           % Node to node connectivity

close(w)          
          
return