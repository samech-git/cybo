function [p,t,f] = fixmesh(p,t,f)

%  FIXMESH: Ensure that triangular mesh data is consistent.
%
%  [p,t,f] = fixmesh(p,t,f);
%
%  p: Nx2 array of nodal XY coordinates, [x1,y1; x2,y2; etc]
%  t: Mx3 array of triangles as indices, [n11,n12,n13; n21,n22,n23; etc]
%  f: (Optional) NxK array of nodal function values. Each column in F 
%     corresponds to a dependent function, F(:,1) = F1(P), F(:,2) = F2(P)
%     etc.
%
% The following checks are performed:
%
%  1. Nodes not refereneced in T are removed.
%  2. Duplicate nodes are removed.
%  3. Triangles are ordered counter-clockwise.
%  4. Triangles with an area less than 100*eps*norm(A,'inf') are removed

% Darren Engwirda - 2007.

if nargin<=3
    if nargin<=2
        gotF = false;
        if nargin<2
            error('Wrong number of inputs');
        end
    else
        gotF = true;
    end
else
    error('Wrong number of inputs');
end
if (gotF&&(nargout>3)) || (~gotF&&(nargout>2))
    error('Wrong number of outputs');
end
if numel(p)~=2*size(p,1)
    error('P must be an Nx2 array');
end
if numel(t)~=3*size(t,1)
    error('T must be an Mx3 array');
end
if (any(t(:))<1) || (max(t(:))>size(p,1))
    error('Invalid T');
end
if gotF && ((size(f,1)~=size(p,1)) || (ndims(f)>2))
    error('F must be an NxK array');
end

% Remove un-used nodes
[i,j,j] = unique(t(:));
if gotF
    f = f(i,:);
end
p = p(i,:);
t = reshape(j,size(t));

% Remove duplicate nodes
[i,i,j] = unique(p,'rows');
if gotF
    f = f(i,:);
end
p = p(i,:);
t = reshape(j(t),size(t));

% Triangle area
d12 = p(t(:,2),:)-p(t(:,1),:);
d13 = p(t(:,3),:)-p(t(:,1),:);
A   = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1));
Ai  = A<0;
Aj  = abs(A)>100*eps*norm(A,'inf');

% Flip node numbering to give a counter-clockwise order
t(Ai,[1,2]) = t(Ai,[2,1]);

% Remove zero area triangles
t = t(Aj,:);
