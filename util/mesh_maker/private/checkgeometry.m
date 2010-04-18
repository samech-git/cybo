function [node,edge,hdata] = checkgeometry(node,edge,hdata)

% CHECKGEOMETRY: Check a geometry input for MESH2D.
%  
%  node  : Nx2 array of XY geometry nodes
%  edge  : Mx2 array of connections between nodes in NODE. (Optional)
%  hdata : Structure array defining size function data. 
%
% The following checks are performed:
%
%  1. Unique edges in EDGE.
%  2. Only nodes referenced in EDGE are kept.
%  3. Unique nodes in NODE.
%  4. No "hanging" nodes and no "T-junctions".
%
% Checks for self-intersecting geometry are NOT done because this can be 
% expensive for big inputs.
%
% HDATA may be re-indexed to maintain consistency.

% Darren Engwirda - 2007.

if nargin<2
   edge = [];
   if nargin<1
      error('Insufficient inputs');
   end
elseif nargin>3
   error('Wrong number of inputs');
end
if nargout>3
   error('Wrong number of outputs');
end
nNode = size(node,1);
if isempty(edge)                                                           % Build cnect if not passed
   edge = [(1:nNode-1)',(2:nNode)'; nNode,1];
end
if numel(node)~=2*nNode                                                    % Check inputs
   error('NODE must be an Nx2 array');
end
if numel(edge)~=2*size(edge,1)
   error('EDGE must be an Mx2 array');
end
if (max(edge(:))>nNode)||(min(edge(:))<=0)                                 % Check geometry indexing
   error('Invalid EDGE');
end
if isfield(hdata,'edgeh') && ~isempty(hdata.edgeh)                         % Check if we've got size data attached to edges
   edgeh = true;
else
   edgeh = false;
end

nEdge = size(edge,1);
[edge,i,j] = unique(sort(edge,2),'rows');                                  % Remove duplicate edges
if edgeh
   hdata.edgeh(:,1) = j(hdata.edgeh(:,1));
end
del = nEdge-size(edge,1);
if del>0
   disp(['WARNING: ',num2str(del),' duplicate edge(s) removed']);
   nEdge = size(edge,1);
end

i = unique(edge(:));                                                       % Remove un-used nodes and re-index
del = nNode-i(end);
if del>0
   node = node(i,:);
   j = zeros(size(node,1),1);
   j(i) = 1;
   j = cumsum(j);
   edge = j(edge);
   
   disp(['WARNING: ',num2str(del),' un-used node(s) removed']);
   nNode = i(end);
end

[i,i,j] = unique(node,'rows');                                             % Remove duplicate nodes and re-index
del = nNode-length(i);
if del>0
   node = node(i,:);
   edge = j(edge);
   
   disp(['WARNING: ',num2str(del),' duplicate node(s) removed']);
   nNode = size(node,1);
end

nEdge = size(edge,1);
S = sparse(edge(:), [1:nEdge,1:nEdge], 1, nNode, nEdge);                   % Sparse node-to-edge connectivity matrix
i = find(sum(S,2)<2);
if ~isempty(i)                                                             % Check for closed geometry loops
   error(['Open geometry contours detected at node(s): ',num2str(i')]);
end
i = find(sum(S,2)>2);                                                      % Check for geometry T-junctions
if ~isempty(i)
   error(['Multiple geometry branches detected at node(s): ',num2str(i')]);
end
