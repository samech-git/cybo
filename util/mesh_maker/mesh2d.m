function [p,t] = mesh2d(node,edge,hdata,options)

%  MESH2D: 2D unstructured triangular mesh generation.
%
% A 2D unstructured triangular mesh is generated based on a piecewise-
% linear geometry input. An iterative method is implemented to optimise 
% mesh quality. General multiply connected domains and element size 
% functions can be specified.
%
% Returns the final coordinates of the nodes p, and their triangulation t
% (with a counter-clockwise node ordering).
%
%  SHORT SYNTAX:
%
%  [p,t] = mesh2d(node);
%
% NODE defines the geometry nodes via an Nx2 array:
%
%  node  = [x1 y1; x2 y2; etc], geometry nodes specified in consectutive
%                               order, such that NODE(2,:) is joined with
%                               NODE(1,:) etc.
%
% An element size function is automatically generated based on the 
% complexity of the geometry. Generally this produces meshes with the 
% fewest number of triangles.
%
%  LONG SYNTAX:
%
%  [p,t] = mesh2d(node,edge,hdata,options);
%
% Blank arguments can be passed using the empty placeholder "[]".
%
% EDGE defines the connectivity between the points in NODE as a list of
% edges:
%
%   edge = [n1 n2; n2 n3; etc], connectivity between nodes to form
%                               geometry edges. If EDGE is specified it is
%                               not required that NODE be consectutive.
%
% HDATA is a structure for user defined element size information. HDATA can 
% include the following fields:
%
%  hdata.hmax  = h0;                   Max allowable global element size.
%  hdata.edgeh = [e1,h1; e2,h2; etc];  Element size on specified geometry 
%                                      edges.
%  hdata.fun   = 'fun' or @fun;        User defined size function.
%  hdata.args  = {arg1, arg2, etc};    Additional arguments for HDATA.FUN.
%
% Calls to user specified functions must accept vectorised input of the 
% form H = FUN(X,Y,ARGS{:}), where X,Y are the xy coordinates where the
% element size will be evaluated and ARGS are optional additional arguments 
% as passed by HDATA.ARGS.
%
% An automatic size function is always generated to ensure that the
% geometry is adequately resolved. The overall size function is the minimum
% of the user specified and automatic functions.
%
% OPTIONS is a structure array that allows some of the "tuning" parameters
% used in the solver to be modified:
%
%   options.mlim   - Specifies the maximum allowable change in edge length
%                    per iteration that must be obtained [default = 0.05
%                    (5%)]. Larger values should accelerate convergence,
%                    but could decrease mesh quality.
%   options.maxit  - Specifies the maximum number of iterations the
%                    algorithm will perform [default = 20].
%   options.dhmax  - Specifies the maximum allowable (relative) gradient in
%                    the size function [default = 0.3].
%   options.output - Displays the mesh and the mesh statistics upon
%                    completion [default = true].
%
% EXAMPLE:
%
%   meshdemo                  % Will run the standard demos
%   mesh_collection(n)        % Will run some additional demos
%
% See also REFINE, SMOOTHMESH, DELAUNAYN


% Mesh2d is a delaunay based algorithm with a "Laplacian-like" smoothing
% operation built into the mesh generation process. 
% 
% An unbalanced quadtree decomposition is used to evaluate the element size 
% distribution required to resolve the geometry. The quadtree is 
% triangulated and used as a backgorund mesh to store the element size 
% data.  
%
% The main method attempts to optimise the node location and mesh topology 
% through an iterative process. In each step a constrained delaunay 
% triangulation is generated with a series of "Laplacian-like" smoothing 
% operations used to improve triangle quality. Nodes are added or removed 
% from the mesh to ensure the required element size distribution is 
% approximated.  
%
% The optimisation process generally returns well shaped meshes with no
% small angles and smooth element size variations. Mesh2d shares some 
% similarities with the Distmesh code: 
%
%   [1] P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB.
%       SIAM Review, Volume 46 (2), pp. 329-345, June 2004
%
%   Darren Engwirda : 2005-07
%   Email           : d_engwirda@hotmail.com
%   Last updated    : 08/07/2007 with MATLAB 7.0
%
% Mesh2d is Copyright (C) 2007 Darren Engwirda. See "copyright.m" for full 
% details.
%
% Please email me any un-meshable geometries, meshing benchmarks or
% suggestions!

tic
wbar = waitbar(0,'Forming geometry faces');

% Error checks
if nargin<4
   options = [];
   if nargin<3
      hdata = [];
      if nargin<2
         edge = [];
         if nargin<1
            error('Wrong number of inputs');
         end
      end
   end
elseif nargin>4
   error('Wrong number of inputs');
end
if nargout>2
   error('Wrong number of outputs');
end

% Get user options
[mlim,maxit,dhmax,output] = getoptions(options);

% Check geometry
[node,edge,hdata] = checkgeometry(node,edge,hdata);                        % Error checking for geometry
edgexy = [node(edge(:,1),:), node(edge(:,2),:)];                           % Geometry as edge endpoints [x1,y1,x2,y2]

% QUADTREE DECOMPOSITION
%  PH    : Background mesh nodes
%  TH    : Background mesh triangles
%  HH    : Size function value at PH
[ph,th,hh] = quadtree(edgexy,hdata,dhmax,wbar);

waitbar(0,wbar,'Initialising mesh');

% INITIALISE MESH
%  P     : Initial nodes
%  T     : Initial triangulation
%  WNDX  : Closest edge for each node as indices into EDGEXY
%  TNDX  : Enclosing triangle for each node as indices into TH
%  FIX   : Indices of FIXED nodes in P
[p,t,wndx,tndx,fix,fixed] = initmesh(ph,th,node,edge,edgexy);

% MAIN LOOP
retri = false;
subit = 3;
for iter = 1:maxit

   % Re-triangulation
   if retri
      t = MyDelaunayn(p);                                                  % Delaunay triangulation via QHULL
      pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;                        % Centroids
      t = t(inpoly(pc,node,edge),:);                                       % Take triangles with internal centroids
   end
   
   [e,bnd] = getedges(t,size(p,1));                                        % Unique edges and boundary nodes
   bnd(fix) = false;                                                       % Don't flag fixed nodes
   bnd = find(bnd);

   nume = size(e,1);
   S = sparse(e(:),[1:nume,1:nume],1,size(p,1),nume);                      % Sparse node-to-edge connectivity matrix

   tndx = mytsearch(ph(:,1),ph(:,2),th,p(:,1),p(:,2),tndx);                % Find enclosing triangle in background mesh for nodes
   h = tinterp(ph,th,hh,p,tndx);                                           % Size function at nodes via linear interpolation
   h = 0.5*(h(e(:,1))+h(e(:,2)));                                          % from the background mesh. Average to edge midpoints.
   
   % Inner smoothing iterations
   L = max(sqrt(sum((p(e(:,1),:)-p(e(:,2),:)).^2,2)),eps);                 % Edge length
   done = false;
   subit = max(subit,iter);                                                % Increment sub-iters with outer iters to aid convergence
   for subiter = 1:subit   
      r = sqrt(L./h);                                                      % Ratio of actual to desired edge length

      pm = 0.5*[r,r].*( p(e(:,1),:)+p(e(:,2),:) );                         % Edge midpoints, size function weighted

      W = max(S*r,eps);                                                    % Size function weighting
      p = (S*pm)./[W,W];                                                   % Weighted Laplacian-like smoothing
      p(fix,:) = fixed;                                                    % Don't update fixed nodes

      [p,wndx] = project2poly(p,bnd,edgexy,wndx);                          % Project bnd nodes onto the closest geometry edge
      
      Lnew = max(sqrt(sum((p(e(:,1),:)-p(e(:,2),:)).^2,2)),eps);           % Edge length
      move = norm((Lnew-L)./L,inf);                                        % Percentage change
      L = Lnew;

      if move<mlim                                                         % Test convergence
         done = true;
         break
      end
   end
   
   msg = ['Generating mesh (Iteration ',num2str(iter),')'];                % Show progress
   waitbar(mlim/max(move,eps),wbar,msg);

   r = L./h;
   if done && max(r)<3                                                     % Main loop convergence
      break
   end

   % Nodal density control
   retri = false;
   if iter<maxit
      test = find(r<=0.5);                                                 % Remove both nodes for edges with r<0.5
      hang = sum(S,2)<2;                                                   % Remove nodes connected to less than 2 edges
      if ~isempty(test) || any(hang)
         prob = false(size(p,1),1);                                        % True for nodes to be removed
         prob(e(test,:)) = true;                                           % Edges with r<0.5
         prob(hang) = true;
         prob(fix) = false;                                                % Make fixed nodes work
         pnew = p(~prob,:);                                                % New node list

         tmp_wndx = wndx(~prob);
         tmp_tndx = tndx(~prob);
         
         j = zeros(size(p,1),1);                                           % Re-index fix to make fixed nodes work
         j(~prob) = 1;
         j = cumsum(j);
         fix = j(fix);
         
         retri = true;
      else
         pnew = p;                                                         % No nodes removed
         tmp_wndx = wndx;
         tmp_tndx = tndx;
      end
      test = find(r>=2);                                                   % Add node at midpoint for edges with r>=2
      if ~isempty(test)
         p = [pnew; 0.5*(p(e(test,1),:)+p(e(test,2),:))];
         wndx = [tmp_wndx; zeros(length(test),1)];
         tndx = [tmp_tndx; zeros(length(test),1)];
         retri = true;
      else
         p = pnew;                                                         % No nodes added
         wndx = tmp_wndx;
         tndx = tmp_tndx;
      end
   end

end
close(wbar)

if iter==maxit
   disp('WARNING: Maximum number of iterations reached. Solution did not converge!')
   disp('Please email the geometry and settings to d_engwirda@hotmail.com')
end

[p,t] = fixmesh(p,t);                                                      % Ensure triangulation unique and
                                                                           % CCW node ordered
tfinal = toc;
if output

   figure('Name','Mesh')
   patch('faces',t,'vertices',p,'facecolor','None','edgecolor','b')
   hold on
   patch('faces',edge,'vertices',node,'facecolor','None','edgecolor','r')
   axis equal off; hold off

   % Mesh measures
   q = quality(p,t);
   
   disp(struct('Iterations',       iter, ...
               'Time',             tfinal, ...
               'Triangles',        size(t,1), ...
               'Nodes',            size(p,1), ...
               'Mean_quality',     mean(q), ...
               'Min_quality',      min(q)) );
               %'Mean_size_ratio',  mean(r), ...
               %'Min_size_ratio',   min(r), ...
               %'Max_size_ratio',   max(r)) );
end


%% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,t,wndx,tndx,fix,fixed] = initmesh(p,t,node,edge,edgexy)

% Initialise the nodes, triangulation and data structures for the mesh 
% based on the quadtree data and geometry.

in = inpoly(p,node,edge);                                                  % True for internal nodes
tin = any(in(t),2);                                                        % True for triangles with at least one internal node
ok = false(size(p,1),1);                                                   % Nodes considered approximately internal
ok(t(tin,:)) = true;

fixed = node;
ndx = tsearch(p(:,1),p(:,2),t,fixed(:,1),fixed(:,2));                      % Find enclosing triangle for fixed nodes

% At this stage some nodes have been accepted that are not inside the 
% geometry. This is done because the quadtree triangulation will overlap 
% the geometry in some cases, so nodes invloved in the overlap are accepted 
% to get a reasonable distribution near the edges.
[p,wndx] = project2poly(p,find(~in&ok),edgexy);     

% Find the closest node in P for each FIXED node and replace P(i,:) with
% FIXED(i,:)
fix = zeros(size(fixed,1),1);
for k = 1:size(ndx,1)
   
   x = fixed(k,1);
   y = fixed(k,2);

   if isnan(ndx(k))
      % Slow search for all p(ok,:)
      [tmp,tmp] = min( (p(ok,1)-x).^2+(p(ok,2)-y).^2 );
      fix(k) = tmp;
   else
      % Search nodes in ndx(k)
      d = inf;
      j = 1;
      while j<=3
         cn = t(ndx(k),j);
         if ok(cn)
            dkj = (p(cn,1)-x)^2+(p(cn,2)-y)^2;
            if dkj<d
               fix(k) = cn;
               d = dkj;
            end
         end
         j = j+1;
      end
   end
   
   if fix(k)==0      
      % Slow search for all p(ok,:)
      [tmp,tmp] = min( (p(ok,1)-x).^2+(p(ok,2)-y).^2 );
      fix(k) = tmp;
   end

end
p(fix,:) = fixed;

% Take internal nodes
p = p(ok,:);

% Re-index to keeps lists consistent
wndx = wndx(ok);
j = zeros(length(ok),1);
j(ok) = 1;
j = cumsum(j);
t = j(t(tin,:));
fix = j(fix);
tndx = zeros(size(p,1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e,bnd] = getedges(t,n)

% Get the unique edges and boundary nodes in a triangulation.

e = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];                                  % Non-unique edges

swap = e(:,2)<e(:,1);                                                      % Ensure e(:,1) contains the lower value
e(swap,:) = e(swap,[2,1]);

e = sortrows(e);
idx = all(diff(e,1)==0,2);                                                 % Find shared edges
idx = [idx;false]|[false;idx];                                             % True for all shared edges
bnde = e(~idx,:);                                                          % Boundary edges
e = e(idx,:);                                                              % Internal edges
e = [bnde; e(1:2:end-1,:)];                                                % Unique edges and bnd edges for tri's

bnd = false(n,1);                                                          % True for boundary nodes
bnd(bnde) = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,ndxnew] = project2poly(p,bnd,edgexy,ndx)

% Project the points in P(BND) onto the closest edge of the polygon defined
% by the edge segments in EDGEXY. NDX is an optional argument defining
% the edge to project onto: P(BND) is projected onto EDGEXY(NDX(BND)).

% Uses (something like?) a double sweep-line approach to reduce the number
% of edges that are required to be tested in order to determine the closest
% edge for each point. On average only size(EDGEXY)/4 comparisons need to
% be made for each point.

if nargin<4 || isempty(ndx)
   ndx = zeros(size(p,1),1);
end
ndxnew = zeros(size(p,1),1);
todo = true(size(bnd));

% Check NDX first
for k = 1:length(bnd)
   cn = bnd(k);
   if ndx(cn)>0
      j = ndx(cn);

      x1 = edgexy(j,1); x2mx1 = edgexy(j,3)-x1;
      y1 = edgexy(j,2); y2my1 = edgexy(j,4)-y1;

      r = ((p(cn,1)-x1)*x2mx1+(p(cn,2)-y1)*y2my1)/(x2mx1^2+y2my1^2);
      if (r>0) && (r<1)
         todo(k) = false;
         p(cn,1) = x1+r*x2mx1;
         p(cn,2) = y1+r*y2my1;
         ndxnew(cn) = j;
      end

   end
end

% Do a full search for points not already projected
if any(todo)

   bnd = bnd(todo);

   % Choose the direction with the biggest range as the "y-coordinate" for the
   % test. This should ensure that the sorting is done along the best
   % direction for long and skinny problems wrt either the x or y axes.
   dxy = max(p)-min(p);
   if dxy(1)>dxy(2)
      % Flip co-ords if x range is bigger
      p = p(:,[2,1]);
      edgexy  = edgexy(:,[2,1,4,3]);
      flip = true;
   else
      flip = false;
   end

   % Ensure edgexy(:,[1,2]) contains the lower y value
   swap = edgexy(:,4)<edgexy(:,2);
   edgexy(swap,:) = edgexy(swap,[3,4,1,2]);

   % Sort edges
   [ilower,ilower] = sort(edgexy(:,2));                                    % Sort edges by lower y value
   edgexy_lower = edgexy(ilower,:);
   [iupper,iupper] = sort(edgexy(:,4));                                    % Sort edges by upper y value
   edgexy_upper = edgexy(iupper,:);

   % Mean edge y value
   ne = size(edgexy,1);
   ymean = 0.5*( sum(sum(edgexy(:,[2,4]))) )/ne;

   % Loop through points
   for k = 1:length(bnd)

      cn = bnd(k);
      x = p(cn,1);
      y = p(cn,2);
      d = inf;

      if y<ymean

         % Loop through edges bottom up
         for j = 1:ne
            y2 = edgexy_lower(j,4);
            if y2>=(y-d)
               y1 = edgexy_lower(j,2);
               if y1<=(y+d)

                  % Calculate the distance along the normal projection from [x,y] to the jth edge
                  x1 = edgexy_lower(j,1); 
                  x2mx1 = edgexy_lower(j,3)-x1;
                  y2my1 = y2-y1;

                  r = ((x-x1)*x2mx1+(y-y1)*y2my1)/(x2mx1^2+y2my1^2);
                  if r>1                                                   % Limit to wall endpoints
                     r = 1;
                  elseif r<0
                     r = 0;
                  end
                  xn = x1+r*x2mx1;
                  yn = y1+r*y2my1;

                  dj = (xn-x)^2+(yn-y)^2;
                  if ( dj<d^2 )
                     d = sqrt(dj);
                     p(cn,1) = xn;
                     p(cn,2) = yn;
                     ndxnew(cn) = ilower(j);
                  end

               else
                  break
               end
            end
         end

      else

         % Loop through edges top down
         for j = ne:-1:1
            y1 = edgexy_upper(j,2);
            if y1<=(y+d)
               y2 = edgexy_upper(j,4);
               if y2>=(y-d)

                  % Calculate the distance along the normal projection from [x,y] to the jth edge
                  x1 = edgexy_upper(j,1); 
                  x2mx1 = edgexy_upper(j,3)-x1;
                  y2my1 = y2-y1;

                  r = ((x-x1)*x2mx1+(y-y1)*y2my1)/(x2mx1^2+y2my1^2);
                  if r>1                                                   % Limit to wall endpoints
                     r = 1;
                  elseif r<0
                     r = 0;
                  end
                  xn = x1+r*x2mx1;
                  yn = y1+r*y2my1;

                  dj = (xn-x)^2+(yn-y)^2;
                  if ( dj<d^2 )
                     d = sqrt(dj);
                     p(cn,1) = xn;
                     p(cn,2) = yn;
                     ndxnew(cn) = iupper(j);
                  end

               else
                  break
               end
            end
         end

      end

   end

   if flip
      p = p(:,[2,1]);
   end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mlim,maxit,dhmax,output] = getoptions(options)

% Extract the user defined options

if ~isempty(options)
   if ~isstruct(options)
      error('OPTIONS must be a structure array');
   end
   if numel(options)~=1
      error('Options cannot be an array of structures');
   end
   fields = fieldnames(options);
   names = {'mlim','maxit','dhmax','output'};
   for k = 1:length(fields)
      if strcmp(fields{k},names)
         error('Invalid field in OPTIONS');
      end
   end
   if isfield(options,'mlim')                                              % Movement tolerance
      mlim = checkposscalar(options.mlim,'options.mlim');
   else
      mlim = 0.05;
   end
   if isfield(options,'maxit')                                             % Maximum iterations
      maxit = round(checkposscalar(options.maxit,'options.maxit'));
   else
      maxit = 20;
   end
   if isfield(options,'dhmax')                                             % Size function gradient limit
      dhmax = checkposscalar(options.dhmax,'options.dhmax');
   else
      dhmax = 0.3;
   end
   if isfield(options,'output')                                            % Output on/off
      output = checklogicalscalar(options.output,'options.output');
   else
      output = true;
   end
else                                                                       % Default values
   mlim = 0.05;
   maxit = 20;
   dhmax = 0.3;
   output = true;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = checkposscalar(var,name)

% Helper function to check if var is a positive scalar.

if var<0 || any(size(var)>1)
   error([name,' must be a positive scalar']);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function var = checklogicalscalar(var,name)

% Helper function to check if var is a logical scalar.

if ~islogical(var) || any(size(var)>1)
   error([name,' must be a logical scalar']);
end
