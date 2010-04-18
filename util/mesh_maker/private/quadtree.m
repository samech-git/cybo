function [p,t,h] = quadtree(edgexy,hdata,dhmax,wbar)

%  QUADTREE: 2D quadtree decomposition.
%
% A quadtree decomposition of the linear geometry input is performed. The 
% bounding box is recursively subdivided until the dimension of each box 
% matches the local geometry feature size. The geometry feature size is 
% based on the minimum distance between linear geometry segments.
%
% A size function is obtained at the quadtree vertices based on the minimum
% neighbouring box dimension at each vertex. This size function is gradient
% limited to produce a smooth function.
%
% The quadtree is triangulated for use as a background mesh to store the
% size function.
%
%   edgexy  : [x1,y1,x2,y2] edge endpoints
%   hdata   : User defined size function structure
%   dhmax   : Maximum allowalble relative gradient in the size function
%
%   p       : Background mesh nodes
%   t       : Background mesh triangles
%   h       : Size function value at p

%   Darren Engwirda : 2007
%   Email           : d_engwirda@hotmail.com
%   Last updated    : 08/07/2007 with MATLAB 7.0

%% LOCAL FEATURE SIZE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(0,wbar,'Estimating local geometric feature size');

% Get size function data
[hmax,edgeh,fun,args] = gethdata(hdata);

% Insert test points along the boundaries at which the LFS can be
% approximated.
wm = 0.5*(edgexy(:,[1,2])+edgexy(:,[3,4]));                                % Use the edge midpoints as a first pass
len = sqrt(sum((edgexy(:,[3,4])-edgexy(:,[1,2])).^2,2));                   % Edge length
L = 2*dist2poly(wm,edgexy,len);                                            % Estimate the LFS at these points by calculating
                                                                           % the distance to the closest edge segment                           
% In cases where edges are separated by less than their length
% we will need to add more points to capture the LFS in these
% regions. This allows us to pick up "long and skinny" geometry
% features
r = 2*len./L;                                                              % Compare L (LFS approximation at wm) to the edge lengths
r = round((r-2)/2);                                                        % Amount of points that need to be added
add = find(r);                                                             % at each edge
if ~isempty(add)
   num = 2*sum(r(add));                                                    % Total number of points to added
   start = size(wm,1)+1;
   wm = [wm; zeros(num,2)];                                                % Alloc space
   next = start;
   for j = 1:length(add)                                                   % Loop through edges to be subdivided
      cw = add(j);                                                         % Current edge
      num = r(cw);
      tmp = (1:num)'/(num+1);                                              % Subdivision increments
      num = next+2*num-1;

      x1 = edgexy(cw,1); x2 = edgexy(cw,3); xm = wm(cw,1);                 % Edge values
      y1 = edgexy(cw,2); y2 = edgexy(cw,4); ym = wm(cw,2);

      wm(next:num,:) = [x1+tmp*(xm-x1), y1+tmp*(ym-y1)                     % Add to list
                        xm+tmp*(x2-xm), ym+tmp*(y2-ym)];
      next = num+1;
   end
   L = [L; dist2poly(wm(start:next-1,:),edgexy)];                          % Estimate LFS at the new points
end

% Compute the required size along the edges for any boundary layer size
% functions and add additional points where necessary.
if ~isempty(edgeh)
   h0 = edgeh(:,2);
   i = edgeh(:,1);
   for j = 1:length(i)
      if L(i(j))>h0
         cw = i(j);
         r = 2*len(cw)/h0(j);
         r = round((r-2)/2);                                               % Number of points to be added
         tmp = (1:r)'/(r+1);

         x1 = edgexy(cw,1); x2 = edgexy(cw,3); xm = wm(cw,1);              % Edge values
         y1 = edgexy(cw,2); y2 = edgexy(cw,4); ym = wm(cw,2);

         wm = [wm; x1+tmp*(xm-x1), y1+tmp*(ym-y1);                         % Add to list
                   xm+tmp*(x2-xm), ym+tmp*(y2-ym)];

         L(cw) = h0(j);                                                    % Update LFS
         L = [L; h0(j)*ones(2*r,1)];
      end
   end
end

% To speed the point location in the quadtree decomposition
% sort the LFS points based on y-value
[i,i] = sort(wm(:,2));
wm = wm(i,:);
L = L(i);
nw = size(wm,1);

%% UNBALANCED QUADTREE DECOMPOSITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(0.4,wbar,'Quadtree decomposition');

xymin = min([edgexy(:,[1,2]); edgexy(:,[3,4])]);                           % Bounding box
xymax = max([edgexy(:,[1,2]); edgexy(:,[3,4])]);

dim = max(xymax-xymin);                                                    % Bbox dimensions
xm = 0.5*(xymin(1)+xymax(1));
ym = 0.5*(xymin(2)+xymax(2));

p = [xm-0.5*dim, ym-0.5*dim                                                % Initial nodes
     xm+0.5*dim, ym-0.5*dim
     xm+0.5*dim, ym+0.5*dim
     xm-0.5*dim, ym+0.5*dim];
b = [1,2,3,4];                                                             % Initial boxes
h = userhfun(p(:,1),p(:,2),fun,args,hmax,xymin,xymax);

pblock = 5*nw;                                                             % Alloc memory in blocks
bblock = pblock;

np = 4;
nb = 1;
test = true;
while true                                                           
   
   vec = find(test(1:nb));                                                 % Boxes to be checked at this step
   if isempty(vec)
      break
   end

   N = np;
   for k = 1:length(vec)                                                   % Loop through boxes to be checked for subdivision
      
      m  = vec(k);                                                         % Current box
      n1 = b(m,1);  n2 = b(m,2);
      n3 = b(m,3);  n4 = b(m,4);
      x1 = p(n1,1); y1 = p(n1,2);
      x2 = p(n2,1); y4 = p(n4,2);

      % Binary search to find first wm with y>=ymin for current box
      if wm(1,2)>=y1
         start = 1;
      elseif wm(nw,2)<=y1
         start = nw;
      else
         lower = 1;
         upper = nw;
         for i = 1:nw
            start = round(0.5*(lower+upper));
            if wm(start,2)<y1
               lower = start;
            elseif wm(start-1,2)<y1
               break;
            else
               upper = start;
            end
         end
      end

      LFS = 4*min([h(n1),h(n2),h(n3),h(n4),2*hmax/4]);
      for i = start:nw                                                     % Loop through points (acending y-value order)
         if wm(i,2)<=y4                                                    % Check box bounds and current min
            if wm(i,2)>=y1 && wm(i,1)>=x1 && wm(i,1)<=x2 && L(i)<LFS
               LFS = L(i);                                                 % New min found - reset
            end
         else                                                              % Due to the sorting
            break;
         end
      end

      if (x2-x1)>LFS                                                       % Split current box
         if (np+5)>=size(p,1)                                              % Alloc memory on demand
            p = [p; zeros(pblock,2)];
            pblock = 2*pblock;
         end
         if (nb+3)>=size(b,1)
            b = [b; zeros(bblock,4)];
            test = [test; true(bblock,1)];
            bblock = 2*bblock;
         end

         xm = x1+0.5*(x2-x1);                                              % Current midpoints
         ym = y1+0.5*(y4-y1);

         p(np+1,:) = [xm,ym];                                              % New nodes
         p(np+2,:) = [xm,y1];
         p(np+3,:) = [x2,ym];
         p(np+4,:) = [xm,y4];
         p(np+5,:) = [x1,ym];

         b(m,:)    = [n1,np+2,np+1,np+5];                                  % New boxes
         b(nb+1,:) = [np+2,n2,np+3,np+1];
         b(nb+2,:) = [np+1,np+3,n3,np+4];
         b(nb+3,:) = [np+5,np+1,np+4,n4];

         nb = nb+3;
         np = np+5;
      else
         test(m) = false;
      end
   end

   h = [h; userhfun(p(N+1:np,1),p(N+1:np,2),fun,args,hmax,xymin,xymax)];
end
p = p(1:np,:);
b = b(1:nb,:);

[p,i,j] = unique(p,'rows');                                                % Delete replicated nodes
h = h(i);
b = reshape(j(b),size(b));

%% FORM SIZE FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
waitbar(0.6,wbar,'Forming element size function')

E = [b(:,[1,2]); b(:,[2,3]); b(:,[3,4]); b(:,[4,1])];                      % Edges
L = sqrt(sum((p(E(:,1),:)-p(E(:,2),:)).^2,2));                             % Edge length

nE = size(E,1);
for k = 1:nE                                                               % Initial h - minimum neighbouring edge length
   if L(k)<h(E(k,1))
      h(E(k,1)) = L(k); 
   end             
   if L(k)<h(E(k,2))
      h(E(k,2)) = L(k); 
   end
end

% Gradient limiting
tol = 1e-4;
while true                                                                 % Loop over the edges of the background mesh ensuring
   h_old = h;                                                              % that dh satisfies the dhmax tolerance
   for k = 1:nE                                                            % Loop over edges
      n1 = E(k,1);
      n2 = E(k,2);
      if h(n1)>h(n2)                                                       % Ensure grad(h)<=dhmax
         dh = (h(n1)-h(n2))/L(k);
         if dh>dhmax
            h(n1) = h(n2) + dhmax*L(k);
         end
      else
         dh = (h(n2)-h(n1))/L(k);
         if dh>dhmax
            h(n2) = h(n1) + dhmax*L(k);
         end
      end
   end
   if norm((h-h_old)./h,inf)<tol                                           % Test convergence
      break
   end
end

waitbar(0.75,wbar,'Triangulating quadtree');

t = MyDelaunayn(p);

waitbar(1,wbar);


%% SUB-FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function L = dist2poly(p,edgexy,lim)

% Find the minimum distance from the points in P to the polygon defined by
% the edges in EDGEXY. LIM is an optional argument that defines an upper
% bound on the distance for each point.

% Uses (something like?) a double sweep-line approach to reduce the number
% of edges that are required to be tested in order to determine the closest
% edge for each point. On average only size(EDGEXY)/4 comparisons need to
% be made for each point.

if nargin<3
   lim = [];
end
np = size(p,1);
ne = size(edgexy,1);
if isempty(lim)
   lim = inf*ones(np,1);
end

% Choose the direction with the biggest range as the "y-coordinate" for the
% test. This should ensure that the sorting is done along the best
% direction for long and skinny problems wrt either the x or y axes.
dxy = max(p)-min(p);
if dxy(1)>dxy(2)
    % Flip co-ords if x range is bigger
    p       = p(:,[2,1]);
    edgexy  = edgexy(:,[2,1,4,3]);
end

% Ensure edgexy(:,[1,2]) contains the lower y value
swap           = edgexy(:,4)<edgexy(:,2);
edgexy(swap,:) = edgexy(swap,[3,4,1,2]);

% Sort edges
[i,i]          = sort(edgexy(:,2));                                        % Sort edges by lower y value
edgexy_lower   = edgexy(i,:);
[i,i]          = sort(edgexy(:,4));                                        % Sort edges by upper y value
edgexy_upper   = edgexy(i,:);

% Mean edge y value
ymean = 0.5*( sum(sum(edgexy(:,[2,4]))) )/ne;

% Alloc output
L = zeros(np,1);

% Loop through points
tol = 100*eps*max(dxy);
for k = 1:np

   x = p(k,1);
   y = p(k,2);
   d = lim(k);

   if y<ymean

      % Loop through edges bottom up
      for j = 1:ne
         ymax = edgexy_lower(j,4);
         if ymax>=(y-d)
            ymin = edgexy_lower(j,2);
            if ymin<=(y+d)

               % Calculate the distance along the normal projection from [x,y]
               % to the jth edge
               x1 = edgexy_lower(j,1); x2mx1 = edgexy_lower(j,3)-x1;
               y1 = ymin;              y2my1 = ymax-y1;

               r = ((x-x1)*x2mx1+(y-y1)*y2my1)/(x2mx1^2+y2my1^2);
               if r>1                                                      % Limit to wall endpoints
                  r = 1;
               elseif r<0
                  r = 0;
               end

               dj = (x1+r*x2mx1-x)^2+(y1+r*y2my1-y)^2;
               if (dj<d^2) && (dj>tol)
                  d = sqrt(dj);
               end

            else
               break
            end
         end
      end

   else

      % Loop through edges top down
      for j = ne:-1:1
         ymin = edgexy_upper(j,2);
         if ymin<=(y+d)
            ymax = edgexy_upper(j,4);
            if ymax>=(y-d)

               % Calculate the distance along the normal projection from [x,y]
               % to the jth edge
               x1 = edgexy_upper(j,1); x2mx1 = edgexy_upper(j,3)-x1;
               y1 = ymin;              y2my1 = ymax-y1;

               r = ((x-x1)*x2mx1+(y-y1)*y2my1)/(x2mx1^2+y2my1^2);
               if r>1                                                      % Limit to wall endpoints
                  r = 1;
               elseif r<0
                  r = 0;
               end

               dj = (x1+r*x2mx1-x)^2+(y1+r*y2my1-y)^2;
               if (dj<d^2) && (dj>tol)
                  d = sqrt(dj);
               end

            else
               break
            end
         end
      end

   end

   L(k) = d;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = userhfun(x,y,fun,args,hmax,xymin,xymax)

% Evaluate user defined size function.

if ~isempty(fun)
   h = feval(fun,x,y,args{:});
   if size(h)~=size(x)
      error('Incorrect user defined size function. SIZE(H) must equal SIZE(X).');
   end
else
   h = inf*ones(size(x));
end
h = min(h,hmax);

% Limit to domain
out = (x>xymax(1))|(x<xymin(1))|(y>xymax(2))|(y<xymin(2));
h(out) = inf;

if any(h<=0)
   error('Incorrect user defined size function. H must be positive.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hmax,edgeh,fun,args] = gethdata(hdata)

% Check the user defined size functions

if ~isempty(hdata)
   if ~isstruct(hdata)
      error('HDATA must be a structure');
   end
   if numel(hdata)~=1
      error('HDATA cannot be an array of structures');
   end
   fields = fieldnames(hdata);
   names = {'hmax','edgeh','fun','args'};
   for k = 1:length(fields)
      if ~any(strcmp(fields{k},names))
         error('Invalid field in HDATA');
      end
   end
   if isfield(hdata,'hmax')
      if (numel(hdata.hmax)~=1) || (hdata.hmax<=0)
         error('HDATA.HMAX must be a positive scalar');
      else
         hmax = hdata.hmax;
      end
   else
      hmax = inf;
   end
   if isfield(hdata,'edgeh')
      if (numel(hdata.edgeh)~=2*size(hdata.edgeh,1)) || any(hdata.edgeh(:)<0)
         error('HDATA.EDGEH must be a positive Kx2 array');
      else
         edgeh = hdata.edgeh;
      end
   else
      edgeh = [];
   end
   if isfield(hdata,'fun')
      if ~ischar(hdata.fun) && ~isa(hdata.fun,'function_handle')
         error('HDATA.FUN must be a function name or a function handle');
      else
         fun = hdata.fun;
      end
   else
      fun = '';
   end
   if isfield(hdata,'args')
      if ~iscell(hdata.args)
         error(['HDATA.ARGS must be a cell array of additional' ...
            'inputs for HDATA.FUN']);
      else
         args = hdata.args;
      end
   else
      args = {};
   end
else
   hmax = inf;
   edgeh = [];
   fun = '';
   args = {};
end
