% Mesh Maker will triangulate the mesh and prompt user for the BC's of the
% domain by use of a GUI

function mesh_maker

% Select the mesh you want to make here
% Use the mesh collection or write your own 
% mesh.  Somehow, get:
%                   -p(npts,2)  = x,y location of nodes
%                   -t(ntri,3) = CCW ordering of nodes of triangles

[p,t] = mesh_collection(13);

% Run the data structure on the mesh and 
% set the boundary edges
set_bounds(p,t);


end