function set_bounds_v2(p,t)

clc;
close all;

% GUI data structure
data = struct('buttons'  ,[], ...
              'init'     ,[], ...   % Initial conditions & flow variables
              'animation',[], ...   % Animation settings 
              'settings' ,[], ...   % Integration settings
              'mesh'     ,[], ...   % Mesh data
              'bc'       ,[], ...   % Boundary conditions
              'flag'     ,[]);      % Edges marked for lift/drag calc
          
          % Run fixmesh to make sure the mesh is CCW and non-duplicate
[p,t] = fixmesh(p,t);

% Setup the mesh based data structures
data.mesh      = data_structure(p,t);
data.bc        = repmat(-1,size(data.mesh.e,1),1);
data.init      = [];
data.settings  = [];
data.animation = [];
data.flag      = false(size(data.mesh.e,1),1);

% Set GUI data
set(figure(1),'UserData',data);
          

% Main window
set(figure(1), ...
             'Name'       ,'CYBO PreProcessor'  , ...
             'Units'      ,'Normalized', ...
             'NumberTitle','Off'       , ...
             'UserData',data); axis off
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
%                                 Frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frmcol = [0.9,0.9,0.825];

% Mesh frame
mfrm = [0.1,0.6,0.7,0.3];
uicontrol('Style'          ,'Frame'      , ...
          'Units'          ,'Normalized' , ...
          'Position'       ,mfrm         , ...
          'BackgroundColor',frmcol);
      
% BC frame
bfrm = [0.1,0.2,0.7,0.3];
uicontrol('Style'          ,'Frame'      , ...
          'Units'          ,'Normalized' , ...
          'Position'       ,bfrm         , ...
          'BackgroundColor',frmcol);    
           
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Frame Headers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mesh
uicontrol('Style'          ,'Text'                        , ...
          'Units'          ,'Normalized'                  , ...
          'Position'       ,[mfrm(1:2)+0.05,mfrm(3:4)-0.1], ...
          'String'         ,'Mesh View'                , ...
          'FontSize'       ,12                            , ...
          'BackgroundColor',frmcol);

% Boundary conditions      
uicontrol('Style'          ,'Text'                        , ...
          'Units'          ,'Normalized'                  , ...
          'Position'       ,[bfrm(1:2)+0.05,bfrm(3:4)-0.1], ...
          'String'         ,'Set Boundary Conditions'         , ...
          'FontSize'       ,12                            , ...
          'BackgroundColor',frmcol);    
               

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%                                Buttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

btnx = 0.1;
btny = 0.05;
   
% Show mesh                       
data.button(5) = uicontrol('Style'   ,'PushButton'                          , ...
                           'Units'   ,'Normalized'                          , ...
                           'Position',[mfrm(1)+.1,mfrm(2)+0.05,btnx,btny], ...
                           'String'  ,'View'                                , ...
                           'Callback',@show_mesh);     
% Median mesh                       
data.button(6) = uicontrol('Style'   ,'PushButton'                          , ...
                           'Units'   ,'Normalized'                          , ...
                           'Position',[mfrm(1)+.5,mfrm(2)+0.05,btnx,btny], ...
                           'String'  ,'Median'                              , ...
                           'Callback',@show_median);                            
                       
% Set boundary conditions           
data.button(7) = uicontrol('Style'   ,'PushButton'                              , ...
                           'Units'   ,'Normalized'                              , ...
                           'Position',[bfrm(1)+.18+0.075,bfrm(2)+0.05,btnx+0.1,btny], ...
                           'String'  ,'Velocity/Pressure'                       , ...
                           'Callback',@set_bc);    

                                                                    
% Save
data.button(15) = uicontrol('Style'   ,'PushButton'        , ...
                            'Units'   ,'Normalized'        , ...
                            'Position',[0.4,0.05,btnx,btny], ...
                            'String'  ,'Save'              , ...
                            'Callback',@save_data);  
                               
                                            
% Pass button handles                        
set(figure(1),'UserData',data);   

                       
return

           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function save_data(varargin)

% Save flow data

data = get(figure(1),'UserData');

%if ~isfield(data.init,'time')
%    errordlg('Simulation has not been run.','Error')
%    return
%end

% Load vars into current workspace
p    = data.mesh.p;
t    = data.mesh.t;
%U    = data.init.U;
%V    = data.init.V;
%S    = data.init.S;
%P    = data.init.P;
%W    = data.init.W;
%time = data.init.time;
%resx = data.init.resx;
%resy = data.init.resy;
epts = data.mesh.e;                 % List of 2 points which connect each edge
e2t = data.mesh.e2t;        % List of two triangles separated by the edge
bc = data.bc;  % This is the bc flag for each edge bc value 

% Write edge data for diamond structure - C.Yu 5/5/10
n1 = epts(1);
n2 = epts(2);
t1 = e2t(1);
t2 = e2t(2);

if (t1 == 0) 
    nt1 = 0;
else 
    for ind=1:3
        if (t1(ind) ~= n1 & t1(ind) ~= n2) nt1 = t1(ind); it1 = ind; end
    end
end
    
if (t2 == 0) 
    nt2 = 0;
else
    for ind=1:3
        if (t2(ind) ~= n1 & t2(ind) ~= n2) nt2 = t2(ind); it2 = ind; end
    end
end

e(:,1) = nt1; 
if (t1 ~= 0) e(:,2) = t1(it1+1);
else         e(:,2) = t2(it2-1);
e(:,3) = nt2;
if (t2 ~= 0) e(:,4) = t2(it2+1);
else         e(:,4) = t1(it1-1);   

e(:,5) = bc;    


% Clear others
clear('data'); clear('varargin');

% Prompt for save dlg
%uisave

mesh_name = input('Save the mesh (*.msh) file as:','s');
if (mesh_name(end-3:end) == '.msh')
else
    mesh_name = [mesh_name,'.msh'];
end

grid_stats = [size(p,1),size(t,1),size(e,1) ];
comment = 'Mesh file for CYBO';
dlmwrite(mesh_name, comment,'');
dlmwrite(mesh_name,grid_stats,'delimiter',' ','-append');
dlmwrite(mesh_name,p,'delimiter',' ','-append');
dlmwrite(mesh_name,t,'delimiter',' ','-append');
dlmwrite(mesh_name,e,'delimiter',' ','-append');

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function load_mesh(varargin)

% Prompt user to load in a mesh file.

uiload

% Get current GUI data
data = get(figure(1),'UserData');

% Error checking
if ~exist('p','var')||~exist('t','var')
    errordlg('Not a valid mesh file.','Error')
    return
end
if (size(p,2)~=2)||(size(t,2)~=3)
    errordlg('Incorrect dimensions for p or t.','Error')
    return
end
if (max(t(:))>size(p,1))||(min(t(:))<=0)
    errordlg('t is not a valid triangulation of the nodes in p.','Error')
    return
end

% Run fixmesh to make sure the mesh is CCW and non-duplicate
[p,t] = fixmesh(p,t);

% Setup the mesh based data structures
data.mesh      = data_structure(p,t);
data.bc        = repmat(-1,size(data.mesh.e,1),1);
data.init      = [];
data.settings  = [];
data.animation = [];
data.flag      = false(size(data.mesh.e,1),1);

% Set GUI data
set(figure(1),'UserData',data);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_mesh(varargin)

% Show the current mesh

data = get(figure(1),'UserData');

if isempty(data.mesh)
    errordlg('No mesh file loaded.','Error')
    return
end

set(figure,'Name','Mesh'); 

patch('faces',data.mesh.t,'vertices',data.mesh.p,'facecolor','w','edgecolor','b'); 

axis equal, axis off

title([num2str(size(data.mesh.p,1)),' Nodes, ', ...
       num2str(size(data.mesh.t,1)),' Triangles'])

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function show_median(varargin)

% Show the median dual mesh used by NS.m

% Get GUI data
data = get(figure(1),'UserData');

if isempty(data.mesh)
    errordlg('No mesh file loaded.','Error')
    return
end

set(figure,'Name','Median Dual Mesh'); 

p    = data.mesh.p;     % Nodes
t    = data.mesh.t;     % Triangulation
pc   = data.mesh.pc;    % Centroids
n2n  = data.mesh.n2n;   % Node to node connectivity
n2e  = data.mesh.n2e;   % Node to edge connectivity
e2t  = data.mesh.e2t;   % Edge to triangle connectivity
numn = size(p,1);

w = waitbar(0,'Building median mesh'); drawnow

xline = repmat(0,2*size(n2n,2)*numn,2);
yline = xline;
next  = 1;
lim   = 0;
for k = 1:numn
    m = 1;
    n = k/numn;
    if (n>lim)||(n==1)
        waitbar(n,w); lim = n+0.1;
    end
    while n2n(k,m)>0    % Loop around neighbours of node k
        % Triangles associated with current edge
        t1 = e2t(n2e(k,m),1);
        t2 = e2t(n2e(k,m),2);
        % Edge midpoint
        mx = 0.5*(p(k,1)+p(n2n(k,m),1));
        my = 0.5*(p(k,2)+p(n2n(k,m),2));
        % Median edge 1 (joins pc1 and pm)
        xline(next,1) = mx;
        xline(next,2) = pc(t1,1);
        yline(next,1) = my;
        yline(next,2) = pc(t1,2);
        next          = next+1;
        % Median edge 2 (joins pc2 and pm)
        if t2>0
            xline(next,1) = mx;
            xline(next,2) = pc(t2,1);
            yline(next,1) = my;
            yline(next,2) = pc(t2,2);
            next          = next+1; 
        end
        % Counter
        m = m+1;
    end
end
close(w)

% Triangles
patch('faces',t,'vertices',p,'facecolor','none','edgecolor','w'); axis equal, axis off, hold on
% Median CV's
patch('xdata',xline','ydata',yline','facecolor','none','edgecolor','b'); 

title([num2str(size(data.mesh.p,1)),' Nodes, ', ...
       num2str(size(data.mesh.t,1)),' Triangles'])

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_bc(varargin)

% Set the boundary conditions

data = get(figure(1),'UserData');

if isempty(data.mesh)
    errordlg('No mesh file loaded.','Error')
    return
end

set(figure, ...
          'Name'        ,'Boundary Conditions', ...
          'DoubleBuffer','On'                 , ...
          'Units'       ,'Normalized');
      
axes('Units'   ,'Normalized', ...
     'Position',[0.25,0.1,0.7,0.8]);      
     
btnx = 0.1;
btny = 0.05;

% Buttons
b1 = uicontrol('Style','PushButton'           , ...
               'Units','Normalized'           , ...
               'Position',[0.05,0.3,btnx,btny], ...
               'String','Select'              , ...
               'Callback',@select);
           
b2 = uicontrol('Style','PushButton'           , ...
               'Units','Normalized'           , ...
               'Position',[0.05,0.1,btnx,btny], ...
               'String','Set'                 , ...
               'Callback',@bc_type);   
           
b3 = uicontrol('Style','PushButton'           , ...
               'Units','Normalized'           , ...
               'Position',[0.05,0.2,btnx,btny], ...
               'String','Clear'               , ...
               'Callback',@clear_sel);  
           
% b4 = uicontrol('Style'   ,'PushButton'        , ...
%                'Units'   ,'Normalized'        , ...
%                'Position',[0.05,0.4,btnx,btny], ...
%                'String'  ,'Help'              , ...
%                'Callback',@help_bc);              
           
% Headers           
uicontrol('Style'          ,'Text'              , ...
          'Units'          ,'Normalized'        , ...
          'Position'       ,[0.025,0.8,0.2,0.05], ...
          'String'         ,'Black = Unassigned', ...
          'BackgroundColor',[0.8,0.8,0.8]);  
uicontrol('Style'          ,'Text'               , ...
          'Units'          ,'Normalized'         , ...
          'Position'       ,[0.025,0.75,0.2,0.05], ...
          'String'         ,'Blue = Freestream'    , ...
          'BackgroundColor',[0.8,0.8,0.8]);  
uicontrol('Style'          ,'Text'              , ...
          'Units'          ,'Normalized'        , ...
          'Position'       ,[0.025,0.7,0.2,0.05], ...
          'String'         ,'Green = Slipwall'  , ...
          'BackgroundColor',[0.8,0.8,0.8]);  
uicontrol('Style'          ,'Text'               , ...
          'Units'          ,'Normalized'         , ...
          'Position'       ,[0.025,0.65,0.2,0.05], ...
          'String'         ,'Yellow = Outflow'  , ...
          'BackgroundColor',[0.8,0.8,0.8]);  

% Boundary edge geometry
e  = data.mesh.e;
be = data.mesh.be;
p  = data.mesh.p;
pe = data.mesh.pe;
pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));

boundary     = false(size(e,1),1);
boundary(be) = true;
unassigned   = data.bc(:,1)==-1 & boundary;
freestream     = data.bc(:,1)== 1;
slipwall     = data.bc(:,1) == 2;
outflow     = data.bc(:,1) == 3;
%velocity     = data.bc(:,1)== 1;
%pressure     = data.bc(:,5) > 0;
%gradient     = xor(data.bc(:,1),data.bc(:,3));

% Plot midpoints
plot(pe(unassigned,1),pe(unassigned,2),'k.', ...
     pe(freestream,1)  ,pe(freestream,2)  ,'b.', ...
     pe(slipwall,1)  ,pe(slipwall,2)  ,'g.', ...
     pe(outflow,1)  ,pe(outflow,2)  ,'y.'), axis equal, axis off, hold
%      pe(velocity,1)  ,pe(velocity,2)  ,'b.', ...
%      pe(pressure,1)  ,pe(pressure,2)  ,'g.', ...
%      pe(gradient,1)  ,pe(gradient,2)  ,'y.'), axis equal, axis off, hold on

% Plot edges
patch('faces',e(unassigned,:),'vertices',data.mesh.p,'facecolor','none','edgecolor','k'); 
patch('faces',e(freestream,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','b'); 
patch('faces',e(slipwall,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','g'); 
patch('faces',e(outflow,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','y'); 

% Plot arrows for velocity type
%if any(velocity)
%   quiver(pe(velocity,1),pe(velocity,2),data.bc(velocity,2),data.bc(velocity,4))
%end


% New GUI data just for the "set bc" window
bcdata = struct('p'         ,p                  , ...
                'pe'        ,pe                 , ...
                'e'         ,e                  , ...
                'be'        ,be                 , ...
                'in'        ,false(size(pm,1),1), ...
                'unassigned',unassigned         , ...
                'velocity'  ,freestream           , ...
                'pressure'  ,slipwall           , ...
                'gradient'  ,outflow);

% Set GUI data
set(gcf,'UserData',bcdata);

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clear_sel(varargin)

% Clear current mouse selection 

% Get GUI data
fig    = gcf; 
ax     = gca;
bcdata = get(fig,'UserData');

% Boundary edge geometry
e  = bcdata.e;
be = bcdata.be;
p  = bcdata.p;
pe = bcdata.pe;
pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));

% Clear selection
in        = false(size(pm,1),1);
bcdata.in = in;

% Plot midpoints
plot(pe(bcdata.unassigned,1),pe(bcdata.unassigned,2),'k.', ...
     pe(bcdata.velocity,1)  ,pe(bcdata.velocity,2)  ,'b.', ...
     pe(bcdata.pressure,1)  ,pe(bcdata.pressure,2)  ,'g.', ...
     pe(bcdata.gradient,1)  ,pe(bcdata.gradient,2)  ,'y.', ...
     pm(in,1)               ,pm(in,2)               ,'r.');

set(fig,'UserData',bcdata); 
 
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function select(varargin)

% Select boundary midpoints using the mouse

% Get GUI data
fig    = gcf; 
ax     = gca;
bcdata = get(fig,'UserData');

% Boundary edge geometry
e  = bcdata.e;
be = bcdata.be;
p  = bcdata.p;
pe = bcdata.pe;
in = bcdata.in;
pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));

% Plot midpoints
plot(pe(bcdata.unassigned,1),pe(bcdata.unassigned,2),'k.', ...
     pe(bcdata.velocity,1)  ,pe(bcdata.velocity,2)  ,'b.', ...
     pe(bcdata.pressure,1)  ,pe(bcdata.pressure,2)  ,'g.', ...
     pe(bcdata.gradient,1)  ,pe(bcdata.gradient,2)  ,'y.', ...
     pm(in,1)               ,pm(in,2)               ,'r.');

title('Right click to confirm') 
 
% Set mouse to "crosshair"
set(fig,'Pointer','crosshair');

while true
    
    % Wait for mouse click
    waitforbuttonpress
    
    % Grab type of mouse click
    btn = get(fig,'SelectionType');
    
    if strcmp(btn,'normal')
        % Draw selection box
        p1 = get(ax,'CurrentPoint');
        rbbox
        p2 = get(ax,'CurrentPoint');
        % xy co-ords within axis
        p1 = p1(1,1:2); p2 = p2(1,1:2);
        % Sorted (left, right, bottom, top)
        x1 = min(p1(1),p2(1)); x2 = max(p1(1),p2(1));
        y1 = min(p1(2),p2(2)); y2 = max(p1(2),p2(2));
        % Find nodes within selection
        in = xor((pm(:,1)>=x1&pm(:,1)<=x2&pm(:,2)>=y1&pm(:,2)<=y2),in);
        % Show selected
        plot(pe(bcdata.unassigned,1),pe(bcdata.unassigned,2),'k.', ...
             pe(bcdata.velocity,1)  ,pe(bcdata.velocity,2)  ,'b.', ...
             pe(bcdata.pressure,1)  ,pe(bcdata.pressure,2)  ,'g.', ...
             pe(bcdata.gradient,1)  ,pe(bcdata.gradient,2)  ,'y.', ...
             pm(in,1)               ,pm(in,2)               ,'r.');
    else
        break
    end

end
bcdata.in = in; title('')

% Pass selection as userdata
set(fig,'UserData',bcdata,'Pointer','arrow');

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bc_type(varargin)

% Set the BC type and value

% Get GUI data from the "set bc" window
fig    = gcf;
bcdata = get(fig,'UserData');

% Boundary edge geometry
e  = bcdata.e;
be = bcdata.be;
p  = bcdata.p;
pe = bcdata.pe;
in = bcdata.in;
pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));

if sum(in)==0
    errordlg('No edges selected.')
    return
end

% Get main GUI data
data = get(figure(1),'Userdata'); figure(fig);

% Prompt with listbox
[i,ok] = listdlg('PromptString' , 'BC type:'            , ...
                 'Name'         , 'Boundary Conditions' , ...
                 'SelectionMode', 'Single'              , ...
                 'ListString'   , {'Freestream','Slipwall','Outflow'});
             
if ok
    if i==1     % FreeStream
        
        data.bc(be(in),1) = repmat( 1 ,sum(in),1);  
    else        % Pressure
        if i==2
            % Assign to main GUI data
            data.bc(be(in),1) = repmat( 2 ,sum(in),1);      % Extrapolated BC
        else
            % Assign to main GUI data
            data.bc(be(in),1) = repmat( 3 ,sum(in),1);      % Fixed pressure BC
        end
    end
end

% Re-evaluate
boundary     = false(size(e,1),1);
boundary(be) = true;
unassigned   = data.bc(:,1)==-1 & boundary;
freestream     = data.bc(:,1)== 1;
slipwall     = data.bc(:,1)== 2;
outflow     = data.bc(:,1)== 3;


% Clear selection
in                = false(size(pm,1),1);
bcdata.in         = in;
bcdata.unassigned = unassigned;
bcdata.freestream   = freestream;
bcdata.slipwall   = slipwall;
bcdata.outflow   = outflow;

% Plot midpoints
hold off

plot(pe(bcdata.unassigned,1),pe(bcdata.unassigned,2),'k.', ...
     pe(bcdata.freestream,1)  ,pe(bcdata.freestream,2)  ,'b.', ...
     pe(bcdata.slipwall,1)  ,pe(bcdata.slipwall,2)  ,'g.', ...
     pe(bcdata.outflow,1)  ,pe(bcdata.outflow,2)  ,'y.', ...
     pm(in,1)               ,pm(in,2)               ,'r.'), axis equal, axis off, hold on

% Plot edges
patch('faces',e(unassigned,:),'vertices',data.mesh.p,'facecolor','none','edgecolor','k'); 
patch('faces',e(freestream,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','b'); 
patch('faces',e(slipwall,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','g');  
patch('faces',e(outflow,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','y');  

% Plot arrows for velocity type
% if any(velocity)
%     quiver(pe(velocity,1),pe(velocity,2),data.bc(velocity,2),data.bc(velocity,4))
% end
 
set(fig,'UserData',bcdata);

% Set main GUI data
set(figure(1),'UserData',data); figure(fig)

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p,t] = fixmesh(p,t)

% Remove duplicate nodes & setup CCW (Counter-clockwise) node ordering.

% Darren Engwirda - 2006. Adapted from Distmesh.

% Remove duplicate nodes
snap        = max(max(p,[],1)-min(p,[],1),[],2)*1024*eps;
[foo,ix,jx] = unique(round(p/snap)*snap,'rows');
p           = p(ix,:);
t           = jx(t);
[pix,ix,jx] = unique(t);
t           = reshape(jx,size(t));
p           = p(pix,:);

% CCW ordering
p1  = p(t(:,1),:); 
d12 = p(t(:,2),:)-p1;
d13 = p(t(:,3),:)-p1;

% Negative area?
flip          = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))<0;
t(flip,[1,2]) = t(flip,[2,1]);

return

