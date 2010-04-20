function Navier2d

% This is the GUI for Navier2d.
%
% Type "navier2d" to start the program.
%
% See the *.pdf tutorial for more information.
%
% The actual Navier-Stokes solver is in tvd_rk2.m and is where most of 
% the work is done, although the mesh based data structures are setup 
% in this function.
%
%
% Darren Engwirda - 2006
%
% Naver2d is Copyright (C) 2005-2006 Darren Engwirda. See "copyright.m" for
% full details.


clc, close all

% All units "normalized"

% GUI data structure
data = struct('buttons'  ,[], ...
              'init'     ,[], ...   % Initial conditions & flow variables
              'animation',[], ...   % Animation settings 
              'settings' ,[], ...   % Integration settings
              'mesh'     ,[], ...   % Mesh data
              'bc'       ,[], ...   % Boundary conditions
              'flag'     ,[]);      % Edges marked for lift/drag calc

% Main window
set(figure(1), ...
             'Name'       ,'Navier2d'  , ...
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
      
% % Integration frame
% ifrm = [0.55,0.6,0.35,0.3];
% uicontrol('Style'          ,'Frame'      , ...
%           'Units'          ,'Normalized' , ...
%           'Position'       ,ifrm         , ...
%           'BackgroundColor',frmcol);    
%       
% % Animation frame
% afrm = [0.55,0.2,0.35,0.3];
% uicontrol('Style'          ,'Frame'      , ...
%           'Units'          ,'Normalized' , ...
%           'Position'       ,afrm         , ...
%           'BackgroundColor',frmcol);          
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Frame Headers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mesh
uicontrol('Style'          ,'Text'                        , ...
          'Units'          ,'Normalized'                  , ...
          'Position'       ,[mfrm(1:2)+0.05,mfrm(3:4)-0.1], ...
          'String'         ,'Load Mesh'                , ...
          'FontSize'       ,12                            , ...
          'BackgroundColor',frmcol);
 
% % Integration      
% uicontrol('Style'          ,'Text'                        , ...
%           'Units'          ,'Normalized'                  , ...
%           'Position'       ,[ifrm(1:2)+0.05,ifrm(3:4)-0.1], ...
%           'String'         ,'Integration Settings'        , ...
%           'FontSize'       ,12                            , ...
%           'BackgroundColor',frmcol);  
      
% Boundary conditions      
uicontrol('Style'          ,'Text'                        , ...
          'Units'          ,'Normalized'                  , ...
          'Position'       ,[bfrm(1:2)+0.05,bfrm(3:4)-0.1], ...
          'String'         ,'Set Boundary Conditions'         , ...
          'FontSize'       ,12                            , ...
          'BackgroundColor',frmcol);    
      
% % Animation Settings      
% uicontrol('Style'          ,'Text'                        , ...
%           'Units'          ,'Normalized'                  , ...
%           'Position'       ,[afrm(1:2)+0.05,afrm(3:4)-0.1], ...
%           'String'         ,'Animation Settings'          , ...
%           'FontSize'       ,12                            , ...
%           'BackgroundColor',frmcol);          

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%                                Buttons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

btnx = 0.1;
btny = 0.05;

% % Integration settings
% data.button(1) = uicontrol('Style'   ,'PushButton'                         , ...
%                            'Units'   ,'Normalized'                         , ...
%                            'Position',[ifrm(1)+0.05,ifrm(2)+0.05,btnx,btny], ...
%                            'String'  ,'Set'                                , ...
%                            'Callback',@set_integration);
% % Main run           
% data.button(2) = uicontrol('Style'   ,'PushButton'         , ...
%                            'Units'   ,'Normalized'         , ...
%                            'Position',[0.45,0.05,btnx,btny], ...
%                            'String'  ,'Run'                , ...
%                            'Callback',@run);           
% Load mesh                       
data.button(3) = uicontrol('Style'   ,'PushButton'                          , ...
                           'Units'   ,'Normalized'                          , ...
                           'Position',[mfrm(1)+0.025,mfrm(2)+0.05,btnx,btny], ...
                           'String'  ,'Load'                                , ...
                           'Callback',@load_mesh);    
% Show mesh                       
data.button(5) = uicontrol('Style'   ,'PushButton'                          , ...
                           'Units'   ,'Normalized'                          , ...
                           'Position',[mfrm(1)+.18+0.125,mfrm(2)+0.05,btnx,btny], ...
                           'String'  ,'View'                                , ...
                           'Callback',@show_mesh);     
% Median mesh                       
data.button(6) = uicontrol('Style'   ,'PushButton'                          , ...
                           'Units'   ,'Normalized'                          , ...
                           'Position',[mfrm(1)+.36+0.225,mfrm(2)+0.05,btnx,btny], ...
                           'String'  ,'Median'                              , ...
                           'Callback',@show_median);                            
                       
% Set boundary conditions           
data.button(7) = uicontrol('Style'   ,'PushButton'                              , ...
                           'Units'   ,'Normalized'                              , ...
                           'Position',[bfrm(1)+.18+0.075,bfrm(2)+0.05,btnx+0.1,btny], ...
                           'String'  ,'Velocity/Pressure'                       , ...
                           'Callback',@set_bc);    
% % Animation select 
% data.button(8) = uicontrol('Style'   ,'PushButton'                          , ...
%                            'Units'   ,'Normalized'                          , ...
%                            'Position',[afrm(1)+0.05,afrm(2)+0.125,btnx,btny], ...
%                            'String'  ,'Set'                                 , ...
%                            'Callback',@set_anim);  
                       
% Main help button           
% data.button(9) = uicontrol('Style'   ,'PushButton'        , ...
%                            'Units'   ,'Normalized'        , ...
%                            'Position',[0.8,0.05,btnx,btny], ...
%                            'String'  ,'Help'              , ...
%                            'Callback',@help_main);                          

% % Initial conditions
% data.button(10) = uicontrol('Style'   ,'PushButton'                          , ...
%                             'Units'   ,'Normalized'                          , ...
%                             'Position',[ifrm(1)+0.05,ifrm(2)+0.125,btnx,btny], ...
%                             'String'  ,'Initial'                             , ...
%                             'Callback',@init);                         
%                        
% % Animation popup                         
% data.button(11) = uicontrol('Style'           ,'Popup'                                    , ...
%                             'Units'           ,'Normalized'                               , ...
%                             'Position'        ,[afrm(1)+0.2,afrm(2)+0.125,btnx+0.025,btny], ...
%                             'String'          ,{'Mesh';'Surf';'Contour'}                  , ...
%                             'BackgroundColor' ,frmcol);   
% % Animation popup                          
% data.button(12) = uicontrol('Style'           ,'Popup'                                   , ...
%                             'Units'           ,'Normalized'                              , ...
%                             'Position'        ,[afrm(1)+0.2,afrm(2)+0.05,btnx+0.025,btny], ...
%                             'String'          ,{'2D view','3D view'}                     , ...
%                             'BackgroundColor' ,frmcol); 
%                                                                      
% Save
data.button(15) = uicontrol('Style'   ,'PushButton'        , ...
                            'Units'   ,'Normalized'        , ...
                            'Position',[0.4,0.05,btnx,btny], ...
                            'String'  ,'Save'              , ...
                            'Callback',@save_data);  
                        
% % Flow type (UVP, tracer, thermal coupled)                         
% data.button(16) = uicontrol('Style'           ,'Popup'                                    , ...
%                             'Units'           ,'Normalized'                               , ...
%                             'Position'        ,[ifrm(1)+0.2,ifrm(2)+0.125,btnx+0.025,btny], ...
%                             'String'          ,{'Normal','Tracer','Thermal'}              , ...
%                             'BackgroundColor' ,frmcol); 

% % Set boundary conditions  (tracer)         
% data.button(17) = uicontrol('Style'  ,'PushButton'                            , ...
%                             'Units'   ,'Normalized'                           , ...
%                             'Position',[bfrm(1)+0.125,bfrm(2)+0.125,btnx,btny], ...
%                             'String'  ,'Tracer'                               , ...
%                             'Callback',@set_bc_tracer);  
                       
% % Select edges for lift/drag calculation         
% data.button(18) = uicontrol('Style'  ,'PushButton'                          , ...
%                             'Units'   ,'Normalized'                         , ...
%                             'Position',[afrm(1)+0.05,afrm(2)+0.05,btnx,btny], ...
%                             'String'  ,'Monitors'                           , ...
%                             'Callback',@residuals);                      
                       
                        
% Pass button handles                        
set(figure(1),'UserData',data);   


% % GNU Copyright
% copyright = {'Navier2d version 2.3, Copyright (C) 2005-2006 Darren Engwirda.'
%              ''
%              'Navier2d comes with ABSOLUTELY NO WARRANTY; for details see the'
%              'GNU license included.  This is free software, and you are welcome to'
%              'redistribute it under certain conditions; see the GNU license for details.'
%              };
% 
% uiwait(msgbox(copyright,'Copyright','none')); 
                       
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
e = data.mesh.e; % List of 2 points which connect each edge
bc = data.bc;  % This is the bc flag for each edge bc value 

e(:,3) = bc;      % Form 3 element array for EDGE (pt1,pt2,bc_type)

% Clear others
clear('data'); clear('varargin');

% Prompt for save dlg
uisave

mesh_name = 'mesh.msh'
grid_stats = [size(p,1),size(t,1),size(e,1) ];
comment = 'Mesh file for CYBO';
dlmwrite(mesh_name, comment,'');
dlmwrite(mesh_name,grid_stats,'delimiter',' ','-append');
dlmwrite(mesh_name,p,'delimiter',' ','-append');
dlmwrite(mesh_name,t,'delimiter',' ','-append');
dlmwrite(mesh_name,e,'delimiter',' ','-append');

return


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function set_integration(varargin)
% 
% % Prompt user for the integration settings
% 
% % Get current GUI data           
% data = get(figure(1),'UserData');   
% 
% prompts = {'Maximum number of steps'
%            'Maximum integration time (seconds)'
%            'CFL Number'
%            'Kinematic viscosity (m^2/s)'
%            'Output frequency (steps)'
%            'Kinematic viscosity (m^2/s) (tracer)'
%            'Temperature/velocity coupling'};
% 
% % Default answers
% if isempty(data.settings)
%     defaults = {'100','Inf','1.0','0.01','25','0.01','1'};
% else
%     for k = 1:length(prompts)
%         defaults{k} = num2str(data.settings(k));
%     end
% end
%      
% % Prompt for user input
% settings = inputdlg(prompts,'Integration Settings',1,defaults);
% 
% % Convert to double
% if ~isempty(settings)
%     data.settings = str2double(settings);
% else
%     return
% end
% 
% % Set new GUI data
% set(figure(1),'UserData',data);
% 
% return
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function init(varargin)
% 
% % Prompt for initial conditions
% 
% % Get GUI data
% data = get(figure(1),'UserData');
% 
% if isempty(data.mesh)
%     errordlg('Mesh file has not been loaded.','Error')
%     return
% end
% 
% % Prompt
% init = inputdlg({'U velocity (m/s)'; 'V velocity (m/s)'; 'S (tracer)'}, ...
%                  'Initial Conditions (Use ".*" and "./" in functions!)',1,{'0'; '0'; '0'});
% 
% % Nodes
% x    = data.mesh.p(:,1);
% y    = data.mesh.p(:,2);
% numn = length(x);
% 
% % Convert to double
% if ~isempty(init)
%     % Evaluate inputs
%     if checkvar(init{1}), U = eval(init{1}); else return, end
%     if checkvar(init{2}), V = eval(init{2}); else return, end
%     if checkvar(init{3}), S = eval(init{3}); else return, end
%     % Deal with scalar inputs
%     if numel(U)==1
%         U = repmat(U,numn,1);
%     end
%     if numel(V)==1
%         V = repmat(V,numn,1);
%     end
%     if numel(S)==1
%         S = repmat(S,numn,1);
%     end
% else
%     return
% end
% 
% % UVP fields
% data.init.U = U;
% data.init.V = V;
% data.init.S = S;
% data.init.P = repmat(0,numn,1);
% 
% % Set new GUI data
% set(figure(1),'UserData',data);
% 
% % Show IC's
% set(figure,'Name','Initial Conditions')
% subplot(1,3,1), trimesh(data.mesh.t,x,y,U), title('U velocity'), axis square
% subplot(1,3,2), trimesh(data.mesh.t,x,y,V), title('V velocity'), axis square
% subplot(1,3,3), trimesh(data.mesh.t,x,y,S), title('Tracer')    , axis square
% 
% return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ok = checkvar(expr)
% 
% % Helper for init
% 
% var = symvar(expr);     % Get variables in expr
% num = numel(var);       % Number of variables
% 
% ok = false;
% if num==0
%     ok = true;
% elseif num<=2
%     for i = 1:num       % Check that expr only contains 'x' and/or 'y'
%         if any(strcmp(var{i},{'x','y'}))
%             ok = true;
%         end
%     end
% end
% if ~ok
%     errordlg('Inputs must contain only constants or functions of [x,y]')
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function run(varargin)
% 
% % Main run button
% 
% % Get GUI data
% data = get(figure(1),'UserData');
% 
% % Check to see that everything has been specified
% if isempty(data.settings)
%     errordlg('Integration settings have not been specified','Error')
%     return
% end
% if isempty(data.mesh)
%     errordlg('Mesh file has not been loaded','Error')
%     return
% end
% if isempty(data.init)
%     errordlg('Initial conditions have not been set')
%     return
% end
% 
% % UVP BC check
% e            = data.mesh.e;
% be           = data.mesh.be;
% boundary     = false(size(e,1),1);
% boundary(be) = true;
% unassigned   = data.bc(:,1)==-1 & boundary;
% 
% % Solver settings
% solver = get(data.button(16),'Value');              % UVP, tracer, thermal
% if solver>1
%     unassigned = (data.bc(:,7)==-1 & boundary) | unassigned;
% end
% data.solver = solver;
% 
% if any(unassigned)
%     errordlg('Boundary conditions have not been specified','Error')
%     return
% end
% 
% % Animation settings
% data.animation(3) = get(data.button(11),'Value');   % Plot type
% data.animation(4) = get(data.button(12),'Value');   % Plot view
% 
% % CALL NAVIER-STOKES SOLVER
% tvd_rk2(data);      % 2 stage TVD Runge-Kutta solver 
% 
% return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function help_main(varargin)
% 
% % Quick help for the main window
% 
% helptext = {'Quick help for Navier2d.'
%             ''
%             ['The main window is broken into four sections, corresponding ', ...
%              'to the four main steps involved in setting up a simulation.']
%             ''
%             'Step 1: Load a mesh file'
%             'Step 2: Specify the boundary conditions'
%             'Step 3: Specify the integration settings and initial conditions'
%             'Step 4: Specify the animation settings'
%             ''
%             'The "Run" button will then start the simulation.'
%             ''
%             ['Mesh files must be created separately. The mesh generator mesh2d ', ...
%              'was developed for this purpose.']   
%            };
% 
% uiwait(msgbox(helptext,'Help','none')); 
% 
% return


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
% function help_bc(varargin)
% 
% % Quick help for the BC window
% 
% helptext = {'Quick help for boundary conditions.'
%             ''
%             ['This console shows the boundary edges of the mesh, with the '  , ...
%              'edge midpoints plotted. Boundary conditions must be specified ', ...
%              'for each edge.']
%             ''
%             ['Edges are first selected using the mouse by pressing the '   , ...
%              '"Select" button. Multiple edges can be selected/de-selected ', ...
%              'using this mode. "Clear" will clear the selection.']
%             ''
%             ['Boundary conditions are then specified for the highlighted edges ', ...
%              'by pressing the "Set" button. This will prompt the user for the ', ...
%              'boundary condition type.']
%             ''
%             'Velocity boundary conditions are used when the velocity is known'
%             'Pressure boundary conditions are used for outlets/outflows'
%             ''
%             ['If a velocity type is chosen the user is then asked to enter the ', ...
%              'xy velocity components.']
%             ''
%             'Velocity conditions are shown in blue, pressure conditions in green.'
%             ''
%             'Boundary conditions can be reassigned.'
%            };
% 
% uiwait(msgbox(helptext,'Help','none')); 
% 
% return
% 

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
        
%         while true
%             % Prompt for symmetry or wall
%             [j,ok] = listdlg('PromptString' , 'BC type: (CTRL+click)', ...
%                              'Name'         , 'Boundary Conditions'  , ...
%                              'ListString'   , {'U','dU/dn = 0','V','dV/dn = 0'});
% 
%             if ~ok
%                 return
%             end  
%             if length(j)~=2 || (j(1)~=1)&&(j(1)~=2) || (j(2)~=3)&&(j(2)~=4)
%                 uiwait(errordlg('Must make selections for U and V','Error'));
%             else
%                 break
%             end
%         end
        
%         % BC type
%         if j(1)==2, utype = 0; else utype = 1; end
%         if j(2)==4, vtype = 0; else vtype = 1; end
%         
%         % Prompt for user input
%         if (utype==1) && (vtype==1)
%             value = inputdlg({'U velocity (m/s)','V velocity (m/s)'},'Velocity components',1,{'0','0'});
%         elseif utype==1
%             value = inputdlg('U velocity (m/s)','Velocity components',1,{'0'});
%         elseif vtype==1
%             value = inputdlg('V velocity (m/s)','Velocity components',1,{'0'});
%         else
%             value = ['0','0'];
%         end
%         
%         % Convert to double
%         if ~isempty(value)
%             value = str2double(value);
%             if (utype==1) && (vtype==0), value = [value,0]; end
%             if (vtype==1) && (utype==0), value = [0,value]; end
%             % Assign to main GUI data
%             data.bc(be(in),1:6) = repmat([utype,value(1),vtype,value(2),0,0],sum(in),1);
%         else
%             return
%         end
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
unassigned   = data.bc(:,1)==0 & boundary;
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


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function set_bc_tracer(varargin)
% 
% % Set the boundary conditions
% 
% data = get(figure(1),'UserData');
% 
% if isempty(data.mesh)
%     errordlg('No mesh file loaded.','Error')
%     return
% end
% 
% set(figure, ...
%           'Name'        ,'Boundary Conditions (tracer)', ...
%           'DoubleBuffer','On'                          , ...
%           'Units'       ,'Normalized');
%       
% axes('Units'   ,'Normalized', ...
%      'Position',[0.25,0.1,0.7,0.8]);      
%      
% btnx = 0.1;
% btny = 0.05;
% 
% % Buttons
% b1 = uicontrol('Style','PushButton'           , ...
%                'Units','Normalized'           , ...
%                'Position',[0.05,0.3,btnx,btny], ...
%                'String','Select'              , ...
%                'Callback',@select);
%            
% b2 = uicontrol('Style','PushButton'           , ...
%                'Units','Normalized'           , ...
%                'Position',[0.05,0.1,btnx,btny], ...
%                'String','Set'                 , ...
%                'Callback',@bc_type_tracer);   
%            
% b3 = uicontrol('Style','PushButton'           , ...
%                'Units','Normalized'           , ...
%                'Position',[0.05,0.2,btnx,btny], ...
%                'String','Clear'               , ...
%                'Callback',@clear_sel);  
%            
% b4 = uicontrol('Style'   ,'PushButton'        , ...
%                'Units'   ,'Normalized'        , ...
%                'Position',[0.05,0.4,btnx,btny], ...
%                'String'  ,'Help'              , ...
%                'Callback',@help_bc);              
%            
% % Headers           
% uicontrol('Style'          ,'Text'              , ...
%           'Units'          ,'Normalized'        , ...
%           'Position'       ,[0.025,0.8,0.2,0.05], ...
%           'String'         ,'Black = Unassigned', ...
%           'BackgroundColor',[0.8,0.8,0.8]);  
% uicontrol('Style'          ,'Text'               , ...
%           'Units'          ,'Normalized'         , ...
%           'Position'       ,[0.025,0.75,0.2,0.05], ...
%           'String'         ,'Blue = Value'       , ...
%           'BackgroundColor',[0.8,0.8,0.8]);  
% uicontrol('Style'          ,'Text'               , ...
%           'Units'          ,'Normalized'         , ...
%           'Position'       ,[0.025,0.7,0.2,0.05], ...
%           'String'         ,'Yellow = Gradient'  , ...
%           'BackgroundColor',[0.8,0.8,0.8]);  
% 
% % Boundary edge geometry
% e  = data.mesh.e;
% be = data.mesh.be;
% p  = data.mesh.p;
% pe = data.mesh.pe;
% pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));
% 
% boundary     = false(size(e,1),1);
% pressure     = boundary;            % "Pressure" type doesn't exist for tracer. Set = false
% boundary(be) = true;
% unassigned   = data.bc(:,7)==-1 & boundary;
% velocity     = data.bc(:,7)== 1;    % Re-use "velocity" as value so that the other
% gradient     = data.bc(:,7)== 0;    % sub-functions will work...  
% 
% % Plot midpoints
% plot(pe(unassigned,1),pe(unassigned,2),'k.', ...
%      pe(velocity,1)  ,pe(velocity,2)  ,'b.', ...
%      pe(gradient,1)  ,pe(gradient,2)  ,'y.'), axis equal, axis off, hold on
% 
% % Plot edges
% patch('faces',e(unassigned,:),'vertices',data.mesh.p,'facecolor','none','edgecolor','k'); 
% patch('faces',e(velocity,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','b'); 
% patch('faces',e(gradient,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','y'); 
% 
% % New GUI data just for the "set bc" window
% bcdata = struct('p'         ,p                  , ...
%                 'pe'        ,pe                 , ...
%                 'e'         ,e                  , ...
%                 'be'        ,be                 , ...
%                 'in'        ,false(size(pm,1),1), ...
%                 'unassigned',unassigned         , ...
%                 'velocity'  ,velocity           , ...
%                 'pressure'  ,pressure           , ...
%                 'gradient'  ,gradient);
% 
% % Set GUI data
% set(gcf,'UserData',bcdata);
% 
% return


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function bc_type_tracer(varargin)
% 
% % Set the BC type and value
% 
% % Get GUI data from the "set bc" window
% fig    = gcf;
% bcdata = get(fig,'UserData');
% 
% % Boundary edge geometry
% e  = bcdata.e;
% be = bcdata.be;
% p  = bcdata.p;
% pe = bcdata.pe;
% in = bcdata.in;
% pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));
% 
% if sum(in)==0
%     errordlg('No edges selected.')
%     return
% end
% 
% % Get main GUI data
% data = get(figure(1),'Userdata'); figure(fig);
% 
% % Prompt with listbox
% [i,ok] = listdlg('PromptString' , 'BC type:'            , ...
%                  'Name'         , 'Boundary Conditions' , ...
%                  'SelectionMode', 'Single'              , ...
%                  'ListString'   , {'S','dS/dn = 0'});
%              
% if ok
%    if i==1      % Value
%        % Prompt for value
%        value = inputdlg('Tracer value','Set tracer value',1,{'0'});
%        if ~isempty(value)
%            % Assign to main GUI data
%            data.bc(be(in),7:8) = repmat([1,str2double(value)],sum(in),1);
%        else
%            return
%        end
%    else         % Gradient
%        % Assign to main GUI data
%        data.bc(be(in),7:8) = repmat([0,0],sum(in),1); 
%    end
% end
% 
% % Re-evaluate
% boundary     = false(size(e,1),1);
% pressure     = boundary;            % "Pressure" type doesn't exist for tracer. Set = false
% boundary(be) = true;
% unassigned   = data.bc(:,7)==-1 & boundary;
% velocity     = data.bc(:,7)== 1;    % Re-use "velocity" as value so that the other
% gradient     = data.bc(:,7)== 0;    % sub-functions will work...  
% 
% % Clear selection
% in                = false(size(pm,1),1);
% bcdata.in         = in;
% bcdata.unassigned = unassigned;
% bcdata.velocity   = velocity;
% bcdata.pressure   = pressure;
% bcdata.gradient   = gradient;
% 
% % Plot midpoints
% hold off
% 
% % Plot midpoints
% plot(pe(unassigned,1),pe(unassigned,2),'k.', ...
%      pe(velocity,1)  ,pe(velocity,2)  ,'b.', ...
%      pe(gradient,1)  ,pe(gradient,2)  ,'y.'), axis equal, axis off, hold on
% 
% % Plot edges
% patch('faces',e(unassigned,:),'vertices',data.mesh.p,'facecolor','none','edgecolor','k'); 
% patch('faces',e(velocity,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','b'); 
% patch('faces',e(gradient,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','y');  
%  
% set(fig,'UserData',bcdata);
% 
% % Set main GUI data
% set(figure(1),'UserData',data); figure(fig)   
% 
% return
% 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function set_anim(varargin)
% 
% % Prompt for animation variables
% 
% % Get main GUI data
% data = get(figure(1),'Userdata');
% 
% list = {'U velocity','V velocity','Total velocity','Pressure','Vorticity','Tracer'};
% 
% % Prompt with listbox
% [i,ok] = listdlg('PromptString', 'Select variables: (CTRL+click)', ...
%                  'Name'        , 'Animation'                     , ...
%                  'ListString'  , list);
%        
% if ok
%     if length(i)>2
%         errordlg('Can''t select more than 2 variables.')
%         return
%     else
%         if length(i)==0
%             data.animation(1) = 0;
%         else
%             data.animation(1) = i(1);
%         end
%         if length(i)==2
%             data.animation(2) = i(2);
%         else
%             data.animation(2) = 0;
%         end
%     end
% end
% 
% % Save GUI data
% set(figure(1),'Userdata',data);
%             
% return
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function residuals(varargin)
% 
% % Select residual type
% 
% data = get(figure(1),'UserData');
% 
% list = {'[U,V] residuals','[x,y] forces'};
% 
% % Prompt with listbox
% [i,ok] = listdlg('PromptString' , 'Select type', ...
%                  'Name'         , 'Monitors'   , ...
%                  'ListString'   , list         , ...
%                  'SelectionMode','Single');
% 
% if ok
%     if i==1     % Residuals
%         % Un-mark any edges for lift/drag
%         data.flag = false(size(data.flag));
%         % Reset main GUI data
%         set(figure(1),'UserData',data);
%     else        % Lift/drag
%         lift_drag;
%     end
% end
% 
% return
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function lift_drag(varargin)
% 
% % Mark edges at which the lift and drag forces will be calculated
% 
% data = get(figure(1),'UserData');
% 
% if isempty(data.mesh)
%     errordlg('No mesh file loaded.','Error')
%     return
% end
% 
% set(figure, ...
%           'Name'        ,'Lift/drag Forces', ...
%           'DoubleBuffer','On'              , ...
%           'Units'       ,'Normalized');
%       
% axes('Units'   ,'Normalized', ...
%      'Position',[0.25,0.1,0.7,0.8]);      
%      
% btnx = 0.1;
% btny = 0.05;
% 
% % Buttons
% b1 = uicontrol('Style','PushButton'           , ...
%                'Units','Normalized'           , ...
%                'Position',[0.05,0.3,btnx,btny], ...
%                'String','Select'              , ...
%                'Callback',@select);
%            
% b2 = uicontrol('Style','PushButton'           , ...
%                'Units','Normalized'           , ...
%                'Position',[0.05,0.1,btnx,btny], ...
%                'String','Mark'                , ...
%                'Callback',@set_lift_drag);   
%            
% b3 = uicontrol('Style','PushButton'           , ...
%                'Units','Normalized'           , ...
%                'Position',[0.05,0.2,btnx,btny], ...
%                'String','Clear'               , ...
%                'Callback',@clear_sel);           
%            
% % Headers           
% uicontrol('Style'          ,'Text'              , ...
%           'Units'          ,'Normalized'        , ...
%           'Position'       ,[0.025,0.8,0.2,0.05], ...
%           'String'         ,'Black = Unassigned', ...
%           'BackgroundColor',[0.8,0.8,0.8]);  
% uicontrol('Style'          ,'Text'               , ...
%           'Units'          ,'Normalized'         , ...
%           'Position'       ,[0.025,0.75,0.2,0.05], ...
%           'String'         ,'Blue = Marked'      , ...
%           'BackgroundColor',[0.8,0.8,0.8]);  
% 
% % Boundary edge geometry
% e  = data.mesh.e;
% be = data.mesh.be;
% p  = data.mesh.p;
% pe = data.mesh.pe;
% pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));
% 
% boundary     = false(size(e,1),1);
% pressure     = boundary;            % "Pressure" type doesn't exist. Set = false
% gradient     = boundary;            % "Gradient" type doesn't exist. Set = false
% boundary(be) = true;
% unassigned   = data.flag==0 & boundary;
% velocity     = data.flag==1;        % Re-use "velocity" as value so that the other
%                                     % sub-functions will work...  
% 
% % Plot midpoints
% plot(pe(unassigned,1),pe(unassigned,2),'k.', ...
%      pe(velocity,1)  ,pe(velocity,2)  ,'b.'), axis equal, axis off, hold on
% 
% % Plot edges
% patch('faces',e(unassigned,:),'vertices',data.mesh.p,'facecolor','none','edgecolor','k'); 
% patch('faces',e(velocity,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','b'); 
% 
% % New GUI data just for the "set bc" window
% bcdata = struct('p'         ,p                  , ...
%                 'pe'        ,pe                 , ...
%                 'e'         ,e                  , ...
%                 'be'        ,be                 , ...
%                 'in'        ,false(size(pm,1),1), ...
%                 'unassigned',unassigned         , ...
%                 'velocity'  ,velocity           , ...
%                 'pressure'  ,pressure           , ...
%                 'gradient'  ,gradient);
% 
% % Set GUI data
% set(gcf,'UserData',bcdata);
% 
% return
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function set_lift_drag(varargin)
% 
% % Record marked edges (Lift/drag calc) in the GUI data
% 
% % Get GUI data from the "set bc" window
% fig    = gcf;
% bcdata = get(fig,'UserData');
% 
% % Boundary edge geometry
% e  = bcdata.e;
% be = bcdata.be;
% p  = bcdata.p;
% pe = bcdata.pe;
% in = bcdata.in;
% pm = 0.5*(p(e(be,1),:)+p(e(be,2),:));
% 
% if sum(in)==0
%     errordlg('No edges selected.')
%     return
% end
% 
% % Get main GUI data
% data = get(figure(1),'Userdata'); figure(fig);
% 
% % Set main GUI data
% data.flag(be) = in;
% 
% boundary     = false(size(e,1),1);
% pressure     = boundary;            % "Pressure" type doesn't exist. Set = false
% gradient     = boundary;            % "Gradient" type doesn't exist. Set = false
% boundary(be) = true;
% unassigned   = data.flag==0 & boundary;
% velocity     = data.flag==1;        % Re-use "velocity" as value so that the other
%                                     % sub-functions will work... 
% % Clear selection
% in                = false(size(pm,1),1);
% bcdata.in         = in;
% bcdata.unassigned = unassigned;
% bcdata.velocity   = velocity;
% bcdata.pressure   = pressure;
% bcdata.gradient   = gradient;
% 
% % Plot midpoints
% hold off
% 
% % Plot midpoints
% plot(pe(unassigned,1),pe(unassigned,2),'k.', ...
%      pe(velocity,1)  ,pe(velocity,2)  ,'b.'), axis equal, axis off, hold on
% 
% % Plot edges
% patch('faces',e(unassigned,:),'vertices',data.mesh.p,'facecolor','none','edgecolor','k'); 
% patch('faces',e(velocity,:)  ,'vertices',data.mesh.p,'facecolor','none','edgecolor','b'); 
%  
% set(fig,'UserData',bcdata);
% 
% % Set main GUI data
% set(figure(1),'UserData',data); figure(fig)
% 
% return


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
