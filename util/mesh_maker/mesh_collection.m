function [p,t] = mesh_collection(num)

%  MESH_COLLECTION: Collection of meshing examples from MESH2D users. 
%
%  mesh_collection(n) will run the nth example.
%
%  1. Simple square domain. Used for "driven cavity" CFD studies.
%
%  2. Rectangular domain with circular hole. Used in thermally coupled CFD
%     studies to examine the flow around a heated pipe.
%
%  3. Rectangular domain with circular hole and user defined size
%     functions. Used in a CFD study to examine vortex shedding about
%     cylinders.
%
%  4. Rectangular domain with 2 circular holes and user defined size
%     functions. Used in a CFD study to examine the unsteady flow between
%     cylinders.
%
%  5. Rectangular domain with square hole and user defined size functions.
%     Used in a CFD study to examine vortex shedding about square prisms.
%
%  6. 3 element airfoil with user defined size functions and boundary layer
%     size functions. Used in a CFD study to examin the lift/drag
%     characteristics.
%
%  7. U shaped domain.
%
%  8. Rectangular domain with step. Used for "backward facing step" CFD
%     studies.
%
%  9. NACA airfoil with boundary layer size functions. Used in a CFD study
%     to examine the lift/drag vs. alpha characteristics.
%
%  10. Wavy channel from Kong Zour. Used in a CFD study to examine unsteady
%      behaviour.
%
%  11. Tray of glass beads from Falk Hebe. Used in a CFD study to examine the flow
%      through past a collection of beads.
%
%  12. Ideally expanded nozzle geometry
%
%  13. NACA airfoil generator
%
% I am always looking for new meshes to add to the collection, if you would
% like to contribute please send me an email with an m-file description of
% the NODE, EDGE, HDATA and OPTIONS used to setup the mesh.
%
% Darren Engwirda    : 2006-2007
% Email              : d_engwirda@hotmail.com


switch(num)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 1

      node  = [0,0; 1,0; 1,1; 0,1];

      hdata.hmax = .1;%0.01;

      [p,t] = mesh2d(node,[],hdata);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 2

      theta = 0:pi/50:(2*pi-pi/50);
      x = cos(theta)/2;
      y = sin(theta)/2;

      node = [ x',y'; -5,-5; 5,-5; 5,5; -5,5];

      n = size(node,1)-4;
      edge = [(1:n-1)' (2:n)'; n,1; n+1,n+2; n+2,n+3; n+3,n+4; n+4,n+1];

      hdata.hmax = 0.15;
      [p,t] = mesh2d(node,edge,hdata);
      
      save sink p t

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 3

      theta = 0:pi/50:(2*pi-pi/50);
      x = cos(theta)/2;
      y = sin(theta)/2;

      node = [ x',y'; -5,-10; 25,-10; 25,10; -5,10];

      n = size(node,1)-4;
      edge = [(1:n-1)' (2:n)'; n,1; n+1,n+2; n+2,n+3; n+3,n+4; n+4,n+1];

      hdata.fun = @const_h;
      hdata.args = {-1,25,-3,3,0.1};
      options.dhmax = 0.2;

      [p,t] = mesh2d(node,edge,hdata,options);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 4

      theta = 0:pi/36:(2*pi-pi/36);
      x = cos(theta)/2;
      y = sin(theta)/2;

      cyL1 = [x' y'+1];
      cyL2 = [x' y'-1];
      box  = [-5,-10; 25,-10; 25,10; -5,10];

      n1 = size(cyL1,1);
      n2 = size(cyL2,1);
      c1 = [(1:n1-1)',(2:n1)'; n1,1];
      c2 = [(1:n2-1)',(2:n2)'; n2,1];
      c3 = [1,2; 2,3; 3,4; 4,1];

      node = [cyL1; cyL2; box];
      edge = [c1; c2+n1; c3+n1+n2];

      hdata.fun = @const_h;
      hdata.args = {-1,25,-4,4,0.2};

      [p,t] = mesh2d(node,edge,hdata);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 5

      node = [0,-10; 20,-10; 20,10; 0,10; 5,-0.5; 6,-0.5; 6,0.5; 5,0.5];
      edge = [1,2; 2,3; 3,4; 4,1; 5,6; 6,7; 7,8; 8,5];

      hdata.fun = @case5;
      hdata.edgeh = [5, 0.05; 6, 0.05; 7, 0.05; 8, 0.05];
      hdata.fun = @const_h;
      hdata.args = {5,20,-3,3,0.1};

      options.dhmax = 0.15;

      [p,t] = mesh2d(node,edge,hdata,options);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 6

      temp = [
         0.027490  0.017991 0.899870  0.017200 0.052890 -0.021420 1.214624 -0.113530
         0.021231  0.013241 0.871260  0.019430 0.048640 -0.021060 1.193291 -0.106968
         0.011552  0.004325 0.835990  0.021190 0.043550 -0.019330 1.150970 -0.093895
         -0.004135 -0.011795 0.802400  0.021820 0.039060 -0.015140 1.123617 -0.085397
         -0.012160 -0.021418 0.766470  0.020920 0.037490 -0.009480 1.097197 -0.077260
         -0.018975 -0.031128 0.753490  0.019960 0.039190 -0.001600 1.058932 -0.065674
         -0.022844 -0.037478 0.735200  0.017900 0.043650  0.006190 1.034620 -0.058528
         -0.027273 -0.046152 0.729540  0.017030 0.048240  0.010880 0.996129 -0.047562
         -0.029451 -0.052343 0.718900  0.015170 0.055920  0.016870 0.967841 -0.039703
         -0.030201 -0.061288 0.709580  0.013210 0.062110  0.020990 0.953663 -0.035852
         -0.028411 -0.069148 0.701600  0.011010 0.072390  0.026910 0.940446 -0.032286
         -0.024596 -0.073614 0.690290  0.007120 0.083730  0.032240 0.926244 -0.028500
         -0.018756 -0.075514 0.682300  0.003430 0.104760  0.040220 0.917211 -0.026096
         -0.017118 -0.075335 0.676980 -0.000170 0.117300  0.044210 0.909167 -0.023871
         -0.016780 -0.076060 0.672990 -0.004160 0.134530  0.047900 0.905108 -0.022393
         -0.025493 -0.079097 0.669990 -0.011380 0.153460  0.050670 0.901391 -0.019944
         -0.035315 -0.082430 0.673320 -0.019060 0.168330  0.052460 0.900339 -0.012133
         -0.042170 -0.084269 0.677310 -0.020330 0.172890  0.052960 0.907722 -0.004509
         -0.049084 -0.085176 0.677310 -0.020960 0.182440  0.054030 0.915378 -0.001953
         -0.055933 -0.084663 0.668660 -0.021720 0.200270  0.055850 0.920255 -0.001185
         -0.059101 -0.083382 0.634730 -0.024850 0.203160  0.056150 0.926453 -0.000706
         -0.062122 -0.081635 0.602130 -0.028010 0.234300  0.058920 0.930333 -0.000745
         -0.066395 -0.076400 0.567860 -0.031570 0.266130  0.061280 0.937713 -0.000877
         -0.067831 -0.070173 0.534930 -0.034700 0.301060  0.063370 0.941848 -0.001212
         -0.067150 -0.063754 0.500670 -0.037820 0.335000  0.064970 0.950351 -0.002146
         -0.066302 -0.060865 0.466730 -0.040350 0.366270  0.066030 0.960400 -0.003718
         -0.063478 -0.055113 0.432470 -0.042250 0.401530  0.066770 0.971630 -0.005932
         -0.059726 -0.049766 0.400200 -0.043450 0.434460  0.067000 0.977522 -0.007257
         -0.053690 -0.043024 0.368260 -0.043910 0.468400  0.066730 0.989221 -0.010143
         -0.046491 -0.036434 0.333670 -0.043750 0.499000  0.066030 0.995344 -0.011797
         -0.034454 -0.026043 0.300070 -0.042880 0.533270  0.064740 1.001430 -0.013554
         -0.019737 -0.014303 0.267800 -0.041420 0.567860  0.062840 1.020150 -0.019516
         -0.006940 -0.004739 0.232870 -0.039160 0.599800  0.060480 1.034272 -0.024486
         0.008486  0.006239 0.202590 -0.036530 0.635400  0.057520 1.063809 -0.035726
         0.013998  0.010012 0.167000 -0.032930 0.668660  0.053560 1.097886 -0.050076
         0.019714  0.013891 0.136730 -0.029770 0.701260  0.049500 1.125638 -0.063296
         0.027025  0.018988 0.101460 -0.026280 0.734530  0.044940 1.155608 -0.079195
         0         0        0.084170 -0.024550 0.765140  0.040450 1.183517 -0.094908
         0         0        0.071120 -0.023250 0.799070  0.035100 1.214740 -0.113210
         0         0        0.067860 -0.022820 0.833330  0.029440 0         0
         0         0        0.061240 -0.022290 0.867930  0.023520 0         0
         0         0        0.054720 -0.021660 0.899870  0.017900 0         0];

      slat = temp(:,1:2);
      wing = temp(:,3:6);
      flap = temp(:,7:8);

      slat = slat(slat(:,1)~=0,:);
      flap = flap(flap(:,1)~=0,:);
      wing = [wing(:,1:2); wing(:,3:4)];
      box  = [-0.75,-1; 2.25,-1; 2.25,1; -0.75,1];
      
      
      
      n1 = size(slat,1);
      n2 = size(wing,1);
      n3 = size(flap,1);
      c1 = [(1:n1-1)',(2:n1)'; n1,1];
      c2 = [(1:n2-1)',(2:n2)'; n2,1];
      c3 = [(1:n3-1)',(2:n3)'; n3,1];
      c4 = [1,2; 2,3; 3,4; 4,1];

      node = [slat; wing; flap; box];
      edge = [c1; c2+n1; c3+n1+n2; c4+n1+n2+n3];

      options.dhmax = 0.15;

      hdata.edgeh(:,1) = [unique(c1(:)); unique(c2(:))+n1; unique(c3(:))+n1+n2];
      hdata.edgeh(:,2) = 0.001;
      hdata.fun = @const_h;
      hdata.args = {-.15,2.25,-.5,.5,0.025};
      [p,t] = mesh2d(node,edge,hdata,options);
      save 3foil p t

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 7

      node  = [0,0; 4,0; 4,1; 2,1; 2,2; 4,2; 4,3; 0,3];

      hdata.hmax = 0.05;

      [p,t] = mesh2d(node,[],hdata);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 8

      node = [-2,1; 1,1; 1,0; 20,0; 20,2; -2,2];

      hdata.hmax = 0.1;
      hdata.fun = @const_h;
      hdata.args = {-1,5,0,2,0.05};

      options.dhmax = 0.1;

      [p,t] = mesh2d(node,[],hdata,options);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 9

      wing = [
         1.00003  0.00126
         0.99730  0.00170
         0.98914  0.00302
         0.97563  0.00518
         0.95693  0.00812
         0.93324  0.01176
         0.90482  0.01602
         0.87197  0.02079
         0.83506  0.02597
         0.79449  0.03145
         0.75070  0.03712
         0.70417  0.04285
         0.65541  0.04854
         0.60496  0.05405
         0.55335  0.05924
         0.50117  0.06397
         0.44897  0.06811
         0.39733  0.07150
         0.34681  0.07402
         0.29796  0.07554
         0.25131  0.07597
         0.20738  0.07524
         0.16604  0.07320
         0.12732  0.06915
         0.09230  0.06265
         0.06203  0.05382
         0.03730  0.04324
         0.01865  0.03176
         0.00628  0.02030
         0.00015  0.00956
         0.00000  0.00000
         0.00533 -0.00792
         0.01557 -0.01401
         0.03029 -0.01870
         0.04915 -0.02248
         0.07195 -0.02586
         0.09868 -0.02922
         0.12954 -0.03282
         0.16483 -0.03660
         0.20483 -0.04016
         0.24869 -0.04283
         0.29531 -0.04446
         0.34418 -0.04510
         0.39476 -0.04482
         0.44650 -0.04371
         0.49883 -0.04188
         0.55117 -0.03945
         0.60296 -0.03655
         0.65360 -0.03327
         0.70257 -0.02975
         0.74930 -0.02607
         0.79330 -0.02235
         0.83407 -0.01866
         0.87118 -0.01512
         0.90420 -0.01180
         0.93279 -0.00880
         0.95661 -0.00621
         0.97543 -0.00410
         0.98901 -0.00254
         0.99722 -0.00158
         0.99997 -0.00126];

      wing = rotate(wing,45.0);
      wall = [-1,-2; 2,-2; 2,2; -1,2];

      nwing = size(wing,1);
      nwall = size(wall,1);

      cwing = [(1:nwing-1)', (2:nwing)'; nwing, 1];
      cwall = [(1:nwall-1)', (2:nwall)'; nwall, 1];

      edge = [cwing; cwall+nwing];
      node = [wing; wall];

      hdata.edgeh(:,1) = (1:size(cwing,1))';
      hdata.edgeh(:,2) = 0.0025;

      options.dhmax = 0.1;

      [p,t] = mesh2d(node,edge,hdata,options);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 10

      % DEFINE GEOMETRY
      H   = 1;                % Height
      L   = 5;                % Length
      n   = 4;                % Cycles
      k   = 2*pi*n*H/L;
      dx  = 0.05;             % Streamwise spatial increment
      x   = 0:dx:L;
      num = length(x);

      % Wavy channel walls
      ytop = 1+0.1*cos(k*x);

      node = [0,0; L,0; x(num:-1:1)' ytop(num:-1:1)'];

      hdata.hmax = 0.05;

      [p,t] = mesh2d(node,[],hdata);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 11

      % A mesh for an industrial process from an engineer in Germany (for a CFD
      % study invloving fluid transport through a tray of glass beads...) which
      % caused eariler versions of mesh2d big problems...


      Radius  = 2;                    % Radius des Kreises
      xNumber = 15;                    % Anzahl in x-Richtung
      yNumber = 8;                    % Anzahl in y-Richtung
      Anzahl  = xNumber*yNumber;      % Anzahl der Kreise
      Abstand = 2*Radius+1;           % Abstand der Kreise
      dtheta  = 0.1;                 % Schrittweite der Winkelaul"osung

      theta = 0:pi*dtheta:(2*pi-pi*dtheta);

      x = zeros(Anzahl,length(theta));
      y = x;
      Radii = 1.5+rand(Anzahl,1);

      % setzen der Kreise
      for i = 1:Anzahl
         x(i,:) = Radii(i)*cos(theta);
         y(i,:) = Radii(i)*sin(theta);
         dy = floor((i-1)/xNumber);
         dx = (i-1)-dy*xNumber;
         node((i-1)*20+1:i*20,1) = [x(i,:)'+dx*Abstand];
         node((i-1)*20+1:i*20,2)  = [y(i,:)'-dy*Abstand];
      end

      % setzen des begrenzenden Rechtecks
      node = [    node(1:Anzahl*20,1) node(1:Anzahl*20,2)
         -3                  -((yNumber-1)*5+3)
         ((xNumber-1)*5+3)    -((yNumber-1)*5+3)
         ((xNumber-1)*5+3)    3
         -3                  3
         ];


      n = size(node,1)-4;

      cnect = [(1:n-1)' (2:n)'
         n        1
         n+1      n+2
         n+2      n+3
         n+3      n+4
         n+4      n+1
         ];

      for i = 1:Anzahl
         cnect(i*20,2) = 1+(i-1)*20;
      end

      options.dhmax = 0.3;

      [p,t] = mesh2d(node,cnect,[],options);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 12
        p1 = [ 0 0
             .5 0
            1 0
            1 1
            0 1
            0 0];
        % Put nozzle case here....
            
        
        % data defining the boundary...
       

      p1 = p1(1:end-1,:);

      mesh2d(p1);

      
    case 13
        
        
%   Example of 2-D wing.. airfoil of length 2 m +/-
iaf.designation='4215';
iaf.designation='0012';
iaf.designation='0008';
iaf.n=64;
iaf.HalfCosineSpacing=1;
iaf.wantFile=1;
iaf.datFilePath='./'; % Current folder
iaf.is_finiteTE=0;

af = naca4gen(iaf);

clear wing
off = iaf.n / 8;
wing(:,1) = af.x(off+1:end-off);
wing(:,2) = af.z(off+1:end-off);

     wing = rotate(wing,1.5);

     wall = [-3,-2; 3,-2; 3,2; -3,2];
      %wall = [-3,-15;30,-15;30,20;-3,20];

      nwing = size(wing,1);
      nwall = size(wall,1);

      cwing = [(1:nwing-1)', (2:nwing)'; nwing, 1];
      cwall = [(1:nwall-1)', (2:nwall)'; nwall, 1];

      edge = [cwing; cwall+nwing];
      node = [wing; wall];

      hdata.edgeh(:,1) = (1:size(cwing,1))';
      hdata.edgeh(:,2) = 0.1;
        
      options.dhmax = 0.1;
      hdata.fun = @const_h;
      hdata.args = {-.2,3,-.25,.3,0.1};

      [p,t] = mesh2d(node,edge,hdata,options);

      %save naca p t
        
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = move(p,xm,ym)

% Move a node set p by [xm,ym]

n = size(p,1);
p = p + [xm*ones(n,1), ym*ones(n,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = rotate(p,A)

% Rotate a node set p by A degrees.

A = A*pi/180;
T = [ cos(A), sin(A)
   -sin(A), cos(A)];
p = (T*p')';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = const_h(x,y,x1,x2,y1,y2,h0)

% User defined size function specifying a constant size of h0 within the
% rectangle bounded by [x1,y1] & [x2,y2].

h = inf*ones(size(x,1),1);

in = (x>=x1)&(x<=x2)&(y>=y1)&(y<=y2);
h(in) = h0;


