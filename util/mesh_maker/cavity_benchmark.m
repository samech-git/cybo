function cavity_benchmark(p,U,V,Re)

% Simple benchmarking function for the driven cavity flow
%
% Uses data from Erturk et. al: www.cavityflow.com
%
%   cavity_benchmark(p,U,V,Re)
%
% Needs:    p  - nodes
%           U  - U velocity at p
%           V  - V velocity at p
%           Re - 1,2,3, comparisons for Reynolds of 1000,2500,5000.
%
% Will plot comparisons of U,V along the xy cavity midlines. It also
% returns the average %error along these lines.
%
% Darren Engwirda - 2005-06.


% Benchmark U velocity profile along x = 0.5
yt = [
     0.00
     0.02
     0.04
     0.06
     0.08
     0.10
     0.12
     0.14
     0.16
     0.18
     0.20
     0.50 
     0.90
     0.91 
     0.92
     0.93
     0.94
     0.95
     0.96
     0.97
     0.98
     0.99
     1.00
     ];

xt = 0.5*ones(size(yt));

% Re = [1000, 2500, 5000]
ut = [
     0.0000,  0.0000,  0.0000
    -0.0757, -0.1517, -0.2223
    -0.1392, -0.2547, -0.3480
    -0.1951, -0.3372, -0.4272    
    -0.2472, -0.3979, -0.4419
    -0.2960, -0.4250, -0.4168
    -0.3381, -0.4200, -0.3876
    -0.3690, -0.3965, -0.3652
    -0.3854, -0.3688, -0.3467
    -0.3869, -0.3439, -0.3285
    -0.3756, -0.3228, -0.3100
    -0.0620, -0.0403, -0.0319
     0.3838,  0.4141,  0.4155
     0.3913,  0.4256,  0.4307
     0.3993,  0.4353,  0.4452
     0.4101,  0.4424,  0.4582
     0.4276,  0.4470,  0.4683
     0.4582,  0.4506,  0.4738
     0.5102,  0.4607,  0.4739
     0.5917,  0.4971,  0.4749
     0.7065,  0.5624,  0.5159
     0.8486,  0.7704,  0.6866
     1.0000,  1.0000,  1.0000
     ];
   
% Interpolate to find total error     
errU = abs(ut(:,Re) - griddata(p(:,1),p(:,2),U,xt,yt)); errU = sum(errU)/length(xt)   
     
% Interpolate to plot
yi = (0:0.02:1)';     
xi = 0.5*ones(size(yi));     
ui = griddata(p(:,1),p(:,2),U,xi,yi);  
    
figure, grid on, axis square, hold on 

plot(ut(:,Re),yt,'bo') 
plot(ui,yi,'r',ui,yi,'r.')

title('U velocity profile')
xlabel('U velocity')
ylabel('Y')
legend('Erturk et. al.','Present study')
   
   
% Benchmark V velocity along y = 0.5
xt = [
     0.000
     0.015
     0.030
     0.045
     0.060
     0.075
     0.090
     0.105
     0.120
     0.135
     0.150
     0.500
     0.850
     0.865
     0.880
     0.895
     0.910
     0.925
     0.940
     0.955
     0.970
     0.985
     1.000
     ];
 
yt = 0.5*ones(size(xt));

vt = [
     0.0000,  0.0000,  0.0000
     0.1019,  0.1607,  0.2160
     0.1792,  0.2633,  0.3263
     0.2349,  0.3238,  0.3868
     0.2746,  0.3649,  0.4258
     0.3041,  0.3950,  0.4426
     0.3273,  0.4142,  0.4403
     0.3460,  0.4217,  0.4260
     0.3605,  0.4187,  0.4070
     0.3705,  0.4078,  0.3878
     0.3756,  0.3918,  0.3699
     0.0258,  0.0160,  0.0117
    -0.4028, -0.3671, -0.3624  
    -0.4407, -0.3843, -0.3806
    -0.4803, -0.4042, -0.3982
    -0.5132, -0.4321, -0.4147
    -0.5263, -0.4741, -0.4318 
    -0.5052, -0.5268, -0.4595
    -0.4417, -0.5603, -0.5139
    -0.3400, -0.5192, -0.5700
    -0.2173, -0.3725, -0.5019
    -0.0973, -0.1675, -0.2441
     0.0000,  0.0000,  0.0000
    ];
     
     
% Interpolate to find total error     
errV = abs(vt(:,Re) - griddata(p(:,1),p(:,2),V,xt,yt)); errV = sum(errV)/length(yt)   
     
% Interpolate to plot
xi = (0:0.02:1)';     
yi = 0.5*ones(size(xi));     
vi = griddata(p(:,1),p(:,2),V,xi,yi);  
    
figure, grid on, axis square, hold on 

plot(xt,vt(:,Re),'bo') 
plot(xi,vi,'r',xi,vi,'r.')     
     
title('V velocity profile')
xlabel('X')
ylabel('V velocity')
legend('Erturk et. al.','Present study')
   