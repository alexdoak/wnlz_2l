% ----------------------------------------------------------------------- %
% This code integrates the Boussinesq system found in **. We spatially 
% discretise with second-order differences  and integrate in time using
% forward Euler
%
% N - Number of points in x
% L - Length of domain
% H - Undistrurbed depth
% h(x) - Fluid depth
% q(x) - Flux = bar(u)*h
% g - Gravity
% t0 - Initial time
% tprint - Time-spacing for data storage
% tfinal - Final time
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Problem setup
% ----------------------------------------------------------------------- %
% ---- Clear workspace, plots
clear all; close all;

% ---- Add path to aux files (the folder this script lives in)
codeDir = fileparts(mfilename('fullpath'));
if isempty(codeDir)
    codeDir = pwd;
end
addpath(genpath(codeDir))


% ---- Set up plot
axis tight manual % this ensures that getframe() returns a consistent size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig=figure(1);
set(fig,'Units','centimeters','Position',[10 8 20 12])
videofilename = 'lockreleasesim.avi'

% ---- Parameters
g = 1;
H = 1;
Lx= 30;
dx = 0.01;%0.0073242187500;%0.007*2; 0.025;
Nx = round(Lx/dx)+1; 
t0 = 0;
dt= 0.01;
tprint = 10*dt;
tfinal = 27;%$2000*tprint;

% ----  Boussinesq approximaition on densities? This will actually be set 
% ----  in IG if loading waves
opts.Boussinesq = true;

% ---- Rescale gravity to give c0=1? or alternatively choose rescale?
opts.RescaleSpeed = true;
opts.RescaleSpeedValue = [];  

% ---- Physical mesh, spectral coefficients
if exist('dt')
    [x,dx,dt,kx] = mesh(Lx,Nx,dt);
else
    [x,dx,dt,kx] = mesh(Lx,Nx);
end
 

% ---- Bathymetry
opts.Bathymetry.on = false;
opts.Bathymetry.stepHeight = 0.1;
opts.Bathymetry.stepLocation = -Lx+150;
opts.Bathymetry.stepWidth = 8; 
opts.Bathymetry.type = 'tanh';

PATH = 'C:\Users\maaad\Documents\MATLAB\wnlz_2l';
addpath(genpath(PATH))

% ----------------------------------------------------------------------- %
% General Initial Guess
% ----------------------------------------------------------------------- %
% --- Parameter values
r1 = 1;
r2 = 1.01;
H = 1;
H1 = .2;
H2 = .8;

% --- Choose g to make c0=1. Start by evaluating c0 assuming g=1, then set
% --- g=1/c0^2, since c0~g^.5
if opts.Boussinesq
    c0 = ((r2-r1)*H1*H2/(H2+H1))^.5;
else
    c0 = ((r2-r1)*H1*H2/(r1*H2+r2*H1))^.5;
end

g = g/c0^2;
c0 = 1;
% ---- Physical mesh, spectral coefficients
[x,dx,dt,kx] = mesh(Lx,Nx,dt);

% ---- Bathymetry
if opts.Bathymetry.on
    switch opts.Bathymetry.type  
        case 'tanh'
            % ---- Tanh step
            stepHeight = opts.Bathymetry.stepHeight;
            stepLocation = opts.Bathymetry.stepLocation;
            stepWidth = opts.Bathymetry.stepWidth;
            b = zeros(Nx,1) + H2 - stepHeight/2*tanh((x-stepLocation)/stepWidth) -stepHeight/2;
        
        case 'CLC'
            % ---- Constant-Linear-Constant
            stepHeight = opts.Bathymetry.stepHeight;
            stepLocation = opts.Bathymetry.stepLocation;
            stepWidth = opts.Bathymetry.stepWidth;
            b =  H2 - stepHeight/2 - stepHeight/(2*stepWidth)*(x-stepLocation);
            b(b>H2) = H2;
            b(b<H2-stepHeight) = H2-stepHeight;

        case 'Gaussian'
            % ---- Gaussian 

    end

else 
    b=H2;
end

% ---- Initial guess
optsSteady = {};
cWave = [];

dval = 0.25;
wval = 0;

gateInterfaceWidth = .05;
gateWidth = 1.3;
d = H2-H1 + dval*(tanh((x-(-Lx+gateWidth))/gateInterfaceWidth)-1);
w = 0*x;



q = [d;w];


% ----------------------------------------------------------------------- %
% Final details
% ----------------------------------------------------------------------- %

% ---- Boundary condition? "wall", "periodic"
opts.BC = 'wall';

% ---- Spatial derivatives: true - finite differences, false - spectral
opts.FinDiff = true;

% ---- Differentiation matrices
switch opts.BC
    case "periodic"
        [Dx,Dxx] = FDmatricesPeriodic(Nx,dx);
    case "wall"
        [Dx,Dxx] = FDmatricesOnesided(Nx,dx);
end

% ---- Make a video live?
opts.VideoMaker = true;
if opts.VideoMaker
    % ================================================================
    % Video initialisation
    % ================================================================
    % ---- Create figure
    fig = figure(1);
    set(fig,'Units','centimeters','Position',[10 8 20 12])
    
    % ---- Create VideoWriter
    v = VideoWriter(videofilename);
    v.FrameRate = 1/0.075;       
    open(v);
    
    % ---- Plot initial condition
    zeta = 1/2*(d-(H2-H1));
    u2 = (H-d)/(2*H).*w;
    u1 = (H+d)/(2*H).*w;
    
    % ---- Create zeta plot
    p1 = plot(x,zeta+H2(1),'k');
    hold on;
    
    % ---- Create title
    title1 = title(['$t=0$']);
    
    % ---- Set limits
    ylim([-0.1*H(1),1.1*H(1)])
    
    % ---- Plot bathymetry
    if opts.Bathymetry.on
        p2 = plot(x,-H2+H2(1),'k','Linewidth',2);
        p3 = plot(x,H2(1)-H2+H,'k','Linewidth',2);
    else    
        p2 = plot([x(1),x(end)],[0,0],'k','Linewidth',2);
        p3 = plot([x(1),x(end)],[H(1),H(1)],'k','Linewidth',2);
    end
    
    % ---- Capture first frame
    drawnow;
    frame = getframe(fig);
    writeVideo(v,frame);    
end


% ----------------------------------------------------------------------- %
% RK4 integration
% ----------------------------------------------------------------------- %
% ---- Store data at print times
d_data = d;
w_data = w;
t_data = [0];

if opts.FinDiff
    F=@(q,t)twolDtFindif(q,Dx,Dxx,H1,H2,r1,r2,g,opts);
else
    F=@(q,t)twolDtSpectral(q,kx,H1,H2,r1,r2,g,Lx,opts);
end

cnt=0;
for t=t0+tprint:tprint:tfinal
    for tint=0:dt:tprint-dt
        [dq1] = F(q,t);
        [dq2] = F(q+dq1*dt/2,t+dt/2);
        [dq3] = F(q+dq2*dt/2,t+dt/2);
        [dq4] = F(q+dq3*dt,t+dt);
        
        dq = (1/6)*(dq1+2*dq2+2*dq3+dq4)*dt;

        q = q + dq;
        d = q(1:Nx);
        w = q(Nx+1:2*Nx);
    end
   
   % ---- Store new values
   d_data=[d_data,d];
   w_data=[w_data,w];
   t_data=[t_data,t];
   t;
   cnt=cnt+1;
   
    if opts.VideoMaker
   
        zeta =1/2*(d-(H2-H1));
      %  u2 = (H-d)/(2*H).*w;
      %  u1 = (H+d)/(2*H).*w;
    
        
        % ---- Update zeta plot
        zeta = 1/2*(d-(H2-H1));
        p1.YData = zeta + H2(1);
        
        % ---- Update title
        title1.String = ['$t=$' num2str(t)];
        
        % --- Render updated figure
        drawnow;
        
        % ---- Capture frame
        frame = getframe(fig);
        
        % ---- Write frame to video
        writeVideo(v,frame);
        

    end


zeta_data =1/2*(d_data-(H2-H1));
end

if opts.VideoMaker
    close(v);
end
