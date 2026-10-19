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
filename = 'Unstable_Maxamp_R=2.2_ep=6.gif';

% ---- Parameters
g = 1;
H = 1;
Lx=200;
dx = 0.01;%0.0073242187500;%0.007*2; 0.025;
Nx = round(Lx/dx)+1; 
t0 = 0;
dt= 0.01;
tprint = 25*dt;
tfinal = 150;

% ----  Boussinesq approximaition on densities? This will actually be set 
% ----  in IG if loading waves
opts.Boussinesq = false;

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

% ---- Add a perturbation to test instability?
opts.Perturb.on = true;
opts.Perturb.epsilon = 10^(-6);
opts.Perturb.k = 40;
opts.Perturb.Width = 4;  
opts.Perturb.Location = -175;

% ---- Filter initial guess ? 
opts.FilerInitialGuess.on = false;
opts.FilerInitialGuess.kcut = 0.2;
opts.FilerInitialGuess.order = 2;
opts.FilerInitialGuess.strength = 20;
    
% ---- Choose files to load
opts.WaveFileName =  {'tabletop_stable'} %{"maxamp_2","smallamp",}%{"H2=0.3_bous_branch"};%,"H2=0.3_bous_branch","H2=0.3_bous_branch"};
opts.WaveFileIndex = {[]};%,5,2};
opts.WaveLocation = {-Lx+25};%,-Lx+28,-Lx+40};
opts.WaveDirection = {"right","left"};%,'r','l'};

% ---- Build IG
[w,d,x,dx,dt,kx,g,c,c0,cWave,Nx,H1,H2,H,r1,r2,opts,optsSteady] = twolIG(x,Lx,Nx,dt,opts);
q = [d;w];

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
    v = VideoWriter(filename);
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

    

%save('overtake_polarity_switch_1.mat','x','d_data','w_data','t_data','zeta_data','opts','optsSteady','H1','H2','c','c0','c1','c2','g','r1','r2')

%save(filename,'x','d_data','w_data','t_data','zeta_data','opts','optsSteady','H1','H2','c','c0','c1','c2','g','r1','r2','dt')
    

   %{
   % ---- Spectral analysis for instability checking
    beta = H.^2-d.^2;
    
   
    betahat = fft(extensionEven(beta));
    betahat = fftshift(abs(betahat));   
    tail_beta = [tail_beta,max(betahat(round(Nx/2)-20:round(Nx/2)+20))];
    bulk_beta = [bulk_beta,max(abs(betahat))];

    dhat = fft(extensionEven(d));
    dhat = fftshift(abs(dhat));
    tail_d = [tail_d,max(dhat(round(Nx/2)-20:round(Nx/2)+20))];
    bulk_d = [bulk_d,max(abs(dhat))];


    what = fft(extensionOdd(w));
    what = fftshift(abs(what));
    tail_w = [tail_w,max(what(round(Nx/2)-20:round(Nx/2)+20))];
    bulk_w = [bulk_w,max(abs(what))];

    betawhat = fft(extensionOdd(beta.*w));
    betawhat = fftshift(abs(betawhat));   
    tail_betaw = [tail_betaw,max(betawhat(round(Nx/2)-20:round(Nx/2)+20))];

    d2what = fft(extensionOdd(d.^2.*w));
    d2what = fftshift(abs(d2what));   
    tail_d2w = [tail_d2w,max(d2what(round(Nx/2)-20:round(Nx/2)+20))];
   %}


  
