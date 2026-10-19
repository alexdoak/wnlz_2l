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
%%% Clear workspace, plots
clear all; close all;


for stepHeight = [0.55]
close all;

%%% Add path to aux files
PATH = 'C:\Users\maaad\Documents\MATLAB\wnlz_2l';
addpath(genpath(PATH))

%%% Set up plot
axis tight manual % this ensures that getframe() returns a consistent size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig=figure(1);
set(fig,'Units','centimeters','Position',[10 8 20 12])

%%% Parameters
g = 1;
H = 1;
%Lx = 400;
Lx = 600;
dx = 0.05;
Nx = round(Lx/dx); 
t0 = 0;
dt=.05;
tprint = 5;
tfinal = 650;
%tfinal = 350;
%%% Boussinesq approximaition on densities? This will actually be set in IG
opts.Boussinesq = true;

%%% Rescale gravity to give c0=1? or alternatively choose rescale?
opts.RescaleSpeed = true;
opts.RescaleSpeedValue = [];  

%%% Physical mesh, spectral coefficients
if exist('dt')
    [x,dx,dt,kx] = mesh(Lx,Nx,dt);
else
    [x,dx,dt,kx] = mesh(Lx,Nx);
end
 
%%% Boundary condition? "wall", "periodic"
opts.BC = 'wall';

%%% Spatial derivatives: "spectral", "findif"
opts.FinDiff = true;


%%% Bathymetry
opts.Bathymetry.on = true;
opts.Bathymetry.stepHeight = stepHeight;
opts.Bathymetry.stepLocation = -Lx+125%300;%+300;%+125;
opts.Bathymetry.stepWidth = 4;%14; 
opts.Bathymetry.type = 'tanh';

PATH = 'C:\Users\maaad\Documents\MATLAB\wnlz_2l';
addpath(genpath(PATH))


%%% Choose the type of simulation;
opts.Simulation = 'Overtake';

opts.WaveOneFileName = 'H2=0.7_bous_branch';%'H2=0.65_branch.mat';
opts.WaveTwoFileName = []; %'H2=0.2_bous_branch';
opts.WaveOneFileIndex = 10;
opts.WaveTwoFileIndex = 20;
opts.WaveOneLocation = -Lx+30;
opts.WaveTwoLocation = -Lx+70;

[w,d,x,dx,dt,kx,g,c,c0,c1,c2,Nx,H1,H2,H,r1,r2,opts,optsSteady] = twolIG(x,Lx,Nx,dt,opts);
q = [d;w]; zeta=1/2*(d-H2+H1);


%%% Boundary condition? "wall", "periodic"
opts.BC = 'wall';

%%% Spatial derivatives: "spectral", "findif"
opts.FinDiff = true;

%%% Differentiation matrices
switch opts.BC
    case "periodic"
        [Dx,Dxx] = FDmatricesPeriodic(Nx,dx);
    case "wall"
        [Dx,Dxx] = FDmatricesOnesided(Nx,dx);
end

%%% Make a video live?
opts.VideoMaker = true;
filename = ['amp=' num2str(round(max(abs(zeta)),3),'%4.2f') '_A=' num2str(opts.Bathymetry.stepHeight,'%4.3f') '_W=' num2str(opts.Bathymetry.stepWidth)  '_H2=' num2str(H2(1)) '_fission.gif']


% ----------------------------------------------------------------------- %
% RK4 integration
% ----------------------------------------------------------------------- %
%%% Store data at print times
d_data = d;
w_data = w;
t_data = [0];

if opts.FinDiff
    F=@(q,t)twolDtFindif(q,Dx,Dxx,H1,H2,r1,r2,g,opts);
else
    F=@(q,t)twolDtSpectral(q,kx,H1,H2,r1,r2,g,opts);
end

tail_beta=[];
bulk_beta=[];

tail_d=[];
bulk_d=[];

tail_w=[];
bulk_w=[];


tail_betaw = [];
tail_d2w = [];

cnt=0;
opts.Boussinesq=false;
for t=t0+tprint:tprint:tfinal
    tic
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
   toc
   %%% Store new values
   d_data=[d_data,d];
   w_data=[w_data,w];
   t_data=[t_data,t];
   t
   cnt=cnt+1;
    


    if opts.VideoMaker
   
       zeta =1/2*(d-(H2-H1));
       u2 = (H-d)/(2*H).*w;
       u1 = (H+d)/(2*H).*w;
    
       figure(1)
       p1=plot(x,zeta+H2(1),'k'); hold on;
       title1=title(['$t=$' num2str(t) ', $\zeta(x)$']);
       ylim([-0.1,1.1])
       if opts.Bathymetry.on
           p2=plot(x,-H2+H2(1),'k','Linewidth',2);
           p3=plot(x,H2(1)-H2+H,'k','Linewidth',2);
       end

      
       % Capture the plot as an image 
       frame = getframe(fig); 
       im = frame2im(frame); 
       [imind,cm] = rgb2ind(im,256); 
       % Write to the GIF File 
       if cnt==1
           imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
       else 
           imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.075); 
       end
    
    
        delete(title1); delete(p1);
        if opts.Bathymetry.on
            delete(p2);
            delete(p3);
        end
    end   
    
    

   %{
   %%% Spectral analysis for instability checking
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


   

zeta_data =1/2*(d_data-(H2-H1));

zeta=zeta_data(:,1);


end


%save('overtake_polarity_switch_1.mat','x','d_data','w_data','t_data','zeta_data','opts','optsSteady','H1','H2','c','c0','c1','c2','g','r1','r2')

epsilonWave = max(abs(zeta))/H2(1)
effectiveWavelength = abs(trapz(x,zeta)/max(abs(zeta)))
epsilonStep = effectiveWavelength/(2*opts.Bathymetry.stepWidth)

save(['amp=' num2str(round(max(abs(zeta_data(:,1))),3),'%4.2f') '_A=' num2str(opts.Bathymetry.stepHeight,'%4.3f') '_W=' num2str(opts.Bathymetry.stepWidth)  '_H2=' num2str(H2(1)) '_fission.mat'], ...
    'x','d_data','w_data','t_data','zeta_data','opts','optsSteady','H1','H2', ...
    'c','c0','c1','c2','g','r1','r2','epsilonWave','epsilonStep','effectiveWavelength');

end