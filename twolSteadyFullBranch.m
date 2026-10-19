% ----------------------------------------------------------------------- %
% This code solve the steady Boussinesq system found in **. We spatially 
% discretise with FFTs
%
% N - Number of points in x
% L - Length of domain
% H - Undistrurbed depth
% h(x) - Fluid depth
% q(x) - Flux = bar(u)*h
% g - Gravity
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Problem setup
% ----------------------------------------------------------------------- %
%%% Clear workspace, plots
clear all; close all;

% ---- Add path to aux files (the folder this script lives in)
codeDir = fileparts(mfilename('fullpath'));
if isempty(codeDir)
    codeDir = pwd;
end
addpath(genpath(codeDir))


%%% Parameters 
N = 2^12;
L = 30;
g = 1;
H = 1;
H1 = 0.99;
H2 = H-H1;
r1 = 1;
r2 = 2.1;


%%% Physical mesh, spectral coefficients
[x,dx,dt,kx] = meshSteady(L,N);

%%% Half a wavelength for x, compute solution over half a domain
x_h = x(1:N/2+1);

%%% Initial guess for h
amp=0.01;

[c0,y,y_x,y_xx,y_xxx,c,amp] = KDV_twolayer(amp,H1,H2,r1,r2,g,x_h);
amp = y(1);
%%% Load solution if you want
Load = 0;
if Load == 1
    filename = 'H2=0.35_bous_branch';
    load(filename)
    %CNT=30;
    %y=y_list{30}; c=c_list(30); 
   
    [x,dx,dt,kx] = meshSteady(L,N);
    x_h = x(1:N/2+1);
end

%%% Boussinesq approx on densities?
optsSteady.Boussinesq = false;

%%% Choices for linearised dispersion term: qxxt, hxtt, uxxt, full
optsSteady.linDisp = 'hxtt';

%%% Choices for differentiation: spectral, findif
optsSteady.differentiation = 'findif';

%%% Interpolate onto new solution
%[x,y,dx,kx,N,x_h] = interpNewMesh(x_h,y,L,250)

%%% Rescale
%g = g/c^2;
%c= c/c;

% Interpo solution into new grid
%L_new = 50*2*2; N_new = 2000*2*2;
%[L,x,x_h,dx,dt,kx,N,y] = interpMesh(x,y,L_new,N_new); 

%%% Unknowns stored in a vector. We take only a wavelength
unvec = [y;c];
M = numel(unvec);
amp = y(end);
area=trapz(x_h,y); 


%%% Choose amplitudes to sweep through
amp_max = 0.5-H2;
if amp_max>0
    amp_start = 0.01;
    amp_space = 0.01;
    amp_end = amp_max - 0.02;
else
    amp_start = -0.01;
    amp_space = -0.01;
    amp_end = amp_max + 0.02;
end


li.y={}; li.amp=[]; li.area=[];
li.c=[]; li.effective_wavelength=[];

CNT=0;

%%% First sweep in amplitude
optsSteady.ParamFix = 'amp';

for amp=amp_start:amp_space:amp_end


    [unvec,y,c] = twolSteadySolve(H,H2,N,M,dx,kx,x,x_h,g,r1,r2,unvec,amp,area,optsSteady);

    
    amp = y(end);
    area=trapz(x_h,y); 

    CNT=CNT+1;
    li.y{CNT} = y;
    li.c(CNT) = c;
    li.amp(CNT) = y(end);
    li.area(CNT) =trapz(x_h,y);
    li.effective_wavelength(CNT) = 2*abs(trapz(x_h,y)/max(y));


    plot(x_h,y+H2,'k'); hold on; plot(x_h,x_h*0,'k','Linewidth',2) ; plot(x_h,x_h./x_h*H,'k','Linewidth',2)
    plot(-x_h,y+H2,'k'); hold on; plot(-x_h,x_h*0,'k','Linewidth',2) ; plot(-x_h,x_h./x_h*H,'k','Linewidth',2)
end




%%% Next sweep in area
optsSteady.ParamFix = 'area';

for area=linspace(area,area*5*2*2,10*2*2)
    
    [unvec,y,c] = twolSteadySolve(H,H2,N,M,dx,kx,x,x_h,g,r1,r2,unvec,amp,area,optsSteady);
    
    amp = y(end);
    area=trapz(x_h,y); 

    CNT=CNT+1;
    li.y{CNT} = y;
    li.c(CNT) = c;
    li.amp(CNT) = y(end);
    li.area(CNT) =trapz(x_h,y);
    li.effective_wavelength(CNT) = 2*abs(trapz(x_h,y)/max(abs(y)));


    plot(x_h,y+H2,'k'); hold on; plot(x_h,x_h*0,'k','Linewidth',2) ; plot(x_h,x_h./x_h*H,'k','Linewidth',2)
    plot(-x_h,y+H2,'k'); hold on; plot(-x_h,x_h*0,'k','Linewidth',2) ; plot(-x_h,x_h./x_h*H,'k','Linewidth',2)

end




save(['H2=' num2str(H2) '_bous_branch.mat'],'x','c','H1','H2','g','r1','r2','y','li','optsSteady')


clear J; clear yx; clear yxx; clear yxxx; clear unvec; clear unvec_store; clear Dx; clear Dxx; clear Dxx; clear Dxxx
