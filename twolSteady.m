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
% ---- Clear workspace, plots
clear all; close all;

% ---- Add path to aux files (the folder this script lives in)
codeDir = fileparts(mfilename('fullpath'));
if isempty(codeDir)
    codeDir = pwd;
end
addpath(genpath(codeDir))

% ---- Parameters 
N = 2^12;
L = 30;
g = 1;
H = 1;
H1 =  10/11;
H2 = H-H1;
r1 = 1;
r2 = 5;


% ---- Physical mesh, spectral coefficients
[x,dx,dt,kx] = meshSteady(L,N);

% ---- Half a wavelength for x, compute solution over half a domain
x_h = x(1:N/2+1);

% ---- Initial guess for h
amp=0.01;

% ---- Boussinesq approx on densities?
optsSteady.Boussinesq = false;

[c0,y,y_x,y_xx,y_xxx,c,amp] = KDV_twolayer(amp,H1,H2,r1,r2,g,x_h,optsSteady.Boussinesq);
amp = y(end);
% ---- Load solution if you want
Load = true;
if Load 
    filename = 'CompareSpinsSteady_SOL1';
    load(filename)
   % CNT=1; N=numel(x); L=-x(1)*2;
   % y=li.y{CNT}; c=li.c(CNT);
    L = -x(1)*2; N=numel(x);
    [x,dx,dt,kx] = meshSteady(L,N);
    x_h = x(1:N/2+1);
end


% ---- Choices for linearised dispersion term: qxxt, hxtt, uxxt, full
optsSteady.linDisp = 'hxtt';

% ---- Choices for differentiation: spectral, findif
optsSteady.differentiation = 'findif';

% ---- Interpolate onto new solution
%[x,y,dx,kx,N,x_h] = interpNewMesh(x_h,y,L,250)

% ---- Rescale
%g = g/c^2;
%c= c/c;

% ---- Interp solution into new grid
%L_new = 50; N_new = 5000;
%[L,x,x_h,dx,dt,kx,N,y] = interpMesh(x,y,L_new,N_new);

% ---- Unknowns stored in a vector. We take only a wavelength
unvec = [y;c];
M = numel(unvec);
amp = y(end);
area=trapz(x_h,y); 

% ---- What parameter to param_fix?
optsSteady.ParamFix = 'amp';
% ---- Solve
li.y={}; li.amp=[]; li.area=[]; li.c=[]; li.effective_wavelength=[];li.r1=[];li.r2=[];

CNT=numel(li.c);

%r1 = 1;
%for r2=1.05:0.05:5
%for amp=0.01:0.01:0.55
%for amp=amp


for amp=[-0.1,-0.11,-0.12,-0.13]
%% 
    %r2=5:-0.1:1.1 %area=linspace(area,area*20,40)
    H1=H-H2;
    % ---- Differentiation matrices 
    [Dx,Dxx,Dxxx] = FDmatricesNeumann(N/2+1,dx);

    % ---- Pre-allocated memory
    J = zeros(M,M); 
    Hi = ones(1,M);
    Hi_D=Hi;
    while max(abs(Hi))>10^(-10)
        
        % Unpack unvec
        y = unvec(1:N/2+1);
        c = unvec(N/2+2);
       
        % Differentiate
        [yx,yxx,yxxx] =  SpatialDifferentiation(y,N,dx,Dx,Dxx,Dxxx,kx,optsSteady.differentiation);
        
        % Equations
        [Hi] = twolSteadyEqnY(y,yx,yxx,yxxx,H,c,g,r1,r2,H1,H2,optsSteady);
        

        switch optsSteady.ParamFix
            case 'amp'
            Hi(M) = y(end)-amp;
            case 'area'
            Hi(M) = trapz(x_h,y)-area;
        end
        % --------------------------------------------------------------- %
        % Jacobian, computed numerically.
        % --------------------------------------------------------------- %
        delta =10^(-11);
        unvec_store=unvec;
        for p=1:M
            unvec=unvec_store;
            unvec(p)=unvec(p)+delta;

            % Unpack unvec
            y = unvec(1:N/2+1);
            c = unvec(N/2+2);
           
            % Differentiate
            [yx,yxx,yxxx] =  SpatialDifferentiation(y,N,dx,Dx,Dxx,Dxxx,kx,optsSteady.differentiation);
            
            % Equations
            [Hi_D] = twolSteadyEqnY(y,yx,yxx,yxxx,H,c,g,r1,r2,H1,H2,optsSteady);
        
            
    
            switch optsSteady.ParamFix
                case 'amp'
                Hi_D(M) = y(end)-amp;
                case 'area'
                Hi_D(M) = trapz(x_h,y)-area;
            end
            
            J(:,p) = (Hi_D-Hi)/delta;
        end
    
    unvec=unvec_store;
    Residuals=max(abs(Hi));
    RL2 = sum(Hi.^2).^(1/2)
    unvec = unvec + J\(-Hi);
    end


    
    amp = y(1);
    area=trapz(x_h,y); 


    plot(x_h,y+H2,'k'); hold on; plot(x_h,x_h*0,'k','Linewidth',2) ; plot(x_h,x_h./x_h*H,'k','Linewidth',2)
    plot(-x_h,y+H2,'k'); hold on; plot(-x_h,x_h*0,'k','Linewidth',2) ; plot(-x_h,x_h./x_h*H,'k','Linewidth',2)
%    plot(x_h,y1_KDV,'k'); hold on; plot(x_h,y1,'r')
    CNT=CNT+1;
    li.y{CNT} = y;
    li.c(CNT) = c;
    li.amp(CNT) = amp;
    li.area(CNT) =trapz(x_h,y);
    li.effective_wavelength(CNT) = 2*abs(trapz(x_h,y)/max(abs(y)))
    li.r2(CNT) = r2;
    li.r1(CNT) = r1;
end

clear J; clear yx; clear yxx; clear yxxx; clear unvec; clear unvec_store; clear Dx; clear Dxx; clear Dxx; clear Dxxx

%save('P=10_maxAmp_SweepR','x','N','dx','H','c','H1','H2','g','r1','r2','y','li','optsSteady');




