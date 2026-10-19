%function [w,d,x,dx,dt,kx,g,c,c0,Nx,H2,H,opts] = twolIG(H1,H2,g,r1,r2,x,Lx,Nx,dt,opts)
function [w,d,x,dx,dt,kx,g,c,c0,cWave,Nx,H1,H2,H,r1,r2,opts,optsSteady] = twolIG(x,Lx,Nx,dt,opts)
% ---------------------------------------------------------------------- %
% Load First Solution
% ---------------------------------------------------------------------- %
load(opts.WaveFileName{1});
if ~isempty(opts.WaveFileIndex{1})
    M = opts.WaveFileIndex{1};
else
    M=numel(li.c);
end
zeta = li.y{M};
x_h=x(1:(numel(x)/2+1))
effective_wavelength = 2*abs(trapz(x_h,y)/max(abs(y)));

if opts.Bathymetry.on
    epsilonStep = effective_wavelength/(2*opts.Bathymetry.stepWidth)
    epsilonWave = max(abs(li.y{M}))/H2
    %epsilonStep = li.effective_wavelength(M)/(2*opts.Bathymetry.stepWidth)
end
if isfield(li,'r1')
    r1=li.r1(M);
    r2=li.r2(M);
end
if isfield(li,'H1')
    H1=li.H1(M);
    H2=li.H2(M);
end


opts.Boussinesq = optsSteady.Boussinesq;

% ---- Useful paramters
H=H1+H2;

if opts.Boussinesq
    c0 = (g*(r2-r1)*H1*H2/(H2+H1))^.5;
else
    c0 = (g*(r2-r1)*H1*H2/(r1*H2+r2*H1))^.5;
end
% Rescale speed
if opts.RescaleSpeed
    if isempty(opts.RescaleSpeedValue)
        c_scale = c0;
    else
        c_scale = opts.RescaleSpeedValue 
    end
g_rescale = g/c_scale^2;
cWave{1} = c/c_scale;
else
c_scale = 1;
g_rescale=g;
cWave{1}=c;
end

% ---- Reflect
%zeta = [zeta(end);zeta(end:-1:2);zeta(1:end)];
zeta = [zeta(1:end);zeta(end-1:-1:2)];
% ---- Extend domain for interpolation
x_solwave{1}=[-200000*H;x;200000*H]; zeta_solwave{1}=[0;zeta;0];

% ---------------------------------------------------------------------- %
% Load additional waves
% ---------------------------------------------------------------------- %
numWaves = numel(opts.WaveFileName);
for p = 2:numWaves
    load(opts.WaveFileName{p});
    if ~isempty(opts.WaveFileIndex{p})
        zeta=li.y{opts.WaveFileIndex{p}};
        c = li.c(opts.WaveFileIndex{p});
    else
        zeta=li.y{end};
        c = li.c(end)
    end

    cWave{p} = c/c_scale;
    % ---- Reflect
    zeta = [zeta(1:end);zeta(end-1:-1:2)];
    % ---- Extend domain for interpolation
    x_solwave{p}=[-2000000*H;x;2000000*H]; zeta_solwave{p}=[0;zeta;0];
    
end


% ---------------------------------------------------------------------- %
% Mesh, Bathymetry
% ---------------------------------------------------------------------- %
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

% ---------------------------------------------------------------------- %
% Build initial guess from waves 
% ---------------------------------------------------------------------- %
for p=1:numWaves
    % ---- Interpolate solution on new grid
    zeta = interp1(x_solwave{p},zeta_solwave{p},x-opts.WaveLocation{p});
    
    % ---- Incorporate bathymetry
    H2=b;
    H=H2+H1;
    
    % ---- Recover (d,w) from zeta, H_i.
    h2 = H2+zeta;
    h1 = H1-zeta;
    q2 = cWave{p}*(h2-H2);
    q1 = cWave{p}*(h1-H1);
    if strcmp(opts.WaveDirection{p},"right")
        u2 = q2./h2;
        u1 = q1./h1;
    else
        u2 = -q2./h2;
        u1 = -q1./h1;
    end
    if p==1
        w = u2-u1;
        d = h2-h1;
    else
        w = w + u2-u1;
        d = d-d(1)+ h2-h1 - (H2-H2(1));
    end
end

% ---------------------------------------------------------------------- %
% Update legacy 
% ---------------------------------------------------------------------- %
if isfield(opts,'linDisp')
    opts = rmfield(opts,'linDisp');
end
if isfield(opts,'differentiation')
    opts = rmfield(opts,'differentiation');
end


% ---------------------------------------------------------------------- %
% Fix gravity rescale (was reloaded to old values when evaluating waves)
% ---------------------------------------------------------------------- %
g=g_rescale;


% ---------------------------------------------------------------------- %
% Perturb the solution?
% ---------------------------------------------------------------------- %
if opts.Perturb.on
    epsilon = opts.Perturb.epsilon
    k = opts.Perturb.k;
    W = opts.Perturb.Width;             % perturbation width
    xc = opts.Perturb.Location
    %[~,I] = max(d);
    %xc = x(I);     % crest
    
    pert = cos(k*(x-xc)) .* exp(-(x-xc).^2/(2*W^2));

    d=d+ epsilon*pert;
end

% ---------------------------------------------------------------------- %
% Filter the initial guess?
% ---------------------------------------------------------------------- %
if isfield(opts,"FilterInitialGuess") && opts.FilerInitialGuess.on
    kcut = opts.FilerInitialGuess.kcut;
    order =  opts.FilerInitialGuess.order;
    strength = opts.FilerInitialGuess.strength;
    
    % ---- Extend d and w 
    dext = [d;d(end-1:-1:2)];
    wext = [w;-w(end-1:-1:2)];

    % ---- FFT
    dhat = real(fft(dext));
    what = imag(fft(wext));

    % ---- Build filter
    Lext = 2*Lx; Next = numel(dext);
    kx = (2*pi/Lext)*[0:Next/2-1 0 -Next/2+1:-1]';
    kNyq = pi*Next/Lext;
    I = abs(kx)<kcut*kNyq;
    filter = exp(-strength*((abs(kx)-kcut*kNyq)./(kNyq-kcut*kNyq)).^order);
    filter(I)=1; 

    % ---- IFFT 
    d2 = real(ifft(filter.*dhat));
    w2 = imag(ifft(filter.*what));
    
    % ---- Take only physical domain
    d = d2(1:Nx);
    w = -w2(1:Nx);
end




%{
% ----------------------------------------------------------------------- %
% STABILITY ANALYSIS
% ----------------------------------------------------------------------- %
figure()
W = w(N/2+1);
D = (2*zeta(N/2+1)+H2-H1)
[omega2]=plotTwolayerStabilityAnalysis(H1,H2,r1,r2,g,W,gamma)
%}


%{
% ---------------------------------------------------------------------- %
% Second wave travelling left
% ---------------------------------------------------------------------- %
% Interpolate solution on new grid
zeta = interp1(x_solwave,zeta_solwave,x-10);

% Recover (d,w) from zeta, H_i.
h2 = H2+zeta;
h1 = H1-zeta;
q2 = c*(h2-H2);
q1 = c*(h1-H1);
u2 = -q2./h2;
u1 = -q1./h1;

w = w + u2-u1;
d = d-d(1)+ h2-h1;
%}


%w = (w-w(1)) + u2-u1;
%d = (d-d(1)) + h2-h1;
%}
