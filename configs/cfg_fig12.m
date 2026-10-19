function p = cfg_fig12()
% Configuration for Figure 12.
%
% Any field left out here will be asked for interactively.
% Derived quantities (H1, c0, criticality, max amplitude) are computed by
% the setup script and should NOT be set here.
%
% Name configuration files with the prefix cfg_ so they never clash with
% function names in the main code.

% ---- Background state
p.H          = 1;          % total depth
p.Lx         = 100;        % domain length, x in [-Lx, 0]
p.r1         = 0.99;          % upper layer density
p.r2         = 1;       % lower layer density
p.Boussinesq = true;       % Boussinesq approximation on densities
p.H2         = 0.35;        % lower layer depth (criticality is at H/2 here)

% ---- Gravity: the script's default, chosen so that the long-wave speed c0 = H
H1    = p.H - p.H2;
if p.Boussinesq
    c0sq  = (p.r2 - p.r1)*H1*p.H2/p.H;    % c0^2 with g = 1 (Boussinesq)
else
    c0sq =  ((p.r2-p.r1)*H1*p.H2/(p.r1*p.H2+p.r2*H1)); % c0^2 with g = 1 (non-Boussinesq)
end
p.g   = p.H^2/c0sq;

% ---- Mesh and time stepping
p.defaultMesh = false;
% To override the defaults instead, set defaultMesh = false and give:
p.dx          = 0.02;
p.dt          = 0.01;
p.tprintSteps = 10;     % time steps between stored outputs
p.tfinal      = 80;


% ---- Solitary waves (one struct per wave)
% Amplitude : same sign as the maximum amplitude, and no larger in size
% Location  : in [-Lx, 0]
% Direction : "right" or "left"
p.Waves(1) = struct('Amplitude', 0.1177, 'Location', -p.Lx+15, 'Direction', "right");

% ---- Bathymetry (Boussinesq case only)
p.Bathymetry.on           = true;
p.Bathymetry.stepHeight   = 0.2;
p.Bathymetry.stepLocation = -65;
p.Bathymetry.stepWidth    = 6;
p.Bathymetry.type =  'tanh';

% ---- Output
p.VideoMaker = true;

% ---- Expert options (omit for defaults)
p.Perturb.on      = false;
p.Perturb.epsilon = 1e-6;
p.Perturb.k       = 40;
p.Perturb.Width   = 4;
p.Perturb.Location = -175;


end

