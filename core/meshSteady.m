function [x,dx,dt,kx] = meshSteady(L,Nx,dt)

% ---- Spacing in x
dx = L/Nx;
x = -L/2:dx:L/2-dx; x=x';

% ---- Spectral mesh 
kx = (2*pi/L)*[0:Nx/2-1 0 -Nx/2+1:-1]'; 

% ---- Here, put CFL condition
if nargin<3
    dt = 0.01;
end





