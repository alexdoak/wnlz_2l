function [x,dx,dt,kx] = mesh(L,Nx,dt)

% ---- Spacing in x
dx = L/(Nx-1);
x  = linspace(-L,0,Nx)';

% ---- Spectral mesh 
Next = 2*Nx-2;
Lext = 2*L;
kx = (2*pi/Lext)*[0:Next/2-1 0 -Next/2+1:-1]'; 
x_ext = -L:dx:L-dx;

% ---- TO DO: put CFL condition
if nargin<3
    dt = dx;
end




