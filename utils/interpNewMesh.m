function [x,y,dx,kx,N,x_h] = interpNewMesh(x_h_old,y_old,L,N)

[x,dx,dt,kx] = meshSteady(L,N);
x_h = x(1:N/2+1);
y = interp1(x_h_old,y_old,x_h,'PCHIP');






