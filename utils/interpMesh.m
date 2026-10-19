function [L,x,x_h,dx,dt,kx,N,y] = interpMesh(x_old,y,L_new,N_new)
N_old = numel(x_old);
[x,dx,dt,kx] = meshSteady(L_new,N_new);
y=interp1(x_old(1:N_old/2+1),y,x(1:N_new/2+1),'linear');
y(isnan(y))=0;
L=L_new; N=N_new;
x_h = x(1:N/2+1);