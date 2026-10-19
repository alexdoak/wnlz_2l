function [Dx,Dxx] = FDmatrices_periodic(N,dx)
% ----------------------------------------------------------------------- %
% First derivative: centered difference, periodic
% ----------------------------------------------------------------------- %

e = ones(N,1);

Dx = spdiags([-e,e],[-1,1],N,N)/(2*dx);

% periodic wrap-around
Dx(1,N) = -1/(2*dx);
Dx(N,1) =  1/(2*dx);


% ----------------------------------------------------------------------- %
% Second derivative: centered difference, periodic
% ----------------------------------------------------------------------- %

Dxx = spdiags([e,-2*e,e],[-1,0,1],N,N)/(dx^2);

% periodic wrap-around
Dxx(1,N) = 1/(dx^2);
Dxx(N,1) = 1/(dx^2);