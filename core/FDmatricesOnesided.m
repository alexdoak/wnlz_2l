function [Dx,Dxx] = FDmatrices_onesided(N,dx)
% ----------------------------------------------------------------------- %
% First derivative: centered difference, periodic
% ----------------------------------------------------------------------- %

e = ones(N,1);

Dx = spdiags([-e,e],[-1,1],N,N)/(2*dx);

% One-sided differences at x_1 and x_N
Dx(1,:) = 0;
Dx(1,1) = -3/(2*dx);
Dx(1,2) = +4/(2*dx);
Dx(1,3) = -1/(2*dx);


Dx(N,:) = 0;
Dx(N,N) = 3/(2*dx);
Dx(N,N-1) = -4/(2*dx);
Dx(N,N-2) = 1/(2*dx);

% ----------------------------------------------------------------------- %
% Second derivative: centered difference, periodic
% ----------------------------------------------------------------------- %

Dxx = spdiags([e,-2*e,e],[-1,0,1],N,N)/(dx^2);

% One-sided differences at x_1 and x_N
Dxx(1,:) = 0;
Dxx(1,1) = +2/(2*dx);
Dxx(1,2) = -5/(2*dx);
Dxx(1,3) = +4/(2*dx);
Dxx(1,4) = -1/(2*dx);

Dxx(N,:) = 0;
Dxx(N,N) = 2/(dx^2);
Dxx(N,N-1) = -5/(dx^2);
Dxx(N,N-2) = +4/(dx^2);
Dxx(N,N-3) = -1/(dx^2);