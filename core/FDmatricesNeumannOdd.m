function [Dx,Dxx,Dxxx] = FDmatrices_neumann_odd(N,dx)
% Sparse derivative matrices for ODD functions
%
% Ghost-point conditions:
%
%   f_0     = -f_2
%   f_{N+1} = -f_{N-1}
%
% corresponding to
%
%   f(-x) = -f(x)
%
% Inputs:
%   N  - number of grid points
%   dx - grid spacing
%
% Outputs:
%   Dx   - first derivative matrix
%   Dxx  - second derivative matrix
%   Dxxx - third derivative matrix

e = ones(N,1);

%--------------------------------------------------
% First derivative
%--------------------------------------------------

Dx = spdiags([-e e],[-1 1],N,N)/(2*dx);

% Left boundary:
% (f2-f0)/(2dx)
% = f2/dx
Dx(1,:) = 0;
Dx(1,2) = 1/dx;

% Right boundary:
% (f(N+1)-f(N-1))/(2dx)
% = -f(N-1)/dx
Dx(N,:) = 0;
Dx(N,N-1) = -1/dx;

%--------------------------------------------------
% Second derivative
%--------------------------------------------------

Dxx = spdiags([e -2*e e],[-1 0 1],N,N)/(dx^2);

% Left boundary:
% (f2-2f1+f0)/dx^2
% = -2f1/dx^2
Dxx(1,:) = 0;
Dxx(1,1) = -2/dx^2;

% Right boundary:
% (f(N+1)-2fN+f(N-1))/dx^2
% = -2fN/dx^2
Dxx(N,:) = 0;
Dxx(N,N) = -2/dx^2;

%--------------------------------------------------
% Third derivative
%--------------------------------------------------

Dxxx = spdiags([-e,2*e,-2*e,e],[-2 -1 1 2], ...
               N,N)/(2*dx^3);

% i = 1
%
% (f(-2)-2f(-1)+2f(1)-f(2))/(2dx^3)
%
% Using odd symmetry:
% f(-1)=-f(1), f(-2)=-f(2)
%
% = 2(f1-f2)/dx^3

Dxxx(1,:) = 0;
Dxxx(1,1) =  2/dx^3;
Dxxx(1,2) = -2/dx^3;

% i = 2

Dxxx(2,:) = 0;
Dxxx(2,1) = -2/(2*dx^3);
Dxxx(2,2) = -1/(2*dx^3);
Dxxx(2,3) =  2/(2*dx^3);
Dxxx(2,4) =  1/(2*dx^3);

% i = N-1

Dxxx(N-1,:) = 0;
Dxxx(N-1,N-3) = -1/(2*dx^3);
Dxxx(N-1,N-2) = -2/(2*dx^3);
Dxxx(N-1,N-1) =  1/(2*dx^3);
Dxxx(N-1,N)   =  2/(2*dx^3);

% i = N

Dxxx(N,:) = 0;
Dxxx(N,N-1) =  2/dx^3;
Dxxx(N,N)   = -2/dx^3;

end