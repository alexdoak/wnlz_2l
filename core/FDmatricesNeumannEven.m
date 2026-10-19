function [Dx,Dxx,Dxxx] = FDmatrices_neumann_even(N,dx)
% Sparse first- and second-derivative matrices
% using second-order central differences and
% ghost-point Neumann BCs:
%   
%   f_0   = f_2
%   f_N+1  = f_N-1
%
% Inputs:
%   N  - number of grid points
%   dx - grid spacing
%
% Outputs:
%   Dx  - first derivative matrix
%   Dxx - second derivative matrix
%   Dxxx - third derivative matrix

e = ones(N,1);

%---- First derivative matrix
Dx = spdiags([-e e],[-1 1],N,N)/(2*dx);

% Boundary rows from ghost-point condition
Dx(1,:) = 0;
Dx(N,:) = 0;

%----- Second derivative matrix
Dxx = spdiags([e -2*e e],[-1 0 1],N,N)/(dx^2);

% Left boundary:
% (f2 - 2f1 + f(-1))/dx^2
% = 2(f2-f1)/dx^2
Dxx(1,:) = 0;
Dxx(1,1) = -2/dx^2;
Dxx(1,2) =  2/dx^2;

% Right boundary:
% (f(N+1)-2fN+f(N-1))/dx^2
% = 2(f(N-1)-fN)/dx^2
Dxx(N,:) = 0;
Dxx(N,N-1) =  2/dx^2;
Dxx(N,N)   = -2/dx^2;

%---- Third derivative matrix
Dxxx = spdiags([-e,2*e,-2*e,e],[-2 -1 1 2], N,N)/(2*dx^3);

% i = 1
Dxxx(1,:) = 0;

% i = 2
Dxxx(2,:) = 0;
Dxxx(2,1) =  2/(2*dx^3);
Dxxx(2,2) = -1/(2*dx^3);
Dxxx(2,3) = -2/(2*dx^3);
Dxxx(2,4) =  1/(2*dx^3);

% i = N-1
Dxxx(N-1,:) = 0;
Dxxx(N-1,N-3) = -1/(2*dx^3);
Dxxx(N-1,N-2) =  2/(2*dx^3);
Dxxx(N-1,N-1) =  1/(2*dx^3);
Dxxx(N-1,N)   = -2/(2*dx^3);

% i = N
Dxxx(N,:) = 0;

end