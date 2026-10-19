clear all
close all

opts.Bathymetry=true
%%% Parameters
H1 = 0.2*4^.4;
H2 = pi;
% b is bathmetry, no bathymetry implies b=H2, bx=0
if opts.Bathymetry
    b = 2.3;
    bx = -1.1;
else
    b=H2; 
    bx=0;
end

r1 = 0.97;
r2 = 1.1;
g=2.3;

u1 = 0.2^.4;
u2 = 2.123213;

zeta = .2;
h2 = b+zeta;
h1 = H1-zeta;

D = H2-H1;
d = h2-h1;
w = u2-u1;

%%% Test some equations

% Fix some choices for b, dx

dx = 9.2;
Hx=dx;

H = H1+b;

a = (H-d)/(2*H);
b = (H+d)/(2*H);
ax = (Hx*d-H*dx)/(2*H^2);
bx = (H*dx-Hx*d)/(2*H^2);


LHS = a*ax-b*bx;
RHS = ax;

LHS = a^2-b^2;
RHS = -d/H;