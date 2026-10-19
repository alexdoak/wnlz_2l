clear all, close all
addpath '/home/adoak/Documents/MATLAB/wnlz/theory'

%%% Densities
g=1 ; r2=1; r1=0.99;  g1=g*(r2-r1);

%%% Dephts
H1 = 0.3; H2=1.2; 
H=H1+H2; D=H2-H1;

%{
%---------  Boussinesq Solitary wave parameters ---------%
%%% Tabletop solitary wave parameters Boussinesq
aB = (H1-H2)/2;
h1B = H1-aB;
h2B = H2+aB;
cB = (g1*h1B*h2B/H)^.5;

%%% Parameters for flow at crest
u1B = cB*(1-H1./h1B)/H;
u2B = cB*(1-H2./h2B)/H;
WB=u2B-u1B;

W = WB;
%}

% ---------- NUMERICAL EXAMPLE
W = 0.060302266260532;
D = 0.1;
H1=0.45;
H2=0.55;
H=H1+H2; D=H2-H1;

%---------  Stability analysis ---------%
%%% Mesh spacing 
dx = 0.04;


%%% Continious stability
k_list = linspace(0,pi/dx,2000);
[omega1,omega2]=twolStabilityAnalysis(D,H,r1,r2,g,W,k_list,0);

%%% Numerical stability
[lambda1,lambda2]=twolStabilityAnalysisNumerical(D,H,r1,r2,g,W,k_list,dx);

dt = 22;
z1 = lambda1*dt;
z2 = lambda2*dt;

R = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24;


R1 = R(z1);
R2 = R(z2);
plot(real(R1),imag(R1),'kx');
hold on
plot(real(R2),imag(R2),'rx');
plot(exp(1j*linspace(-0.01,2*pi,20000)),'k','Linewidth',2)