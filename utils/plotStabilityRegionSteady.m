clear all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig = figure();
%set(fig,'Units','centimeters','Position',[2 2 12 7]);
FONTSIZE = 10;

% ----------------------------------------------------------------------- %
% Plot steady solution
% ----------------------------------------------------------------------- %

filename = 'P=10_maxAmp_SweepR'
load(filename)


for CNT=25:numel(li.c);
figure()
y=li.y{CNT}; c=li.c(CNT);
r1=li.r1(CNT); r2=li.r2(CNT)

h2 = H2+y;
h1 = H1-y;
q2 = c*(h2-H2);
q1 = c*(h1-H1);
u2 = q2./h2;
u1 = q1./h1;


plot(h2-h1,u2-u1,'b')
hold on;

% ----------------------------------------------------------------------- %
% Stability region
% ----------------------------------------------------------------------- %

H=H1+H2;
%%% Sweep D
D=-H:0.01:H;
h2 = (H+D)/2
h1 = (H-D)/2;
g1 = (r2-r1)*g;

Wc = (g1*h1.*h2.*H.^2./(r1.*h2.^3+r2.*h1.^3)).^.5;
colour='k';
plot(D,Wc,'k')
hold on
plot(D,-Wc,'k')

% Right hand wall
D=linspace(-H,H,200);
h2 = (H+D)/2
h1 = (H-D)/2;
RHS = (g*(r2-r1)/(r2*r1).*(r1.*h2+r2.*h1)).^(1/2);
plot(D,RHS,'k--')
plot(D,-RHS,'k--')


plot([-H,-H],[-RHS(1),+RHS(1)],'k--')
plot([+H,+H],[-RHS(end),+RHS(end)],'k--')
xlim([-1.05,1.05])
xlabel('$h_2-h_1$')
ylabel('$u_2-u_1$')%,FONTSIZE)
title(['r2=' num2str(r2)])
end