function [] = plotStabilityRegion(H1,H2,r1,r2,g)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig = figure();
set(fig,'Units','centimeters','Position',[1,1,5,5]);
FONTSIZE = 10;


H=H1+H2;
%%% Sweep D
D=-H:0.01:H;
h2 = (H+D)/2
h1 = (H-D)/2;
g1 = (r2-r1)*g;

Wc = (g1*h1.*h2.*H.^2./(r1.*h2.^3+r2.*h1.^3)).^.5;
colour='k';
plot(D,Wc,colour)
hold on
plot(D,-Wc,colour)

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
xlabel('$d$')
ylabel('$w$')%,FONTSIZE)