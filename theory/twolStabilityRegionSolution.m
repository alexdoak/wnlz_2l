
%clear all; close all
load P=0.1_A=0.5_R=5.mat
%load P=0.1_A=limiting_R=5.mat

% Recover (d,w) from zeta, H_i.
h2 = H2+y;
h1 = H1-y;
q2 = c*(h2-H2);
q1 = c*(h1-H1);
u2 = q2./h2;
u1 = q1./h1;


w = u2-u1;
d = h2-h1;


%%% Non-dimensional parameters
g1 = g*(r2-r1);
H=1;


%%% Colour for plots
COLOUR = ['k','k','k','k']
cnt=0;

colour=COLOUR(1);

%%% Sweep D
D=-H:0.01:H;
h2 = (H+D)/2
h1 = (H-D)/2;


Wc = (g1*h1.*h2.*H.^2./(r1.*h2.^3+r2.*h1.^3)).^.5;
plot(D,Wc,colour)
hold on
plot(D,-Wc,colour)

plot(d,w,'b')

xlabel('D')
ylabel('W')
%{
Wc_hyp = ((1-r*D)/(1-r^2)).^.5

plot(D,Wc_hyp,colour,'Linestyle','--')
hold on
plot(D,-Wc_hyp,colour,'Linestyle','--')
plot([D(1),D(1)],[Wc_hyp(1),-Wc_hyp(1)],colour,'Linestyle','--')
plot([D(end),D(end)],[Wc_hyp(end),-Wc_hyp(end)],colour,'Linestyle','--')

xlim([-1.1,1.1])
ylim([-max(Wc_hyp)-.1,max(Wc_hyp)+.1])
title(['r=' num2str(r)])
%}













