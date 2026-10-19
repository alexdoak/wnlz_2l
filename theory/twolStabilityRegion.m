
clear all; close all
load P=0.1_A=0.4_R=5.mat

%%% Non-dimensional parameters
g1 = 1;
H=1;

%%% Colour for plots
COLOUR = ['k','k','k','k']
cnt=0;
for r=[0,.5,.9,.99]
cnt=cnt+1;
colour=COLOUR(cnt);
%%% Fix densities
r2=1+r;
r1=1-r; 

%%% Sweep D
D=-H:0.01:H;
h2 = (H+D)/2
h1 = (H-D)/2;


Wc = (g1*h1.*h2.*H.^2./(r1.*h2.^3+r2.*h1.^3)).^.5;
subplot(2,2,cnt)
plot(D,Wc,colour)
hold on
plot(D,-Wc,colour)


Wc_hyp = ((1-r*D)/(1-r^2)).^.5

plot(D,Wc_hyp,colour,'Linestyle','--')
hold on
plot(D,-Wc_hyp,colour,'Linestyle','--')
plot([D(1),D(1)],[Wc_hyp(1),-Wc_hyp(1)],colour,'Linestyle','--')
plot([D(end),D(end)],[Wc_hyp(end),-Wc_hyp(end)],colour,'Linestyle','--')

xlim([-1.1,1.1])
ylim([-max(Wc_hyp)-.1,max(Wc_hyp)+.1])
title(['r=' num2str(r)])
xlabel('D')
ylabel('W')
end













