clear all, close all
addpath '/home/adoak/Documents/MATLAB/wnlz/theory'

g=1 ; r2=1.2; r1=1;  


R2 = linspace(1.01,5,40);
H_RATIO = linspace(0.01:0.015,5);


%for r2=R2
%    for H_ration=H_RATION
%        H1/
        


g1=g*(r2-r1);






W_list=[];
Wc_list = [];
WB_list=[];
WcB_list=[];
H1_list=[];
H2_list=[];

r1=1;
r2=5;
H2=1/0.15;
H1=1.2;
H=H1+H2;

%for H1=linspace(0.0001,H-0.001,1000)
H2=H-H1;

%---------  BOUSSINESQ ---------%
%%% Tabletop solitary wave parameters Boussinesq
aB = (H1-H2)/2;
h1B = H1-aB;
h2B = H2+aB;
cB = (g1*h1B*h2B/H)^.5;

cB2 = (1/4*g1*H)^.5


%%% Parameters for flow at crest
u1B = cB*(1-H1./h1B)/H;
u2B = cB*(1-H2./h2B)/H;
WB=u2B-u1B;
%W = c*A*(1/h2+1/h1); 
WcB = (h1B*h2B*(h1B+h2B)/(h1B^2+h2B^2-h1B*h2B)*(r2-r1)*g)^.5;

WB_list = [WB_list,WB];
WcB_list = [WcB_list,WcB];

%---------  NONBOUSSINESQ ---------%
%%% Tabletop solitary wave parameters non-Boussinesq
K = (r2^.5+r1^.5);
a = (r2^.5*H1-r1^.5*H2)/ K;
h1 = H1-a;
h2 = H2+a;
c0 = ((g*(r2-r1)*h1*h2)./(r1*H1+r2*H2))^.5;
c = (g*H*(r2^.5-r1.^.5)/(r2^.5+r1^.5))^.5;

%%% Parameters for flow at crest
u1 = c*(1-H1./h1)/H;
u2 = c*(1-H2./h2)/H;
W=u2-u1;
%W = c*A*(1/h2+1/h1);
k=10^(5);
Wc = (g*(r2-r1)*H^2*(r2*h1+r1*h2+(h1*h2/3)*(h2*r1+h1*r2)*k.^2) ...
   ./(r1*r2*H^2+k^2/3*(r2*h1^3+r1*h2^3)*(h2*r1+h1*r2)))^.5;;

Wc2 = ((g*(r2-r1)*H^2*(h1*h2/3)*(h2*r1+h1*r2))./(1/3*(r2*h1^3+r1*h2^3)*(h2*r1+h1*r2))).^.5

Wc3= (g*(r2-r1)*h1*h2*H^2*(h2*r1+h1*r2)/((r2*h1^3+r1*h2^3)*(h2*r1+h1*r2))).^.5

W_list = [W_list,W];
Wc_list = [Wc_list,Wc];



H1_list = [H1_list,H1];
H2_list = [H2_list,H2];
%end
subplot(1,2,1)
plot(H1_list,abs(W_list));
hold on
plot(H1_list,Wc_list);
plot(H1_list,H1_list*0,'k--');
title(['non-Boussinesq, rho=' num2str(r1)])
xlabel('H_1')
ylabel('|W|')


subplot(1,2,2)
plot(H1_list,abs(WB_list));
hold on
plot(H1_list,WcB_list);
plot(H1_list,H1_list*0,'k--');
title('Boussinesq')
xlabel('H_1')
ylabel('|W|')



%{
k = linspace(0.000001,200,200)
Wck = (g1*H*(1+h1*h2/3*k.^2)./(1+k.^2*(h2^2+h1^2-h1*h2)/3)).^.5;

A=(1+h1*h2/3*k.^2);
B=(1+k.^2*(h2^2+h1^2-h1*h2));
%}







%{
k_list = 10^(-8);
H=1;
H2=0.3;
H1=H-H2;
D=H2-H1;
W=7;
gamma=0;
g=1;
r2=1;
r1=0.99;

[omega1,omega2]=twolStabilityAnalysis(D,H,r1,r2,g,W,k_list,gamma)
%}