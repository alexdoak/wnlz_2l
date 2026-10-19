
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig = figure();
set(fig,'Units','centimeters','Position',[2 2 12 7]);
FONTSIZE = 10;



load('P=10_maxAmp_SweepR');
p=1;
y=li.y{p}; c=li.c(p);

N_new = round(numel(x)*(0.7324218750000/.7));
[L,x,x_h,dx,dt,kx,N,y] = interpMesh(x,y,L,N_new)

N = numel(x);
dx = x(2)-x(1);
x_h=x(1:N/2+1);
%%% Rescale gravity to make speed 1
%g = g/c^2;
%c=c/c;

g1=g*(r2-r1);


%%% Useful parameters
H=H1+H2; D=H2-H1;
rm = r2-r1;
rp = r2+r1;

%%% Boussinesq parameter
opts.Boussinesq = false;



%%% Get w,d for whole wave
h2 = H2+y;
h1 = H1-y;
q2 = c*(h2-H2);
q1 = c*(h1-H1);
u2 = q2./h2;
u1 = q1./h1;

x
w = u2-u1;
d = h2-h1;


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
%{
% ---------- NUMERICAL EXAMPLE
W = 0.060302266260532;
D = 0.1;
H1=0.45;
H2=0.55;
H=H1+H2; D=H2-H1;

%---------  Stability analysis ---------%
%%% Mesh spacing 
dx = 0.04;
%}
Nk = 200;
dt_list=[];
lam1_imag = [];
lam1_real = [];
lam2_imag = [];
lam2_real = [];

Wc_list = [];
Wck_list = zeros(N/2,Nk);
for p = numel(w)-500:numel(w)

W=w(p);
D=d(p);

%%% Continious stability
k_list = linspace(0,pi/dx,Nk);
[omega1,omega2]=twolStabilityAnalysis(D,H,r1,r2,g,W,k_list,opts);

%%% Numerical stability
[lambda1,lambda2,a2,a1,a0]=twolStabilityAnalysisNumerical(D,H,r1,r2,g,W,k_list,dx,opts);

%lambda1=omega1*1j;
%lambda2=omega2*1j;


if max(real(lambda1)) <  10^(-10) & max(real(lambda2)) <10^(-10)
dt = 2*sqrt(2)/max(abs([lambda1;lambda2]));
z1 = lambda1*dt;
z2 = lambda2*dt;

dt_list=[dt_list,dt];

else 

dt_list=[dt_list,0];

end


lam1_real = [lam1_real,max(abs(real(lambda1)))];
lam1_imag = [lam1_imag,max(abs(imag(lambda1)))];
lam2_real = [lam2_real,max(abs(real(lambda2)))];
lam2_imag = [lam2_imag,max(abs(imag(lambda2)))];

%{
figure()
plot(k_list,real(lambda1),'b'); hold on;
plot(k_list,imag(lambda1),'r')
plot(k_list,real(lambda2),'b--');
plot(k_list,imag(lambda2),'r--')
legend(['1-real';'1-imag';'2-real';'2-imag'])
xlabel('k')
%}


h1 = (H-D)/2;
h2 = (H+D)/2;
Wck = (g*rm*(3*(r1*h2+r2*h1) + k_list.^2*h1*h2*(r1*h1+r2*h2)) ...
     ./(3*r1*r2+k_list.^2/H*(r1*h2^3+r2*h1^3).*(r1*h1+r2*h2))).^.5;

Wc = (g*rm*h1*h2*H^2/(r1*h2^3+r2*h1^3)).^.5;

Wc_list=[Wc_list,Wc];
Wck_list(p,:) = Wck;

growth_num = max(real(lambda1),real(lambda2));
growth =  max(real(omega1),real(omega2));

[gmax,idx] = max(growth_num);

if gmax>10^(-10)
k_fastest(p) = k_list(idx);
else
    k_fastest(p)=nan;
end
gmax_list(p) = gmax;


plot(k_list,growth_num); hold on;
plot(k_list,growth,'--'); hold on;
end

R = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24;
dt=dx;
G1 = R(lambda1*dt);
G2 = R(lambda1*dt);

growth_RK4_1 = log(abs(G1))/dt;
growth_RK4_2 = log(abs(G2))/dt;


plot(k_list, real(lambda1), 'DisplayName', 'FD continuous time','Linewidth',4,'Color','k')
hold on
plot(k_list, growth_RK4_1, 'DisplayName', 'FD + RK4','Linewidth',4,'Color','k')
plot(k_list, real(omega1), 'DisplayName', 'True','Linewidth',4,'Color','k')
yline(0,'k--')


%{
plot(x_h,k_fastest,'k')
hold on
plot(x_h,k_fastest,'kx')
plot(x_h,gmax_list,'r')
plot(x_h,gmax_list,'rx')

figure()
plot(x_h,w,'k','Linewidth',3)
hold on
plot(x_h,Wc_list,'r','Linewidth',3)
plot(x_h,Wck_list,'--');
plot(-x_h,w,'k','Linewidth',3)
plot(-x_h,Wc_list,'r','Linewidth',3)
plot(-x_h,Wck_list,'--');
title('W shortwave stability condition')
legend('W','critical W')
xlabel('x')
ylabel('W(x;k)')

figure()
plot(x_h,lam1_real,'b'); hold on; plot(x_h,lam1_imag,'r'); plot(x_h,lam2_real,'b--'); plot(x_h,lam2_imag,'r--')
xlabel('x')
legend(['1-real';'1-imag';'2-real';'2-imag'])


dt_list=[dt_list,dt_list(end-1:-1:2)];
y=[y;y(end-1:-1:2)];
figure()
subplot(1,2,1)
plot(x,dt_list,'k');
hold on;
xlabel('$x$');
ylabel('$\Delta_{t,max}$','Rotation',0)

subplot(1,2,2)
plot(x,y+H2,'k'); hold on;
plot(x,x*0,'k','Linewidth',2)
plot(x,x./x,'k','Linewidth',2)
xlabel('$x$');
ylabel('$z$','Rotation',0)

%}

%{
R = @(z) 1 + z + z.^2/2 + z.^3/6 + z.^4/24;


R1 = R(z1);
R2 = R(z2);
plot(real(R1),imag(R1),'kx');
hold on
plot(real(R2),imag(R2),'rx');
plot(exp(1j*linspace(-0.01,2*pi,20000)),'k','Linewidth',2)
%}