clear all; close all;

% Depth parameters
H1=0.65; 
H2=0.15; 
H=H1+H2;

D = H2-H1;

% Density paramters
r1=0.99; r2=1;  g=1;
g1 = g*(r2-r1); 

% Shear parameters
W = .06;
%W = (H1*H2*(H1+H2)/(H1^2+H2^2-H1*H2) *g1)^.5;


% MCC parameter: gamma=1 is on, gamma=0 is off!
gamma = 0;
A1u = 0;
A2u = 0;
A1h = 1-A1u;
A2h = 1-A2u;

% -------------------------------------------------------------
% Linear operator L
% -------------------------------------------------------------

alpha = D*W/H; beta=H1*H2/3;

Ldd = alpha; 
Ldw = -2*H1*H2/H;
Lwd = 1/2*(W^2/H - g1);
Lww = Ldd;

L = [Ldd,Ldw;
     Lwd,Lww];



%LHS1=Ldd*Lww-Lwd*Ldw;
%RHS1=alpha^2+3*beta/H*(W^2/H-g1);


% -------------------------------------------------------------
% K matrix
% -------------------------------------------------------------
C = -W/(2*H)*(D*(H2*A2h + H1*A1h)/(3) + H2^2*A2u - H1^2*A1u );
Kd = C*Ldd
Kw = C*Ldw;

K = [0,0; Kd,Kw];

% -------------------------------------------------------------
%  solve
% -------------------------------------------------------------

% Wavenumber
%k_list=logspace(0,8,100);
k_list=linspace(0,20,100);
cm=[]; cp=[];
for k = k_list

%%% Operator for w_t
A=1+k^2/12*(H^2-D^2);

a = H1/H; b=H2/H;
beta = 1/3*(H1*H2*(H2*A2h+H1*A1h)/H + H2^2*a*A2u+H1^2*b*A1u);
A =1+ beta*k^2;
%%% LHS matrix
P = [1,0;
     0,A];

  
syms omega
omega = double(solve(det(P*omega-k*L+ k^3*K) == 0, omega));


%{
A1 = gamma*D*W*beta/H;
A2 = -2*alpha;
A3 = W^2*beta/4*gamma;
A4 =  (alpha^2+3*beta/H*(W^2/H-g1 ));


a2 = 1+beta*k^2;
a1 = A1*k^3 + A2*k;
a0 = A3*k^4 + A4*k^2;
omega2 = roots([a2,a1,a0]);
%}
subplot(1,2,1)
plot(k,real(omega),'kx'); hold on;
plot(k,imag(omega),'rx'); hold on;
ylabel('lambda')
xlabel('log(k)')
subplot(1,2,2)
plot(log(k),real(omega/k),'kx'); hold on;
ylabel('c')
xlabel('log(k)')

end


LHS = (Kw-beta*Ldd);
RHS = W*H1*H2/H^2 * (H2^2*A2u - H1^2*A1u - 1/3*D*(H2*A2u+H1*A1u))

omega_BOUSA_lim2 = [0;  -1/beta*(Kw-beta*Ldd)]*k;
omega_BOUSA_lim = [0;  -1/beta* (W*H1*H2/H^2 * (H2^2*A2u - H1^2*A1u - 1/3*D*(H2*A2u+H1*A1u)))]*k;


omega_MCC_lim = [W/(2*H)*(-D+2*(-H1*H2)^.5); W/(2*H)*(-D-2*(-H1*H2)^.5)]*k;



omega_BOUS = (2*alpha*k + [-1;1]*(4*alpha^2*k^2 - 4*(1+beta*k^2)*(alpha^2+3*beta/H*(W^2/H-g1 ))*(k^2)).^.5)/(2*(1+beta*k.^2))
omega_BOUS = (2*alpha*k + [-1;1] ...
   * (4*alpha^2*k.^2 - 4*(1+beta*k.^2).*(alpha^2 + 3*beta/H*(W^2/H-g1))*k.^2)^.5)...
   ./(2*(1+beta*k.^2))

%{

LHS = a1^2-4*a0*a2
RHS1 = (4*alpha^2*k^2- 4*(1+beta*k.^2).*(alpha^2 + 3*beta/H*(W^2/H-g1))*k.^2)
RHS2 = -4*k^2*(3*beta/H*(W^2/H-g1) + k^2*(beta*alpha^2 + 3*beta^2/H*(W^2/H-g1 ))           )
%
%}

omega_BOUS_lim = 1/beta*(H1*H2/(3*H)*(H1*H2*g1-W^2/H*(H1^2+H2^2-H1*H2)))^.5;


%{
LHS = beta*alpha^2+3*beta^2/H*(W^2/H-g1)
RHS1 = beta*(D^2*W^2/H^2 + 3*beta/H^2*W^2 - 3*beta/H*g1)
RHS2 = H1*H2/(3*H)*(W^2/H*(D^2+3*beta)-H1*H2*g1)
%}


H=1;
H1=0.01:0.01:0.99;
H2=H-H1;

D=H2-H1;
F = H1./H2 .* (3+D./H)./(3-D./H)

%{

%%%% MCC two-layer, using three-layer code and symmetry
H3=H1;  H22=H2*2;

t1 = k.^2*H1^2/3 + 1; t2 = k.^2*H22^2/3 + 1; t3 = k.^2*H3^2/3 + 1; 
t4 = k.^2*H22^2/12 + 1;

D1=0.01; D2=0.01;
r1=1; r2=1; r3=1;

a4 = r1*r2*H3*t1.*t2 + r1*r3*H22*t1.*t3 + r2^2*H1*H22*H3*k.^2.*t4 + r2*r3*H1*t2.*t3;
a2 = r1*D2*H3*H22*t1 + r2*(D1+D2)*H1*H3*t2 + r3*D1*H1*H22*t3;
a0 = H1*H22*H3*D1*D2;



for p=1:numel(k)
r2=roots([a4(p),-a2(p).*g,a0.*g^2]);
cp = [cp,max(r2).^(0.5)];
cm = [cm,min(r2).^(0.5)];
end



end

plot(log(k_list),cm,'r')
hold on
plot(log(k_list),cp,'k')

c_MCC2 = (H1*H2*g1./(H*(1+1/3*H1*H2*k_list.^2))).^.5;
plot(log(k_list),c_MCC2,'g')
%}