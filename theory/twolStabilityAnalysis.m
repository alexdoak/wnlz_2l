function [omega1,omega2]=twolStabilityAnalysis(D,H,r1,r2,g,W,k_list,opts)

% -------------------------------------------------------------
% Useful terms
% -------------------------------------------------------------
rm = r2-r1;
rp = r2+r1;

H1 = (H-D)/2;
H2 = (H+D)/2;

H = H1+H2;
g1 = (r2-r1)*g;
% -------------------------------------------------------------
% Linear operator L, K matrix, L matrix
% -------------------------------------------------------------
if opts.Boussinesq
alpha = D*W/H; beta=H1*H2/3;

Ldd = D*W/H; 
Ldw = -2*H1*H2/H;
Lwd = 1/2*(W^2/H - g1);
Lww = Ldd;


Kd = -D*W/6*Ldd;
Kw = -D*W/6*Ldw;

else
alpha = D*W/H; beta=H1*H2/3;

A0 = 1 + rm*(2*D)/(rp*H-rm*D);
Gamma0 = rm/(2*(rp*H-rm*D));

Ldd = D*W/H; 
Ldw = -2*H1*H2/H;
Lwd = 1/2*(W^2/H - g1) - Gamma0*(g*rm*D+g*H*(2-rp)-2*W^2*D/H);
Lww = Ldd - 2*Gamma0*(H^2-D^2)*W/H;


Kd = -D*W*A0/6*Ldd;
Kw = -D*W*A0/6*Ldw;

end



L = [Ldd,Ldw;
 Lwd,Lww];
K = [0,0; Kd,Kw];



% -------------------------------------------------------------
%  solve
% -------------------------------------------------------------

% Wavenumber
Nk = numel(k_list);
% Store roots
omega1 = zeros(Nk,1);
omega2 = zeros(Nk,1);

i=0;
for k = k_list
i=i+1;

%%% Operator for w_t
if opts.Boussinesq
A=1+k^2/12*(H^2-D^2);
else
A=1+k^2/12*(H^2-D^2)*A0;
end

%%% LHS matrix
P = [1,0;
     0,A];
%{
syms omega
omega = double(solve(det(P*omega-k*L+ k^3*K) == 0, omega));
r=omega;
%}



if opts.Boussinesq
    gamma=0;
    A1 = gamma*D*W*beta/H;
    A2 = -2*alpha;
    A3 = W^2*beta/4*gamma;
    A4 =  (alpha^2+3*beta/H*(W^2/H-g1 ));
    
    
    a2 = 1+beta*k^2;
    a1 = A1*k^3 + A2*k;
    a0 = A3*k^4 + A4*k^2;
    r = roots([a2,a1,a0]);
else

    a2 = 1+beta*A0*k^2;
    a1 = 1j*k*(Ldd+Lww-k^2*(Kw-beta*A0*Ldd));
    a0 = -k^2*(Ldd*Lww-Ldw*Lwd - k^2*(Ldd*Kw-Ldw*Kd));
    r = roots([a2,a1,a0]);

    
end

% Store roots
if k==k_list(1)
    omega1(i) = r(1);
    omega2(i) = r(2);
else
 % Two possible assignments
    cost1 = abs(r(1)-omega1(i-1)) + abs(r(2)-omega2(i-1));
    cost2 = abs(r(2)-omega1(i-1)) + abs(r(1)-omega2(i-1));

    if cost1 <= cost2
        omega1(i) = r(1);
        omega2(i) = r(2);
    else
        omega1(i) = r(2);
        omega2(i) = r(1);
    end
end

end



%{

U1 = -H2/H*W;
U2 = H1/H*W;
theta1 = 1+1/3*k^2*H1^2;
theta2 = 1+1/3*k^2*H2^2;

b2 = r1*H2*(theta1) + r2*H1*(theta2);
b1 = -2*r1*H2*U1 - 2*r2*H1*U2;
b0 = r1*H2*(U1^2+g*H1) + r2*H1*(U2^2-g*H2);

R = roots([b2,b1,b0])




LHS = Ldd+Lww
RHS = 2*(alpha-Gamma0*(H^2-D^2)*W/H)

LHS2 = Ldd*Lww-Ldw*Lwd
RHS2 = alpha^2 + H1*H2/H*(W^2/H - g*rm)-2*H1*H2*Gamma0*(4*alpha*W/H + 1/H*(g*rm*D+g*H*(2-rp)-2*W^2*D/H))
RHS22 = alpha^2 + H1*H2/H*(W^2/H - g*rm)-2*H1*H2*Gamma0*(2*alpha*W/H + 1/H*(g*rm*D+g*H*(2-rp)))
RHS222 = D^2*W^2/H^2+ H1*H2/H*(W^2/H - g*rm-2*Gamma0*(2*D/H*W^2 + 1*(g*rm*D+g*H*(2-rp))))



d1 = 1j*k; 
d2 = -k^2;
LHS = a1^2-4*a0*a2
RHS = 4*d1^2*(Gamma0*4*H1*H2/H-D/H)^2*W^2 - 4*a2*d1^2*(D^2/H^2*W^2 + H1*H2/H*(W^2/H - g*rm-2*Gamma0*(2*D*W^2/H + g*(rm*D+H*(2-rp)))))
RHSS = 4*d1^2*(((Gamma0*4*H1*H2-D)^2 - (1-H1*H2/3*A0*d2)*(D^2+H1*H2-4*Gamma0*D*H1*H2) )*W^2/H^2  + g*H1*H2/H*(1-H1*H2/3*A0*d2)*(rm + 2*Gamma0*(rm*D+H*(2-rp))  ))


RHSSS = 4*d1^2*((4*Gamma0*H1*H2-D)^2 - (1-H1*H2/3*A0*d2)*(D^2+H1*H2-4*Gamma0*D*H1*H2))*W^2/H^2 ...
       +4*d1^2*g*H1*H2/H*(1-H1*H2/3*A0*d2)*(rm + 2*Gamma0*(rm*D+H*(2-rp)))


S=r1*H2+r2*H1
LHS =rm^2*H1*H2 - rm*D*S - S^2
RHS = -r1*r2*H^2
RHS2 = -H1*H2*2*r1*r2 - H1^2*(r2^2+2*r1*r2+r1^2) - H2^2*(r1^2+2*r1*r2+r1^2)
RHS22 = -r1*r2*(D^2+4*H1*H2)


LHS =  (D^2+H1*H2)*S - rm*D*H1*H2
RHS = r1*H2^3+r2*H1^3

LHS = 2*S+rm*D+H*(2-rp)
RHS = 2*H

%}

%{
omega_MCC_lim = [W/(2*H)*(-D+2*(-H1*H2)^.5); W/(2*H)*(-D-2*(-H1*H2)^.5)]*k;




omega_BOUS = (2*alpha*k + [-1;1]*(4*alpha^2*k^2 - 4*(1+beta*k^2)*(alpha^2+3*beta/H*(W^2/H-g1 ))*(k^2)).^.5)/(2*(1+beta*k.^2))
omega_BOUS = (2*alpha*k + [-1;1] ...
   * (4*alpha^2*k.^2 - 4*(1+beta*k.^2).*(alpha^2 + 3*beta/H*(W^2/H-g1))*k.^2)^.5)...
   ./(2*(1+beta*k.^2))


LHS = a1^2-4*a0*a2
RHS1 = (4*alpha^2*k^2- 4*(1+beta*k.^2).*(alpha^2 + 3*beta/H*(W^2/H-g1))*k.^2)
RHS2 = -4*k^2*(3*beta/H*(W^2/H-g1) + k^2*(beta*alpha^2 + 3*beta^2/H*(W^2/H-g1 ))           )


omega_BOUS_lim = 1/beta*(H1*H2/(3*H)*(H1*H2*g1-W^2/H*(H1^2+H2^2-H1*H2)))^.5;



LHS = beta*alpha^2+3*beta^2/H*(W^2/H-g1)
RHS1 = beta*(D^2*W^2/H^2 + 3*beta/H^2*W^2 - 3*beta/H*g1)
RHS2 = H1*H2/(3*H)*(W^2/H*(D^2+3*beta)-H1*H2*g1)



H=1;
H1=0.01:0.01:0.99;
H2=H-H1;

D=H2-H1;
F = H1./H2 .* (3+D./H)./(3-D./H)
%}
%{

%%%% MCC two-layer, using three-layer code and symmetry
H3=H1;  H22=H2*2;

t1 = k.^2*H1^2/3 + 1; t2 = k.^2*H22^2/3 + 1; t3 = k.^2*H3^2/3 + 1; 
t4 = k.^2*H22^2/12 + 1;

D1=0.01; D2=0.01;
omega1=1; omega2=1; r3=1;

a4 = omega1*omega2*H3*t1.*t2 + omega1*r3*H22*t1.*t3 + omega2^2*H1*H22*H3*k.^2.*t4 + omega2*r3*H1*t2.*t3;
a2 = omega1*D2*H3*H22*t1 + omega2*(D1+D2)*H1*H3*t2 + r3*D1*H1*H22*t3;
a0 = H1*H22*H3*D1*D2;



for p=1:numel(k)
omega2=roots([a4(p),-a2(p).*g,a0.*g^2]);
cp = [cp,max(omega2).^(0.5)];
cm = [cm,min(omega2).^(0.5)];
end



end

plot(log(k_list),cm,'r')
hold on
plot(log(k_list),cp,'k')

c_MCC2 = (H1*H2*g1./(H*(1+1/3*H1*H2*k_list.^2))).^.5;
plot(log(k_list),c_MCC2,'g')
%}