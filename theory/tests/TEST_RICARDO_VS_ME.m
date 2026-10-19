
clear all

% -------------------------------------------------------------
% parameters
% -------------------------------------------------------------

k_list=0.823;
k=k_list;
opts.Boussinesq = false
g=1;
r1 = 0.9;
r2 = 1;
H1=1.3;
H2=0.7;
H = H1+H2;
D = H2-H1;



% -------------------------------------------------------------
% Useful terms
% -------------------------------------------------------------
rm = r2-r1;
rp = r2+r1;

H1 = (H-D)/2;
H2 = (H+D)/2;

H = H1+H2;
g1 = (r2-r1)*g;

rm = r2-r1;
rp = r2+r1;



d2 = -k^2
T=H1*r1+H2*r2;
S=H1*r2+H2*r1;
Wc = (g*rm*(3*(r1*H2+r2*H1)-d2*H1*H2*T)/((3*r1*r2-d2/H^2*T*(r1*H2^3+r2*H1^3))))^.5+0.01;
W = Wc;


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

syms omega
omega = double(solve(det(P*omega-k*L+ k^3*K) == 0, omega));
r_solve=omega;
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




%%
%----- DISPERSION RELATION
%\begin{equation}\label{disp_rel}
%\rho_1 H_2 \left[ (c-U_1)^2+gH_1 + \frac{1}{3} k^2 H_1^2 c^2\right] + \rho_2 H_1 \left[ (c-U_2)^2-gH_2 + \frac{1}{3} k^2 H_2^2 c^2\right]=0.
%\end{equation}


U1 = -H2/H*W;
U2 = H1/H*W;
theta1 = 1+1/3*k^2*H1^2;
theta2 = 1+1/3*k^2*H2^2;

b2 = r1*H2*(theta1) + r2*H1*(theta2);
b1 = -2*r1*H2*U1 - 2*r2*H1*U2;
b0 = r1*H2*(U1^2+g*H1) + r2*H1*(U2^2-g*H2);

R = roots([b2,b1,b0]);


sort(1j*r/k)
sort(R)

