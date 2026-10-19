function [lambda1,lambda2,a2_list,a1_list,a0_list]=twolStabilityAnalysisNumerical(D,H,r1,r2,g,W,k_list,dx,opts)

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
lambda1 = zeros(Nk,1);
lambda2 = zeros(Nk,1);

a2_list = [];
a1_list = [];
a0_list = [];
i=0;
for k = k_list
i=i+1;

%%% Mesh spacing terms
delta1 = 1j*sin(k*dx)/(dx);
delta2 = -4*(sin(.5*k*dx)).^2/(dx^2);


if opts.Boussinesq
%%% Polynomial coefficients
A2 = -2*alpha;
A4 =  (alpha^2+3*beta/H*(W^2/H-g1 ));


a2 = 1-beta*delta2;
a1 =  A2*delta1;
a0 =  A4*delta1^2;
r = roots([a2,a1,a0]);


else
a2 = 1-beta*A0*delta2;
a1 = -delta1*(Ldd+Lww+delta2*(Kw-beta*A0*Ldd));
a0 = delta1^2*(Ldd*Lww-Ldw*Lwd +delta2*(Ldd*Kw-Ldw*Kd));
r = roots([a2,a1,a0]);

end


% Store roots
if k==k_list(1)
    lambda1(i) = r(1);
    lambda2(i) = r(2);
else
 % Two possible assignments
    cost1 = abs(r(1)-lambda1(i-1)) + abs(r(2)-lambda2(i-1));
    cost2 = abs(r(2)-lambda1(i-1)) + abs(r(1)-lambda2(i-1));

    if cost1 <= cost2
        lambda1(i) = r(1);
        lambda2(i) = r(2);
    else
        lambda1(i) = r(2);
        lambda2(i) = r(1);
    end
end

a2_list=[a2_list,a2];
a1_list=[a1_list,a1];
a0_list=[a0_list,a0];
end





