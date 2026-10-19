function [dq] = twolDtFindif(q,Dx,Dxx,H1,H2,r1,r2,g,opts)
% ----------------------------------------------------------------------- %
% Function to recover dh/dt, and dv/dt, using h_i and v_i at time t_i 
% ----------------------------------------------------------------------- %
% ----- unpack q
N = numel(q)/2;
d = q(1:N);
w = q(N+1:2*N);

% ----- Useful terms
dx = Dx*d;
wx = Dx*w;
b = H2;
bx = Dx*b;

H=H1+H2;
N = numel(w);
rm = r2-r1;
rp = r2+r1;
D = H2-H1;

h1 = (H+d)/2;
h2 = (H-d)/2;

% ---- Check h1 and h2 are not too small (touch up/down singularity)
if min(h1) < 1e-8 || min(h2) < 1e-8
    error('twolDtFindif:DryCell', ...
        'Layer depth approaching zero. min(h1) = %.3e, min(h2) = %.3e', ...
        min(h1), min(h2));
end

% ----- Conservation of mass
dd = -1/2 * Dx*((H.^2-d.^2).*w./H);

% ----- Conservation of momentum
if ~opts.Boussinesq
    factor = (1+rm*(D+d)./(rp*H-rm*d));
    RHS = 1/(2*H)*Dx*(w.^2.*d) - 1/2*g*rm*dx ...
        + 1/(12*H).*factor.*Dxx*(w.*d.*Dx*((H^2-d.^2).*w)) ...
        - rm./(2*(rp*H-rm*d)).* ...
        Dx*((H^2-d.^2).*w.^2/H + g*H*(2-rp)*d + g*rm*d.^2/2);

     A = speye(N) - (1/12)*spdiags(factor,0,N,N)*Dxx*spdiags(H^2-d.^2,0,N,N);


     % factor2 = (rp*H+rm*D)./(rp*H-rm*d);
     % RHS2 = 1/(2*H)*Dx*(w.^2.*d) ...
     %    + 1/(12*H).*factor2.*Dxx*(w.*d.*Dx*((H^2-d.^2).*w)) ...
     %    - rm./(2*(rp*H-rm*d)).* ...
     %    Dx*((H^2-d.^2).*w.^2/H + 2*g*H*d);
     % 
     % TEST1=- 1/2*g*rm*dx  - rm./(2*(rp*H-rm*d)).* ...
     %    Dx*((H^2-d.^2).*w.^2/H + g*H*(2-rp)*d + g*rm*d.^2/2);
     % TEST2 =  - rm./(2*(rp*H-rm*d)).* ...
     %    Dx*((H^2-d.^2).*w.^2/H + 2*g*H*d )%- g*rm*d.^2);

else
 
    if ~opts.Bathymetry.on
        RHS = 1./(2*H)*Dx*(w.^2.*d) - 1/2*g*rm*dx ...
            + 1./(12*H)*Dxx*(w.*d.*Dx*((H.^2-d.^2).*w));

         A = speye(N) - (1/12)*Dxx*spdiags(H^2-d.^2,0,N,N);
    else
        RHS = 1./(2*H).*Dx*(w.^2.*d) ...
            - 1/2*g*rm*(dx-bx) ...
            - (H./6).*Dxx*(d.*w.*dd./H) ...
            - (bx./12).*Dx*(d.*w.*dd./H) ...
            - (b.*bx./4).*Dx*(w.*dd./H);

        A = speye(N) ...
           - spdiags(H./12,0,N,N)*Dxx*spdiags((H.^2-d.^2)./H,0,N,N) ...
           - spdiags(bx./24,0,N,N)*Dx*spdiags((H.^2-d.^2)./H,0,N,N) ...
           - spdiags(b.*bx./4,0,N,N)*Dx*spdiags((H-d)./H,0,N,N);
    end

%$        A = speye(N) - (H/12).*Dxx*spdiags((H^2-d.^2)./H,0,N,N) ...
  %          - (bx/24).*Dx*spdiags((H^2-d.^2)./H,0,N,N) ...
  %          - (b.*bx/4).*Dx*spdiags((H-d)./H,0,N,N);
end



switch opts.BC
    case "wall"
        A(1,:)   = 0;
        A(1,1)   = 1;
        
        A(end,:) = 0;
        A(end,end)=1;

        RHS(1)=0; RHS(end)=0;
end


dw = A\RHS;

dq = [dd;dw];


%{
MAT = spdiag(
%}

%{
%-------- DEBUGGING ----------
%-------- Check the steady state version of the equations ---------
RHS_CON_MASS =  (-1/(2*H) * Dx*((H^2-d.^2).*w));
LHS_CON_MASS = -c*dx ;
plot(RHS_CON_MASS); hold on
plot(LHS_CON_MASS,'x')

RHS_CON_MOM = RHS;
LHS_CON_MOM = -c*A*wx;
figure()
plot(RHS_CON_MOM); hold on
plot(LHS_CON_MOM)
%}