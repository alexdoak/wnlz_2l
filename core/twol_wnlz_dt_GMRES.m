function [dd,dw] = twol_wnlz_dt_GMRES(d,w,Dx,Dxx,H1,H2,r1,r2,g,BC)
% ----------------------------------------------------------------------- %
% Function to recover dh/dt, and dv/dt, using h_i and v_i at time t_i 
% ----------------------------------------------------------------------- %

% ----- Useful terms
dx = Dx*d;
wx = Dx*w;
H=H1+H2;
N = numel(w);

% ----- Conservation of mass
dd = -1/(2*H) * Dx*((H^2-d.^2).*w);

% ----- Conservation of momentum
RHS = 1/(2*H)*Dx*(w.^2.*d) - 1/2*g*(r2-r1)*dx ...
    + 1/(12*H)*Dxx*(w.*d.*Dx*((H^2-d.^2).*w));


switch BC
    case "wall"
        RHS(1)=0; RHS(end)=0;
end


[dw,flag,relres,iter] = gmres(@(v) twol_applyL(v,Dxx,H,d,BC), RHS);


%{
%-------- DEBUGGING ----------
%-------- Check the steady state version of the equations ---------
RHS_CON_MASS =  (-1/(2*H) * Dx*((H^2-d.^2).*w));
LHS_CON_MASS = -c*dx ;
plot(RHS_CON_MASS); hold on
plot(LHS_CON_MASS)

RHS_CON_MOM = RHS;
LHS_CON_MOM = -c*A*wx;
figure()
plot(RHS_CON_MOM); hold on
plot(LHS_CON_MOM)
%}