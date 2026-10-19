function [dq] = twolDtSpectral(q,kx,H1,H2,r1,r2,g,Lx,opts)
% ----------------------------------------------------------------------- %
% Function to recover dh/dt, and dv/dt, using h_i and v_i at time t_i 
% ----------------------------------------------------------------------- %
% ----- unpack q
N = numel(q)/2;
d = q(1:N);
w = q(N+1:2*N);

% ----- Useful terms
H=H1+H2;
N = numel(w);
rm = r2-r1;
rp = r2+r1;
D = H2-H1;

% ----- Filter varibles
%d = spectralFilter(d,'even');
%w = spectralFilter(w,'odd');

% ---- Differential operators
Dx = @(u,parity) DxSpectral(u,kx,parity,1);
Dxx =  @(u,parity) DxSpectral(u,kx,parity,2);

% ---- Product of two functions, using double padding
Lext = 2*Lx;
mpoo = @(f,g,f_order,g_order) productPadWall(f,g,f_order,g_order,'odd','odd',Lext);
mpoe = @(f,g,f_order,g_order) productPadWall(f,g,f_order,g_order,'odd','even',Lext);
mpee = @(f,g,f_order,g_order) productPadWall(f,g,f_order,g_order,'even','even',Lext);

% ----- Nonlinear differentiation operators (i.e. with dealiasing)
cutfactor = 2/3;
DxNL = @(u,parity) DxSpectral(u,kx,parity,1,cutfactor);
DxxNL =  @(u,parity) DxSpectral(u,kx,parity,2,cutfactor);

% ----- Useful expressions
dx = Dx(d,'even');
wx = Dx(w,'odd');

d2 = mpee(d,d,0,0);
w2 = mpoo(w,w,0,0);


% ----- Conservation of mass
d2w = mpoe(w,d2,0,0);
dd = -1/(2*H) * (Dx(H^2*w,'odd') - Dx(d2w,'odd'));

% ----- Conservation of momentum
if ~opts.Boussinesq
  
    %%% NEED TO CODE STILL LEAVE EMPTY FOR NOW
    warning('Code missing! Boussinesq only please')
else
   % RHS2 = 1/(2*H)*DxNL(w.^2.*d,'even') - 1/2*g*rm*Dx(d,'even') ...
   %     + 1/(12*H)*DxxNL(w.*d.*DxNL((H^2-d.^2).*w,'odd'),'odd');
   
   
    term1 = mpee(w2,d,0,0);
    term1 = 1/(2*H)*Dx(term1,'even');

    term2 = - 1/2*g*rm*Dx(d,'even');
    
    temp3 = mpoe(w,H^2-d2,0,0);
    temp3 = Dx(temp3,'odd');
    temp3 = mpee(temp3,d,0,0);
    temp3 = mpoe(w,temp3,0,0);
    term3 =  1/(12*H)*Dxx(temp3,'odd');

    RHS = term1 + term2 + term3; 
    
%{
    term1_old = 1/(2*H)*Dx(w.^2.*d,'even');
    term2_old = - 1/2*g*rm*Dx(d,'even');
    term3_old = + 1./(12*H)*Dxx(w.*d.*Dx((H.^2-d.^2).*w,'odd'),'odd');

    RHS_old = term1_old+term2_old+term3_old;
%}

         

end


% ------ Solve Lw_t = RHS, where we filter the L spectrum by killing
% ------ high-frequency terms in beta
beta =  (H^2-d2);

%{
betahat = fft(extensionEven(beta));
N = length(betahat);
k = [0:N/2-1 0 -N/2+1:-1]';

alpha = 36;
p = 8;

sigma = exp(-alpha*(abs(k)/max(abs(k))).^p);

betahat = betahat .* sigma;

beta = real(ifft(betahat));
beta = beta(1:numel(d));
%}

% ------ Pre-condition GMRES with undisturbed depth operator P=L(D,W)
beta0 = (H^2-D.^2);
Mfun = @(v) twolApplyPinv(v,kx,beta0);
[dw,flag,relres,iter] = gmres( ...
      @(v) twolApplyL(v,beta,Dxx,H,d,Lext,opts), ...
      RHS,[],1e-10,50,Mfun);


if flag~=0
    warning('gmres failed to converge to prescribed tol')
end

dq = [dd;dw];


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