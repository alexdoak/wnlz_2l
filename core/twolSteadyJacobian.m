function J = twolSteadyJacobian(Hi,N,unvec,x_h,dx,delta,y,yx,yxx,yxxx,H,c,g,r1,r2,H1,H2,Dx,Dxx,Dxxx,kx,optsSteady,area,amp)
% ----------------------------------------------------------------------- % 
% To compute the jacobian of the central differences, consider what happens
% as you add delta to the left, central, or right node of the stencil
% ----------------------------------------------------------------------- %
% ---- Pre-allocate J
Nh = numel(x_h);
M = numel(unvec);
J = zeros(M,M);
% ---- Perform on central meshpoints
k = [3:Nh-3];

% ---- Right node (k+1)
yx_d  = yx  + delta/(2*dx);   
yxx_d = yxx + delta/(dx^2);   
          
[Hi_D] = twolSteadyEqnY(y,yx_d,yxx_d,yxxx,H,c,g,r1,r2,H1,H2,optsSteady);
J(k + ((k-1)+1)*M) = (Hi_D(k)-Hi(k))/delta;

% ---- Left node (k-1)
yx_d  = yx  - delta/(2*dx);  
yxx_d = yxx + delta/(dx^2);   

[Hi_D] = twolSteadyEqnY(y,yx_d,yxx_d,yxxx,H,c,g,r1,r2,H1,H2,optsSteady);

J(k + ((k-1)-1)*M)   = (Hi_D(k)-Hi(k))/delta;


% ---- central node (k) 
yxx_d = yxx -2*delta/(dx^2);   

[Hi_D] = twolSteadyEqnY(y+delta,yx_d,yxx_d,yxxx,H,c,g,r1,r2,H1,H2,optsSteady);

J(k + ((k-1)+0)*M)   = (Hi_D(k)-Hi(k))/delta;

% ----------------------------------------------------------------------- % 
% Edge meshpoints
% ----------------------------------------------------------------------- %
unvec_store=unvec;
for p=[1:3,Nh-3:M]
    
    unvec=unvec_store;
    unvec(p)=unvec(p)+delta;
    
    % Unpack unvec
    y = unvec(1:N/2+1);
    c = unvec(N/2+2);
    
    % Differentiate
    [yx,yxx,yxxx] =  SpatialDifferentiation(y,N,dx,Dx,Dxx,Dxxx,kx,optsSteady.differentiation);
    
    % Equations
    [Hi_D] = twolSteadyEqnY(y,yx,yxx,yxxx,H,c,g,r1,r2,H1,H2,optsSteady);
    
    switch optsSteady.ParamFix
        case 'amp'
        Hi_D(M) = y(end)-amp;
        case 'area'
        Hi_D(M) = trapz(x_h,y)-area;
    end
    
    J(:,p) = (Hi_D-Hi)/delta;
end
unvec=unvec_store;



% ----------------------------------------------------------------------- % 
% Parameter equation
% ----------------------------------------------------------------------- %
J(M,:) = 0;
switch optsSteady.ParamFix
    case 'amp'
        J(M,Nh) = 1;
    case 'area'
        J(M,2:Nh-1) = dx;
        J(M,1) = dx/2;
        J(M,Nh) = dx/2;
        
        %Hi_D(M) = trapz(x_h,y)-area;
end








end











