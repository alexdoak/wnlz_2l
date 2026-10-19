function [unvec,y,c] = twolSteadySolve(H,H2,N,M,dx,kx,x,x_h,g,r1,r2,unvec,amp,area,optsSteady,verbose)


if nargin<16
    verbose = true;
end

H1=H-H2;
%%% Differentiation matrices 
[Dx,Dxx,Dxxx] = FDmatricesNeumann(N/2+1,dx);

%%% Pre-allocated memory
J = zeros(M,M); 
Hi = ones(1,M);
Hi_D=Hi;

while max(abs(Hi))>10^(-10)%*H^2
    
    % Unpack unvec
    y = unvec(1:N/2+1);
    c = unvec(N/2+2);
   
    % Differentiate
    [yx,yxx,yxxx] =  SpatialDifferentiation(y,N,dx,Dx,Dxx,Dxxx,kx,optsSteady.differentiation);
    
    % Equations
    [Hi] = twolSteadyEqnY(y,yx,yxx,yxxx,H,c,g,r1,r2,H1,H2,optsSteady);
    

    switch optsSteady.ParamFix
        case 'amp'
        Hi(M) = y(end)-amp;
        case 'area'
        Hi(M) = trapz(x_h,y)-area;
    end
    % --------------------------------------------------------------- %
    %%% JACOBIAN, computed numerically.
    % --------------------------------------------------------------- %

    
    switch optsSteady.differentiation
        case 'findif'
            delta =10^(-11);
            J=twolSteadyJacobian(Hi,N,unvec,x_h,dx,delta,y,yx,yxx,yxxx,H,c,g,r1,r2,H1,H2,Dx,Dxx,Dxxx,kx,optsSteady,area,amp);

        case 'spectral'
            delta =10^(-11);
            unvec_store=unvec;
            for p=1:M
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
        end
    
    Residuals=max(abs(Hi));
    if verbose
        RL2 = sum(Hi.^2).^(1/2)
    end

    J=sparse(J);
    unvec = unvec + J\(-Hi);
end