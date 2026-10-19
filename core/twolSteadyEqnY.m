function [Hi] =twolSteadyEqnY(y,yx,yxx,yxxx,H,c,g,r1,r2,H1,H2,opts)

h1 = H1-y;
h2 = H2+y;

switch opts.linDisp
    
    % ------------------------------------------------------------------ %
    case 'qxxt'
        Hi = 0;

    % ------------------------------------------------------------------ %
    case 'hxtt'
        switch opts.Boussinesq
            case false
            Hi = 1/3*(r1*H1+r2*H2)*c^2*yxx - g*(r1-r2)*y - ...
                 c^2/2*(r1*H1^2./h1.^2 - r2*H2^2./h2.^2 - (r1-r2));
            case true
            Hi =  1/3*(H1+H2)*c^2*yxx - g*(r1-r2)*y - ...
                 c^2/2*(H1^2./h1.^2 - H2^2./h2.^2 );
        end
    
    % ------------------------------------------------------------------ %
    case 'uxxt'
        Hi = 0;
    
    % ------------------------------------------------------------------ %
    case 'full'
        switch opts.Boussinesq
            case 'off'
            Hi = 1/3*(r1*H1^2./h1+r2*H2^2./h2)*c^2.*yxx ... 
                 + 1/6*(r1*H1^2./h1.^2-r2*H2^2./h2.^2)*c^2.*yx.^2 ...
                 - g*(r1-r2)*y ...
                 - c^2/2*(r1*H1^2./h1.^2 - r2*H2^2./h2.^2 - (r1-r2));

             %%% As written in Choi+Camassa (1999)
             %%% Hi2 = 3*y.^2.*(r1*c^2*h2 + r2*c^2*h1 - g*(r2-r1).*h1.*h2)./(r1*c^2*H1^2*h2 + r2*c^2*H2^2*h1) - yx.^2;
            case 'on'
            Hi = 0;
        end
        

end



