function [c0,y,y_x,y_xx,y_xxx,c,amp] = KDV_twolayer(amp,H1,H2,r1,r2,g,x,Boussinesq)

if Boussinesq
    c0 = (g*(r2-r1)*H1*H2/(H1+H2))^.5;
    c0=abs(c0);
    c1 = 3*c0/(2*H1*H2) * (H1^2-H2^2)/(H2+H1);
    c2 = c0*H1*H2/6;
else
    c0 = (g*(r2-r1)*H1*H2/(r1*H1+r2*H2))^.5;
    c0=abs(c0);
    c1 = 3*c0/(2*H1*H2) * (r2/r1) * (r2*H1^2-r1*H2^2)/(r1*H2+r2*H1);
    c2 = c0*H1*H2/6* (r2/r1) * (r1*H1+r2*H2)/(r1*H2+r2*H1);
end




if c1*c2<0
    amp=-abs(amp);
else
    amp=abs(amp);
end


La=(amp*c1/(12*c2))^(1/2);
y = amp*cosh(x*La).^(-2); 
y_x = amp*La*(-2)* cosh(x*La).^(-3) .*  sinh(x*La);
y_xx =  amp*La^2*(6)* cosh(x*La).^(-4) .*  sinh(x*La).^2 + amp*La^2*(-2)* cosh(x*La).^(-2);
y_xxx =  amp*La^3*(-24)* cosh(x*La).^(-5) .*  sinh(x*La).^3 +   amp*La^3*(12)* cosh(x*La).^(-3) .*  sinh(x*La) +  amp*La^3*(4)* cosh(x*La).^(-3) .* sinh(x*La);


c=c0+c1/3*amp;   

KDV=(c0-c)*y_x + c1*y.*y_x + c2*y_xxx;

