clear all
H1 = .4234;
H2 = .1239123
H=H1+H2;
D=H2-H1;



g=1.2

zeta = .3;

h2=H2+zeta;
h1=H1-zeta;


d = h2-h1;

w=pi;
u1 = -h2/H*w;
u2 = h1/H*w;


LHS = h2*H2-h1*H1;
RHS = 1/2*(H*D+H*d)


r2=1.434;
r1 = 0.1323;
rm = r2-r1;
rp = r2+r1;



LHS2 = 1+(rm*(D+d))/(rp*H-rm*d)
RHS2 = (rp*H+rm*D)/(rp*H-rm*d)
RHS22 = 2*(r1*H1+r2*H2)/(rp*H-rm*d)




LHS = h1/r1 + h2/r2;
RHS = (rp*H-rm*d)/(2*r1*r2)


LHS = u1^2*h1+u2^2*h2 
RHS = (H^2-d^2)*w^2/(4*H)


LHS = + g/2*h2^2 - g/2*h1^2
RHS =  + g*H*d/2
