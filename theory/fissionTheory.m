% Choose heights
H1=0.05;
H2=1-H1;

% Find F(N)
N=1:10;
F = N.*(N+1)/2;

% Sweep A inbetween A=[0,H2]
A = 0:0.01:H2-0.01;

% Rescale according to H2, as in the paper
Ah = (H2-A)/H2;
H1h = H1/H2;

% RHS
G = Ah.^(-5/4).*(Ah-H1h)./(1-H1h).*((1+H1h)./(Ah+H1h)).^(3/4)

plot(A,G); hold on;

for N=1:10
    plot([0,A(end)],[1,1]*(N.*(N+1)/2),'r')
end

wavePolarityInfinity = sign(H1-H2)
wavePolarityStep = sign(H1-(H2-A))
A_PolarityChange = H2-H1;
if A_PolarityChange>0
    plot([A_PolarityChange,A_PolarityChange],[0,10],'k','Linewidth',2)
end



dG_dAh = ((-H1h+1)./(Ah+H1h)).^(3/4) .* (-4*Ah.^2 + 7*Ah*H1h + 5*H1h^2) ...
    ./(4*Ah.^(9/4).*(H1h-1).*(Ah+H1h))


