clear all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig = figure();
%set(fig,'Units','centimeters','Position',[2 2 12 7]);
FONTSIZE = 10;


filename = 'stable_Maxamp_R=1.8';
load(filename)
[~,M] = size(zeta_data);

for p = M:-1:M-20
    zeta = abs(zeta_data(:,p));
    if max(zeta)>10 | sum(isnan(zeta))>0
        zeta_data = zeta_data(:,1:end-1);
        d_data = d_data(:,1:end-1);
        w_data = w_data(:,1:end-1);
        t_data = t_data(1:end-1);
    end    

end
[~,M] = size(zeta_data);
indspace = round(M/20);


subplot(1,2,1)
%xwindow=[x(1),x(end)];
xwindow = [-20,20]
ywindow=[-0.05,1.05];
t_final = t_data(end);
%indspace=100;
[IND]=plotWaterfall(filename,indspace,xwindow,ywindow,true,t_final);


H=H1+H2;


subplot(1,2,2)
H=H1+H2;
%%% Sweep D
D=-H:0.01:H;
h2 = (H+D)/2
h1 = (H-D)/2;
g1 = (r2-r1)*g;

IND=1:M;

Wc = (g1*h1.*h2.*H.^2./(r1.*h2.^3+r2.*h1.^3)).^.5;
colour='k';
plot(D,Wc,colour)
hold on
plot(D,-Wc,colour)


plot(d_data(:,IND),w_data(:,IND),'b')



% Right hand wall
D=linspace(-H,H,200);
h2 = (H+D)/2
h1 = (H-D)/2;
RHS = (g*(r2-r1)/(r2*r1).*(r1.*h2+r2.*h1)).^(1/2);
plot(D,RHS,'k--')
plot(D,-RHS,'k--')


plot([-H,-H],[-RHS(1),+RHS(1)],'k--')
plot([+H,+H],[-RHS(end),+RHS(end)],'k--')
xlim([-1.05,1.05])
xlabel('$h_2-h_1$')
ylabel('$u_2-u_1$')%,FONTSIZE)