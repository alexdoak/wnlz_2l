function [IND] = plotWaterfall(filename,indspace,xwindow,ywindow,movewindow,rescalespeed,firstcol,backwards,caz_cel,IND)

if nargin<10
    IND=[];
end

if nargin<9
    caz = 3.520000000000001e+02, 
    cel  =35;
else
    caz=caz_cel(1); cel=caz_cel(2);
end

if nargin<8
    backwards = true;
    'backwards default'
end

if nargin<7
    firstcol = false;
end
if nargin<6
    rescalespeed = false;
    'rescale speed default off'
end

red  = [0.75 0.10 0.10];
blue = [0.00 0.35 0.75];


%%% Load data and clean 
load(filename)
c1=cWave{1};



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


%%% Rescale the speed
if rescalespeed
    c0 = (g*(r2-r1)*H1*H2/(r1*H2+r2*H1))^.5;
    t_data=t_data*c0;
    c1=c1/c0;
end



%%% Select t data to use, and extract zeta(x,t) from data
if isempty(IND)
    if backwards
        IND = numel(t_data):-indspace:1; %IND=IND(end:-1:1)
    else
        IND=1:indspace:numel(t_data);
    end
end

IND=sort(IND);
t=t_data(IND);
Z = 1/2*(d_data(:,IND)'-H2+H1);

%%% Build mesh grid using x,t
[X,T] = meshgrid(x,t);

for n = 1:length(t)
    ti=t(n);
    if movewindow
        xshift = c1*ti+opts.WaveLocation{1};
    else
        xshift=0;
    end
    if n==length(t) & firstcol
        plot3(x-xshift,t(n)*ones(size(x)),Z(n,:),'Color',blue);
    elseif n==1 & firstcol
        plot3(x-xshift,t(n)*ones(size(x)),Z(n,:),'Color',red);
    else
        plot3(x-xshift,t(n)*ones(size(x)),Z(n,:),'k');
    end
        hold on;
    xlim(xwindow)
    zlim(ywindow)
end

if exist('caz')
    view(caz,cel)
end
%{

view([-37.5 30])

%xlim([-20,20])

xlabel('$X$')
ylabel('$t$')
zlabel('$\zeta$','Rotation',0)


subplot(1,2,2);
colour='k';

%%% Sweep D
D=-H:0.01:H;
h2 = (H+D)/2
h1 = (H-D)/2;
g1 = (r2-r1)*g;

%g1=g;
Wc = (g1*h1.*h2.*H.^2./(r1.*h2.^3+r2.*h1.^3)).^.5;
plot(D,Wc,colour)
hold on
plot(D,-Wc,colour)

plot(d_data(:,1),w_data(:,1),'b')

%Wc_hyp = ((1-r*D)/(1-r^2)).^.5

%plot(D,Wc_hyp,colour,'Linestyle','--')
%hold on
%plot(D,-Wc_hyp,colour,'Linestyle','--')
%plot([D(1),D(1)],[Wc_hyp(1),-Wc_hyp(1)],colour,'Linestyle','--')
%plot([D(end),D(end)],[Wc_hyp(end),-Wc_hyp(end)],colour,'Linestyle','--')

xlim([-1.1,1.1])
%ylim([-max(Wc_hyp)-.1,max(Wc_hyp)+.1])
%title(['r=' num2str(r)])
xlabel('D')
ylabel('W','Rotation',0)

%}