clearvars
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')

%% Filenames of waves to plot
filename1 = "P=10_R=0.2_branch";
filename2 = "P=10_R=0.2_branch";
filename3 = "P=10_R=0.2_branch";
filename4 = "P=10_R=0.2_branch";

%% Figure setup

fig = figure();
set(fig,'Units','centimeters','Position',[2 2 13 6]);

ax1=gca;

FONTSIZE = 10;

LW = 1.0;
LW_initial = 1.8;
LW_critical = 1.2;

red  = [0.75 0.10 0.10];
blue = [0.00 0.35 0.75];

% Tiled layout

%{
t = tiledlayout(2,4);
t.TileSpacing = 'compact';
t.Padding = 'compact';
%}

%% Plot the stability region


% ---- Define r1,h1
r1 = 1;
h1 = 1;

% ---- Sweep r2,h2
r2 = linspace(1.01,8,1000);
h2 = linspace(0.01,5,1000);
[R2,H2] = meshgrid(r2,h2);

% ---- Algebraic curve
F = @(r2,h2) r1*h2.^2 + r2*h1.^2 ...
    - (h1.^2 + 4*h1.*h2 + h2.^2).*(r1*r2).^0.5;

Z = F(R2,H2);

hold(ax1,'on');

% Maximum amplitude
h1 = 1;
r1 = 1;

Amax = @(r2,h2) ...
    (r2.^0.5.*h1 - r1.^0.5.*h2) ...
    ./(r2.^0.5+r1.^0.5);

Am = Amax(R2,H2);

% Plot amplitude
pcolor(ax1,R2,H2,abs(Am)./(h1+H2));
shading(ax1,'interp');

%colormap(ax1,parula);
colormap(ax1,turbo);
cb = colorbar(ax1);
%cb.Label.String = '$\frac{|a_m|}{H}$';
cb.Label.Rotation = 0;
cb.Label.Interpreter = 'latex';

hold(ax1,'on');

% Stability boundary
contour(ax1,R2,H2,Z,[0 0], ...
    'k','LineWidth',1.5);

text(ax1,4,1.15,'\textbf{STABLE}', ...
    'Interpreter','latex', ...
    'FontSize',11, ...
    'HorizontalAlignment','center', ...
    'Color','white');

text(ax1,6,0.15,'\textbf{UNSTABLE}', ...
    'Interpreter','latex', ...
    'FontSize',11, ...
    'HorizontalAlignment','center', ...
    'Color','black');

%text(ax1,4.0,0.48,'$F=0$', ...
%    'Interpreter','latex', ...
%    'FontSize',FONTSIZE, ...
%    'Rotation',10, ...
%    'Color','black');


%xlabel(ax1,'$\rho_2/\rho_1$','FontSize',FONTSIZE,'Rotation',0);
%ylabel(ax1,'$H_2/H_1$','FontSize',FONTSIZE,'Rotation',0);

xlim(ax1,[1 max(r2)]);
ylim(ax1,[0 2]);

set(ax1,'FontSize',FONTSIZE);

ax1.Position = [0.1 0.15 0.72 0.8];

%% Save the picture
% ---- Add an invisible object occupying the full desired area

% Invisible points defining the desired export boundary
plot([0 1], [0 1], ...
    'Color', 'none', ...
    'HandleVisibility', 'off');

hold off
exportgraphics(gcf,'stability_region_smaller.jpg','Resolution',300,'Padding','figure')



%{
for filename=[filename1,filename2,filename3,filename4]
    load(filename);
    plot(r2./r1,H2/H1,'kx')
end
%}




%{
% ----------------------------------------------------------------------- %
%% Plot solutions
ax2 = nexttile(t,3,[1,2])
filename = filename1;
plotWave(filename);

ax3 = nexttile(t,7,[1,2])
filename = filename2;
plotWave(filename);

%}
% ----------------------------------------------------------------------- %
%% plotWave function for plotting the solutions
function plotWave(filename,xrange)
    if nargin<2
        xrange=[-15,15];
    end
    
    % Load wave
    load(filename);
    % Reflect y
    y=[y;y(end-1:-1:2)];
    % Plot profile
    plot(x,y+H2,'k','LineWidth',1.2);
    hold on;
    % Plot walls
    H=H1+H2;
    plot([-100,100],[0,0],'k','LineWidth',2);
    plot([-100,100],[H,H],'k','LineWidth',2);
    % Extend wave 
    plot([-100,x(1)],[H2,H2],'k','LineWidth',1.2);
    plot([x(end),100],[H2,H2],'k','LineWidth',1.2);

    % ylim
    ylim([-0.02,1.02])
    xlim(xrange);
end


