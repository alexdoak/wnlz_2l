% ----------------------------------------------------------------------- %
% This code integrates the Boussinesq system found in **. We spatially 
% discretise with second-order differences  and integrate in time using
% forward Euler
%
% N - Number of points in x
% L - Length of domain
% H - Undistrurbed depth
% h(x) - Fluid depth
% q(x) - Flux = bar(u)*h
% g - Gravity
% t0 - Initial time
% tprint - Time-spacing for data storage
% tfinal - Final time
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
% Problem setup
% ----------------------------------------------------------------------- %
%%% Clear workspace, plots
clear all; close all;

%%% Add path to aux files
addpath core

%%% Set up plot
axis tight manual % this ensures that getframe() returns a consistent size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig=figure(1);
set(fig,'Units','centimeters','Position',[10 8 20 12])

load 'run2_twol.mat';

zeta_data = 1/2*(d_data-(H2-H1));; hold on;


writerObj = VideoWriter('myVideo3.avi','Motion JPEG AVI');
writerObj.Quality = 100;
writerObj.FrameRate = 10;
open(writerObj);


for cnt=1:100;numel(t_data)
 
   zeta = zeta_data(:,cnt);
   t = t_data(cnt);

  
   %%% Plot
   figure(1); 
   p1=plot(x,H2+zeta,'k');
   title1=title(['$t=$' num2str(t) ]);
   xlabel('$x$')
   ylabel('$z$')
   xlim([min(x),max(x)])
   ylim([-0.01,1.01])
 
   

   
    
    drawnow limitrate
    writeVideo(writerObj, getframe(fig));


    filename = 'wnl_2layer_run3.gif';
    % Capture the plot as an image 
    frame = getframe(fig); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if cnt==1
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.025); 
    end
    
    
    delete(title1); delete(p1); 

   
   
   
    


end

