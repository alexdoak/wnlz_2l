function [] = makeVideo(filename,indspace,framerate)

if nargin<2
    indspace = 1;
end
if nargin<3
    framerate= 1/0.075;    
end


load(filename)

H=H1+H2;
t_data=t_data(1:indspace:end);
zeta_data=zeta_data(:,1:indspace:end);
d_data=d_data(:,1:indspace:end);
w_data=w_data(:,1:indspace:end);



axis tight manual % this ensures that getframe() returns a consistent size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
fig=figure(1);
set(fig,'Units','centimeters','Position',[10 8 20 12])
filename = 'video.avi';



% ---- Create VideoWriter
v = VideoWriter(filename);
v.FrameRate = framerate; 
open(v);

% ---- Plot initial condition
zeta = zeta_data(:,1);

% ---- Create zeta plot
p1 = plot(x,zeta+H2(1),'k');
hold on;

% ---- Create title
title1 = title(['$t=0$']);

% ---- Set limits
ylim([-0.1*H(1),1.1*H(1)])

% ---- Plot bathymetry
if opts.Bathymetry.on
    p2 = plot(x,-H2+H2(1),'k','Linewidth',2);
    p3 = plot(x,H2(1)-H2+H,'k','Linewidth',2);
else    
    p2 = plot([x(1),x(end)],[0,0],'k','Linewidth',2);
    p3 = plot([x(1),x(end)],[H(1),H(1)],'k','Linewidth',2);
end

% ---- Capture first frame
drawnow;
frame = getframe(fig);
writeVideo(v,frame);    


for n=1:numel(t_data)
  
   t=t_data(n);
   d=d_data(:,n);
   w=w_data(:,n);
   
    zeta =1/2*(d-(H2-H1));
     
    
    % ---- Update zeta plot
    zeta = 1/2*(d-(H2-H1));
    p1.YData = zeta + H2(1);
    
    % ---- Update title
    title1.String = ['$t=$' num2str(t)];
    
    % --- Render updated figure
    drawnow;
    
    % ---- Capture frame
    frame = getframe(fig);
    
    % ---- Write frame to video
    writeVideo(v,frame);
    

end


if opts.VideoMaker
    close(v);
end

    