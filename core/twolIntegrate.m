function [t_data,d_data,w_data,zeta_data] = twolIntegrate(filename,q,H1,H2,H,r1,r2,g,x,Nx,Lx,Dx,Dxx,t0,dt,tprint,tfinal,opts)

d = q(1:Nx);
w = q(Nx+1:2*Nx);

if opts.VideoMaker
    % ================================================================
    % Video initialisation
    % ================================================================
    % ---- Create figure
    fig = figure(1);
    set(fig,'Units','centimeters','Position',[10 8 20 12])
    
    % ---- Create VideoWriter
    v = VideoWriter(filename);
    v.FrameRate = 1/0.075;       
    open(v);
    
    % ---- Plot initial condition
    zeta = 1/2*(d-(H2-H1));
    u2 = (H-d)/(2*H).*w;
    u1 = (H+d)/(2*H).*w;
    
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
end


% ----------------------------------------------------------------------- %
% RK4 integration
% ----------------------------------------------------------------------- %
% ---- Store data at print times
d_data = d;
w_data = w;
t_data = [0];

if opts.FinDiff
    F=@(q,t)twolDtFindif(q,Dx,Dxx,H1,H2,r1,r2,g,opts);
else
    F=@(q,t)twolDtSpectral(q,kx,H1,H2,r1,r2,g,Lx,opts);
end

cnt=0;
for t=t0+tprint:tprint:tfinal
    for tint=0:dt:tprint-dt
        [dq1] = F(q,t);
        [dq2] = F(q+dq1*dt/2,t+dt/2);
        [dq3] = F(q+dq2*dt/2,t+dt/2);
        [dq4] = F(q+dq3*dt,t+dt);
        
        dq = (1/6)*(dq1+2*dq2+2*dq3+dq4)*dt;

        q = q + dq;
        d = q(1:Nx);
        w = q(Nx+1:2*Nx);
    end
   
   % ---- Store new values
   d_data=[d_data,d];
   w_data=[w_data,w];
   t_data=[t_data,t];
   t;
   cnt=cnt+1;
   
    if opts.VideoMaker
   
        zeta =1/2*(d-(H2-H1));
      %  u2 = (H-d)/(2*H).*w;
      %  u1 = (H+d)/(2*H).*w;
    
        
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


zeta_data =1/2*(d_data-(H2-H1));
end

if opts.VideoMaker
    close(v);
end
