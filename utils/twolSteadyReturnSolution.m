function [] = twolSteadyReturnSolution(Amplitude,filename,H1,H2,r1,r2,g,H,optsSteady)
%%% Parameters 
N = 2^12;
L = 60*H;
[x,dx,dt,kx] = meshSteady(L,N);

%%% Half a wavelength for x, compute solution over half a domain
x_h = x(1:N/2+1);

%%% Initial guess for h
amp=0.01*H;

[c0,y,y_x,y_xx,y_xxx,c,amp] = KDV_twolayer(amp,H1,H2,r1,r2,g,x_h,optsSteady.Boussinesq);
amp = y(1);

%%% Choices for linearised dispersion term: qxxt, hxtt, uxxt, full
optsSteady.linDisp = 'hxtt';

%%% Choices for differentiation: spectral, findif
optsSteady.differentiation = 'findif';

%%% Unknowns stored in a vector. We take only a wavelength
unvec = [y;c];
M = numel(unvec);
amp = y(end);
area=trapz(x_h,y); 

%%% Choose amplitudes to sweep through
if optsSteady.Boussinesq 
    amp_max = (H/2-H2);
else
    amp_max = (r2.^0.5.*H1 - r1.^0.5.*H2)./(r2.^0.5+r1.^0.5);
end

if isfloat(Amplitude)
    if amp_max>0
        amp_start = 0.01*H;
        amp_space = 0.01*H;
        amp_end = Amplitude;
    else
        amp_start = -0.01*H;
        amp_space = -0.01*H;
        amp_end = Amplitude;
    end
else
    if amp_max>0
        amp_start = 0.01*H;
        amp_space = 0.01*H;
        amp_end = amp_max-0.02*H;
    else
        amp_start = -0.01*H;
        amp_space = -0.01*H;
        amp_end = amp_max+0.02*H;
    end
end

li.y={}; li.amp=[]; li.area=[];
li.c=[]; li.effective_wavelength=[];

CNT=0;

%%% First sweep in amplitude
optsSteady.ParamFix = 'amp';

for amp=[amp_start:amp_space:amp_end, amp_end]

    fprintf("Solving Amplitude = %g \n" , amp)
    [unvec,y,c] = twolSteadySolve(H,H2,N,M,dx,kx,x,x_h,g,r1,r2,unvec,amp,area,optsSteady,false);

    
    amp = y(end);
    area=trapz(x_h,y); 

    CNT=CNT+1;
    li.y{CNT} = y;
    li.c(CNT) = c;
    li.amp(CNT) = y(end);
    li.area(CNT) =trapz(x_h,y);
    li.effective_wavelength(CNT) = 2*abs(trapz(x_h,y)/max(y));


    %plot(x_h,y+H2,'k'); hold on; plot(x_h,x_h*0,'k','Linewidth',2) ; plot(x_h,x_h./x_h*H,'k','Linewidth',2)
    %plot(-x_h,y+H2,'k'); hold on; plot(-x_h,x_h*0,'k','Linewidth',2) ; plot(-x_h,x_h./x_h*H,'k','Linewidth',2)
end


% ---- Next sweep in area if we are going for maximum amplitude
if ~isfloat(Amplitude)
    optsSteady.ParamFix = 'area';

    for area=linspace(area,area*5,10)
        
        fprintf("Solving Area = %g \n" , area)

        [unvec,y,c] = twolSteadySolve(H,H2,N,M,dx,kx,x,x_h,g,r1,r2,unvec,amp,area,optsSteady,false);
        
        amp = y(end);
        area=trapz(x_h,y); 
    
        CNT=CNT+1;
        li.y{CNT} = y;
        li.c(CNT) = c;
        li.amp(CNT) = y(end);
        li.area(CNT) =trapz(x_h,y);
        li.effective_wavelength(CNT) = 2*abs(trapz(x_h,y)/max(abs(y)));


%        plot(x_h,y+H2,'k'); hold on; plot(x_h,x_h*0,'k','Linewidth',2) ; plot(x_h,x_h./x_h*H,'k','Linewidth',2)
%        plot(-x_h,y+H2,'k'); hold on; plot(-x_h,x_h*0,'k','Linewidth',2) ; plot(-x_h,x_h./x_h*H,'k','Linewidth',2)

    end
end


% ---- Plot and save solution
figure()
plot(x_h,y+H2,'k'); hold on; plot(-x_h,y+H2,'k');
title(filename)



clear li
li.c=[c];
li.y{1}=y;

save(filename,'x','N','H','c','H1','H2','g','r1','r2','y','li','L','x_h','optsSteady');
clear J; clear yx; clear yxx; clear yxxx; clear unvec; clear unvec_store; clear Dx; clear Dxx; clear Dxx; clear Dxxx
   
end