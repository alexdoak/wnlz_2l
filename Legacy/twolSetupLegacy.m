clear all
close all
% ---- Add path to aux files
PATH = 'C:\Users\maaad\Documents\MATLAB\wnlz_2l';
addpath(genpath(PATH))

clc
format short

% ----------------------------------------------------------------------- %
%      Introduction
% ----------------------------------------------------------------------- %
fprintf( "********************************************************* \n" + ...
        "This code accompanies the paper (reference paper here). \n" + ...
        "It solves the two-layer strongly nonlinear Boussinesq-type  \n" + ...
        "system. It can include any number of solitary waves, as \n" +...
        "well as bathymetry, either a tanh step or a Gaussian. \n" + ...
        "If you wish, you can explore the file 'twol.m' to setup  \n " +...
        "parameters faster. \n\n" +...
        "Otherwise, continue through the prompts to set up a simulation \n" +...
         "********************************************************* \n");

input("Press any key to continue \n\n")
% ----------------------------------------------------------------------- %
%      Choose parameters
% ----------------------------------------------------------------------- %
fprintf( "********************************************************* \n" + ...
        "First, we will choose some parameter values of the  \n" + ...
        "background state. You will be asked for the total depth (H), \n" + ...
        "and domain length (Lx). These form the computati-onal domain. \n" + ...
        "Upon choosing densities, a choice of a lower \n" +...
        "layer depth will be requested. Since this code explores solitary   \n" +...
        "wave interaction, please choose a H2 sufficiently far from  \n" + ...
        "the criticality condition for which no solitary waves exist \n" +...
        " (the criticality condition will be given) \n" + ...
         "********************************************************* \n");


% ---- Physical parameters
H = chooseFloatValue("Value for total depth, H? ",0,[],1);
 
Lx = chooseFloatValue("Value of domain length, Lx, x=[-Lx,0]? ", 0, [],200*H);

r1 = chooseFloatValue("Value for upper layer density, r1? ", 0,[],1);
r2 = chooseFloatValue("Value for lower layer density, r2? ", r1,[],1.01);

% ---- Boussinesq approximation
optsSteady.Boussinesq = chooseBooleanValue("Boussinesq approximation on densities? (y/n) ", true);
opts.Boussinesq = optsSteady.Boussinesq;

% ---- H2, with criticality avoided
if opts.Boussinesq
    H2_criticality = 0.5*H;
else
    H2_criticality = H*(r2)^.5/(r1^.5 + r2^.5);

end

fprintf("\n************************************************************\n" + ...
    "NOTE:\nFor these parameter values, criticality (no solitary waves) occurs\n " + ...
    "for H2=%g. Please keep this in mind when choosing H2!\n" +...
    "************************************************************\n\n",H2_criticality)
H2 = chooseFloatValue("Value of lower layer depth, H2? ",0,H);
H1 = H-H2;


% ---- Long-wave speed (assuming g=1) 
g=1;
if opts.Boussinesq
    c0 = (g*(r2-r1)*H1*H2/(H2+H1))^.5;
else
    c0 = (g*(r2-r1)*H1*H2/(r1*H2+r2*H1))^.5;
end

% ---- Choose gravity, and rescale c0 accordingly 
opts.RescaleSpeed = false;
g = chooseFloatValue('Value for gravity g? The default is to choose g such that\nthe long-wave speed c0 = H. (recommended)', 0,[],H^2/c0^2);
if opts.Boussinesq
    c0 = (g*(r2-r1)*H1*H2/(H2+H1))^.5;
else
    c0 = (g*(r2-r1)*H1*H2/(r1*H2+r2*H1))^.5;
end

% ---- Mesh spacing and time step (depends on velocity scale)
defaultMesh = chooseBooleanValue("Would you like to use default spatial and temporal spacing (recommended)?",true);
if defaultMesh
    dx = 0.02*H;
    Nx = round(Lx/dx)+1; 
    t0 = 0;
    dt = dx/c0;
    tprint = 20*dt;
    tfinal = Lx;
else
    dx = chooseFloatValue("Value for mesh spacing? ", 0,[],0.02*H);
    Nx = round(Lx/dx)+1; 
    fprintf('PLEASE NOTE, if Lx/dx is not an integer, dx will be minimally modified to make this condition true \n\n')
    t0 = 0;
    dt = chooseFloatValue("Value for time step dt? Default is to choose dt such that dt=dx*c0 ", 0,[],dx*c0);
    tprint = chooseFloatValue("How many time steps to pass before storing data? Default is 20", 0, [], 20) * dt;
    tfinal = chooseFloatValue("What to set as the final time?", 0, [], Lx);
end

% ----------------------------------------------------------------------- %
%      Solitary wave parameters
% ----------------------------------------------------------------------- %
% ---- Maximum amplitude
if optsSteady.Boussinesq 
    amp_max = (H/2-H2);
else
    amp_max = (r2.^0.5.*H1 - r1.^0.5.*H2)./(r2.^0.5+r1.^0.5);
end

% ---- Polarity
if amp_max>0
    polarity='ELEVATION';
else
    polarity='DEPRESSION';
end
fprintf("********************************************************* \n");
fprintf('Given these parameter choices, solitary waves are of %s. \n', polarity);
fprintf('The maximum AMPLITUDE of a solitary wave is given by %g \n', amp_max);
fprintf("Solitary waves broaden as they approach the maximum value.\n" + ...
       "You will be given the option of adding any number of solitary waves.\n" + ...
       "You must specify the chosen amplitude, starting position, and direction.\n" + ...
       "Choosing amplitude=max does not give a wavefront;\n" +...
       "instead, a close-to-limiting tabletop solitary wave is computed.\n\n" + ...
       "WARNING:\nSteady-state solutions are recovered numerically via continuation.\n" + ...
       "In rare instances, this fails. You will be presented with a plot of the\n" + ...
       "recovered solution before proceeding. If the method has failed, you\n" + ...
       "can use the steady-state solver to recover the solution manually, and\n" + ... 
       "use the code twol.m to load in that solution.\n\n" +...
       "WARNING:\nWhen choosing the location of the wave, consider the domain length\n" +...
       "and the width of the computed steady solution, such that the wave \n" + ...
       "fits inside the domain\n");
fprintf("********************************************************* \n");
% ----------------------------------------------------------------------- %
%      Choose waves
% ----------------------------------------------------------------------- %
numWaves = chooseIntegerValue("How many solitary waves should be included in the simulation?",0,[],[]);

for p=1:numWaves

    if p>1
        
        for n=1:p-1
            fprintf('Wave %d: amplitude = %s, direction = %s, location = %g\n', ...
            n, string(opts.WaveAmplitude{n}), ...
            opts.WaveDirection{n}, opts.WaveLocation{n});
        end
    end
    fprintf("\n\n********************************************************* \n");
    fprintf(['\n                  Wave ' num2str(p) ': \n'])

    opts.WaveFileName{p} = sprintf('soln%d.mat', p);

    fprintf(['\nREMINDER: The maximum amplitude of a solitary wave is %g \n'], amp_max);
    fprintf("********************************************************* \n\n");
    % ---- Wave amplitude
    opts.WaveAmplitude{p} = chooseAmplitude(['\nWhat amplitude would you like for wave ' num2str(p) '?'], amp_max);

    % ---- Solve for wave
    fprintf("--- \n--- Please wait, recovering solution numerically.  \n")
    twolSteadyReturnSolution(opts.WaveAmplitude{p},opts.WaveFileName{p},H1,H2,r1,r2,g,H,optsSteady);
    load(opts.WaveFileName{p});
    fprintf('--- Wave recovered, saving as %s \n---\n', opts.WaveFileName{p});
    
    % ---- wave location and direction
    opts.WaveLocation{p} = chooseFloatValue(['\nWhere would you like wave ' num2str(p) ' placed? Inbetween [' num2str(-Lx) ',0] '], -Lx,0);
    opts.WaveDirection{p} = chooseWaveDirection(['\nWhich direction is wave' num2str(p) ' travelling? (r/l)'], "right");
    opts.WaveFileIndex{p} = []; 
end

% ----------------------------------------------------------------------- %
%      Bathymetry
% ----------------------------------------------------------------------- %
if ~opts.Boussinesq
    fprintf("To include Bathymetry, we currently only have the code for the Boussinesq case. Bathymetry switched off.\n")
    opts.Bathymetry.on = false;
else
    opts.Bathymetry.on = chooseBooleanValue("\nInclude a bottom step? (y/n) ", true);
    if opts.Bathymetry.on
        opts.Bathymetry.stepHeight = chooseFloatValue("\nHeight of step?",0,H2,H2/2);
        opts.Bathymetry.stepLocation = chooseFloatValue(['\nCenter of step? Inbetween [' num2str(-Lx) ',0] '] ,-Lx,0,-Lx/2);
        opts.Bathymetry.stepWidth = chooseFloatValue("\nWidth of step?",0,[],8*H);
        opts.Bathymetry.type = 'tanh';
end
end

% ----------------------------------------------------------------------- %
%      Final setting up
% ----------------------------------------------------------------------- %
% ---- Add a perturbation to test instability?
opts.Perturb.on = false;
opts.Perturb.epsilon = 10^(-6);
opts.Perturb.k = 40;
opts.Perturb.Width = 4;  

% ---- Filter initial guess ? 
opts.FilerInitialGuess.on = false;
opts.FilerInitialGuess.kcut = 0.2;
opts.FilerInitialGuess.order = 2;
opts.FilerInitialGuess.strength = 20;

% ---- Make a video live?
opts.VideoMaker = false;
opts.VideoMaker = chooseBooleanValue("Render video during simulation? NOTE: This may slow the computation down" ,true)


% ---- Build IG
dx=Lx/Nx;
[x,dx,dt,kx] = mesh(Lx,Nx,dt);
[w,d,x,dx,dt,kx,g,c,c0,cWave,Nx,H1,H2,H,r1,r2,opts,optsSteady] = twolIG(x,Lx,Nx,dt,opts);
q = [d;w];

% ---- Boundary condition? "wall", "periodic"
opts.BC = 'wall';

% ---- Spatial derivatives: true - finite differences, false - spectral
opts.FinDiff = true;

% ---- Differentiation matrices
switch opts.BC
    case "periodic"
        [Dx,Dxx] = FDmatricesPeriodic(Nx,dx);
    case "wall"
        [Dx,Dxx] = FDmatricesOnesided(Nx,dx);
end


% ----------------------------------------------------------------------- %
% Figure Setup
% ----------------------------------------------------------------------- %
axis tight manual % this ensures that getframe() returns a consistent size
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextInterpreter','latex')
close all

filename = 'my_video.avi';

[t_data,d_data,w_data,zeta_data] = twolIntegrate(filename,q,H1,H2,H,r1,r2,g,x,Nx,Lx,Dx,Dxx,t0,dt,tprint,tfinal,opts);
