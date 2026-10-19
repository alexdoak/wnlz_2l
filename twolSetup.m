clear
close all
% ---- Add path to aux files (the folder this script lives in)
codeDir = fileparts(mfilename('fullpath'));
if isempty(codeDir)
    codeDir = pwd;
end
addpath(genpath(codeDir))

clc
format short

% ---- Validators used when taking values from a configuration
isPos     = @(x) isnumeric(x) && isscalar(x) && x > 0;
isNonNeg  = @(x) isnumeric(x) && isscalar(x) && x >= 0;
isBool    = @(b) isscalar(b) && (islogical(b) || (isnumeric(b) && any(b == [0 1])));
isBetween = @(x, lo, hi) isnumeric(x) && isscalar(x) && x > lo && x < hi;   % strict, like the choose* bounds

% ----------------------------------------------------------------------- %
%      Introduction
% ----------------------------------------------------------------------- %
fprintf( "********************************************************* \n" + ...
        "This code accompanies the paper (reference paper here). \n" + ...
        "It solves the two-layer strongly nonlinear Boussinesq-type  \n" + ...
        "system. It can include any number of solitary waves, as \n" +...
        "well as bathymetry, either a tanh step or a Gaussian. \n" + ...
        "If you wish, you can explore the file 'twol.m' to setup  \n" +...
        "parameters faster. \n\n" +...
        "Otherwise, continue through the prompts to set up a simulation \n" +...
         "********************************************************* \n");

input("Press Enter to continue \n\n", 's');

% ----------------------------------------------------------------------- %
%      Optional pre-made configuration
% ----------------------------------------------------------------------- %
cfgDir    = fullfile(codeDir, 'configs');
cfgFiles  = [dir(fullfile(cfgDir, '*.m')); dir(fullfile(cfgDir, '*.json'))];

% ---- Sort by the first number in each name (fig9 before fig10), then alphabetically
[~, idx] = sort(lower({cfgFiles.name}));
cfgFiles = cfgFiles(idx);
figNum = cellfun(@(s) str2double(regexp(s, '\d+', 'match', 'once')), {cfgFiles.name});
[~, idx] = sort(figNum);              % stable, so 7a stays before 7b
cfgFiles = cfgFiles(idx);

cfg       = struct();
cfgLoaded = false;

if ~isempty(cfgFiles)
    fprintf( "********************************************************* \n" + ...
            "Pre-made configurations are available. Loading one fills in \n" + ...
            "the parameters below; anything it does not specify will \n" + ...
            "still be asked for. \n\n");
    fprintf("  0) none (interactive setup)\n");
    for k = 1:numel(cfgFiles)
        fprintf("  %d) %s\n", k, cfgFiles(k).name);
    end
    fprintf("********************************************************* \n");

    % chooseIntegerValue's bounds are strict, so check the range here instead
    choice = -1;
    while choice < 0 || choice > numel(cfgFiles)
        choice = chooseIntegerValue("Which configuration would you like to load?", [], [], 0);
        if choice < 0 || choice > numel(cfgFiles)
            fprintf("Please enter a number from 0 to %d.\n", numel(cfgFiles));
        end
    end
    if choice > 0
        cfg = loadConfig(fullfile(cfgDir, cfgFiles(choice).name));
        cfgLoaded = true;
        fprintf("\nLoaded configuration %s \n\n", cfgFiles(choice).name);
    end
end

% ----------------------------------------------------------------------- %
%      Choose parameters
% ----------------------------------------------------------------------- %
fprintf( "********************************************************* \n" + ...
        "First, we will choose some parameter values of the  \n" + ...
        "background state. You will be asked for the total depth (H), \n" + ...
        "and domain length (Lx). These form the computational domain. \n" + ...
        "Upon choosing densities, a choice of a lower \n" +...
        "layer depth will be requested. Since this code explores solitary   \n" +...
        "wave interaction, please choose a H2 sufficiently far from  \n" + ...
        "the criticality condition for which no solitary waves exist \n" +...
        " (the criticality condition will be given) \n" + ...
         "********************************************************* \n");


% ---- Physical parameters
H  = getParam(cfg, 'H', ...
        @() chooseFloatValue("Value for total depth, H? ",0,[],1), isPos);

Lx = getParam(cfg, 'Lx', ...
        @() chooseFloatValue("Value of domain length, Lx, x=[-Lx,0]? ", 0, [],200*H), isPos);

r1 = getParam(cfg, 'r1', ...
        @() chooseFloatValue("Value for upper layer density, r1? ", 0,[],1), isPos);
r2 = getParam(cfg, 'r2', ...
        @() chooseFloatValue("Value for lower layer density, r2? ", r1,[],1.01), ...
        @(x) isPos(x) && x > r1);

% ---- Boussinesq approximation
optsSteady.Boussinesq = getParam(cfg, 'Boussinesq', ...
        @() chooseBooleanValue("Boussinesq approximation on densities? (y/n) ", true), isBool);
optsSteady.Boussinesq = logical(optsSteady.Boussinesq);
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
H2 = getParam(cfg, 'H2', ...
        @() chooseFloatValue("Value of lower layer depth, H2? ",0,H), ...
        @(x) isBetween(x, 0, H) && x > 0 && x < H);
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
g = getParam(cfg, 'g', ...
        @() chooseFloatValue('Value for gravity g? The default is to choose g such that\nthe long-wave speed c0 = H. (recommended)', 0,[],H^2/c0^2), ...
        isPos);
if opts.Boussinesq
    c0 = (g*(r2-r1)*H1*H2/(H2+H1))^.5;
else
    c0 = (g*(r2-r1)*H1*H2/(r1*H2+r2*H1))^.5;
end

% ---- Mesh spacing and time step (depends on velocity scale)
defaultMesh = getParam(cfg, 'defaultMesh', ...
        @() chooseBooleanValue("Would you like to use default spatial and temporal spacing (recommended)?",true), ...
        isBool);
if defaultMesh
    dx = 0.02*H;
    Nx = round(Lx/dx)+1;
    t0 = 0;
    dt = dx/c0;
    tprintSteps = 20;
    tprint = tprintSteps*dt;
    tfinal = Lx;
else
    dx = getParam(cfg, 'dx', ...
            @() chooseFloatValue("Value for mesh spacing? ", 0,[],0.02*H), isPos);
    Nx = round(Lx/dx)+1;
    fprintf('PLEASE NOTE, if Lx/dx is not an integer, dx will be minimally modified to make this condition true \n\n')
    t0 = 0;
    dt = getParam(cfg, 'dt', ...
            @() chooseFloatValue("Value for time step dt? Default is to choose dt such that dt=dx/c0 ", 0,[],dx/c0), isPos);
    tprintSteps = getParam(cfg, 'tprintSteps', ...
            @() chooseFloatValue("How many time steps to pass before storing data? Default is 20", 0, [], 20), isPos);
    tprint = tprintSteps * dt;
    tfinal = getParam(cfg, 'tfinal', ...
            @() chooseFloatValue("What to set as the final time?", 0, [], Lx), isPos);
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
% ---- Number of waves: taken from the configuration's wave list, if given
if isfield(cfg, 'Waves')
    cfg.numWaves = numel(cfg.Waves);
end
numWaves = getParam(cfg, 'numWaves', ...
        @() chooseIntegerValue("How many solitary waves should be included in the simulation?",0,[],[]), ...
        @(n) isNonNeg(n) && n == round(n));

% ---- Validators for wave entries
isValidAmp = @(a) ~isnumeric(a) || ...   % non-numeric (e.g. "max") passed on as-is
        (isscalar(a) && sign(a) == sign(amp_max) && abs(a) <= abs(amp_max));
isValidDir = @(d) (ischar(d) || isstring(d)) && any(strcmpi(string(d), ["right", "left"]));

for k=1:numWaves

    % ---- Configuration entries for this wave, if any
    waveCfg = struct();
    if isfield(cfg, 'Waves') && k <= numel(cfg.Waves)
        if iscell(cfg.Waves)          % jsondecode gives a cell if waves differ in fields
            waveCfg = cfg.Waves{k};
        else
            waveCfg = cfg.Waves(k);
        end
    end

    if k>1
        for n=1:k-1
            fprintf('Wave %d: amplitude = %s, direction = %s, location = %g\n', ...
            n, string(opts.WaveAmplitude{n}), ...
            opts.WaveDirection{n}, opts.WaveLocation{n});
        end
    end
    fprintf("\n\n********************************************************* \n");
    fprintf(['\n                  Wave ' num2str(k) ': \n'])

    opts.WaveFileName{k} = sprintf('soln%d.mat', k);

    fprintf(['\nREMINDER: The maximum amplitude of a solitary wave is %g \n'], amp_max);
    fprintf("********************************************************* \n\n");
    % ---- Wave amplitude
    opts.WaveAmplitude{k} = getParam(waveCfg, 'Amplitude', ...
            @() chooseAmplitude(['\nWhat amplitude would you like for wave ' num2str(k) '?'], amp_max), ...
            isValidAmp);

    % ---- Solve for wave
    fprintf("--- \n--- Please wait, recovering solution numerically.  \n")
    twolSteadyReturnSolution(opts.WaveAmplitude{k},opts.WaveFileName{k},H1,H2,r1,r2,g,H,optsSteady);
    % Loading into a struct keeps the solution file's variables from
    % overwriting anything in this workspace (e.g. H, cfg).
    waveSoln = load(opts.WaveFileName{k}); %#ok<NASGU>
    fprintf('--- Wave recovered, saving as %s \n---\n', opts.WaveFileName{k});

    % ---- wave location and direction
    opts.WaveLocation{k} = getParam(waveCfg, 'Location', ...
            @() chooseFloatValue(['\nWhere would you like wave ' num2str(k) ' placed? Inbetween [' num2str(-Lx) ',0] '], -Lx,0), ...
            @(x) isBetween(x, -Lx, 0));
    opts.WaveDirection{k} = lower(string(getParam(waveCfg, 'Direction', ...
            @() chooseWaveDirection(['\nWhich direction is wave' num2str(k) ' travelling? (r/l)'], "right"), ...
            isValidDir)));
    opts.WaveFileIndex{k} = [];
end

% ----------------------------------------------------------------------- %
%      Bathymetry
% ----------------------------------------------------------------------- %
if ~opts.Boussinesq
    fprintf("To include Bathymetry, we currently only have the code for the Boussinesq case. Bathymetry switched off.\n")
    opts.Bathymetry.on = false;
else
    opts.Bathymetry.on = getParam(cfg, 'Bathymetry.on', ...
            @() chooseBooleanValue("\nInclude a bottom step? (y/n) ", true), isBool);
    if opts.Bathymetry.on
        opts.Bathymetry.stepHeight = getParam(cfg, 'Bathymetry.stepHeight', ...
                @() chooseFloatValue("\nHeight of step?",0,H2,H2/2), ...
                @(x) isBetween(x, 0, H2));
        opts.Bathymetry.stepLocation = getParam(cfg, 'Bathymetry.stepLocation', ...
                @() chooseFloatValue(['\nCenter of step? Inbetween [' num2str(-Lx) ',0] '] ,-Lx,0,-Lx/2), ...
                @(x) isBetween(x, -Lx, 0));
        opts.Bathymetry.stepWidth = getParam(cfg, 'Bathymetry.stepWidth', ...
                @() chooseFloatValue("\nWidth of step?",0,[],8*H), isPos);
        opts.Bathymetry.type = 'tanh';
    end
end


% ----------------------------------------------------------------------- %
%      Final setting up
% ----------------------------------------------------------------------- %
% ---- Perturbation to test instability (expert option: config file only)
opts.Perturb = struct('on', false, 'epsilon', 1e-6, 'k', 40, 'Width', 4, 'Location', opts.WaveLocation{1});
if isfield(cfg, 'Perturb')
    fn = fieldnames(cfg.Perturb);
    for i = 1:numel(fn)
        if isfield(opts.Perturb, fn{i})
            opts.Perturb.(fn{i}) = cfg.Perturb.(fn{i});
        else
            warning('Unknown field Perturb.%s in configuration, ignored.', fn{i});
        end
    end
    if opts.Perturb.on
        fprintf("  Perturbation ON (from configuration): epsilon = %g, k = %g, Width = %g\n", ...
            opts.Perturb.epsilon, opts.Perturb.k, opts.Perturb.Width);
    end
end

% ---- Filter initial guess ?
opts.FilerInitialGuess.on = false;
opts.FilerInitialGuess.kcut = 0.2;
opts.FilerInitialGuess.order = 2;
opts.FilerInitialGuess.strength = 20;

% ---- Make a video live?
opts.VideoMaker = getParam(cfg, 'VideoMaker', ...
        @() chooseBooleanValue("Render video during simulation? NOTE: This may slow the computation down" ,true), ...
        isBool);

% ----------------------------------------------------------------------- %
%      Offer to save these choices as a configuration
% ----------------------------------------------------------------------- %
if ~cfgLoaded
    cfgOut = struct();
    cfgOut.H          = H;
    cfgOut.Lx         = Lx;
    cfgOut.r1         = r1;
    cfgOut.r2         = r2;
    cfgOut.Boussinesq = opts.Boussinesq;
    cfgOut.H2         = H2;
    cfgOut.g          = g;
    cfgOut.defaultMesh = defaultMesh;
    if ~defaultMesh
        cfgOut.dx          = dx;
        cfgOut.dt          = dt;
        cfgOut.tprintSteps = tprintSteps;
        cfgOut.tfinal      = tfinal;
    end
    if numWaves > 0
        cfgOut.Waves = struct('Amplitude', opts.WaveAmplitude, ...
                              'Location',  opts.WaveLocation, ...
                              'Direction', opts.WaveDirection);
    else
        cfgOut.Waves = [];
    end
    cfgOut.Bathymetry = opts.Bathymetry;
    cfgOut.VideoMaker = opts.VideoMaker;

    saveName = strtrim(input("\nSave these choices as a configuration? " + ...
        "Enter a name, or press Enter to skip: ", 's'));
    if ~isempty(saveName)
        [~, saveName] = fileparts(saveName);          % drop any extension given
        if ~isfolder(cfgDir)
            mkdir(cfgDir);
        end
        saveFile = fullfile(cfgDir, [saveName '.json']);
        fid = fopen(saveFile, 'w');
        fprintf(fid, '%s', jsonencode(cfgOut, 'PrettyPrint', true));
        fclose(fid);
        fprintf("Configuration saved to %s \n\n", saveFile);
    end
end

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