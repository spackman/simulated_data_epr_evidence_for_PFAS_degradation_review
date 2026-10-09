%% Hu et al. 2023
% Hydrochar coupled with iodide for efficient photodegradation of 
% perfluorooctanoic acid and perfluorooctane sulfonic acid under 
% ultraviolet light

clc; clear;

%% Reading in the data
% Figure 4 b) (H-BC/KI+H2O2, green trace)

csvPath = fullfile(pwd, 'Hu_2023\exp_fig_4b_green_Hu_2023.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Correct the baseline 
% (use the average values in these ranges for linear baseline)
Bvals = [346 347 ; 353.5 355];

% Plot the diagnostic plots
SpinTrapUtils.plot_full_baseline_diagnostics(BRaw, spcRaw, Bvals)
[spcCorr, baseline] = SpinTrapUtils.baseline_correct(BRaw, spcRaw, Bvals);
% Scale the baseline corrected spectra to (-1,1)
spcCorr = spcCorr ./ max(abs(spcCorr));

% Uniformly interpolate the data for fitting
nPoints = 10e3;
SpinTrapUtils.plot_uniform(BRaw, spcCorr, nPoints);
[BExp, spcExp] = SpinTrapUtils.cast_to_uniform(BRaw, spcCorr, nPoints);

%% Plot with original fit (guess MW frequency)

% the MW frequency is not given in the original text
% an appropriate frequency was guessed to center the spectra
Exp.mwFreq = 9.834; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp.nPoints = nPoints;

% the original text only identifies the DMPO-OH adduct signal
SysOH = SpinTrapUtils.DMPO_OH(); % use the default DMPO data

[~, spcOH] = garlic(SysOH, Exp);

golaySmooth = 51; scaleToUnit = true;
SpinTrapUtils.plot_residuals(BExp, spcExp, spcOH, golaySmooth, scaleToUnit);

%% Plot with updated fit

Exp.mwFreq = 9.834; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp.nPoints = nPoints;
Exp.ModAmp = 0.1; % mT

% try to fit with a DMPO-OH starting point
SysOHFit1 = SpinTrapUtils.DMPO_OH(); 

% Parameters to vary
SysOHFit1Vary.g = 0.01;
SysOHFit1Vary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(SysOHFit1Vary.g);
SysOHFit1Vary.lwpp = 0.5;

% Call esfit - fitted system defined below
%esfit(spcExp, @garlic, {SysOHFit1, Exp}, {SysOHFit1Vary});

% Fitted parameters
SysOHFit1.g = 2.00501;
SysOHFit1.A = [42.9259 41.5672]; % MHz
SysOHFit1.lwpp = 0.141644;

% Compare the original A values to the fitted ones
origSys = SpinTrapUtils.DMPO_OH();

fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FIT #1 WITH 1 COMPONENT")
fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO-OH A values [G]: \t[%s]\n', ...
    num2str(origSys.A / SpinTrapUtils.gauss_to_MHz(origSys.g)));

fprintf('Fitted DMPO-OH A values [G]: \t[%s]\n', ...
    num2str(SysOHFit1.A / SpinTrapUtils.gauss_to_MHz(SysOHFit1.g)));

% Simulate with garlic
[~, spcSysOHFit1] = garlic(SysOHFit1, Exp);

% Plot residuals
golaySmooth = 51; 
scaleToUnit = false;
SpinTrapUtils.plot_residuals(BExp, spcExp, spcSysOHFit1, golaySmooth, scaleToUnit);

% There are clearly signals not explained by just the DMPO-OH adduct
% Re-fitting with the DMPO-OOH signal as well is shown below
% Small signals at the tails can be fit with a DMPO-R adduct

%% Fitting to three species (DMPO-OH, DMPO-OOH, DMPO-R)

% We now fit with three species
SysOHFit2 = SpinTrapUtils.DMPO_OH(2.005, 0.118);
SysOOHFit2 = SpinTrapUtils.DMPO_OOH(2.0051, 0.118);
SysR1Fit2 = SpinTrapUtils.DMPO_R(2.005, 0.118);
SysOHFit2.weight = 0.55;
SysOOHFit2.weight = 0.40;
SysR1Fit2.weight = 0.05;


% Parameters to vary
SysOHFit2Vary.g = 0.01; SysOOHFit2Vary.g = 0.01;  
SysR1Fit2Vary.g = 0.01; 

SysOHFit2Vary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(SysOHFit2.g);
SysOOHFit2Vary.A = [2 2 2] * SpinTrapUtils.gauss_to_MHz(SysOOHFit2.g);
SysR1Fit2Vary.A = [3 3] * SpinTrapUtils.gauss_to_MHz(SysR1Fit2Vary.g);

SysOHFit2Vary.lwpp = 0.1; SysOOHFit2Vary.lwpp = 0.1; 
SysR1Fit2Vary.lwpp = 0.1; 
SysOHFit2Vary.weight = 0.2; SysOOHFit2Vary.weight = 0.2;
SysR1Fit2Vary.weight = 0.1;


% Call esfit - fitted system defined below
%esfit(spcExp, @garlic, {{SysOHFit2, SysOOHFit2, SysR1Fit2}, Exp},...
%    {{SysOHFit2Vary, SysOOHFit2Vary, SysR1Fit2Vary}});

% Fitted Parameters
% DMPO-OH
SysOHFit2.g = 2.00503;
SysOHFit2.A = [42.5579 42.3743]; % MHz
SysOHFit3.Nucs = '14N,1H';
SysOHFit2.lwpp = 0.121145;
SysOHFit2.weight = 0.505985;

% DMPO-OOH
SysOOHFit2.g = 2.00508;
SysOOHFit2.A = [40.7008 27.967 3.94831]; % MHz
SysOOHFit2.Nucs = '14N,1H,1H';
SysOOHFit2.lwpp = 0.127773;
SysOOHFit2.weight = 0.482817;

% DMPO-R
SysRFit2.g = 2.00495;
SysRFit2.A = [44.899 63.4172]; % MHz
SysRFit2.Nucs = '14N,1H';
SysRFit2.lwpp = 0.133859;
SysRFit2.weight = 0.063277;

% Collect the fitted spectra
[~, spcOHFit2] = garlic(SysOHFit2, Exp);
[~, spcOOHFit2] = garlic(SysOOHFit2, Exp);
[~, spcRFit2] = garlic(SysRFit2, Exp);


spcFit2 =  spcOHFit2 + spcOOHFit2 + spcRFit2;
golaySmooth = 21; scaleToUnit = false;
residuals = SpinTrapUtils.plot_residuals(BExp, spcExp, spcFit2, golaySmooth, scaleToUnit);


% ---------------------------------------------------
% Compare the original A values to the fitted ones
% ---------------------------------------------------

origSysOH = SpinTrapUtils.DMPO_OH();
origSysOOH = SpinTrapUtils.DMPO_OOH();
origSysR = SpinTrapUtils.DMPO_R();

fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FIT #2 WITH 3 COMPONENTS")
fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(origSysOH.A / SpinTrapUtils.gauss_to_MHz(origSysOH.g)));

fprintf('Fitted  DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(SysOHFit2.A / SpinTrapUtils.gauss_to_MHz(SysOHFit2.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_OOH A values [G]: \t[%s]\n', ...
    num2str(origSysOOH.A / SpinTrapUtils.gauss_to_MHz(origSysOOH.g)));

fprintf('Fitted  DMPO_OOH A values [G]: \t[%s]\n', ...
    num2str(SysOOHFit2.A / SpinTrapUtils.gauss_to_MHz(SysOOHFit2.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_R A values [G]: \t[%s]\n', ...
    num2str(origSysR.A / SpinTrapUtils.gauss_to_MHz(origSysR.g)));

fprintf('Fitted  DMPO_R A values [G]: \t[%s]\n', ...
    num2str(SysRFit2.A / SpinTrapUtils.gauss_to_MHz(SysRFit2.g)));

fprintf("\n------------------------------------------------------------------\n")


%% Fitting to three additional species (DMPO-OH, DMPO-OOH, DMPO-IO3, DMPO-R)

% DMPO-OH
SysOHFit3.g = 2.00504;
SysOHFit3.A = [42.6201 42.2827]; % MHz
SysOHFit3.Nucs = '14N,1H';
SysOHFit3.lwpp = 0.11769;
SysOHFit3.weight = 0.450937;

% DMPO-OOH
SysOOHFit3.g = 2.00513;
SysOOHFit3.A = [40.7206 27.9673 3.94249]; % MHz
SysOOHFit3.Nucs = '14N,1H,1H';
SysOOHFit3.lwpp = 0.122233;
SysOOHFit3.weight = 0.4184;

% DMPO-IO3
% SysR1Fit3 = SpinTrapUtils.DMPO_IO3(2.005, 0.1);
% SysR1Fit3.weight = 0.01;
SysR1Fit3.g = 2.00465;
SysR1Fit3.A = [41.3301 35.5541]; % MHz
SysR1Fit3.Nucs = '14N,1H';
SysR1Fit3.lwpp = 0.0852144;
SysR1Fit3.weight = 0.016336;


% DMPO-R
% SysR2Fit3 = SpinTrapUtils.DMPO_R(2.005, 0.1);
% SysR2Fit3.weight = 0.01;
SysR2Fit3.g = 2.00492;
SysR2Fit3.A = [45.1939 64.0644]; % MHz
SysR2Fit3.Nucs = '14N,1H';
SysR2Fit3.lwpp = 0.15249;
SysR2Fit3.weight = 0.0594965;


% Parameters to vary
SysOHFit3Vary.g = 0.0001; SysOOHFit3Vary.g = 0.0001;  
SysR1Fit3Vary.g = 0.001; SysR2Fit3Vary.g = 0.00001;

SysOHFit3Vary.A = [0.1 0.1] * SpinTrapUtils.gauss_to_MHz(SysOHFit3.g);
SysOOHFit3Vary.A = [0.1 0.1 0.1] * SpinTrapUtils.gauss_to_MHz(SysOOHFit3.g);
SysR1Fit3Vary.A = [1 1] * SpinTrapUtils.gauss_to_MHz(SysR1Fit3.g);
SysR2Fit3Vary.A = [1 1] * SpinTrapUtils.gauss_to_MHz(SysR2Fit3.g);

SysOHFit3Vary.lwpp = 0.1; SysOOHFit3Vary.lwpp = 0.1; 
SysR1Fit3Vary.lwpp = 0.1; SysR2Fit3Vary.lwpp = 0.1;

SysOHFit3Vary.weight = 0.2; SysOOHFit3Vary.weight = 0.2;
SysR1Fit3Vary.weight = 0.1; SysR2Fit3Vary.weight = 0.1; 


% Call esfit - fitted system defined below
% esfit(spcExp, @garlic, {{SysOHFit3, SysOOHFit3, SysR1Fit3}, Exp},...
%     {{SysOHFit3Vary, SysOOHFit3Vary, SysR1Fit3Vary, SysR2Fit3Vary}});

% Collect the fitted spectra
[~, spcOHFit3] = garlic(SysOHFit3, Exp);
[~, spcOOHFit3] = garlic(SysOOHFit3, Exp);
[~, spcR1Fit3] = garlic(SysR1Fit3, Exp);
[~, spcR2Fit3] = garlic(SysR2Fit3, Exp);


spcFit3 =  spcOHFit3 + spcOOHFit3 + spcR1Fit3 + spcR2Fit3;
golaySmooth = 21; scaleToUnit = false;
residuals = SpinTrapUtils.plot_residuals(BExp, spcExp, spcFit3, golaySmooth, scaleToUnit);


% ---------------------------------------------------
% Compare the original A values to the fitted ones
% ---------------------------------------------------

origSysOH = SpinTrapUtils.DMPO_OH();
origSysOOH = SpinTrapUtils.DMPO_OOH();
origSysIO3 = SpinTrapUtils.DMPO_IO3();
origSysR = SpinTrapUtils.DMPO_R();

fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FIT #3 WITH 4 COMPONENTS")
fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(origSysOH.A / SpinTrapUtils.gauss_to_MHz(origSysOH.g)));

fprintf('Fitted DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(SysOHFit3.A / SpinTrapUtils.gauss_to_MHz(SysOHFit3.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original  DMPO_OOH A values [G]: \t[%s]\n', ...
    num2str(origSysOOH.A / SpinTrapUtils.gauss_to_MHz(origSysOOH.g)));

fprintf('Fitted  DMPO_OOH A values [G]: \t[%s]\n', ...
    num2str(SysOOHFit3.A / SpinTrapUtils.gauss_to_MHz(SysOOHFit3.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_IO3 A values [G]: \t[%s]\n', ...
    num2str(origSysIO3.A / SpinTrapUtils.gauss_to_MHz(origSysIO3.g)));

fprintf('Fitted DMPO_IO3 A values [G]: \t[%s]\n', ...
    num2str(SysR1Fit3.A / SpinTrapUtils.gauss_to_MHz(SysR1Fit3.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original  DMPO_R A values [G]: \t[%s]\n', ...
    num2str(origSysR.A / SpinTrapUtils.gauss_to_MHz(origSysR.g)));

fprintf('Fitted  DMPO_R A values [G]: \t[%s]\n', ...
    num2str(SysR2Fit3.A / SpinTrapUtils.gauss_to_MHz(SysR2Fit3.g)));

fprintf("\n------------------------------------------------------------------\n")

%% Fitting to four additional species (DMPO-OH, DMPO-OOH, DMPO-IO3, DMPO-IO4, DMPO-R)

% DMPO-OH
SysOHFit4.g = 2.00504;
SysOHFit4.A = [42.6162 42.2938]; % MHz
SysOHFit4.Nucs = '14N,1H';
SysOHFit4.lwpp = 0.107548;
SysOHFit4.weight = 0.469166;

% DMPO-OOH
SysOOHFit4.g = 2.00513;
SysOOHFit4.A = [40.7547 27.9369 3.9064]; % MHz
SysOOHFit4.Nucs = '14N,1H,1H';
SysOOHFit4.lwpp = 0.110961;
SysOOHFit4.weight = 0.430544;

% DMPO-IO3
% SysR1Fit4 = SpinTrapUtils.DMPO_IO3(2.005, 0.1);
% SysR1Fit4.weight = 0.01;
SysR1Fit4.g = 2.0046;
SysR1Fit4.A = [41.1834 35.8082]; % MHz
SysR1Fit4.Nucs = '14N,1H';
SysR1Fit4.lwpp = 0.0638492;
SysR1Fit4.weight = 0.0122268;

% DMPO-IO4
% SysR2Fit4 = SpinTrapUtils.DMPO_IO4(2.005, 0.1);
% SysR2Fit4.weight = 0.01;
SysR2Fit4.g = 2.005;
SysR2Fit4.A = [39.7679 22.7209]; % MHz
SysR2Fit4.Nucs = '14N,1H';
SysR2Fit4.lwpp = 0.113509;
SysR2Fit4.weight = 0.00136817;



% DMPO-R
% SysR3Fit4 = SpinTrapUtils.DMPO_R(2.005, 0.1);
% SysR3Fit4.weight = 0.01;
SysR3Fit4.g = 2.00492;
SysR3Fit4.A = [45.3732 63.914]; % MHz
SysR3Fit4.Nucs = '14N,1H';
SysR3Fit4.lwpp = 0.151466;
SysR3Fit4.weight = 0.065292;


% Parameters to vary
SysOHFit4Vary.g = 0.0001; SysOOHFit4Vary.g = 0.0001;  
SysR1Fit4Vary.g = 0.001; SysR2Fit4Vary.g = 0.00001; SysR3Fit4Vary.g = 0.00001;

SysOHFit4Vary.A = [0.1 0.1] * SpinTrapUtils.gauss_to_MHz(SysOHFit4.g);
SysOOHFit4Vary.A = [0.1 0.1 0.1] * SpinTrapUtils.gauss_to_MHz(SysOOHFit4.g);
SysR1Fit4Vary.A = [1 1] * SpinTrapUtils.gauss_to_MHz(SysR1Fit4.g);
SysR2Fit4Vary.A = [1 1] * SpinTrapUtils.gauss_to_MHz(SysR2Fit4.g);
SysR3Fit4Vary.A = [1 1] * SpinTrapUtils.gauss_to_MHz(SysR3Fit4.g);


SysOHFit4Vary.lwpp = 0.1; SysOOHFit4Vary.lwpp = 0.1; 
SysR1Fit4Vary.lwpp = 0.1; SysR2Fit4Vary.lwpp = 0.1;
SysR3Fit4Vary.lwpp = 0.1;

SysOHFit4Vary.weight = 0.2; SysOOHFit4Vary.weight = 0.2;
SysR1Fit4Vary.weight = 0.1; SysR2Fit4Vary.weight = 0.1; 
SysR3Fit4Vary.weight = 0.1; 


% Call esfit - fitted system defined below
% esfit(spcExp, @garlic, {{SysOHFit4, SysOOHFit4, SysR1Fit4, SysR2Fit4, SysR3Fit4}, Exp},...
%     {{SysOHFit4Vary, SysOOHFit4Vary, SysR1Fit4Vary, SysR2Fit4Vary, SysR3Fit4Vary}});

% Collect the fitted spectra
[~, spcOHFit4] = garlic(SysOHFit4, Exp);
[~, spcOOHFit4] = garlic(SysOOHFit4, Exp);
[~, spcR1Fit4] = garlic(SysR1Fit4, Exp);
[~, spcR2Fit4] = garlic(SysR2Fit4, Exp);
[~, spcR3Fit4] = garlic(SysR3Fit4, Exp);



spcFit4 =  spcOHFit4 + spcOOHFit4 + spcR1Fit4 + spcR2Fit4 +spcR3Fit4;
golaySmooth = 21; scaleToUnit = false;
residuals = SpinTrapUtils.plot_residuals(BExp, spcExp, spcFit4, golaySmooth, scaleToUnit);


% ---------------------------------------------------
% Compare the original A values to the fitted ones
% ---------------------------------------------------

origSysOH = SpinTrapUtils.DMPO_OH();
origSysOOH = SpinTrapUtils.DMPO_OOH();
origSysIO3 = SpinTrapUtils.DMPO_IO3();
origSysIO4 = SpinTrapUtils.DMPO_IO4();
origSysR = SpinTrapUtils.DMPO_R();

fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FIT #4 WITH 5 COMPONENTS")
fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(origSysOH.A / SpinTrapUtils.gauss_to_MHz(origSysOH.g)));

fprintf('Fitted DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(SysOHFit3.A / SpinTrapUtils.gauss_to_MHz(SysOHFit3.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original  DMPO_OOH A values [G]: \t[%s]\n', ...
    num2str(origSysOOH.A / SpinTrapUtils.gauss_to_MHz(origSysOOH.g)));

fprintf('Fitted  DMPO_OOH A values [G]: \t[%s]\n', ...
    num2str(SysOOHFit3.A / SpinTrapUtils.gauss_to_MHz(SysOOHFit3.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_IO3 A values [G]: \t[%s]\n', ...
    num2str(origSysIO3.A / SpinTrapUtils.gauss_to_MHz(origSysIO3.g)));

fprintf('Fitted DMPO_IO3 A values [G]: \t[%s]\n', ...
    num2str(SysR1Fit3.A / SpinTrapUtils.gauss_to_MHz(SysR1Fit3.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original DMPO_IO4 A values [G]: \t[%s]\n', ...
    num2str(origSysIO4.A / SpinTrapUtils.gauss_to_MHz(origSysIO4.g)));

fprintf('Fitted DMPO_IO4 A values [G]: \t[%s]\n', ...
    num2str(SysR2Fit4.A / SpinTrapUtils.gauss_to_MHz(SysR2Fit4.g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Original  DMPO_R A values [G]: \t[%s]\n', ...
    num2str(origSysR.A / SpinTrapUtils.gauss_to_MHz(origSysR.g)));

fprintf('Fitted  DMPO_R A values [G]: \t[%s]\n', ...
    num2str(SysR3Fit4.A / SpinTrapUtils.gauss_to_MHz(SysR3Fit4.g)));

fprintf("\n------------------------------------------------------------------\n")


%% Save the original and updated fit

dataFolder = fullfile(pwd, 'Hu_2023');  
figHandles = findall(0, 'Type', 'figure'); 
[~, idx] = sort([figHandles.Number]);
figHandles = figHandles(idx);
for k = 1:length(figHandles)
    fig = figHandles(k);
    filename = fullfile(dataFolder, sprintf('fig%d.png', k));
    saveas(fig, filename);           
end


fitcsvPath = fullfile(dataFolder, 'fitted_data_Hu_2023.csv');
SpinTrapUtils.write_spectra(fitcsvPath, BExp,...
    {spcExp(:), spcOHFit4(:), spcOOHFit4(:), spcR1Fit4, spcR3Fit4}, ...
    {"exp", "DMPO_OH", "DMPO_OOH", "DMPO_IO3", "DMPO_R"})


