%% Zhan et al. 2023
% Synergetic degradation of perfluorooctanoic acid (PFOA) in soil using 
% electrical resistance heating induced persulfate activation


clc; clear;

%% Reading in the data
% Figure6 (pink)

csvPath = fullfile(pwd, 'Zhan_2023\exp_fig6_pink_Zhan_2023.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Correct the baseline 
% (use the average values in these ranges for linear baseline)
Bvals = [342.58; 358.3];

% Plot the diagnostic plots
SpinTrapUtils.plot_full_baseline_diagnostics(BRaw, spcRaw, Bvals);
[spcCorr, baseline] = SpinTrapUtils.baseline_correct(BRaw, spcRaw, Bvals);
% Scale the baseline corrected spectra to (-1,1)
spcCorr = spcCorr ./ max(abs(spcCorr));

% Uniformly interpolate the data for fitting
nPoints = 10e3;
SpinTrapUtils.plot_uniform(BRaw, spcCorr, nPoints);
xlim([345,356]);
[BExp, spcExp] = SpinTrapUtils.cast_to_uniform(BRaw, spcCorr, nPoints);

%% Plot with original fit

% the MW frequency is not given in the original text so a guess value was
% selected

Exp.mwFreq = 9.85; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp.nPoints = nPoints;

% the original figure highlights the DMPO-OH adduct and the DMPO-O2 adduct
SysSTDB_OH = SpinTrapUtils.DMPO_OH();
SysSTDB_SO4 = SpinTrapUtils.DMPO_SO4();
SysSTDB_DIMER = SpinTrapUtils.DMPO_DIMER2(2.0055);


SysSTDB_OH.lwpp = 0.126; SysSTDB_SO4.lwpp = 0.126; SysSTDB_DIMER.lwpp = 0.126;


% simulation using reported species from STDB 

[~, spcSTDB_OH] = garlic(SysSTDB_OH, Exp);
[~, spcSTDB_SO4] = garlic(SysSTDB_SO4, Exp);
[~, spcSTDB_DIMER] = garlic(SysSTDB_DIMER, Exp);



figure('Position', [100, 100, 800, 600]); 
plot(BExp, spcExp.*5000, 'DisplayName', 'raw data'); hold on;
plot(BExp, 0.25*spcSTDB_OH+8000, 'DisplayName', 'DMPO-OH', 'LineWidth',1.5); hold on;
plot(BExp, 0.7*spcSTDB_SO4-7000, 'DisplayName', 'DMPO-SO$_4$','LineWidth',1.5); hold on;

% Note that neither DMPO-OH or DMPO-SO4 adducts explain the central 
% resonance feature that is incorrectly attributed to DMPO-SO4 in the
% original analysis

%% Fitting the spectra with DMPO-dimer artifact signal

% Reset guess parameters to help the fit converge
% update the initial guess g values and reset the HFC 
Sys_OH = SpinTrapUtils.DMPO_OH();
Sys_SO4 = SpinTrapUtils.DMPO_SO4();
Sys_DIMER = SpinTrapUtils.DMPO_DIMER2(2.005);

Sys_OH.weight = 0.5;
Sys_SO4.weight = 0.3;
Sys_DIMER.weight = 0.2;

% Parameters to vary
Sys_OH_Vary.g = 0.01; 
Sys_SO4_Vary.g = 0.01;
Sys_DIMER_Vary.g = 0.01; 

% allow deviation of up to 0.5G from reported HFC
Sys_OH_Vary.A = [1 1] * SpinTrapUtils.gauss_to_MHz(Sys_OH.g); 
Sys_SO4_Vary.A = [1 1 1] * SpinTrapUtils.gauss_to_MHz(Sys_SO4.g);
Sys_DIMER_Vary.A = [1] * SpinTrapUtils.gauss_to_MHz(Sys_DIMER.g);

Sys_OH_Vary.lwpp = 0.5; 
Sys_SO4_Vary.lwpp = 0.5;
Sys_DIMER_Vary.lwpp = 0.5; 

Sys_OH_Vary.weight = 0.3; 
Sys_SO4_Vary.weight = 0.3; 
Sys_DIMER_Vary.weight = 0.3; 


% Call esfit - fitted system defined below
%esfit(spcExp, @garlic, {{Sys_OH, Sys_SO4, Sys_DIMER}, Exp},...
%{{Sys_OH_Vary, Sys_SO4_Vary, Sys_DIMER_Vary}});

%% Resulting Fit

% Fitted Parameters (@ 9.87 GHz)
Sys_OH_fit = SpinTrapUtils.DMPO_OH();
Sys_OH_fit.g = 2.0055;
Sys_OH_fit.A = [42.07 41.6483]; % MHz
Sys_OH_fit.lwpp = 0.26232; % note this lwpp is much larger than typical
Sys_OH_fit.weight = 0.53559;

Sys_SO4_fit = SpinTrapUtils.DMPO_SO4();
Sys_SO4_fit.g = 2.0062;
Sys_SO4_fit.A = [38.8279 27.4558 4.247]; % MHz
Sys_SO4_fit.lwpp = 0.2996; % note this lwpp is much larger than typical
Sys_SO4_fit.weight = 0.06172;

Sys_DIMER_fit = SpinTrapUtils.DMPO_DIMER2();
Sys_DIMER_fit.g = 2.00549;
Sys_DIMER_fit.A = [41.2699]; % MHz
Sys_DIMER_fit.lwpp = 0.26566; % note this lwpp is much larger than typical
Sys_DIMER_fit.weight = 0.1679;


% Collect the fitted spectra
[~, spcOH_fit] = garlic(Sys_OH_fit, Exp);
[~, spcSO4_fit] = garlic(Sys_SO4_fit, Exp);
[~, spcDIMER_fit] = garlic(Sys_DIMER_fit, Exp);

% Plot the fit and residuals
spc_fit =  spcOH_fit + spcSO4_fit + spcDIMER_fit;
golaySmooth = 51; scaleToUnit = true;
SpinTrapUtils.plot_residuals(BExp, spcExp, spc_fit, golaySmooth, scaleToUnit);

% NOTE:

% the fitted linewidths are significantly larger than typical which
% suggests that the spectra is over-modulated, however the modulation
% amplitude is not reported. Lack of fit at higher fields may also suggest
% additional species not accounted for.

%%
% ---------------------------------------------------
% Compare the fitted A values [G] 
% ---------------------------------------------------
fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FIT")
fprintf("\n------------------------------------------------------------------\n")

% DMPO-OH Adduct
fprintf('Fitted  DMPO-OH A values [G]:  \t[%s]\n', ...
    num2str(Sys_OH_fit.A / SpinTrapUtils.gauss_to_MHz(Sys_OH_fit.g)));
fprintf('Database DMPO-OH A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_OH().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_OH().g)));
fprintf("\n------------------------------------------------------------------\n")

% DMPO-SO4 Adduct
fprintf('Fitted  DMPO-SO4 A values [G]: \t[%s]\n', ...
    num2str(Sys_SO4_fit.A / SpinTrapUtils.gauss_to_MHz(Sys_SO4_fit.g)));
fprintf('Database DMPO-SO4 A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_SO4().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_SO4().g)));
fprintf("\n------------------------------------------------------------------\n")

% DMPO-DIMER Species
fprintf('Fitted  DMPO-DIMER A values [G]:  \t[%s]\n', ...
    num2str(Sys_DIMER_fit.A / SpinTrapUtils.gauss_to_MHz(Sys_DIMER_fit.g)));
fprintf('Database DMPO-DIMER A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_DIMER2().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_DIMER2().g)));
fprintf("\n------------------------------------------------------------------\n")


%% Save the original and updated fits

dataFolder = fullfile(pwd, 'Zhan_2023');  

figHandles = findall(0, 'Type', 'figure'); 
[~, idx] = sort([figHandles.Number]);
figHandles = figHandles(idx);
for k = 1:length(figHandles)
    fig = figHandles(k);
    filename = fullfile(dataFolder, sprintf('fig%d.png', k));
    saveas(fig, filename);           
end

fitcsvPath = fullfile(dataFolder, 'sim_data_Zhan_2023.csv');
SpinTrapUtils.write_spectra(fitcsvPath, BExp,...
    {spcExp(:), spcOH_fit(:), spcSO4_fit(:),  spcDIMER_fit(:), spcSTDB_OH(:), spcSTDB_SO4(:), spcSTDB_DIMER(:),},...
    {"exp", "DMPO_OH_fit", "DMPO_SO4_fit", "DMPO_DIMER_fit", "DMPO_OH", "DMPO_SO4", "DMPO_DIMER"})

