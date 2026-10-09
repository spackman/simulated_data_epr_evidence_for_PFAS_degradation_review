%% Qi et al. 2024
% Novel Defluorination Pathways of Perfluoroether Compounds (GenX): 
% α-Fe2O3 Nanoparticle Layer Retains Higher Concentrations of 
% Effective Hydrated Electrons

clc; clear;

%% Reading in the data
% Figure 4 a) (green)

csvPath = fullfile(pwd, 'Qi_2024\exp_fig_4a_green_Qi_2024.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Correct the baseline 
% (use the average values in these ranges for linear baseline)
Bvals = [347.3 347.9; 354 354.5];

% Plot the diagnostic plots
SpinTrapUtils.plot_full_baseline_diagnostics(BRaw, spcRaw, Bvals);
[spcCorr, baseline] = SpinTrapUtils.baseline_correct(BRaw, spcRaw, Bvals);
% Scale the baseline corrected spectra to (-1,1)
spcCorr = spcCorr ./ max(abs(spcCorr));

% Uniformly interpolate the data for fitting
nPoints = 10e3;
SpinTrapUtils.plot_uniform(BRaw, spcCorr, nPoints);
xlim([346,356]);
[BExp, spcExp] = SpinTrapUtils.cast_to_uniform(BRaw, spcCorr, nPoints);

%% Plot with original fit

% the MW frequency is given in the original text as 9.6 GHz
% this value is inconsistent with the stnadard g values reported in the
% spin trap database for the adducts identified

Exp.mwFreq = 9.6; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp.nPoints = nPoints;

% the original figure highlights the DMPO-SO3 adduct and the aqueous
% electron (DMPO-H adduct)
SysSTDB_H = SpinTrapUtils.DMPO_H();
SysSTDB_SO3 = SpinTrapUtils.DMPO_SO3();
SysSTDB_H.weight = 0.5;
SysSTDB_SO3.weight = 0.5;

% Update the simulation range to accomodate these species at 9.6 GHz
Exp1 = Exp;
Exp1.Range = [336 356];

% we can update the hyperfine parameters to exactly match the text
% note that the spin trap database hyperfine values in gauss are
% DMPO-H: [16.3 22.6]
% DMPO-SO3: [14.49 16.17] (see SpinTrapUtils file)
SysSTDB_H.A = [16.5 22.6] * SpinTrapUtils.gauss_to_MHz(SysSTDB_H.g);
SysSTDB_SO3.A = [14.7 15.9] * SpinTrapUtils.gauss_to_MHz(SysSTDB_SO3.g);
SysSTDB_H.lwpp = 0.126; SysSTDB_SO3.lwpp = 0.126;

% simulation using reported hyperfine and g values from STDB 


% <save this for final plots>
[~, spcSTDB_H] = garlic(SysSTDB_H, Exp1);
[BOrig, spcSTDB_SO3] = garlic(SysSTDB_SO3, Exp1);



%% Fitting the spectra with reported MW freq allowing g to vary and HFC

% Reset guess parameters to help the fit converge
% update the initial guess g values and reset the HFC 
Sys_H_tightfit = SpinTrapUtils.DMPO_H();
Sys_SO3_tightfit = SpinTrapUtils.DMPO_SO3();
Sys_H_tightfit.weight = 0.5;
Sys_SO3_tightfit.weight = 0.5;
% (guess a smaller g value based on centering the spectra at 9.6 GHz)
Sys_H_tightfit.g = 1.954;
Sys_SO3_tightfit.g = 1.954;

% Parameters to vary
Sys_H_tightfit_Vary.g = 0.01; Sys_SO3_tightfit_Vary.g = 0.01;
% allow deviation of up to 0.1G from reported HFC
Sys_H_tightfit_Vary.A = [0.1 0.1] * SpinTrapUtils.gauss_to_MHz(Sys_H_tightfit.g); 
Sys_SO3_tightfit_Vary.A = [0.1 0.1] * SpinTrapUtils.gauss_to_MHz(Sys_SO3_tightfit.g);
Sys_H_tightfit_Vary.lwpp = 0.5; Sys_SO3_tightfit_Vary.lwpp = 0.5;
Sys_H_tightfit_Vary.weight = 0.2; Sys_SO3_tightfit_Vary.weight = 0.2; 

% Call esfit - fitted system defined below
% esfit(spcExp, @garlic, {{Sys_H_tightfit, Sys_SO3_tightfit}, Exp},...
% {{Sys_H_tightfit_Vary, Sys_SO3_tightfit_Vary}});

% Fitted Parameters (@ 9.6 GHz)
Sys_H_tightfit.g = 1.95434;
Sys_H_tightfit.A = [45.48006 63.2761]; % MHz
Sys_H_tightfit.lwpp = 0.151461;
Sys_H_tightfit.weight = 0.634093;

Sys_SO3_tightfit.g = 1.95445;
Sys_SO3_tightfit.A = [40.5482 45.4133]; % MHz
Sys_SO3_tightfit.lwpp = 0.432305; % note this lwpp is much larger than the DMPO-H adduct
Sys_SO3_tightfit.weight = 0.468696;

% Collect the fitted spectra
[~, spcH_tightfit] = garlic(Sys_H_tightfit, Exp);
[~, spcSO3_tightfit] = garlic(Sys_SO3_tightfit, Exp);

% Plot the fit and residuals
spc_tightfit =  spcH_tightfit + spcSO3_tightfit;
golaySmooth = 51; scaleToUnit = true;
SpinTrapUtils.plot_residuals(BExp, spcExp, spc_tightfit, golaySmooth, scaleToUnit);

% ---------------------------------------------------
% Compare the original A values to the fitted ones
% ---------------------------------------------------
fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FIT #1 WITH MW freq = 9.6, g= 1.95, VaryA = 0.1 G")
fprintf("\n------------------------------------------------------------------\n")

% note g value for original spectra not reported, using the default
fprintf('Reported DMPO_H A values [G]: \t\t\t[%s]\n', ...
    num2str([16.5 22.6]));

fprintf('Fitted  DMPO_H A values [G]: \t\t\t[%s]\n', ...
    num2str(Sys_H_tightfit.A / SpinTrapUtils.gauss_to_MHz(Sys_H_tightfit.g)));

fprintf('Spin Trap Database DMPO_H A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_H().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_H().g)));

fprintf("\n------------------------------------------------------------------\n")

% note g value for original spectra not reported, using the default
fprintf('Reported DMPO_SO3 A values [G]: \t\t[%s]\n', ...
    num2str([14.7 15.9]));

fprintf('Fitted  DMPO_SO3 A values [G]: \t\t\t[%s]\n', ...
    num2str(Sys_SO3_tightfit.A / SpinTrapUtils.gauss_to_MHz(Sys_SO3_tightfit.g)));

fprintf('Spin Trap Database DMPO_SO3 A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_SO3().A/ SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_SO3().g)));

fprintf("\n------------------------------------------------------------------\n")


%% Compare fitted spectra with reported spectra

figure;

t = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

% --- Top subplot: standard spectra ---
ax1 = nexttile;
plot(ax1, BOrig, spcSTDB_H * 0.70, 'b', 'LineWidth', 2, 'DisplayName', 'DMPO-H (as reported HFC, g=2.0055)'); hold on;
plot(ax1, BOrig, spcSTDB_SO3 * 0.44, 'r', 'LineWidth', 2, 'DisplayName', 'DMPO-SO3 (as reported HFC, g=2.0054)');
ylabel(ax1, 'Signal (a.u.)');
title(ax1, 'Reported spectra fit (9.6 GHz)');
legend(ax1, 'Location','best');
grid(ax1,'on'); box(ax1,'on');

% --- Bottom subplot: fitted spectra ---
ax2 = nexttile;
plot(ax2, BExp, spcExp*4e3, 'LineWidth', 1.5, 'DisplayName', 'Raw', 'Color', [0 0 0 0.5]);  hold on;
plot(ax2, BExp, spcH_tightfit, 'b', 'LineWidth', 3, 'DisplayName', 'DMPO-H, g=1.954'); hold on;
plot(ax2, BExp, spcSO3_tightfit, 'r', 'LineWidth', 3, 'DisplayName', 'DMPO-SO3, g= 1.954');
plot(ax2, BExp, spcH_tightfit+spcSO3_tightfit, 'c', 'LineWidth', 2, 'DisplayName', 'Fit', 'Linestyle',"--");

xlabel(ax2, 'Magnetic field (mT)');
ylabel(ax2, 'Signal (a.u.)');
title(ax2, 'Updated Spectra fit (9.6 GHz)');
legend(ax2, 'Location','best');
grid(ax2,'on'); box(ax2,'on');


% --- Link x-axes to share the same limits ---
linkaxes([ax1, ax2], 'x');
xlim([336, 356])


% NOTE:

% The fitted g values do not make sense - it seems likely that the quoted
% MW frequency is incorrect. 
% The fitted linewidths are very different too- broadening of the SO3
% adduct might be occuring to try to fit nearby peaks of different species.

% Note also that using the standard HFC values for DMPO-SO3 [14.49, 16.17]
% would reveal splitting in the two central peaks at a linewidth of ~ 1.2 G
% which is not seen in the quoted spectra due to the similar values of N
% and H HFC (15.9 - 14.7 = 1.2 G). The result merges these peaks fully.

% Re-fitting with a MW frequency of 9.8 GHz is shown below

%% Re-fit the data with appropriate MW frequency

Exp2.mwFreq = 9.85; % GHz
Exp2.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp2.nPoints = nPoints;

% We will now fit using the DMPO-OH adduct as a starting place
% instead of the DMPO-SO3 adduct (four peak signal)
Sys_OH_9_85 = SpinTrapUtils.DMPO_OH();
Sys_H_9_85 = SpinTrapUtils.DMPO_H();
Sys_OH_9_85.weight = 0.5;
Sys_H_9_85.weight = 0.5;


% Parameters to vary
Sys_OH_9_85Vary.g = 0.01; Sys_H_9_85Vary.g = 0.01;
Sys_OH_9_85Vary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(Sys_OH_9_85Vary.g);
Sys_H_9_85Vary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(Sys_H_9_85Vary.g);
Sys_OH_9_85Vary.lwpp = 0.5; Sys_H_9_85Vary.lwpp = 0.5;
Sys_OH_9_85Vary.weight = 0.2; Sys_H_9_85Vary.weight = 0.2; 

% Call esfit - fitted system defined below
% esfit(spcExp, @garlic, {{Sys_OH_9_85, Sys_H_9_85}, Exp2}, ...
% {{Sys_OH_9_85Vary, Sys_H_9_85Vary}});


% Fitted Parameters
Sys_OH_9_85.g = 2.00527;
Sys_OH_9_85.A = [42.102 41.0836]; % MHz
Sys_OH_9_85.lwpp = 0.120986;
Sys_OH_9_85.weight = 0.3000;

Sys_H_9_85.g = 2.00524;
Sys_H_9_85.A = [46.501 63.627]; % MHz
Sys_H_9_85.lwpp = 0.124514;
Sys_H_9_85.weight = 0.7000;

% Collect the fitted spectra
[~, spcOH_9_85] = garlic(Sys_OH_9_85, Exp2);
[~, spcH_9_85] = garlic(Sys_H_9_85, Exp2);

spc_9_85Fit =  spcOH_9_85 + spcH_9_85;
golaySmooth = 21; scaleToUnit = false;
residuals = SpinTrapUtils.plot_residuals(BExp, spcExp, spc_9_85Fit, golaySmooth, scaleToUnit);

% ---------------------------------------------------
% Compare the original A values to the fitted ones
% ---------------------------------------------------
fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FINAL FIT #2 WITH MW freq = 9.85, g = 2.005")
fprintf("\n------------------------------------------------------------------\n")

% note g value for original spectra not reported, using the default
fprintf('Reported DMPO_H A values [G]: \t\t\t[%s]\n', ...
    num2str([16.5 22.6]));

fprintf('Fitted  DMPO_H A values [G]: \t\t\t[%s]\n', ...
    num2str(Sys_H_9_85.A / SpinTrapUtils.gauss_to_MHz(Sys_H_9_85.g)));

fprintf('Spin Trap Database DMPO_H A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_H().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_H().g)));

fprintf("\n------------------------------------------------------------------\n")

% note g value for original spectra not reported, using the default
fprintf('Reported DMPO_OH A values [G]: \t\t[%s]\n', ...
    num2str([14.7 15.9]));

fprintf('Fitted  DMPO_OH A values [G]: \t\t\t[%s]\n', ...
    num2str(Sys_OH_9_85.A / SpinTrapUtils.gauss_to_MHz(Sys_OH_9_85.g)));

fprintf('Spin Trap Database DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_OH().A/ SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_OH().g)));

fprintf("\n------------------------------------------------------------------\n")


%% Checking residuals
% This new fit resulted in reasonable g values and accounts for most peaks,
% however the residuals indicate several unexplained resonances
% Integrating the residuals shows this clearly:

integral = SpinTrapUtils.integrate_spectrum(BExp, residuals(:));
figure;
plot(BExp, residuals, 'Color', [0 0 0 0.5], 'DisplayName', 'Raw residuals'); hold on;
plot(BExp, integral*10 , 'r',  'LineWidth', 2,  'DisplayName', 'Integrated residuals (10x)'); hold on;
xlim([346, 356]);
xlabel('Magnetic field (mT)');
ylabel('Signal (a.u.)');
title('Residuals from updated fit (9.8 GHz)');
legend('Location','best');
grid('on'); box('on');

% these peaks seem to align well with the HFC range for carbon centered
% radicals (6 peak signal) - a new species can be added to the spin system 
% and the data can be re-fitted to account for this


%% Fitting with 3 components (DMPO-OH, DMPO-H and third species)

% Define a spin system matching the residual peaks
Sys.g = 2.0047;
Sys.A = [16 23] * SpinTrapUtils.gauss_to_MHz(Sys.g);
Sys.Nucs = '14N,1H';
Sys.lwpp = 0.1275;
Sys.weight = 0.1;

% Start with DMPO-OH and DMPO-H spectra from fit above
SysOH = SpinTrapUtils.DMPO_OH(2.0052, 0.1275);
SysH = SpinTrapUtils.DMPO_H(2.0052, 0.1275);

SysOH.weight = 0.3;
SysH.weight = 0.6;


% Parameters to vary (start with DMPO-OH and DMPO-H, then the new species)
SysOHVary.g = 0.01; SysHVary.g = 0.01;
SysOHVary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(SysOHVary.g);
SysHVary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(SysHVary.g);
SysOHVary.weight = 0.2; SysHVary.weight = 0.2; 

% new species
SysVary.g = 0.01;
SysVary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(Sys.g);
SysVary.weight = 0.2; 

% Call esfit - fitted system defined below
% esfit(spcExp, @garlic, {{Sys, SysOH, SysH}, Exp2},...
%     {{SysVary, SysOHVary, SysHVary}});


% Fitted Parameters
% DMPO-OH
SysOH.g = 2.00527;
SysOH.A = [41.8768 42.0596]; % MHz
SysOH.lwpp = 0.1275;
SysOH.weight = 0.256872;

% DMPO-H
SysH.g = 2.00524;
SysH.A = [46.4905 63.5578]; % MHz
SysH.lwpp = 0.1275;
SysH.weight = 0.611937;

% Third species (DMPO-R?)
SysR.g = 2.00523;
SysR.Nucs = Sys.Nucs;
SysR.A = [45.2288 64.7035]; % MHz
SysR.lwpp = 0.1275;
SysR.weight = 0.0909232;

[~, spcFitOH] = garlic(SysOH, Exp2);
[~, spcFitH] = garlic(SysH, Exp2);
[~, spcFitR] = garlic(SysR, Exp2);


spcFitFinal =  spcFitOH + spcFitH + spcFitR;
golaySmooth = 101; scaleToUnit = false;
residuals = SpinTrapUtils.plot_residuals(BExp, spcExp, spcFitFinal, golaySmooth, scaleToUnit);

% ---------------------------------------------------
% Compare the original A values to the fitted ones
% ---------------------------------------------------
fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FINAL FIT WITH MW freq = 9.85")
fprintf("\n------------------------------------------------------------------\n")

% note g value for original spectra not reported, using the default
fprintf('Reported DMPO_H A values [G]: \t\t\t[%s]\n', ...
    num2str([16.5 22.6]));

fprintf('Fitted  DMPO_H A values [G]: \t\t\t[%s]\n', ...
    num2str(SysH.A / SpinTrapUtils.gauss_to_MHz(SysH.g)));

fprintf('Spin Trap Database DMPO_H A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_H().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_H().g)));

fprintf("\n------------------------------------------------------------------\n")

% note g value for original spectra not reported, using the default
fprintf('Reported DMPO_OH A values [G]: \t\t[%s]\n', ...
    num2str([14.7 15.9]));

fprintf('Fitted  DMPO_OH A values [G]: \t\t\t[%s]\n', ...
    num2str(SysOH.A / SpinTrapUtils.gauss_to_MHz(SysOH.g)));

fprintf('Spin Trap Database DMPO_OH A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_OH().A/ SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_OH().g)));

fprintf("\n------------------------------------------------------------------\n")

fprintf('Fitted  DMPO_R A values [G]: \t\t\t[%s]\n', ...
    num2str(SysR.A / SpinTrapUtils.gauss_to_MHz(SysR.g)));

fprintf("\n------------------------------------------------------------------\n")

%% Checking residuals from three component fit
% This new fit easily accounts for the residual peaks that appeared from
% the first two component fit, however the residuals from the three
% component fit still appear to have some intensity around the g = 2 center
% of the spectra - Re-integrating these residuals might reveal additional
% structure that could identify problems in the three component fit

integral = SpinTrapUtils.integrate_spectrum(BExp, residuals(:));
figure;
plot(BExp, residuals, 'Color', [0 0 0 0.5], 'DisplayName', 'Raw residuals'); hold on;
plot(BExp, integral*10 , 'r',  'LineWidth', 2,  'DisplayName', 'Integrated residuals (10x)'); hold on;
xlim([346, 356]);
xlabel('Magnetic field (mT)');
ylabel('Signal (a.u.)');
title('Residuals three component fit (9.8 GHz)');
legend('Location','best');
grid('on'); box('on');

% The integrated residuals don't seem to have much structure that would
% warrant another fit. The high frequency 'noise' might be a numeric
% artifact from the data digitization process - digitizing data that is
% nearly vertical (like the sharp peak transitions in the spectra) tends to
% accumulate additional error

%% Save the original and updated fits

dataFolder = fullfile(pwd, 'Qi_2024');  
figHandles = findall(0, 'Type', 'figure'); 
[~, idx] = sort([figHandles.Number]);
figHandles = figHandles(idx);
for k = 1:length(figHandles)
    fig = figHandles(k);
    filename = fullfile(dataFolder, sprintf('fig%d.png', k));
    saveas(fig, filename);           
end


fitcsvPath = fullfile(dataFolder, 'fitted_data_Qi_2024.csv');
SpinTrapUtils.write_spectra(fitcsvPath, BExp,...
    {BOrig(:), spcExp(:), spcSTDB_H(:), spcSTDB_SO3(:),  spcOH_9_85(:),  spcH_9_85(:),...
     spcFitOH(:), spcFitH(:), spcFitR(:)},...
    {"BOrig", "exp", "reported_DMPO_H", "reported_DMPO_SO3",...
    "initial_fit_DMPO_OH", "inital_fit_DMPO_H",...
    "final_fit_DMPO_OH","final_fit_DMPO_H", "final_fit_DMPO_R"})

