%% Zhang et al. 2023
% Photochemical degradation of perfluorooctanoic acid under UV irradiation
% in the presence of Fe (III)-saturated monmorillonite


clc; clear;

%% Reading in the data
% Figure 5 b) (red)

csvPath = fullfile(pwd, 'Zhang_2023\exp_fig_5b_red_zhang_2023.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Correct the baseline 
% (use the average values in these ranges for linear baseline)
Bvals = [345 346; 353.5 356];

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

% the MW frequency is given in the original text as 9.870 GHz
% this seems inconsistent with the data using the standard g values, 9.82
% is a much better fit for the DMPO-OH adduct
Exp.mwFreq = 9.87; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp.nPoints = nPoints;

% the original figure highlights the DMPO-OH adduct and the DMPO-O2 adduct
SysSTDB_OH = SpinTrapUtils.DMPO_OH();
SysSTDB_OO = SpinTrapUtils.DMPO_OOH();
SysSTDB_OH.weight = 0.6;
SysSTDB_OO.weight = 0.4;
SysSTDB_OH.lwpp = 0.126; SysSTDB_OO.lwpp = 0.126;

% simulation using reported species from STDB 

% <save this for final plots>
[~, spcSTDB_OH] = garlic(SysSTDB_OH, Exp);
[~, spcSTDB_OO] = garlic(SysSTDB_OO, Exp);
sum =  spcSTDB_OH + spcSTDB_OO;

golaySmooth = 51; scaleToUnit = true;
SpinTrapUtils.plot_residuals(BExp, spcExp, sum, golaySmooth, scaleToUnit);

%% Fitting the spectra with reported MW freq allowing g to vary and HFC

% Reset guess parameters to help the fit converge
% update the initial guess g values and reset the HFC 
Sys_OH_tightfit = SpinTrapUtils.DMPO_OH();
Sys_OO_tightfit = SpinTrapUtils.DMPO_OOH();
Sys_OH_tightfit.weight = 0.6;
Sys_OO_tightfit.weight = 0.4;
% (guess a larger g value based on centering the spectra at 9.87 GHz)
Sys_OH_tightfit.g = 2.015;
Sys_OO_tightfit.g = 2.015;

% Parameters to vary
Sys_OH_tightfit_Vary.g = 0.01; Sys_OO_tightfit_Vary.g = 0.01;
% allow deviation of up to 0.5G from reported HFC
Sys_OH_tightfit_Vary.A = [0.5 0.5] * SpinTrapUtils.gauss_to_MHz(Sys_OH_tightfit.g); 
Sys_OO_tightfit_Vary.A = [0.5 0.5 0.5] * SpinTrapUtils.gauss_to_MHz(Sys_OO_tightfit.g);
Sys_OH_tightfit_Vary.lwpp = 0.5; Sys_OO_tightfit_Vary.lwpp = 0.5;
Sys_OH_tightfit_Vary.weight = 0.3; Sys_OO_tightfit_Vary.weight = 0.3; 

% Call esfit - fitted system defined below
%esfit(spcExp, @garlic, {{Sys_OH_tightfit, Sys_OO_tightfit}, Exp},...
%{{Sys_OH_tightfit_Vary, Sys_OO_tightfit_Vary}});

%%

% Fitted Parameters (@ 9.87 GHz)
Sys_OH_tightfit.g = 2.015;
Sys_OH_tightfit.A = [42.3115 41.2809]; % MHz
Sys_OH_tightfit.lwpp = 0.106648;
Sys_OH_tightfit.weight = 0.709953;

Sys_OO_tightfit.g = 2.0162;
Sys_OO_tightfit.A = [40.3611 32.5433 3.67544]; % MHz
Sys_OO_tightfit.lwpp = 0.114976; 
Sys_OO_tightfit.weight = 0.257236;

% Collect the fitted spectra
[~, spcOH_tightfit] = garlic(Sys_OH_tightfit, Exp);
[~, spcOO_tightfit] = garlic(Sys_OO_tightfit, Exp);

% Plot the fit and residuals
spc_tightfit =  spcOH_tightfit + spcOO_tightfit;
golaySmooth = 51; scaleToUnit = true;
SpinTrapUtils.plot_residuals(BExp, spcExp, spc_tightfit, golaySmooth, scaleToUnit);

% NOTE:

% The fit is very poor and superoxide does not appear to significantly
% contribute.

% Re-fitting with a MW frequency of 9.82 GHz is shown below using a
% DMPO-OCH3 radical instead of DMPO-OO radical following Makino (1990) 
% Cautionary note for DMPO spin trapping in the presence of iron ion
% https://doi.org/10.1016/0006-291X(90)91556-8

%% Re-fit the data with appropriate MW frequency

Exp2.mwFreq = 9.825; % GHz
Exp2.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp2.nPoints = nPoints;

% We will now fit using the DMPO-OH adduct and DMPO-R instead of DMPO-OO
% The DMPOX species is also considered to account for features near the
% center of the spectra, though this may also be the result of an artifact
% signal
Sys_OH_9_82 = SpinTrapUtils.DMPO_OH();
Sys_OO_9_82 = SpinTrapUtils.DMPO_OOH();
Sys_OX_9_82 = SpinTrapUtils.DMPOX();
Sys_R_9_82 = SpinTrapUtils.DMPO_R(2.0045);
Sys_OH_9_82.weight = 0.6;
Sys_OO_9_82.weight = 0.05;
Sys_OX_9_82.weight = 0.05;
Sys_R_9_82.weight = 0.5;


% Parameters to vary
Sys_OH_9_82Vary.g = 0.01; Sys_OO_9_82Vary.g = 0.01; 
Sys_OX_9_82Vary.g = 0.01; Sys_R_9_82Vary.g = 0.01;

Sys_OH_9_82Vary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(Sys_OH_9_82Vary.g);
Sys_OO_9_82Vary.A = [2 2 2] * SpinTrapUtils.gauss_to_MHz(Sys_OO_9_82Vary.g);
Sys_OX_9_82Vary.A = [2 2] * SpinTrapUtils.gauss_to_MHz(Sys_OX_9_82Vary.g);
Sys_R_9_82Vary.A = [3 3] * SpinTrapUtils.gauss_to_MHz(Sys_R_9_82Vary.g);

Sys_OH_9_82Vary.lwpp = 0.5; Sys_OO_9_82Vary.lwpp = 0.5; 
Sys_OX_9_82Vary.lwpp = 0.5; Sys_R_9_82Vary.lwpp = 0.5;
Sys_OH_9_82Vary.weight = 0.2; Sys_OO_9_82Vary.weight = 0.2; 
Sys_OX_9_82Vary.weight = 0.2;Sys_R_9_82Vary.weight = 0.2; 

% Call esfit - fitted system defined below
esfit(spcExp, @garlic, {{Sys_OH_9_82, Sys_OO_9_82, Sys_OX_9_82 Sys_R_9_82}, Exp2}, ...
{{Sys_OH_9_82Vary, Sys_OO_9_82Vary, Sys_OX_9_82Vary, Sys_R_9_82Vary}});

%%
% Fitting with the OH, O2, DMPOX and DMPO-R species seems to give a
% resonable fit for most features, but the peaks at the edges of the
% spectra are not well fit. This may suggest the presence of more than one
% six-line species, but without further information we cannot improve the
% fit.

% System 1: DMPO-OH adduct
SysOH_9_82_fit.g = 2.00582;
SysOH_9_82_fit.A = [41.822 41.8236]; % MHz
SysOH_9_82_fit.Nucs = "14N,1H";
SysOH_9_82_fit.lwpp = 0.1187;
SysOH_9_82_fit.weight = 0.55704;

% System 2: DMPO-OOH / Superoxide adduct
SysOO_9_82_fit.g = 2.00655;
SysOO_9_82_fit.A = [40.14 32.8438 3.51028]; % MHz
SysOO_9_82_fit.Nucs = "14N,1H,1H";
SysOO_9_82_fit.lwpp = 0.125279;
SysOO_9_82_fit.weight = 0.030792;

% System 3: Oxidized / Secondary species
SysOX_9_82_fit.g = 2.00589;
SysOX_9_82_fit.A = [20.2065 11.7769]; % MHz
SysOX_9_82_fit.Nucs = "14N,1H";
SysOX_9_82_fit.n = [1 2];
SysOX_9_82_fit.lwpp = 0.176985;
SysOX_9_82_fit.weight = 0.0639963;

% System 4: Carbon-centered / Radical adduct
SysR_9_82_fit.g = 2.00575;
SysR_9_82_fit.A = [44.8865 63.398]; % MHz
SysR_9_82_fit.Nucs = "14N,1H";
SysR_9_82_fit.lwpp = 0.17119;
SysR_9_82_fit.weight = 0.59117;

% Simulate component spectra
[~, spcOH_9_82_fit] = garlic(SysOH_9_82_fit, Exp2);
[~, spcOO_9_82_fit] = garlic(SysOO_9_82_fit, Exp2);
[~, spcOX_9_82_fit] = garlic(SysOX_9_82_fit, Exp2);
[~, spcR_9_82_fit]  = garlic(SysR_9_82_fit, Exp2);

% Sum individual components into full fit
spc_9_82_fit = spcOH_9_82_fit + spcOO_9_82_fit + spcOX_9_82_fit + spcR_9_82_fit;

% Plot fit and residuals
golaySmooth = 51; 
scaleToUnit = true;
SpinTrapUtils.plot_residuals(BExp, spcExp, spc_9_82_fit, golaySmooth, scaleToUnit);

% ---------------------------------------------------
% Compare the fitted A values [G] across all 4 systems
% ---------------------------------------------------
fprintf("\n------------------------------------------------------------------\n")
fprintf(" RESULTS FROM FIT (9.82 GHz)")
fprintf("\n------------------------------------------------------------------\n")

% DMPO-OH Adduct
fprintf('Fitted  DMPO-OH A values [G]:  \t[%s]\n', ...
    num2str(SysOH_9_82_fit.A / SpinTrapUtils.gauss_to_MHz(SysOH_9_82_fit.g)));
fprintf('Database DMPO-OH A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_OH().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_OH().g)));
fprintf("\n------------------------------------------------------------------\n")

% DMPO-OOH Adduct
fprintf('Fitted  DMPO-OOH A values [G]: \t[%s]\n', ...
    num2str(SysOO_9_82_fit.A / SpinTrapUtils.gauss_to_MHz(SysOO_9_82_fit.g)));
fprintf('Database DMPO-OOH A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_OOH().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPO_OOH().g)));
fprintf("\n------------------------------------------------------------------\n")

% DMPO-OX Species
fprintf('Fitted  DMPO-OX A values [G]:  \t[%s]\n', ...
    num2str(SysOX_9_82_fit.A / SpinTrapUtils.gauss_to_MHz(SysOX_9_82_fit.g)));
fprintf('Database DMPO-OX A values [G]: \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPOX().A / SpinTrapUtils.gauss_to_MHz(SpinTrapUtils.DMPOX().g)));
fprintf("\n------------------------------------------------------------------\n")

% DMPO-R Adduct
fprintf('Fitted  DMPO-R A values [G]:   \t[%s]\n', ...
    num2str(SysR_9_82_fit.A / SpinTrapUtils.gauss_to_MHz(SysR_9_82_fit.g)));
fprintf('Database DMPO-R A values [G]:  \t[%s]\n', ...
    num2str(SpinTrapUtils.DMPO_R().A / SpinTrapUtils.gauss_to_MHz(SysR_9_82_fit.g)));
fprintf("\n------------------------------------------------------------------\n")



% Instead we elect to save the simulated spectra for potential interferring
% species for direct comparison rather than reporting the fit.

%%
% Simulate DMPO-OH, DMPO_OO and DMPO-CH2OH

Exp.mwFreq = 9.825; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp.nPoints = nPoints;

Sys_OH = SpinTrapUtils.DMPO_OH();
Sys_OO = SpinTrapUtils.DMPO_OOH();
Sys_CH2OH = SpinTrapUtils.DMPO_CH2OH();
Sys_EP1 = SpinTrapUtils.DMPO_EPOXIDE1(2.0055); % artifacct
Sys_EP2 = SpinTrapUtils.DMPO_EPOXIDE2(2.0055); % artifact
Sys_DIMER = SpinTrapUtils.DMPO_DIMER(2.0055); % artifact
Sys_R = SpinTrapUtils.DMPO_R(2.0055); % arbitrary C radical



Sys_OH.lwpp = 0.126;
Sys_OO.lwpp = 0.126;
Sys_CH2OH.lwpp = 0.126;
Sys_EP1.lwpp = 0.126;
Sys_EP2.lwpp = 0.126;
Sys_DIMER.lwpp = 0.126;
Sys_R.lwpp = 0.126;

% simulation using reported species from STDB 

% <save this for final plots>

[~, spc_OH] = garlic(Sys_OH, Exp);
[~, spc_OO] = garlic(Sys_OO, Exp);
[~, spc_CH2OH] = garlic(Sys_CH2OH, Exp);
[~, spc_EP1] = garlic(Sys_EP1, Exp);
[~, spc_EP2] = garlic(Sys_EP2, Exp);
[~, spc_DIMER] = garlic(Sys_DIMER, Exp);
[~, spc_R] = garlic(Sys_R, Exp);



%% Save the original and updated fits

dataFolder = fullfile(pwd, 'Zhang_2023');  
fitcsvPath = fullfile(dataFolder, 'sim_data_Zhang_2023.csv');
SpinTrapUtils.write_spectra(fitcsvPath, BExp,...
    {spcExp(:), spc_OH(:), spc_OO(:),  spc_CH2OH(:), spc_EP1(:), spc_EP2(:), spc_DIMER(:), spc_R(:)},...
    {"exp", "DMPO_OH", "DMPO_OO", "DMPO_CH2OH", "DMPO_EP1", "DMPO_EP2", "DMPO_DIMER","DMPO_R"})

