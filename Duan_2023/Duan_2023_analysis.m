%% Duan et al. 2023
% Insight into the key role of oxygen dopants over ball-milled boron nitride 
% for efficient degradation of PFOS alternative 6:2 fluorotelomer sulfonic acid 

clc; clear;

%% Reading in the data
% Figure 3 (red)

csvPath = fullfile(pwd, 'Duan_2023\exp_fig3_red_Duan_2023.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Correct the baseline 
% (use the average values in these ranges for linear baseline)
Bvals = [349.5 349.7; 352.5 353.2];

% Plot the diagnostic plots
SpinTrapUtils.plot_full_baseline_diagnostics(BRaw, spcRaw, Bvals);
[spcCorr, baseline] = SpinTrapUtils.baseline_correct(BRaw, spcRaw, Bvals);
% Scale the baseline corrected spectra to (-1,1)
spcCorr = spcCorr ./ max(abs(spcCorr));

% Uniformly interpolate the data for fitting
nPoints = 10e3;
SpinTrapUtils.plot_uniform(BRaw, spcCorr, nPoints);
xlim([349.5,353.2]);
[BExp, spcExp] = SpinTrapUtils.cast_to_uniform(BRaw, spcCorr, nPoints);

%% Plot with reference spectra

% Oxygen Vacancy
% Identification of an Oxygen Defect in Hexagonal Boron Nitride
% Li and Gali 2022
SysVO.Nucs = "14N,14N";
SysVO.A = [82.78 86.22 151.73; 82.78 86.22 151.73];
SysVO.g = 2;
SysVO.lw = [0.8 3.3];

% Three Boron Center (TBC)
% Identification and thermal stability of point defects
% in neutron-irradiated hexagonal boron nitride (hBN)
% Toledo and Krambrock 2021
SysTBC.S = 1/2;
SysTBC.g = 1.9963;
SysTBC.lw = [0.8861 5.734];


% use the quoted MW frequency
Exp.mwFreq = 9.85; % GHz
Exp.Range = [300 400]; % mT 
Exp.nPoints = nPoints;

% reference simulations
[Bsim, spcVO] = garlic(SysVO, Exp);
[Bsim, spcTBC] = garlic(SysTBC, Exp);


%% Fitting the spectra with reported MW freq allowing g to vary and HFC

% Reset guess parameters to help the fit converge
% update the initial guess g values and reset the HFC 
SysTBCfit.S = 1/2;
SysTBCfit.g = 2.003;
SysTBCfit.lw = [0.5, 1]; % use Voight fit

% Parameters to vary
SysTBCfit_Vary.g = 0.01; 
SysTBCfit_Vary.lw = [0.5, 1];

Exp.mwFreq = 9.85; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))];
% Call esfit - fitted system defined below
esfit(spcExp, @garlic, {{SysTBCfit}, Exp},{{SysTBCfit_Vary}});

%%

% Fitted Parameters (@ 9.85 GHz)
SysTBCfitted.S = 1/2;
SysTBCfitted.g = 2.00346;
SysTBCfitted.lw = [0.00965876, 0.612186];

% Collect the fitted spectra
[~, spcTBCfit] = garlic(SysTBCfitted, Exp);

% Plot the fit and residuals
golaySmooth = 51; scaleToUnit = true;
SpinTrapUtils.plot_residuals(BExp, spcExp, spcTBCfit, golaySmooth, scaleToUnit);



%% Save the original and updated fits

dataFolder = fullfile(pwd, 'Duan_2023');  

fitcsvPath = fullfile(dataFolder, 'sim_data_Duan_2023.csv');
SpinTrapUtils.write_spectra(fitcsvPath, BExp,...
    {spcExp(:), Bsim, spcVO(:), spcTBC(:),  spcTBCfit(:)},...
    {"exp", "Bsim", "VO", "TBC", "TBCfit"})

