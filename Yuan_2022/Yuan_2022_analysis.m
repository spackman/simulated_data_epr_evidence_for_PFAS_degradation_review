%% Yuan et al. 2022
% Nitrate promoted defluorination of perfluorooctanoic acid in UV/sulfite 
% system: Coupling hydrated electron/reactive nitrogen species-mediated 
% reduction and oxidation

clc; clear;

% For this paper, both the reported spectra and the reported spectra from
% the reference cited for the hyperfine parameters are plotted. The
% original spectra from Yuan (2022) will not be baseline corrected, but the
% referenced spectra from Lian (2021) will be baseline corrected so that
% the fit comparison is charitable.

% The referenced spectra is from Lian et. al (2022)
% Photochemical Transformation and Catalytic Activity of Dissolved
% Black Nitrogen Released from Environmental Black Carbon

%% Reading in the data
% Figure 3 c) Yuan (2022)

csvPath = fullfile(pwd, 'Yuan_2022\exp_fig_3c_Yuan_2022.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Uniformly interpolate the data for fitting
nPoints = 10e3;
SpinTrapUtils.plot_uniform(BRaw, spcRaw, nPoints);
[BExpYuan, spcExpYuan] = SpinTrapUtils.cast_to_uniform(BRaw, spcRaw, nPoints);


%% Reading in the data
% Figure 5 d) red Lian (2021)

csvPath = fullfile(pwd, 'Yuan_2022\exp_fig_5d_red_Lian_2021.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Correct the baseline 
% (use the average values in these ranges for linear baseline)
Bvals = [346 348; 354 356];

% Plot the diagnostic plots
SpinTrapUtils.plot_full_baseline_diagnostics(BRaw, spcRaw, Bvals)
[spcCorr, baseline] = SpinTrapUtils.baseline_correct(BRaw, spcRaw, Bvals);
% Scale the baseline corrected spectra to (-1,1)
spcCorr = spcCorr ./ max(abs(spcCorr));

% Uniformly interpolate the data for fitting
nPoints = 10e3;
SpinTrapUtils.plot_uniform(BRaw, spcCorr, nPoints);
[BExpLian, spcExpLian] = SpinTrapUtils.cast_to_uniform(BRaw, spcCorr, nPoints);



%% Fit Lian data

Exp.mwFreq = 9.85; % GHz
Exp.ModAmp = 0.1; % mT
Exp.Range = [min(BExpLian(:)) max(BExpLian(:))]; % mT 
Exp.nPoints = nPoints;

SysLian = SpinTrapUtils.DETC_NO(2.005, 0.6); % use the default DETC NO data

% Parameters to vary
SysLianVary.g = 0.01;
SysLianVary.A = 2 * SpinTrapUtils.gauss_to_MHz(SysLian.g);
SysLianVary.lwpp = 0.5;

% Call esfit - fitted system defined below
% esfit(spcExpLian, @garlic, {SysLian, Exp}, {SysLianVary});

% Fitted params
SysLian.g = 2.00571;
SysLian.A = 30.3148; % MHz
SysLian.lwpp = 0.634678;

[~, spcFitLian_9_85] = garlic(SysLian, Exp);

fprintf('Fitted A value [G]: \t[%s]\n', ...
    num2str(SysLian.A / SpinTrapUtils.gauss_to_MHz(SysLian.g)));

fprintf('Reported A value [G]: \t[%s]\n', ...
    num2str(16.59));

%% Plot Yuan spectra and Lian spectra with fit

figure;
plot(BExpYuan, spcExpYuan, 'Color', [0.27 0.51 0.71], 'DisplayName', 'Yuan'); hold on;
plot(BExpLian, spcExpLian*0.25, 'Color', [1.00 0.50 0.31], 'DisplayName', 'Lian'); hold on;
plot(BExpLian, spcFitLian_9_85*0.04, 'k', 'DisplayName', 'Fitted'); hold on;

legend show

%% Re-simulate Lian data approximating the Yuan MW frequency

% Shift the Lian experimental data to match a MW frequency of 9.66 GHz
BExpLian_9_66 = SpinTrapUtils.change_mw_freq(9.85, 9.66, BExpLian);

Exp2.mwFreq =  9.66; % GHz
Exp2.ModAmp = 0.1; % mT
Exp2.Range = [min(BExpLian_9_66(:)) max(BExpLian_9_66(:))]; % mT 
Exp2.nPoints = nPoints;
[~, spcFitLian_9_66] = garlic(SysLian, Exp2);


figure;
plot(BExpYuan, spcExpYuan, 'Color', [0.27 0.51 0.71], 'DisplayName', 'Yuan'); hold on;
plot(BExpLian_9_66, spcExpLian*0.25, 'Color', [1.00 0.50 0.31], 'DisplayName', 'Lian (9.66 GHz)'); hold on;
plot(BExpLian_9_66, spcFitLian_9_66*0.04, 'k', 'DisplayName', 'Fitted (9.66 GHz)'); hold on;

legend show

%% Simulate and plot the quoted parameters

SysReported.g = 2.0054;
SysReported.Nucs = '14N';
SysReported.A = 16.59 * SpinTrapUtils.gauss_to_MHz(SysReported.g);
SysReported.lwpp = SysLian.lwpp;

[~, spcReported] = garlic(SysReported, Exp2);

figure;
plot(BExpYuan, spcExpYuan, 'Color', [0.27 0.51 0.71], 'DisplayName', 'Yuan'); hold on;
plot(BExpLian_9_66, spcExpLian*0.25, 'Color', [1.00 0.50 0.31], 'DisplayName', 'Lian (9.66 GHz)'); hold on;
plot(BExpLian_9_66, spcFitLian_9_66*0.04, 'k', 'DisplayName', 'Fitted (9.66 GHz)'); hold on;
plot(BExpLian_9_66, spcReported*0.04, 'red', 'DisplayName', 'Reported'); hold on;

legend show


%% Save the original and updated fit

dataFolder = fullfile(pwd, 'Yuan_2022');  
figHandles = findall(0, 'Type', 'figure'); 
[~, idx] = sort([figHandles.Number]);
figHandles = figHandles(idx);
for k = 1:length(figHandles)
    fig = figHandles(k);
    filename = fullfile(dataFolder, sprintf('fig%d.png', k));
    saveas(fig, filename);           
end


fitcsvPath = fullfile(dataFolder, 'fitted_data_Yuan_2022.csv');
    SpinTrapUtils.write_spectra(fitcsvPath, BExpYuan, ...
    {spcExpYuan(:), BExpLian(:), spcExpLian(:), BExpLian_9_66(:),...
    spcFitLian_9_66(:), spcReported(:)}, {"spcExpYuan", ...
    "BExpLian", "spcExpLian", "BExpLian_9_66", "spcFitLian_9_66", "spcReported"});