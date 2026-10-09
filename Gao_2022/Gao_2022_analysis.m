%% Gao et al. 2022
% Theoretical and experimental insights into the mechanisms of C6/C6 PFPiA 
% degradation by dielectric barrier discharge plasma

clc; clear;

%% Reading in the data
% Figure 2 a) (purple)

csvPath = fullfile(pwd, 'Gao_2022\exp_fig_2a_purple_Gao_2022.csv');

% Read the raw data
[BRaw, spcRaw] = SpinTrapUtils.read_csv(csvPath);

%% Data preparation

% Correct the baseline 
% (use the average values in these ranges for linear baseline)
Bvals = [345 347; 353.1 355];

% Plot the diagnostic plots
SpinTrapUtils.plot_full_baseline_diagnostics(BRaw, spcRaw, Bvals);
[spcCorr, baseline] = SpinTrapUtils.baseline_correct(BRaw, spcRaw, Bvals);
% Scale the baseline corrected spectra to (-1,1)
spcCorr = spcCorr ./ max(abs(spcCorr));

% Uniformly interpolate the data for fitting
nPoints = 10e3;
SpinTrapUtils.plot_uniform(BRaw, spcCorr, nPoints);
[BExp, spcExp] = SpinTrapUtils.cast_to_uniform(BRaw, spcCorr, nPoints);

%% Plot with original fit

% the MW frequency is not given in the original text
% an appropriate frequency was guessed to center the spectra
% The original annotation includes only the DMPO-OO and DMPO-OH adducts
% The field regions annotated as the DMPO-OO adduct are inconsistent with
% simulated spectra

Exp.mwFreq = 9.8225; % GHz
Exp.Range = [min(BExp(:)) max(BExp(:))]; % mT 
Exp.nPoints = nPoints;

Sys1 = SpinTrapUtils.DMPO_OH(2.005);
Sys2 = SpinTrapUtils.DMPO_OOH(2.005);

[~, spc_OH] = garlic(Sys1, Exp);
[~, spc_OOH] = garlic(Sys2, Exp);


figure('Position', [100, 100, 800, 600]); 
plot(BExp, spcExp.*5000, 'DisplayName', 'raw data'); hold on;
plot(BExp, 0.25*spc_OH+8000, 'DisplayName', 'DMPO-OH', 'LineWidth',1.5); hold on;
plot(BExp, 0.25*spc_OOH-8000, 'DisplayName', 'DMPO-O$_2^-$','LineWidth',1.5); hold on;

% Define the approximate field-ranges of the original annotation
xRanges = [
    347.0, 347.5;
    348.3, 348.6;
    351.2, 351.5;
    352.6, 352.8  
];

yLimits = ylim;
yMin = yLimits(1);
yMax = yLimits(2);
grayColor = [0.5, 0.5, 0.5];

for i = 1:size(xRanges, 1)
    xBox = [xRanges(i,1), xRanges(i,2), xRanges(i,2), xRanges(i,1)];
    yBox = [yMin, yMin, yMax, yMax];
    
    if i == 1
        fill(xBox, yBox, grayColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
             'DisplayName', 'original DMPO-O$_2^-$ annotation regions');
    else
        fill(xBox, yBox, grayColor, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
    end
end

legend('Location', 'southwest', 'Interpreter', 'latex', 'FontName','sans-serif');
ylabel("Intensity [a.u.]")
xlabel("Magnetic Field [mT]")
hold off;



%% Save the original and updated fits

dataFolder = fullfile(pwd, 'Gao_2022');  
figHandles = findall(0, 'Type', 'figure'); 
[~, idx] = sort([figHandles.Number]);
figHandles = figHandles(idx);
for k = 1:length(figHandles)
    fig = figHandles(k);
    filename = fullfile(dataFolder, sprintf('fig%d.png', k));
    saveas(fig, filename);           
end


fitcsvPath = fullfile(dataFolder, 'sim_data_Gao_2022.csv');
SpinTrapUtils.write_spectra(fitcsvPath, BExp,...
    {spcExp(:), spc_OH(:), spc_OOH(:)},...
    {"exp", "DMPO_OH", "DMPO_OO"})