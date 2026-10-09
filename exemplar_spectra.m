%% Exemplar Adduct Plots
% Generates exemplar EPR spectra and writes them to CSV

clear; clc;

% --- Base exemplar list ---
exemplars = { ...
    'DMPO_OH', 'DMPO_OH_SHF', 'DMPOX', 'DMPO_OO', 'DMPO_SO4', 'DMPO_H', ...
    'DMPO_SO3', 'DMPO_OOH', 'DMPO_CO2', 'DMPO_CO2H' 'DMPO_CHOHCH3', 'DMPO_Cl', ...
    'DMPO_CH3', 'DMPO_CF3', 'DMPO_CH2CH2CH3', 'DMPO_CF2CF2CF2CF2CF2CF3', 'DMPO_F', ...
    'DETC_NO','PBN_Cl', 'PBN_F', 'DMPO_EPOXIDE2', 'DMPO_DIMER', 'DMPO_ANION', 'DMPO_DIMER2',...
    'TEMPO', 'DMPO_ONO', 'DMPO_OCO2', 'DMPO_IO3', 'DMPO_IO4'
};

gValue  = 2.0023;
lwpp    = 0.15;
MW_freq = 9.8; % GHz

spcAll  = {};
labels  = {};
Bmaster = [];

%% --- Main exemplar loop ---
for k = 1:numel(exemplars)

    spinFunc = exemplars{k};

    % Create spin system
    Sys = SpinTrapUtils.(spinFunc)(gValue, lwpp);

    % Simulate spectrum
    [B, spc] = SpinTrapUtils.get_spectra(Sys, MW_freq);

    % Store B once
    if isempty(Bmaster)
        Bmaster = B(:);
    end

    % Normalize
    spc = spc / max(abs(spc));

    % Store
    spcAll{end+1} = spc(:);
    labels{end+1} = spinFunc;

    % Plot
    hFig = figure('WindowStyle','docked');
    SpinTrapUtils.plot_spectra(B, spc, hFig);
    title(spinFunc, 'Interpreter','none');

end

%% --- Additional reduced-broadening spectra ---
lwpp_narrow = 0.05;
extraSpecies = { ...
    'DMPO_OH', 'DMPO_OH_SHF', 'DMPOX', 'DMPO_OO', 'DMPO_SO4', 'DMPO_H', ...
    'DMPO_SO3', 'DMPO_OOH', 'DMPO_CO2', 'DMPO_CO2H' 'DMPO_CHOHCH3', 'DMPO_Cl', ...
    'DMPO_CH3', 'DMPO_CF3', 'DMPO_CH2CH2CH3', 'DMPO_CF2CF2CF2CF2CF2CF3', 'DMPO_F', ...
    'DETC_NO', 'PBN_Cl', 'PBN_F', 'DMPO_EPOXIDE2', 'DMPO_DIMER', 'DMPO_ANION', 'DMPO_DIMER2', ...
    'TEMPO', 'DMPO_ONO', 'DMPO_OCO2', 'DMPO_IO3', 'DMPO_IO4'
};

for k = 1:numel(extraSpecies)

    spinFunc = extraSpecies{k};

    Sys = SpinTrapUtils.(spinFunc)(gValue, lwpp_narrow);
    [B, spc] = SpinTrapUtils.get_spectra(Sys, MW_freq);

    spc = spc / max(abs(spc));

    % Append
    spcAll{end+1} = spc(:);
    labels{end+1} = [spinFunc '_narrowLW'];

    % Plot
    hFig = figure('WindowStyle','docked');
    SpinTrapUtils.plot_spectra(B, spc, hFig);
    title([spinFunc ' (lwpp = 0.06)'], 'Interpreter','none');

end

%% --- Write all spectra to CSV ---
csvFile = 'exemplar_spectra.csv';
SpinTrapUtils.write_spectra(csvFile, Bmaster, spcAll, labels);
