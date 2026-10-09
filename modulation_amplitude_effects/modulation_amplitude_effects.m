% Modulation Amplitude Effects

gValue  = 2.0023;
Sys = SpinTrapUtils.DMPO_CF3(gValue);
Sys.lwpp = 0.1;

% 0.1 G modulation amplitude
Exp1.ModAmp = 0.01;
Exp1.mwFreq = 9.8; % GHz
Exp1.Range = [344.75 354.5]; % mT 
Exp1.nPoints = 2048;

[B, spc1] = garlic(Sys, Exp1);

% 1 G modulation amplitude
Exp2.ModAmp = 0.1;
Exp2.mwFreq = 9.8; % GHz
Exp2.Range = [344.75 354.5]; % mT 
Exp2.nPoints = 2048;

[B, spc2] = garlic(Sys, Exp2);


% 2 G modulation amplitude
Exp3.ModAmp = 0.2;
Exp3.mwFreq = 9.8; % GHz
Exp3.Range = [344.75 354.5]; % mT 
Exp3.nPoints = 2048;

[B, spc3] = garlic(Sys, Exp3);


% 5 G modulation amplitude
Exp4.ModAmp = 0.5;
Exp4.mwFreq = 9.8; % GHz
Exp4.Range = [344.75 354.5]; % mT 
Exp4.nPoints = 2048;

[B, spc4] = garlic(Sys, Exp4);

%% Plot the data
spc1n = spc1 ./ max(spc1);
spc2n = spc2 ./ max(spc2);
spc3n = spc3 ./ max(spc3);
spc4n = spc4 ./ max(spc4);

figure; hold on;

offset = 1.2;  % vertical spacing between spectra

plot(B, spc1n + 0*offset, 'LineWidth', 1.2)
plot(B, spc2n + 1*offset, 'LineWidth', 1.2)
plot(B, spc3n + 2*offset, 'LineWidth', 1.2)
plot(B, spc4n + 3*offset, 'LineWidth', 1.2)

xlabel('Magnetic Field (mT)')
ylabel('Normalized Intensity (offset)')
title('Effect of Modulation Amplitude on DMPO–CF_3 Spectrum')

yticks(0:offset:3*offset)
yticklabels({'0.1 G','1 G','2 G','5 G'})

box on
set(gca,'FontSize',11)

hold off;



%% Export the data

data = [ ...
    B(:), ...
    spc1n(:), ...
    spc2n(:), ...
    spc3n(:), ...
    spc4n(:) ...
];

T = table( ...
    data(:,1), data(:,2), data(:,3), data(:,4), data(:,5), ...
    'VariableNames', { ...
        'B_mT', ...
        'ModAmp_0p1G_norm', ...
        'ModAmp_1G_norm', ...
        'ModAmp_2G_norm', ...
        'ModAmp_5G_norm' ...
    } ...
);

% Write to CSV
writetable(T, 'DMPO_CF3_Modulation_NormalizedSpectra.csv');
