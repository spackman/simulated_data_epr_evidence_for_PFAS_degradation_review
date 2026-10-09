%% Spin Trap EPR Simulations - Baseline Spectra Utilities

% These functions define baseline values for common adducts and utilities
% to plot them and save the data



classdef SpinTrapUtils
    % spin_trap_utils - Static methods for spin-trap EPR simulations
    % 
    % Example usage:
    %   Sys = spin_trap_utils.DMPO_OH();
    %   [B, spc] = spin_trap_utils.get_spectra(Sys);
    %   spin_trap_utils.plot_spectra(B, spc);
    
    methods(Static)
        % ------------------------
        % Helper: Gauss to MHz conversion
        % ------------------------
        function MHz = gauss_to_MHz(g)
            mu_B = 9.2740100657E-24; % J/T
            h = 6.62607015E-34;      % J/Hz
            tesla_to_gauss = 1e4;    % G/T
            Hz_to_MHz = 1e6;         % MHz/Hz
            MHz = g * (mu_B / tesla_to_gauss) / (h * Hz_to_MHz);
        end

        function newfield = change_mw_freq(startfreq, endfreq, field_mT)
            mu_B = 9.2740100657E-24; % J/T
            h = 6.62607015E-34;      % J/Hz
            tesla_to_mT = 1e3;       % mT/T
            field_T = field_mT ./ tesla_to_mT;
            g = (startfreq .* h) ./ (mu_B .* field_T);
            newfield = (endfreq .* h) ./ (mu_B .* g) .* tesla_to_mT;
        end


        %% Exemplar Adduct Parameters
        % reference values are retrieved from 
        % https://www.niehs.nih.gov/research/resources/databases/spintrap
        % unless otherwise noted

        % ------------------------
        % DMPO Spin Adducts
        % ------------------------
        function Sys = DMPO_OH_SHF(g,lwpp)
            % Long-Range Proton Hyperfine Coupling in Alicyclic Nitroxide Radicals by
            % Resolution-Enhanced Electron Paramagnetic Resonance
            % Mossoba, Makino and Riesz 1984
            if nargin < 1, g = 2.005; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,1H,1H,1H,1H,1H';
            Sys.n = [1 1 1 6 2 1 1];
            Sys.A = [14.9, 14.9, 0.227, 0.224, 0.135, 0.229, 0.370]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_OH(g,lwpp)
            if nargin < 1, g = 2.0055; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [14.9, 14.9]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPOX(g,lwpp)
            % Oxidation of spin trap 5,5-dimethyl-1-pyrroline-1-oxide in an 
            % electron paramagnetic resonance study of the reaction of 
            % methemoglobin with hydrogen peroxide
            % Mao 1994
            if nargin < 1, g = 2.005; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.n = [1 2];
            Sys.A = [7.2, 4.2]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_SO4(g,lwpp)
            if nargin < 1, g = 2.0062; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,1H';
            Sys.A = [13.84, 10.1, 1.35]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_H(g,lwpp)
            if nargin < 1, g = 2.0055; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.n = [1 2];
            Sys.A = [16.3,22.6]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_SO3(g,lwpp)
            if nargin < 1, g = 2.0054; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [14.49, 16.17]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_OOH(g,lwpp)
            if nargin < 1, g = 2.0061; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,1H';
            Sys.A = [14.3, 11.7, 1.25]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_IO3(g,lwpp)
            % Accurate identification of radicals by in-situ electron paramagnetic
            % resonance in ultraviolet-based homogenous advanced oxidation processes
            % Chen, Duan, Du, Sun, Lai, Liu 2022
            if nargin < 1, g = 2.006; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [14.2, 12.9]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_IO4(g,lwpp)
            % Accurate identification of radicals by in-situ electron paramagnetic
            % resonance in ultraviolet-based homogenous advanced oxidation processes
            % Chen, Duan, Du, Sun, Lai, Liu 2022
            if nargin < 1, g = 2.006; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [14.0, 8.0]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_ONO(g,lwpp)
            % Reactive Nitrogen Species Reactivities with Nitrones: 
            % Theoretical and Experimental Studies
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [13.8, 11.7]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_CO2(g,lwpp)
            if nargin < 1, g = 2.0058; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [15.8, 19.1]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_CO2H(g,lwpp)
            % Long-Range Proton Hyperfine Coupling in Alicyclic Nitroxide Radicals by
            % Resolution-Enhanced Electron Paramagnetic Resonance
            % Mossoba, Makino and Riesz 1984
            if nargin < 1, g = 2.0058; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,1H,1H,1H,1H';
            Sys.n = [1 1 6 2 1 1];
            Sys.A = [15.6, 18.7, 0.236, 0.130, 0.243, 0.275]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_OCO2(g,lwpp)
            % Theoretical and Experimental Studies of the Spin Trapping of 
            % Inorganic Radicals by DMPO. 2. Carbonate Radical Anion
            % Villamena 2006
            if nargin < 1, g = 2.00679; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,1H';
            Sys.A = [14.32, 10.68, 1.37]*SpinTrapUtils.gauss_to_MHz(g);
        end


           % CAREFUL - they cite a different paper, the parameters here are
           % from Janzen instead for DMPO-CH2OH!!
        function Sys = DMPO_CH2OH(g,lwpp)
            % Cautionary note for DMPO spin trapping in the presence 
            % of iron ion
            % Makino 1990
            % https://doi.org/10.1016/0006-291X(90)91556-8
            if nargin < 1, g = 2.0047; end % guess value
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [14.66, 20.67]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_CHOHCH3(g,lwpp)
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [16.05 23.2]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_Cl(g,lwpp)
            % Unveiling the mystery of chlorine radical in water by a 
            % simple methanol-based spin trapping
            % https://doi.org/10.1016/j.watres.2025.124100
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,35Cl';
            Sys.A = [11.34, 5.89, 8.14]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_CH3(g,lwpp)
            % Long-Range Proton Hyperfine Coupling in Alicyclic Nitroxide Radicals by
            % Resolution-Enhanced Electron Paramagnetic Resonance
            % Mossoba, Makino and Riesz 1984
            if nargin < 1, g = 2.0054; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,1H,1H,1H,1H,1H';
            Sys.n = [1 1 3 6 2 1 1];
            Sys.A = [16.1, 23.0, 0.473, 0.237, 0.14, 0.238, 0.302]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_CF3(g,lwpp)
            % Janzen (1969); https://doi.org/10.1016/0022-2364(73)90194-7
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,19F';
            Sys.n = [1 1 3];
            Sys.A = [13.22, 15.54, 1.01]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_CH2CH2CH3(g,lwpp)
            % DMPO-Alkyl Radical Spin Trapping: An In Situ Radiolysis Steady-State
            % ESR Study
            % Taniguchi and Madden 2000
            if nargin < 1, g = 2.00536; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [16.32, 23.46]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_CF2CF2CF2CF2CF2CF3(g,lwpp)
            % Julia (2000); https://doi.org/10.1021/jo9918605
            if nargin < 1, g = 2.0062; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [13.5, 14.5]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_F(g,lwpp)
            % Eberson (1997); https://doi.org/10.1039/A607509F
            if nargin < 1, g = 2.0061; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H,19F,1H,1H';
            Sys.n = [1 1 1 2 1];
            Sys.A = [11.8, 6.16, 54.6,2.0, 1.39]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_EPOXIDE1(g,lwpp)
            % An Artifact in the ESR Spectrum Obtained by Spin Trapping with DMPO
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [16.1, 24.5]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_EPOXIDE2(g,lwpp)
            % Electron Paramagnetic Resonance Tracks Condition-Sensitive Water Radical Cation
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [15.9, 22.9]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_DIMER(g,lwpp)
            % How Do Metal Oxides Mislead Spin-Trapping Electron Paramagnetic Resonance Analysis?
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [15.48 22.53]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_DIMER2(g,lwpp)
            % Stability of 5,5-dimethyl-1-pyrroline-N-oxide as a spin-trap 
            % for quantification of hydroxyl radicals in processes based on Fenton reaction
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N';

            Sys.A = [15.1]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_ANION(g,lwpp)
            % Oxidation of the Spin Trap 5,5-Dimethyl-1-pyrroline 
            % N-Oxide by Singlet Oxygen in Aqueous Solution
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N';
            Sys.A = [26.14]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = DMPO_R(g,lwpp)
            % general starting place for carbon centered radicals or
            % artifact species
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,1H';
            Sys.A = [16 22.6]*SpinTrapUtils.gauss_to_MHz(g);
        end

        % ------------------------
        % Non-DMPO Spin Traps
        % ------------------------
        function Sys = DETC_NO(g,lwpp)
            if nargin < 1, g = 2.004; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N';
            Sys.A = 12.7*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = PTIO(g,lwpp)
            % 2-phenyl-4,4,5,5-tetramethylimidazoline-1-oxyl 3-oxide
            % https://doi.org/10.1021/bi00054a013
            if nargin < 1, g = 2.004; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N';
            Sys.A = 12.7*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = PTIO_NO(g,lwpp)
            % 2-phenyl-4,4,5,5-tetramethylimidazoline-1-oxyl 3-oxide
            % https://doi.org/10.1021/bi00054a013
            if nargin < 1, g = 2.004; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N';
            Sys.A = 12.7*SpinTrapUtils.gauss_to_MHz(g);
        end


        % ------------------------
        % PBN Spin Traps
        % ------------------------
        function Sys = PBN_Cl(g,lwpp)
            % Song (2023); https://doi.org/10.1021/acs.est.3c02025
            % Janzen (1970; https://doi.org/10.1021/j100709a033
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,35Cl,1H';
            Sys.A = [12.11, 6.15, 0.75]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = PBN_F(g,lwpp)
            % Eberson (1992); https://doi.org/10.1039/P29920001807
            % note- product of inverted trapping
            if nargin < 1, g = 2.0023; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N,19F,1H';
            Sys.A = [12.7 45.7 10.5]*SpinTrapUtils.gauss_to_MHz(g);
        end

        function Sys = TEMPO(g, lwpp)
            % Herb (2018); https://doi.org/10.1016/j.jmr.2018.02.006
            % 45.4 MHz in MeOH
            if nargin < 1, g = 2.006; end
            if nargin < 2, lwpp = 0.1; end
            Sys.g = g; Sys.lwpp = lwpp; Sys.Nucs = '14N';
            Sys.A = 45.4;
        end


        % ------------------------
        % Simulation
        % ------------------------
        function [B, spc] = get_spectra(Sys, MW_freq, Range)
            if nargin < 2 || isempty(MW_freq), MW_freq = 9.8; end
            if nargin < 3 || isempty(Range),   Range = [344.75 354.5]; end
            Exp.mwFreq   = MW_freq;
            Exp.Range    = Range;
            Exp.nPoints  = 2048;
            Exp.Harmonic = 1; % 1st derivative CW-EPR
            [B, spc] = garlic(Sys, Exp);
        end

        function hFig = plot_spectra(B, spc, hFig)
            if nargin < 3 || isempty(hFig) || ~ishandle(hFig)
                hFig = figure; hold on;
            else
                figure(hFig); hold on;
            end
            plot(B, spc, 'LineWidth', 1.5);
            xlabel('Magnetic field (mT)'); ylabel('Intensity (a.u.)');
            box on; grid on;
        end

        % ------------------------
        % Data Cleaning
        % ------------------------
        function [Buni, spcUni] = cast_to_uniform(B, spc, nPoints)
            % Ensure numeric column vectors
            B = B(:);
            spc = spc(:);
        
            if nargin < 3 || isempty(nPoints)
                nPoints = numel(B);
            end
        
            % Sort B
            [Bsorted, idx] = sort(B);
            spcSorted = spc(idx);
        
            % --- Remove duplicates in B ---
            [Bunique, ia, ic] = unique(Bsorted);  % Bunique = sorted unique B
            spcUnique = accumarray(ic, spcSorted, [], @mean);  % average spc for duplicates
        
            % Uniformly spaced B
            Buni = linspace(min(Bunique), max(Bunique), nPoints).';
        
            % Interpolate
            spcUni = interp1(Bunique, spcUnique, Buni, 'linear', 'extrap');
        end

        function [spcCorr, baseline, BPoints, spcPoints] = baseline_correct(B, spc, Bvals)
            % Baseline correction
            % B       : field axis (Nx1)
            % spc     : spectrum (Nx1)
            % Bvals   : vector of baseline field values (original behavior)
            %           OR 2x2 array [lowMin lowMax; highMin highMax] for range averaging
            %
            % Returns:
            % spcCorr  : baseline-corrected spectrum
            % baseline : computed baseline
            % BPoints  : B-values used for fitting the baseline
            % spcPoints: corresponding spectrum values used for fitting
        
            % Ensure column vectors
            B = B(:);
            spc = spc(:);
        
            % Detect range-based input
            useRanges = (size(Bvals,1) == 2 && size(Bvals,2) == 2);
        
            if useRanges
                % Find indices within low/high ranges
                lowIdx  = find(B >= Bvals(1,1) & B <= Bvals(1,2));
                highIdx = find(B >= Bvals(2,1) & B <= Bvals(2,2));
        
                if isempty(lowIdx) || isempty(highIdx)
                    error('No points found in the specified baseline ranges.');
                end
        
                % Compute average B and spc in each range
                Bavg = [mean(B(lowIdx)); mean(B(highIdx))];
                spcAvg = [mean(spc(lowIdx)); mean(spc(highIdx))];
        
                BPoints = Bavg;
                spcPoints = spcAvg;

                nBase = 2;  % linear baseline
        
            else
                % Original Bvals behavior
                Bvals = Bvals(:);
        
                % Find nearest indices in B
                idx = zeros(numel(Bvals),1);
                for k = 1:numel(Bvals)
                    [~, idx(k)] = min(abs(B - Bvals(k)));
                end
                idx = unique(idx);  % remove duplicates
                nBase = numel(idx);
        
                if nBase < 2
                    error('At least two baseline points are required.');
                end
        
                % Sort by B
                [BPoints, order] = sort(B(idx));
                spcPoints = spc(idx(order));
            end
        
            % Compute baseline
            if nBase == 2
                % Linear baseline
                p = polyfit(BPoints, spcPoints, 1);
                baseline = polyval(p, B);
            else
                % Spline baseline
                if nBase < 4
                    warning('Few baseline points selected; spline may be unstable.');
                end
                baseline = spline(BPoints, spcPoints, B);
            end
        
            % Corrected spectrum
            spcCorr = spc - baseline;
        
            % Safety check
            assert(numel(baseline) == numel(spc), 'Baseline and spectrum must have the same length.');
        end
        
        function intSpc = integrate_spectrum(B, spc)
            intSpc = cumtrapz(B(:), spc(:));
        end
        
        function plot_full_baseline_diagnostics(B, spc, Bvals)
            % Full diagnostic figure for baseline correction and integration
            B = B(:);
            spc = spc(:);
        
            if nargin < 3 || isempty(Bvals)
                Bvals = [min(B), max(B)];
            end
        
            % Baseline correction
            [spcCorr, baseline, BPoints, spcPoints] = SpinTrapUtils.baseline_correct(B, spc, Bvals);
        
            % Integration
            intRaw = SpinTrapUtils.integrate_spectrum(B, spc);
            intCorr = SpinTrapUtils.integrate_spectrum(B, spcCorr);
        
            figure;
        
            % --- Raw + baseline ---
            subplot(3,1,1)
            plot(B, spc, 'k', 'DisplayName', 'Raw'); hold on;
            plot(B, baseline, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
            plot(BPoints, spcPoints, 'bo', 'MarkerFaceColor','b', 'DisplayName', 'Baseline points');
            ylabel('Signal');
            title('Raw spectrum and baseline');
            legend('Location','best');
            grid on;
        
            % --- Corrected spectrum ---
            subplot(3,1,2)
            plot(B, spcCorr, 'b', 'LineWidth',1.2);
            ylabel('Signal');
            title('Baseline-corrected spectrum');
            grid on;
        
            % --- Integrated spectra (dual axes) ---
            subplot(3,1,3)
            yyaxis left
            plot(B, intRaw, 'k', 'DisplayName', 'Raw integrated');
            ylim(max(abs(intRaw)) * [-1 1]);
            ylabel('Raw integrated (a.u.)');
        
            yyaxis right
            plot(B, intCorr, 'b', 'DisplayName', 'Corrected integrated');
            ylim(max(abs(intCorr)) * [-1 1]);
            ylabel('Corrected integrated (a.u.)');
        
            xlabel('Magnetic field (mT)');
            title('Integrated spectra (independent scaling)');
            grid on;
        
            legend({'Raw integrated','Corrected integrated'}, 'Location','best');
        end

        function plot_uniform(B, spc, nPoints)

            % Ensure column vectors
            B = B(:);
            spc = spc(:);
        
            [Buni, spcUni] = SpinTrapUtils.cast_to_uniform(B, spc, nPoints);
        
            figure;
            hold on;
        
            plot(B, spc, 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'Original');
            plot(Buni, spcUni, 'b', 'LineWidth', 1.5, 'DisplayName', 'Uniformly sampled');
        
            xlabel('Magnetic field (mT)');
            ylabel('Signal');
            title('Original and Uniformly Sampled Spectrum');
            legend('Location', 'best');
            grid on;
            hold off;
        end

        function residuals = plot_residuals(B, spcExp, spcFit, smooth, scaleToUnit, useOffset)
            % Plot experimental vs fit and residuals with least-squares scaling
            %
            % B          : magnetic field axis
            % spcExp     : experimental spectrum
            % spcFit     : simulated or fitted spectrum
            % smooth     : window size for Savitzky-Golay smoothing (optional, default 501)
            % scaleToUnit: logical, true to scale spcFit to [-1,1] (optional, default false)
            % useOffset  : logical, true to include intercept in scaling (optional, default true)
        
            if nargin < 6 || isempty(useOffset), useOffset = true; end
            if nargin < 5 || isempty(scaleToUnit), scaleToUnit = false; end
            if nargin < 4 || isempty(smooth), smooth = 501; end
            if mod(smooth,2)==0, smooth = smooth+1; end

        
            % Ensure column vectors
            B = B(:); spcExp = spcExp(:); spcFit = spcFit(:);
        
            % --- Ensure same length ---
            if numel(spcExp) ~= numel(spcFit)
                error('Experimental and fit spectra must have the same number of points.');
            end
        
            % --- Least-squares scaling ---
            if useOffset
                % Fit: spcExp ≈ a*spcFit + b
                X = [spcFit, ones(size(spcFit))];   % design matrix with intercept
                coeff = X \ spcExp;                 % solves least-squares
                a = coeff(1);
                b = coeff(2);
                spcFitScaled = a*spcFit + b;
            else
                % Fit: spcExp ≈ a*spcFit (no offset)
                a = (spcFit' * spcExp) / (spcFit' * spcFit);
                spcFitScaled = a * spcFit;
            end
        
            % --- Optional scaling to [-1,1] ---
            if scaleToUnit
                spcFitScaled = 2*(spcFitScaled - min(spcFitScaled)) / (max(spcFitScaled) - min(spcFitScaled)) - 1;
                spcExpUnit = 2*(spcExp - min(spcExp)) / (max(spcExp) - min(spcExp)) - 1;
            else
                spcExpUnit = spcExp;
            end
        
            % --- Residuals ---
            residuals = spcExpUnit - spcFitScaled;
        
            % Use movmean if sgolayfilt is unavailable
            if exist('sgolayfilt','file')
                resSmooth = sgolayfilt(residuals, 3, smooth);
            else
                resSmooth = movmean(residuals, smooth);
            end
        
            % --- Figure with two subplots ---
            figure;
        
            % Top: Experimental vs Fit
            subplot(2,1,1);
            plot(B, spcExpUnit, 'k', 'LineWidth', 1.0, 'DisplayName', 'Experimental'); hold on;
            plot(B, spcFitScaled, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Fit');
            xlabel('Magnetic field (mT)');
            ylabel('Signal (a.u.)');
            title('Experimental vs Fit');
            legend('Location','best');
            grid on; box on;
            xlim([min(B) max(B)])

        
            % Bottom: Residuals
            subplot(2,1,2);
            plot(B, residuals, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0); hold on;
            plot(B, resSmooth, 'k', 'LineWidth', 2.0);
            xlabel('Magnetic field (mT)');
            ylabel('Residual intensity (a.u.)');
            title('Residuals (raw and smoothed)');
            legend({'Raw residuals','Smoothed residuals'}, 'Location','best');
            grid on; box on;

            xlim([min(B) max(B)])
            ylim([min(residuals) max(residuals)])

        end



        % ------------------------
        % Data Read / Write
        % ------------------------
        
        function [B, spcRaw] = read_csv(filename)
            % Assume the csv file has two columns:
            % 'field' - magnetic field in gauss
            % 'signal' - measured signal
            gauss_to_mT = 0.1; % mT/G
            T = readtable(filename);
            B = T.field(:).*gauss_to_mT; % mT
            spcRaw = T.signal(:);
        end

        function write_spectra(filename, B, spcCell, labels)
            %WRITE_SPECTRA Write multiple EPR spectra to a single CSV file
            %
            %   B        - magnetic field vector (mT)
            %   spcCell  - cell array of spectra vectors
            %   labels   - optional cell array of column headers
        
            if ~iscell(spcCell)
                error('spcCell must be a cell array');
            end
        
            nSpec = numel(spcCell);
        
            % Default labels if not provided
            if nargin < 4 || isempty(labels)
                labels = arrayfun(@(k) sprintf('Spec_%d', k), 1:nSpec, 'UniformOutput', false);
            end
        
            if numel(labels) ~= nSpec
                error('Number of labels must match number of spectra.');
            end
        
            % Ensure column vector for B
            B = B(:);
        
            % Prepare data matrix
            data = zeros(numel(B), nSpec + 1);
            data(:,1) = B;
        
            for k = 1:nSpec
                spc = spcCell{k}(:);
        
                if numel(spc) ~= numel(B)
                    error('Spectrum %d length mismatch with B', k);
                end
       
                data(:, k + 1) = spc;
            end
        
            % Write table
            varNames = [{'B_mT'}, labels];        
            varNames = cellfun(@char, varNames, 'UniformOutput', false); 
            T = array2table(data, 'VariableNames', varNames);
            writetable(T, filename);

        
        end


    end
end
