%% Tensile stress–strain plots (separate figures) + tail trimming + 0.2% offset yield
% Uses files in the current working folder:
%   "first test 1045 full repot.txt"
%   "second test Al full report.txt"
%
% Expected columns (after some header lines):
% Time (s) | Load (N) | Machine Extension (mm) | Stress (MPa) | Elongation (mm) | Percentage Elongation
%
% Output:
%   - Two separate figures (one per material), showing:
%       Data (trimmed), 0.2% offset line, Peak stress point, 0.2% offset yield point
%   - Prints per-material results to Command Window
%   - Writes per-material CSV results to working folder
%   - Saves per-material PNG figures (300 dpi)

clear; clc; close all;

files  = { "first test 1045 full repot.txt", "second test Al full report.txt" };
labels = { "SAE-1045", "Aluminum 6061" };

% ---------------- User-tunable parameters ----------------
trimEpsilonTol    = 1e-6;   % tolerance for removing locally decreasing strain points
initialMaxStrain  = 0.002;  % region used to estimate slope E (e.g., 0.2% strain)
minPointsForFit   = 50;     % minimum points for E estimation
smoothWindow      = 0;      % 0 = no smoothing; else e.g., 21 (odd) for movmean
offset            = 0.002;  % 0.2% engineering strain offset
% ---------------------------------------------------------

for i = 1:numel(files)
    fn = files{i};
    if ~isfile(fn)
        error('File not found in working folder: %s', fn);
    end

    % Robust import: auto-skip header, tolerate multiple delimiters
    opts = detectImportOptions(fn, ...
        'FileType','text', ...
        'Delimiter', {'\t',' ',';',' ,'}, ...
        'MultipleDelimsAsOne', true);

    varNames = string(opts.VariableNames);
    colStress = find(contains(lower(varNames), "stress"), 1);
    colPercEl = find(contains(lower(varNames), "percentage") & contains(lower(varNames), "elong"), 1);

    % Fallback to fixed positions: Stress=4, Percentage Elongation=6
    if isempty(colStress); colStress = 4; end
    if isempty(colPercEl); colPercEl = 6; end

    opts.SelectedVariableNames = opts.VariableNames([colStress, colPercEl]);
    T = readtable(fn, opts);

    stress = double(T{:,1});  % [MPa]
    percEl = double(T{:,2});  % [%]

    % Clean
    ok = isfinite(stress) & isfinite(percEl);
    stress = stress(ok);
    percEl = percEl(ok);

    % Engineering strain
    strain = percEl / 100;

    % Optional smoothing
    if smoothWindow > 1
        strain_sm = movmean(strain, smoothWindow);
        stress_sm = movmean(stress, smoothWindow);
    else
        strain_sm = strain;
        stress_sm = stress;
    end

    % Shift strain to start at 0 (handles small negatives)
    strain_sm = strain_sm - strain_sm(1);

    % ---------- Trim tail where strain decreases ----------
    % Keep up to the point where maximum strain is achieved
    cumMax = cummax(strain_sm);
    [~, idxTrim] = max(cumMax);

    strain_use = strain_sm(1:idxTrim);
    stress_use = stress_sm(1:idxTrim);

    % Remove any remaining local decreases beyond tolerance
    dEps = diff(strain_use);
    keep = [true; dEps >= -trimEpsilonTol];
    strain_use = strain_use(keep);
    stress_use = stress_use(keep);

    % ---------- Peak stress ----------
    [sigmaMax, idxMax] = max(stress_use);
    epsAtMax = strain_use(idxMax);

    % ---------- Estimate E for 0.2% offset line (not plotted as "initial fit") ----------
    idxInit = find(strain_use <= initialMaxStrain);
    if numel(idxInit) < minPointsForFit
        idxInit = 1:min(minPointsForFit, numel(strain_use));
    end
    p = polyfit(strain_use(idxInit), stress_use(idxInit), 1); % sigma = E*eps + b
    E_est = p(1);  % [MPa]
    b_est = p(2);

    % ---------- 0.2% Offset Yield (proof stress) ----------
    sigma_offset_line = E_est .* (strain_use - offset) + b_est;
    diffY = stress_use - sigma_offset_line;

    idxY = find(diffY(1:end-1) < 0 & diffY(2:end) >= 0, 1, 'first');
    if isempty(idxY)
        [~, idxY] = min(abs(diffY));
        epsY = strain_use(idxY);
        sigmaY = stress_use(idxY);
    else
        % Linear interpolation for crossing
        x1 = strain_use(idxY);     x2 = strain_use(idxY+1);
        y1 = diffY(idxY);          y2 = diffY(idxY+1);

        epsY = x1 - y1 * (x2 - x1) / (y2 - y1);
        sigmaY = interp1(strain_use([idxY idxY+1]), stress_use([idxY idxY+1]), epsY, 'linear');
    end

    % ---------- Area under curve ----------
    areaUnder = trapz(strain_use, stress_use); % [MPa] = [MJ/m^3]

    % ---------- Plot (separate figure) ----------
    figure('Color','w');
    hold on; grid on;

    % Data
    plot(strain_use, stress_use, 'LineWidth', 1.5, 'DisplayName', 'Data (trimmed)');

    % 0.2% offset line
    strainOffLine = linspace(0, max(strain_use), 300)';
    stressOffLine = E_est .* (strainOffLine - offset) + b_est;
    plot(strainOffLine, stressOffLine, ':', 'LineWidth', 1.2, 'DisplayName', '0.2% offset line');

    % Peak stress marker
    plot(epsAtMax, sigmaMax, 's', 'MarkerSize', 7, 'LineWidth', 1.5, 'DisplayName', 'Peak stress');

    % 0.2% offset yield marker
    plot(epsY, sigmaY, 'd', 'MarkerSize', 8, 'LineWidth', 1.6, 'DisplayName', '0.2% offset yield');

    xlabel('Engineering Strain, \epsilon [–]');
    ylabel('Engineering Stress, \sigma [MPa]');
    title(sprintf('Stress–Strain (Engineering) — %s', labels{i}));

    legend('Location','best');
    xlim([0, max(strain_use)*1.02]);
    ylim([0, max(stress_use)*1.05]);

    % Save figure
    figName = sprintf('stress_strain_%s.png', regexprep(labels{i}, '\s+', '_'));
    exportgraphics(gcf, figName, 'Resolution', 300);

    hold off;

    % ---------- Output results (separate from plot) ----------
    results = table( ...
        string(labels{i}), ...
        E_est/1000, ...          % [GPa]
        epsY, sigmaY, ...
        epsAtMax, sigmaMax, ...
        areaUnder, ...
        'VariableNames', { ...
            'Material', 'E_GPa', ...
            'eps_yield_02', 'sigma_yield_02_MPa', ...
            'eps_peak', 'sigma_peak_MPa', ...
            'area_MPa' ...
        });

    disp(results);

    outName = sprintf('results_%s.csv', regexprep(labels{i}, '\s+', '_'));
    writetable(results, outName);

    fprintf('\n=== %s ===\n', labels{i});
    fprintf('E_est = %.3f [GPa]\n', E_est/1000);
    fprintf('0.2%% offset yield: eps=%.6f, sigma0.2=%.2f [MPa]\n', epsY, sigmaY);
    fprintf('Peak: eps=%.6f, sigma_max=%.2f [MPa]\n', epsAtMax, sigmaMax);
    fprintf('Area = %.2f [MPa] (~%.2f [MJ/m^3])\n', areaUnder, areaUnder);

    % Store curves for comparison plot
    strain_all{i} = strain_use;
    stress_all{i} = stress_use;

end

%% Comparison plot: Steel vs Aluminum (clean curves only)

figure('Color','w');
hold on; grid on;

plot(strain_all{1}, stress_all{1}, 'LineWidth', 1.6);
plot(strain_all{2}, stress_all{2}, 'LineWidth', 1.6);

xlabel('Engineering Strain, \epsilon [–]');
ylabel('Engineering Stress, \sigma [MPa]');
title('Stress–Strain Comparison');

legend({'SAE-1045', 'Aluminum 6061'}, 'Location', 'best');

xlim([0, max([strain_all{1}; strain_all{2}]) * 1.02]);
ylim([0, max([stress_all{1}; stress_all{2}]) * 1.05]);

% Save figure
exportgraphics(gcf, 'stress_strain_comparison.png', 'Resolution', 300);

hold off;

