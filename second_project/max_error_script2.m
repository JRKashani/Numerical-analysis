% Driver Script for Heat Equation Solver Testing (with Diff Analysis)
% ===================================================================
filename = 'Heat_Calc_Results_200_matrices.xlsx';
if isfile(filename)
    delete(filename);
end
rng('shuffle');

% Storage for summary tab
num_runs = 10;
all_times = zeros(num_runs, 1);
all_avg_diffs = zeros(num_runs, 1);

for k = 1:num_runs
    % 1. Generate Random Temperatures
    limits = [-30, 100];
    random_pair = limits(1) + (limits(2) - limits(1)) * rand(1, 2);
    sorted_temps = sort(random_pair);
    t_low = sorted_temps(1);
    t_high = sorted_temps(2);
    
    if t_low == t_high, t_high = t_high + 0.1; end
    
    % 2. Run and Time the Solver
    try
        timer_val = tic;                 
        % CRITICAL: heat_calc must return [final_mat, initial_guess_mat]
        [temp_mat, ini_mat] = heat_calc(t_low, t_high); 
        run_time = toc(timer_val);       
    catch ME
        warning('Run %d failed: %s', k, ME.message);
        continue;
    end
    
    % 3. Calculate Metrics
    % Laplacian Error (Convergence Check)
    LapT = temp_mat(1:end-2, 2:end-1) + ...
           temp_mat(3:end,   2:end-1) + ...
           temp_mat(2:end-1, 1:end-2) + ...
           temp_mat(2:end-1, 3:end) + ...
           (-4) * temp_mat(2:end-1, 2:end-1);
    max_error = max(abs(LapT(:)));
    
    % Difference Matrix (Final - Initial)
    diff_mat = temp_mat - ini_mat;
    
    % Average Absolute Difference (Scalar)
    % Measures how much the solution "moved" from the guess
    avg_diff_scalar = mean(abs(diff_mat(:)));
    
    % Store for final summary
    all_times(k) = run_time;
    all_avg_diffs(k) = avg_diff_scalar;
    
    % 4. Prepare Data for Export
    % We stack blocks of data vertically
    
    % Block A: Final Result Matrix
    block_final = num2cell(temp_mat);
    width = size(temp_mat, 2);
    
    % Block B: Run Parameters & Scalar Metrics
    spacer = cell(1, width);
    header_row = cell(1, width);
    value_row  = cell(1, width);
    
    header_row(1:5) = {'T_Low', 'T_High', 'Max_Error', 'Time_Sec', 'Avg_Abs_Diff'};
    value_row(1:5)  = {t_low,    t_high,   max_error,   run_time,   avg_diff_scalar};
    
    % Block C: Difference Matrix Label
    label_row = cell(1, width);
    label_row{1} = 'Difference Matrix (Final - Initial Guess)';
    
    % Block D: Difference Matrix Data
    block_diff = num2cell(diff_mat);
    
    % Combine all blocks
    final_output = [
        block_final; 
        spacer; 
        header_row; 
        value_row; 
        spacer; 
        label_row; 
        block_diff
    ];
    
    sheet_name = sprintf('Run_%d', k);
    writecell(final_output, filename, 'Sheet', sheet_name);
    
    fprintf('Run %d: Time=%.4fs, Error=%.4e, AvgDiff=%.4f\n', ...
        k, run_time, max_error, avg_diff_scalar);
end

% 5. Create Summary Tab (Sheet 11)
summary_headers = {'Metric', 'Value'};
summary_data = {
    'Average Run Time (s)', mean(all_times);
    'Average of Avg Differences', mean(all_avg_diffs)
};
fprintf('all Runs stats: Time=%.4fs, AvgDiff=%.4f\n', ...
        mean(all_times), mean(all_avg_diffs));

writecell([summary_headers; summary_data], filename, 'Sheet', 'Summary_Stats');

disp('All simulations finished. Summary saved to Sheet "Summary_Stats".');