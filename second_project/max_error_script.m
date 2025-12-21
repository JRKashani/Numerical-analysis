% Driver Script for Heat Equation Solver Testing (with Timing)
% ============================================================

filename = 'Heat_Calc_Results.xlsx';

if isfile(filename)
    delete(filename);
end

rng('shuffle');

for k = 1:10
    % 1. Generate Random Temperatures
    limits = [-30, 100];
    random_pair = limits(1) + (limits(2) - limits(1)) * rand(1, 2);
    sorted_temps = sort(random_pair);
    t_low = sorted_temps(1);
    t_high = sorted_temps(2);
    
    if t_low == t_high, t_high = t_high + 0.1; end
    
    % 2. Run and Time the Solver
    % 'tic' starts the stopwatch, 'toc' returns the elapsed time
    try
        timer_val = tic;                 % START TIMER
        temp_mat = heat_calc(t_low, t_high);
        run_time = toc(timer_val);       % STOP TIMER
    catch ME
        warning('Run %d failed: %s', k, ME.message);
        continue;
    end
    
    % 3. Calculate Error
    LapT = temp_mat(1:end-2, 2:end-1) + ...
           temp_mat(3:end,   2:end-1) + ...
           temp_mat(2:end-1, 1:end-2) + ...
           temp_mat(2:end-1, 3:end) + ...
           (-4) * temp_mat(2:end-1, 2:end-1);
       
    max_error = max(abs(LapT(:)));
    
    % 4. Prepare Data for Export
    data_cell = num2cell(temp_mat);
    width = size(temp_mat, 2);
    
    % Prepare footer rows (ensure they match the matrix width)
    spacer_row = cell(1, width); 
    header_row = cell(1, width);
    value_row  = cell(1, width);
    
    % Add the Time column to the headers and values
    header_row(1:4) = {'T_Low', 'T_High', 'Max_Error', 'Time_Sec'};
    value_row(1:4)  = {t_low,    t_high,   max_error,   run_time};
    
    % Combine and Write
    final_output = [data_cell; spacer_row; header_row; value_row];
    
    sheet_name = sprintf('Run_%d', k);
    writecell(final_output, filename, 'Sheet', sheet_name);
    
    fprintf('Run %d: Time=%.4fs, Error=%.4e\n', k, run_time, max_error);
end

disp('All simulations finished.');