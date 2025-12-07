function validate_inputs(Biot_number, root_index_m)
% Inputs:
%   Biot_number:  Dimensionless number representing heat transfer ratio
%   root_index_m: The index of the positive root to find (1, 2, 3...)

    % 1. Validate Biot Number Range
    % Source constraint: 0.1 <= Bi <= 60 [cite: 252]
    if Biot_number < 0.1 || Biot_number > 60
        error(['Input Error: Biot number (Bi) must be between 0.1 ' ...
            'and 60. You provided: %.4f'], Biot_number);
    end

    % 2. Validate Root Index (m)
    % Source constraint: m is a natural number [cite: 255]
    % Check if m is positive and if it is an integer
    if root_index_m <= 0 || mod(root_index_m, 1) ~= 0
        error(['Input Error: The root index (m) must be a natural' ...
            ' number (positive integer). You provided: %f'], root_index_m);
    end
end