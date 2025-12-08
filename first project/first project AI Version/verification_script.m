% FINAL_VALIDATION_SCRIPT
% 1. Verifies specific cases required by the project PDF (Bi=24.8).
% 2. Runs the performance benchmark (100,000 iterations).

clc; clear;

fprintf('--- PART 1: Specific Case Validation ---\n');

% Case 1: Bi = 24.8, m = 1 (Expected: ~5.305)
Bi_test = 24.8;
m1 = 1;
root1 = SolveBesselAI(Bi_test, m1);
fprintf('Case 1 (Bi=%.1f, m=%d): x = %.5f\n', Bi_test, m1, root1);

% Case 2: Bi = 24.8, m = 2
m2 = 2;
root2 = SolveBesselAI(Bi_test, m2);
fprintf('Case 2 (Bi=%.1f, m=%d): x = %.5f\n', Bi_test, m2, root2);

% Case 3: Arbitrary Test (Bi = 0.5, m = 3)
Bi_small = 0.5;
m3 = 3;
root3 = SolveBesselAI(Bi_small, m3);
fprintf('Case 3 (Bi=%.1f, m=%d):  x = %.5f\n', Bi_small, m3, root3);

fprintf('\n--- Additional Arbitrary Large Tests ---\n');

% Case 4: Medium Root (m=50)
% Expected location: ~ 50 * pi = 157.07
Bi_4 = 10.0;
m_4 = 50;
root4 = SolveBesselAI(Bi_4, m_4);
fprintf('Case 4 (Bi=%.1f, m=%d):   x = %.5f (Expected ~157)\n', ...
        Bi_4, m_4, root4);

% Case 5: Large Root (m=300)
% Expected location: ~ 300 * pi = 942.47
% Note: This requires the dynamic search limit fix!
Bi_5 = 0.1; 
m_5 = 300;
root5 = SolveBesselAI(Bi_5, m_5);
fprintf('Case 5 (Bi=%.1f, m=%d):  x = %.5f (Expected ~942)\n', ...
        Bi_5, m_5, root5);

% Case 6: Very Large Root (m=1000)
% Expected location: ~ 1000 * pi = 3141.59
% Note: This effectively stress-tests the scanning loop (~31,000 steps)
Bi_6 = 60.0;
m_6 = 1000;
root6 = SolveBesselAI(Bi_6, m_6);
fprintf('Case 6 (Bi=%.1f, m=%d): x = %.5f (Expected ~3141)\n', ...
        Bi_6, m_6, root6);


fprintf('\n--- PART 2: Performance Benchmark (AI Code) ---\n');
fprintf('Running 100,000 iterations for Bi=%.1f, m=%d...\n', Bi_test, m1);


tic; % Start timer
for k = 1:100000
    % Suppress output for speed, just run the logic
    dummy_var = SolveBesselAI(Bi_test, m1);
end
runtime_ai = toc; % Stop timer

fprintf('Total time: %.4f seconds\n', runtime_ai);
fprintf('Average time per call: %.2f microseconds\n', (runtime_ai / 100000) * 1e6);