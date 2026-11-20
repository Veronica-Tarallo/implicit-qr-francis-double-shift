%% MAIN FOR IMPLICIT QR ALGORITHM (Double Shift)
% This script runs several tests to validate the QR algorithm 
% and generates plots for performance analysis (O(n^3) cost and convergence)

clear; clc; close all;
% Set seed for reproducibility
%rng("threefry")     % fino a n<=599
%rng("shuffle")      % fino a n<=488
rng("twister")      % fino a n<=399

%% TEST 1: Convergence Visualization (Plot 1)
fprintf('--- TEST 1: Sub-diagonal Convergence ---\n');
n1 = 20;
A1 = randn(n1); % Generic random matrix (will have complex eigenvalues)
%A1=A1'+A1;   % If you want to test a matrix with only real eigenvalues

% Call the function with showPlot = true
[eigs1, check1, iter1] = implicitQRAlg(A1, 100, 1e-6, true);

if check1 == 0
    fprintf('Test 1 successfully completed.\n');
    fprintf('Average number of iterations per deflation: %.2f\n', mean(iter1));
    fprintf('As predicted by theory, it takes on average O(1) iterations to extract an eigenvalue.\n');
else
    fprintf('Stagnation reached.\n');
end

%% TEST 2: Accuracy and Numerical Stability
fprintf('\n--- TEST 2: Results Validation ---\n');
n2 = 30;
A2 = randn(n2);
%A2=A2'+A2;

[my_eigs, check2, ~] = implicitQRAlg(A2, 100, 1e-6, false);
true_eigs = eig(A2);

% Sort and compare to find the maximum absolute error (Infinity Norm)
diff_vector = abs(sort(my_eigs) - sort(true_eigs));
max_err = max(diff_vector);

fprintf('Matrix %dx%d. Maximum absolute error compared to eig(A): %e\n', n2, n2, max_err);
if max_err < 1e-4
    fprintf('The algorithm is stable and highly accurate.\n');
end

%% TEST 3: Computational Cost Analysis O(n^3) (Plot 2)
fprintf('\n--- TEST 3: Computational Cost ---\n');
% Testing the algorithm on increasing dimensions
n_vals = [10, 20, 30, 45, 65, 90];
times = zeros(length(n_vals), 1);

fprintf('Running tests...\n');
for i = 1:length(n_vals)
    n_test = n_vals(i);
    A_test = randn(n_test);
    %A_test = A_test'+A_test;
    
    tic;
    [~, check3, ~] = implicitQRAlg(A_test, 100, 1e-6, false);
    times(i) = toc;
    
    fprintf('  n = %d \t Time: %.4f sec\n', n_test, times(i));
end

% Log-Log plot to show polynomial growth
figure('Name', 'Computational Cost', 'NumberTitle', 'off');
loglog(n_vals, times, '-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on; grid on;

% Add an O(n^3) reference line for visual comparison
% Calibrate the constant based on the first point to align the lines
const_c = times(1) / (n_vals(1)^3); 
ref_n3 = const_c * (n_vals.^3);
loglog(n_vals, ref_n3, '--r', 'LineWidth', 1.5);

title('QR Algorithm Computational Cost (Log-Log Scale)');
xlabel('Matrix dimension n');
ylabel('Execution time (seconds)');
legend('Actual time', 'Theoretical O(n^3) reference', 'Location', 'northwest');

%% TEST 4: Algorithm Limitations (Stagnation)
fprintf('\n--- TEST 4: Algorithm Limitations (Stagnation on large matrices) ---\n');
% Forza il seed a un caso in cui sappiamo che fallisce per n = 400 
%rng("threefry")     % fino a n<=599
%rng("shuffle")      % fino a n<=488
rng("twister")      % fino a n<=399
n3 = 400;
A3 = randn(n3);

fprintf('Testing matrix %dx%d to demonstrate stagnation...\n', n3, n3);
[~, check_fail, ~] = implicitQRAlg(A3, 100, 1e-6, false);

if check_fail == 1
    fprintf('As expected, the algorithm reached maxIter without converging.\n');
    fprintf('This demonstrates the "stagnation" phenomenon, which requires\n');
    fprintf('"exceptional shifts" (e.g., LAPACK implementation) to be resolved.\n');
else
    fprintf('The algorithm converged (try changing the RNG seed).\n');
end

fprintf('\nTests completed. Check the generated Figures.\n');