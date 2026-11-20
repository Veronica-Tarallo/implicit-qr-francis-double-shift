function [eigs, controllo, iter_history] = implicitQRAlg(A0, maxIter, tol, showPlot)
% Computes eigenvalues of a square matrix using the implicitly shifted 
% QR algorithm (Francis Double Shift, m=2)
%
% INPUT
%   A0          nxn square matrix
%   maxIter     (optional) max number of iterations per eigenvalue (default 100)
%   tol         (optional) tolerance for subdiagonal elements (default 1e-6)
%   showPlot    (optional) boolean, true to plot deflation process (default false)
%
% OUTPUT
%   eigs            vector containing the computed eigenvalues
%   controllo       0 if successful, 1 if maxIter was reached (stagnation)
%   iter_history    array with the number of iterations used for each deflation

    % --- OPTIONAL INPUTS ---
    if nargin < 4
        showPlot = false;
    end
    if nargin < 3
        tol = 1e-6;
    end
    if nargin < 2
        maxIter = 100;
    end

    % --- INITIALIZATION ---
    A = hess(A0);    % reduce A to upper Hessenberg form [cost: O(n^3)]
    [r,n] = size(A);
    if r~=n
        disp('Error: matrix A is not square!');
        eigs = []; controllo = -1; iter_history = [];
        return;
    end
    
    eigs = zeros(n,1);  % initialize the vector to store the eigenvalues
    nActive = n;        
    controllo = 0;      % control variable (1 if the algorithm is not able to find the eigenvalues)
    iter_history = [];  % initialize the iterations' history

    % Initialize figure for plot if requested
    if showPlot
        figure('Name', 'QR Convergence', 'NumberTitle', 'off');
        hold on; grid on;
        title('Convergence of the sub-diagonal elements (log scale)');
        xlabel('Iterations per cycle');
        ylabel('Magnitude of the deflated element');
    end

    % --- CORE OF THE ALGORITHM ---
    while(nActive > 2) 
        numIter = 0;
        
        % Initialize both possible points for the deflation
        subdiag1 = zeros(maxIter+1,1);
        subdiag2 = zeros(maxIter+1,1);
        
        subdiag1(1) = abs(A(nActive, nActive-1));
        if nActive > 2
            subdiag2(1) = abs(A(nActive-1, nActive-2));
        end

        while abs(A(nActive, nActive-1))>tol && abs(A(nActive-1, nActive-2))>tol && numIter<=maxIter
            numIter = numIter+1;
            
            % --- Double Shift (m=2) ---
            H2 = A(nActive-1:nActive, nActive-1:nActive);
            s = H2(1,1) + H2(2,2);                  % trace
            t = H2(1,1)*H2(2,2) - H2(1,2)*H2(2,1);  % determinant
            
            % Vector x = (A^2-sA+tI)*e1
            idx = min(3, nActive);
            xTilde = zeros(idx, 1);
            xTilde(1) = A(1,1)^2 + A(1,2)*A(2,1) - s*A(1,1) + t;
            xTilde(2) = A(2,1)*(A(1,1) + A(2,2) - s);
            if idx == 3
                xTilde(3) = A(2,1)*A(3,2);
            end
            xTilde = xTilde / norm(xTilde);
                              
            % Householder and initial transformation
            [~, H] = householder_vector(xTilde);
            Q0_tilde = H';
            A(1:idx,:) = Q0_tilde'*A(1:idx,:);
            A(:,1:idx) = A(:,1:idx)*Q0_tilde;
        
            % Bulge chasing       
            for j=1:nActive-2
                idx_start = j+1;
                idx_end = min(j+3, nActive); 
        
                v_chase = A(idx_start:idx_end,j);
                [~, Qj_tilde] = householder_vector(v_chase);
        
                A(idx_start:idx_end,j:nActive) = Qj_tilde'*A(idx_start:idx_end,j:nActive);
                A(1:nActive, idx_start:idx_end) = A(1:nActive, idx_start:idx_end)*Qj_tilde;
            end
            
            % Update for the plot
            subdiag1(numIter+1) = abs(A(nActive, nActive-1));
            if nActive > 2
                subdiag2(numIter+1) = abs(A(nActive-1, nActive-2));
            end
        end

        % --- DEFLATION AND PLOT ---
        if abs(A(nActive, nActive-1)) <= tol
            % Real eigenvalue found (1x1 block)
            if showPlot && numIter > 0
                % Plot with a dashed line
                semilogy(0:numIter, subdiag1(1:numIter+1), '--o', 'LineWidth', 1.5, 'MarkerSize', 4);
                drawnow;
            end
            
            eigs(nActive) = A(nActive, nActive);
            nActive = nActive - 1;
            iter_history = [iter_history; numIter]; 
            
        elseif nActive > 1 && abs(A(nActive-1, nActive-2)) <= tol
            % 2x2 block found (eigenvalues can be real or complex)
            if showPlot && numIter > 0
                % Plot with a solid line
                semilogy(0:numIter, subdiag2(1:numIter+1), '-o', 'LineWidth', 1.5, 'MarkerSize', 4);
                drawnow;
            end
            
            ee = eig(A(nActive-1:nActive, nActive-1:nActive));
            eigs(nActive-1:nActive) = ee;
            nActive = nActive - 2;
            iter_history = [iter_history; numIter]; 
            
        elseif numIter > maxIter
            controllo = 1; 
            break;
        end
        
        if nActive > 2
            A = A(1:nActive, 1:nActive);
        end
    end

    % Last block
    if (controllo == 0 && nActive > 0)
        ee = eig(A(1:nActive, 1:nActive));
        eigs(1:nActive) = ee;
        % No iterations are needed for the last block, hence we do not put iter_history here
    end
end