% MATLAB script to estimate the fractal (box-counting) dimension of the graph
% of the fractional derivative D^alpha W(x) from values.csv using the box-counting method.
% Refined: Full wide epsilon range (dx/1000 to dx/2, 40 points) to capture all scales,
% then automatically select the best linear segment (window of 10-20 points with max R^2)
% for the fit, focusing on the self-similar regime to reduce bias from linear tails.

clear; clc; close all;

% Load data from values.csv (columns: x, Dw)
data = readmatrix('values.csv');
x = data(:, 1);
Dw = data(:, 2);

% Ensure x is sorted (should be, but just in case)
[x, sortIdx] = sort(x);
Dw = Dw(sortIdx);

% Parameters for box sizes (wide range for full log-log curve)
dx = max(x) - min(x);
numSizes = 40;           % More points for detailed log-log
minEps = dx / 1000;      % Very fine min epsilon
maxEps = dx / 2;         % Coarse max epsilon
epsilons = logspace(log10(minEps), log10(maxEps), numSizes);

% Compute N(epsilon) for each box size
logInvEps = log(1 ./ epsilons);
logN = zeros(size(epsilons));
for i = 1:length(epsilons)
    epsilon = epsilons(i);
    nBoxesX = ceil(dx / epsilon);
    N_val = 0;
    for j = 1:nBoxesX
        left = min(x) + (j-1) * epsilon;
        right = min(left + epsilon, max(x));
        mask = (x >= left) & (x < right);
        if any(mask)
            yIn = Dw(mask);
            deltaY = max(yIn) - min(yIn);
            nYBoxes = floor(deltaY / epsilon) + 1;
            N_val = N_val + nYBoxes;
        end
    end
    logN(i) = log(N_val);
end

% Find best linear segment: search windows of size 10-20 for max R^2
bestD = 0;
bestR2 = 0;
bestStart = 1;
bestWS = 0;
for ws = 10:20
    for start = 1:(length(logN) - ws + 1)
        subInv = logInvEps(start:start+ws-1);
        subN = logN(start:start+ws-1);
        pSub = polyfit(subInv, subN, 1);
        DSub = pSub(1);
        residualsSub = subN - polyval(pSub, subInv);
        SS_resSub = sum(residualsSub.^2);
        SS_totSub = sum((subN - mean(subN)).^2);
        R2Sub = 1 - (SS_resSub / SS_totSub);
        if R2Sub > bestR2
            bestR2 = R2Sub;
            bestD = DSub;
            bestStart = start;
            bestWS = ws;
        end
    end
end

% Display results
fprintf('Estimated fractal dimension: %.4f\n', bestD);
fprintf('R-squared of fit: %.4f (best window: start %d, size %d)\n', bestR2, bestStart, bestWS);

% Plot full log-log with highlighted best fit
f = figure;
plot(logInvEps, logN, 'ro', 'MarkerSize', 4, 'LineWidth', 1); hold on;
% Highlight best window
plot(logInvEps(bestStart:bestStart+bestWS-1), logN(bestStart:bestStart+bestWS-1), ...
     'go', 'MarkerSize', 6, 'LineWidth', 2, 'MarkerFaceColor', 'g');
% Fit line on best window
pBest = polyfit(logInvEps(bestStart:bestStart+bestWS-1), ...
                logN(bestStart:bestStart+bestWS-1), 1);
bestInv = linspace(min(logInvEps(bestStart:bestStart+bestWS-1)), ...
                   max(logInvEps(bestStart:bestStart+bestWS-1)), 100);
plot(bestInv, polyval(pBest, bestInv), 'g-', 'LineWidth', 3);
xlabel('log(1/\epsilon)');
ylabel('log N(\epsilon)');
title('Box-Counting Dimension Estimate (Best Linear Segment Highlighted)');
legend('Full Data', 'Best Window', sprintf('Fit (D = %.3f)', bestD), 'Location', 'best');
grid on;

% Set figure size to 2400x1800 pixels
f.Position(3:4) = [2400 1800];

% Save plot with exact pixel size using exportgraphics (requires MATLAB R2020a or later)
exportgraphics(f, 'dimension_fit.png', 'Units', 'pixels', 'Width', 2400, 'Height', 1800);

disp('Script completed. Dimension estimate printed above; plot saved as dimension_fit.png (2400x1800 pixels).');