% Fractal Surface Generation via Spectral Synthesis (fBm)
clear; clc;

% --- Parameters ---
N = 4096;           % Resolution (N x N grid)
H = 0.2;            % Hurst exponent (0 < H < 1). Higher = Smoother.
L = 100;            % Physical length of the side

% --- Frequency Domain Setup ---
% Create a grid of frequencies
k_vec = [(0:N/2), (-N/2+1:-1)];
[Kx, Ky] = meshgrid(k_vec, k_vec);
K_sq = Kx.^2 + Ky.^2;
K_sq(1,1) = 1; % Avoid division by zero at the DC component

% --- Generate the Spectral Density ---
% The power spectrum of fBm follows P(f) ~ 1 / f^(2H + 2)
% The amplitude follows A(f) ~ 1 / f^(H + 1)
amplitude = K_sq.^ (-(H + 1) / 2);

% Create random phases
phases = exp(2 * pi * 1i * rand(N, N));

% Apply the power law to the random phases
F_surface = amplitude .* phases;

% Force the surface to be real-valued (Hermitian symmetry)
% (Optional for visual results, as we take the real part anyway)
surface = real(ifft2(F_surface));

% --- Normalization and Visualization ---
[X, Y] = meshgrid(linspace(0, L, N), linspace(0, L, N));

figure('Color', 'w', 'Position', [100, 100, 1000, 800]);
s = surf(X, Y, surface, 'EdgeColor', 'none');

% Create a custom 'Terrain' colormap (Browns, Greens, Whites)
custom_terrain = [
    0.0, 0.2, 0.5;  % Deep Water
    0.0, 0.5, 0.8;  % Shallow Water
    0.9, 0.8, 0.6;  % Sand/Beach
    0.2, 0.5, 0.2;  % Grass/Forest
    0.4, 0.3, 0.2;  % Dirt/Mountain
    1.0, 1.0, 1.0   % Snow Cap
];
colormap(interp1(linspace(0,1,6), custom_terrain, linspace(0,1,256)));

% Aesthetics
lighting gouraud;
camlight headlight;
material dull;      % Reduces "plastic" look
axis tight;
view(45, 45);
title(['Fractal Surface (fBm) with H = ', num2str(H)], 'FontSize', 14);
zlabel('Height');
colorbar;