%% Main_Vortex_Analysis
% Wavelet-based vortex signature detection for periodic free-surface data.
% The pipeline targets MATLAB R2022b+ with Wavelet Toolbox and Image
% Processing Toolbox support. The physical defaults are aligned with the
% specification in this repository:
%   * Dimples are near-circular (eccentricity < 0.85).
%   * Scars are elongated (eccentricity > 0.85).
%   * The viscous layer depth is L_nu = 0.258 * lambda_T.
%   * The lambda_2 threshold is lambda2_th = -6 * (uTilde / lambda_T)^2.
%   * The lifetime filter is tau_min = 0.166 * T_inf.
%
% The code assumes the primary input is a periodic stack etaStack of size
% [Ny x Nx x Nt], with Nx = Ny = 256 and Nt = 12500 in the production use
% case. Spatial periodicity is respected by FFT-based convolution inside
% detectFeatures.

clearvars;
clc;

%% Configuration
config = struct();
config.scales = [1 2 4];
config.indentationSign = -1;  % Negative eta corresponds to a surface indentation.
config.waveletFamily = 'db4';
config.waveletLevel = 2;
config.waveletThresholdMode = 's';
config.waveletThresholdMultiplier = 3.0;
config.responseThresholdStd = 1.5;
config.minArea = 9;
config.connectivity = 8;
config.eccentricityThreshold = 0.85;
config.overlapThreshold = 0.15;
config.dt = 1.0 / 200.0;
config.Tinf = 1.0;
config.tauMin = 0.166 * config.Tinf;
config.lambdaT = 24.0;
config.uTilde = 0.12;
config.Lnu = 0.258 * config.lambdaT;
config.lambda2Threshold = -6.0 * (config.uTilde / config.lambdaT) ^ 2;
config.gaussianDecaySigma = 0.5 * config.lambdaT;
config.returnResponse = false;
config.maxFramesForDemo = 96;

%% Load or synthesize input data
load("C:\Users\sverrsr\Documents\SYNC\flow-statistics\surfElev_RE2500_WEINF.mat");

etaStack = permute(surfElev, [2 3 1]);

fprintf('Loaded etaStack with size [%d x %d x %d].\n', size(etaStack, 1), size(etaStack, 2), size(etaStack, 3));
fprintf('Using L_nu = %.4f and lambda2 threshold = %.6e.\n', config.Lnu, config.lambda2Threshold);

%% Run the pipeline
[featureFrames, detectDiagnostics] = detectFeatures(etaStack, config);
[tracks, trackedFrames, trackingDiagnostics] = trackFeatures(featureFrames, config);
validation = validateVortices(trackedFrames, tracks, inputData, config);

%% Summaries and persistence-aware metrics
summary = struct();
summary.configuration = config;
summary.detection = detectDiagnostics;
summary.tracking = trackingDiagnostics;
summary.validation = validation;
summary.numFrames = numel(trackedFrames);
summary.numTracks = numel(tracks);
summary.numPersistentTracks = nnz([tracks.isPersistent]);
summary.scarNCC = validation.metrics.scarAreaNCC;
summary.dimpleNCC = validation.metrics.dimpleAreaNCC;
summary.scarLagFrames = validation.metrics.referenceLagFrames;
summary.scarLagSeconds = summary.scarLagFrames * config.dt;

save('wavelet_vortex_analysis_results.mat', 'summary', 'trackedFrames', 'tracks', '-v7.3');

fprintf('Persistent tracks: %d / %d\n', summary.numPersistentTracks, summary.numTracks);
fprintf('Scar NCC (area vs beta^2): %.4f\n', summary.scarNCC);
fprintf('Dimple NCC (area vs beta^2): %.4f\n', summary.dimpleNCC);

