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
inputData = loadInputData(config);
etaStack = inputData.etaStack;

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

%% Local helpers
function inputData = loadInputData(config)
%LOADINPUTDATA Load external periodic data or synthesize a demo dataset.
% The preferred external input is a MAT file containing etaStack and any of
% the optional fields beta2Stack, beta2ByDepth, zLevels, lambda2Stack, or
% velocityGradientStack.

    candidateFiles = {'vortex_input.mat', 'input/vortex_input.mat', 'data/vortex_input.mat'};
    dataFile = '';
    for idx = 1:numel(candidateFiles)
        if exist(candidateFiles{idx}, 'file') == 2
            dataFile = candidateFiles{idx};
            break;
        end
    end

    if ~isempty(dataFile)
        inputData = load(dataFile);
        if ~isfield(inputData, 'etaStack')
            error('Input MAT file must contain etaStack.');
        end
        return;
    end

    warning(['No vortex_input.mat file was found. Generating a compact synthetic ' ...
             'periodic dataset so the full MATLAB pipeline remains runnable.']);

    Nt = config.maxFramesForDemo;
    Ny = 256;
    Nx = 256;
    yCoordinates = linspace(0, 2 * pi, Ny + 1);
    xCoordinates = linspace(0, 2 * pi, Nx + 1);
    yCoordinates(end) = [];
    xCoordinates(end) = [];
    [yGrid, xGrid] = ndgrid(yCoordinates, xCoordinates);

    etaStack = zeros(Ny, Nx, Nt, 'single');
    beta2Stack = zeros(Ny, Nx, Nt, 'single');
    zLevels = linspace(0, 2 * config.lambdaT, 16);
    beta2ByDepth = zeros(1, 1, numel(zLevels), Nt, 'single');

    for t = 1:Nt
        phase = 2 * pi * (t - 1) / Nt;
        dimple1 = exp(-((wrapPeriodic(xGrid - (0.8 + 0.6 * cos(phase)))).^2 + ...
                        (wrapPeriodic(yGrid - (1.5 + 0.4 * sin(phase)))).^2) / (2 * 0.11^2));
        dimple2 = exp(-((wrapPeriodic(xGrid - (4.6 - 0.3 * sin(1.7 * phase)))).^2 + ...
                        (wrapPeriodic(yGrid - (4.0 + 0.2 * cos(phase)))).^2) / (2 * 0.14^2));
        scar = exp(-(wrapPeriodic((xGrid - 3.2) * cos(0.35) + (yGrid - 2.8) * sin(0.35))).^2 / (2 * 0.52^2) ...
                   -(wrapPeriodic(-(xGrid - 3.2) * sin(0.35) + (yGrid - 2.8) * cos(0.35))).^2 / (2 * 0.09^2));
        background = 0.03 * sin(4 * xGrid + phase) .* cos(3 * yGrid - 0.3 * phase);

        etaFrame = background - 0.15 * dimple1 - 0.12 * dimple2 - 0.09 * scar;
        etaStack(:, :, t) = single(etaFrame);

        beta2Frame = 0.7 * dimple1 + 0.4 * dimple2 + 0.8 * scar + 0.2 * abs(background);
        beta2Stack(:, :, t) = single(beta2Frame .^ 2);

        depthEnvelope = exp(-0.5 * ((zLevels - config.Lnu) / (0.35 * config.lambdaT)).^2);
        beta2ByDepth(1, 1, :, t) = single((0.2 + 0.8 * depthEnvelope) * mean(beta2Frame(:) .^ 2));
    end

    inputData = struct();
    inputData.etaStack = etaStack;
    inputData.beta2Stack = beta2Stack;
    inputData.beta2ByDepth = beta2ByDepth;
    inputData.zLevels = zLevels;
end

function wrapped = wrapPeriodic(values)
%WRAPPERIODIC Wrap values to the [-pi, pi) interval for periodic offsets.

    wrapped = mod(values + pi, 2 * pi) - pi;
end
