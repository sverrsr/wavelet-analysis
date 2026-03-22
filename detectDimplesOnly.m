function results = detectDimplesOnly(surfElev, varargin)
%DETECTDIMPLESONLY Detect persistent dimple signatures without lambda_2 validation.
% This helper is the dimple-only branch of the wavelet workflow. It accepts
% a surface-elevation stack such as surfElev(256,256,12500), detects only
% near-circular indentation signatures, applies the overlap-based lifetime
% filter, and saves the requested outputs:
%   1) dimpleBinaryMap : Ny x Nx x Nt logical map of persistent dimples.
%   2) dimpleAreaByFrame : 1 x Nt total dimple area per frame.
%   3) dimpleCountByFrame : 1 x Nt number of dimples per frame.
%
% Usage examples
%   results = detectDimplesOnly(surfElev);
%   results = detectDimplesOnly(surfElev, 'OutputFile', 'dimple_results.mat');
%   results = detectDimplesOnly(surfElev, 'dt', 1/200, 'Tinf', 1.0);
%
% Name-value options
%   'OutputFile'            : MAT-file path for saved outputs.
%   'dt'                    : Frame spacing used by the lifetime filter.
%   'Tinf'                  : Integral time scale used for tau_min.
%   'Scales'                : Mexican Hat scales, default [1 2 4].
%   'MinArea'               : Minimum connected-component area in pixels.
%   'ResponseThresholdStd'  : Mean + N*std threshold on peak response.
%   'EccentricityThreshold' : Dimple criterion, default 0.85.
%   'OverlapThreshold'      : Track-linking overlap threshold.
%   'IndentationSign'       : -1 when negative elevation means dimple.
%
% MATLAB requirements
%   - MATLAB R2022b or later.
%   - Wavelet Toolbox.
%   - Image Processing Toolbox.
%
% Physical notes
%   - The Mexican Hat response isolates bowl-shaped indentations.
%   - Only components with eccentricity < 0.85 are kept as dimples.
%   - The lambda_2 criterion is intentionally omitted here.
%   - The lifetime filter keeps coherent dimples longer than
%     tau_min = 0.166 * Tinf.

    parser = inputParser;
    parser.FunctionName = mfilename;
    addRequired(parser, 'surfElev', @(x) isnumeric(x) && ndims(x) == 3);
    addParameter(parser, 'OutputFile', 'dimple_detection_results.mat', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'dt', 1 / 200, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'Tinf', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'Scales', [1 2 4], @(x) isnumeric(x) && isvector(x) && all(x > 0));
    addParameter(parser, 'MinArea', 9, @(x) isnumeric(x) && isscalar(x) && x >= 1);
    addParameter(parser, 'ResponseThresholdStd', 1.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'EccentricityThreshold', 0.85, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
    addParameter(parser, 'OverlapThreshold', 0.15, @(x) isnumeric(x) && isscalar(x) && x >= 0 && x <= 1);
    addParameter(parser, 'IndentationSign', -1, @(x) isnumeric(x) && isscalar(x) && x ~= 0);
    parse(parser, surfElev, varargin{:});
    options = parser.Results;

    etaStack = double(surfElev);
    [Ny, Nx, Nt] = size(etaStack);

    config = struct();
    config.scales = reshape(options.Scales, 1, []);
    config.indentationSign = options.IndentationSign;
    config.waveletFamily = 'db4';
    config.waveletLevel = 2;
    config.waveletThresholdMode = 's';
    config.waveletThresholdMultiplier = 3.0;
    config.responseThresholdStd = options.ResponseThresholdStd;
    config.minArea = options.MinArea;
    config.connectivity = 8;
    config.eccentricityThreshold = options.EccentricityThreshold;
    config.overlapThreshold = options.OverlapThreshold;
    config.dt = options.dt;
    config.Tinf = options.Tinf;
    config.tauMin = 0.166 * options.Tinf;
    config.returnResponse = false;

    [featureFrames, detectionDiagnostics] = detectFeatures(etaStack, config);
    [tracks, trackedFrames, trackingDiagnostics] = trackFeatures(featureFrames, config);

    dimpleBinaryMap = false(Ny, Nx, Nt);
    dimpleAreaByFrame = zeros(1, Nt);
    dimpleCountByFrame = zeros(1, Nt);
    frameLog = [];
    trackLog = [];
    areaLog = [];
    eccentricityLog = [];
    centroidXLog = [];
    centroidYLog = [];
    scaleLog = [];
    persistentLog = [];

    for frameIndex = 1:Nt
        frameFeatures = trackedFrames(frameIndex).features;
        for featureIndex = 1:numel(frameFeatures)
            feature = frameFeatures(featureIndex);
            if ~strcmp(feature.classLabel, 'dimple') || ~feature.isPersistent
                continue;
            end

            dimpleBinaryMap(:, :, frameIndex) = dimpleBinaryMap(:, :, frameIndex) | pixelIdxToMask(feature.pixelIdxList, Ny, Nx);
            dimpleAreaByFrame(frameIndex) = dimpleAreaByFrame(frameIndex) + feature.area;
            dimpleCountByFrame(frameIndex) = dimpleCountByFrame(frameIndex) + 1;

            frameLog(end + 1, 1) = frameIndex; %#ok<AGROW>
            trackLog(end + 1, 1) = double(feature.trackId); %#ok<AGROW>
            areaLog(end + 1, 1) = double(feature.area); %#ok<AGROW>
            eccentricityLog(end + 1, 1) = double(feature.eccentricity); %#ok<AGROW>
            centroidXLog(end + 1, 1) = double(feature.centroid(1)); %#ok<AGROW>
            centroidYLog(end + 1, 1) = double(feature.centroid(2)); %#ok<AGROW>
            scaleLog(end + 1, 1) = double(feature.scale); %#ok<AGROW>
            persistentLog(end + 1, 1) = double(feature.isPersistent); %#ok<AGROW>
        end
    end

    dimpleComponentTable = table(frameLog, trackLog, areaLog, eccentricityLog, centroidXLog, centroidYLog, scaleLog, persistentLog, ...
        'VariableNames', {'Frame', 'TrackId', 'Area', 'Eccentricity', 'CentroidX', 'CentroidY', 'Scale', 'IsPersistent'});

    results = struct();
    results.dimpleBinaryMap = dimpleBinaryMap;
    results.dimpleAreaByFrame = dimpleAreaByFrame;
    results.dimpleCountByFrame = dimpleCountByFrame;
    results.trackedFrames = trackedFrames;
    results.tracks = tracks(strcmp({tracks.classLabel}, 'dimple'));
    results.componentTable = dimpleComponentTable;
    results.configuration = config;
    results.detectionDiagnostics = detectionDiagnostics;
    results.trackingDiagnostics = trackingDiagnostics;
    results.inputSize = [Ny, Nx, Nt];

    outputFile = char(options.OutputFile);
    save(outputFile, '-struct', 'results', '-v7.3');

    fprintf('Saved dimple-only outputs to %s\n', outputFile);
    fprintf('Persistent dimples detected across %d frames.\n', Nt);
end

function mask = pixelIdxToMask(pixelIdxList, Ny, Nx)
%PIXELIDXTOMASK Convert a linear-index list into a logical mask.

    mask = false(Ny, Nx);
    mask(pixelIdxList) = true;
end
