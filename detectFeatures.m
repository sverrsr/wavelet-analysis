function [featureFrames, diagnostics] = detectFeatures(etaStack, config)
%DETECTFEATURES Detect dimple and scar candidates from periodic eta data.
% The detector uses Wavelet Toolbox denoising (wavedec2 + wthresh +
% waverec2) before a vectorized FFT-based convolution with a 2-D Mexican Hat
% kernel. The FFT route preserves the 2*pi periodicity of the interface.
%
% Inputs
%   etaStack : [Ny x Nx x Nt] periodic surface elevation stack.
%   config   : structure with fields documented in Main_Vortex_Analysis.
%
% Outputs
%   featureFrames : per-frame structure array with connected-component
%                   measurements and classification labels.
%   diagnostics   : summary statistics for wavelet responses.

    arguments
        etaStack (:,:,:) {mustBeNumeric}
        config struct
    end

    [Ny, Nx, Nt] = size(etaStack);
    scales = config.scales(:).';
    numScales = numel(scales);

    kernels = cell(1, numScales);
    fftKernels = cell(1, numScales);
    for scaleIndex = 1:numScales
        kernels{scaleIndex} = buildMexicanHatKernel([Ny, Nx], scales(scaleIndex));
        fftKernels{scaleIndex} = fft2(ifftshift(kernels{scaleIndex}));
    end

    emptyFeature = struct( ...
        'frame', [], ...
        'componentIndex', [], ...
        'classLabel', '', ...
        'area', [], ...
        'eccentricity', [], ...
        'centroid', [], ...
        'orientation', [], ...
        'majorAxisLength', [], ...
        'minorAxisLength', [], ...
        'pixelIdxList', [], ...
        'boundingBox', [], ...
        'meanResponse', [], ...
        'peakResponse', [], ...
        'scale', [], ...
        'trackId', [], ...
        'isPersistent', false, ...
        'geodesicCenter', []);

    frameTemplate = struct( ...
        'frame', [], ...
        'candidateMask', false(Ny, Nx), ...
        'responseThreshold', [], ...
        'responseMean', [], ...
        'responseStd', [], ...
        'features', emptyFeature([]));

    featureFrames = repmat(frameTemplate, Nt, 1);
    diagnostics = struct();
    diagnostics.responseThresholds = zeros(Nt, 1);
    diagnostics.numCandidates = zeros(Nt, 1);
    diagnostics.numDimples = zeros(Nt, 1);
    diagnostics.numScars = zeros(Nt, 1);

    for frameIndex = 1:Nt
        etaFrame = double(etaStack(:, :, frameIndex));
        denoisedFrame = denoiseFrameWavelet(etaFrame, config);

        responses = zeros(Ny, Nx, numScales);
        frameSpectrum = fft2(config.indentationSign * denoisedFrame);
        for scaleIndex = 1:numScales
            responses(:, :, scaleIndex) = real(ifft2(frameSpectrum .* fftKernels{scaleIndex}));
        end

        [peakResponse, peakScaleIndices] = max(responses, [], 3);
        responseMean = mean(peakResponse(:));
        responseStd = std(peakResponse(:));
        responseThreshold = responseMean + config.responseThresholdStd * responseStd;

        candidateMask = peakResponse > responseThreshold;
        candidateMask = bwareaopen(candidateMask, config.minArea, config.connectivity);
        candidateMask = imfill(candidateMask, 'holes');

        componentLabels = bwlabel(candidateMask, config.connectivity);
        componentProps = regionprops(componentLabels, peakResponse, ...
            'Area', 'Eccentricity', 'Centroid', 'Orientation', ...
            'MajorAxisLength', 'MinorAxisLength', 'PixelIdxList', ...
            'BoundingBox', 'MeanIntensity', 'MaxIntensity');

        features = repmat(emptyFeature, numel(componentProps), 1);
        for componentIndex = 1:numel(componentProps)
            featureScale = mode(peakScaleIndices(componentProps(componentIndex).PixelIdxList));
            classLabel = 'scar';
            if componentProps(componentIndex).Eccentricity < config.eccentricityThreshold
                classLabel = 'dimple';
            end

            features(componentIndex) = struct( ...
                'frame', frameIndex, ...
                'componentIndex', componentIndex, ...
                'classLabel', classLabel, ...
                'area', componentProps(componentIndex).Area, ...
                'eccentricity', componentProps(componentIndex).Eccentricity, ...
                'centroid', componentProps(componentIndex).Centroid, ...
                'orientation', componentProps(componentIndex).Orientation, ...
                'majorAxisLength', componentProps(componentIndex).MajorAxisLength, ...
                'minorAxisLength', componentProps(componentIndex).MinorAxisLength, ...
                'pixelIdxList', componentProps(componentIndex).PixelIdxList, ...
                'boundingBox', componentProps(componentIndex).BoundingBox, ...
                'meanResponse', componentProps(componentIndex).MeanIntensity, ...
                'peakResponse', componentProps(componentIndex).MaxIntensity, ...
                'scale', scales(featureScale), ...
                'trackId', [], ...
                'isPersistent', false, ...
                'geodesicCenter', []);
        end

        featureFrames(frameIndex).frame = frameIndex;
        featureFrames(frameIndex).candidateMask = candidateMask;
        featureFrames(frameIndex).responseThreshold = responseThreshold;
        featureFrames(frameIndex).responseMean = responseMean;
        featureFrames(frameIndex).responseStd = responseStd;
        featureFrames(frameIndex).features = features;

        diagnostics.responseThresholds(frameIndex) = responseThreshold;
        diagnostics.numCandidates(frameIndex) = numel(features);
        diagnostics.numDimples(frameIndex) = nnz(strcmp({features.classLabel}, 'dimple'));
        diagnostics.numScars(frameIndex) = nnz(strcmp({features.classLabel}, 'scar'));
    end
end

function denoisedFrame = denoiseFrameWavelet(etaFrame, config)
%DENOISEFRAMEWAVELET Wavelet shrinkage using the requested Wavelet Toolbox functions.

    [coeffs, bookkeeping] = wavedec2(etaFrame, config.waveletLevel, config.waveletFamily);
    [horizontalDetail, verticalDetail, diagonalDetail] = detcoef2('all', coeffs, bookkeeping, 1);
    allDetailCoeffs = [horizontalDetail(:); verticalDetail(:); diagonalDetail(:)];
    sigmaEstimate = median(abs(allDetailCoeffs)) / 0.6745;
    threshold = config.waveletThresholdMultiplier * sigmaEstimate;
    thresholdedCoeffs = wthresh(coeffs, config.waveletThresholdMode, threshold);
    denoisedFrame = waverec2(thresholdedCoeffs, bookkeeping, config.waveletFamily);
end

function kernel = buildMexicanHatKernel(frameSize, scale)
%BUILDMEXICANHATKERNEL Construct a zero-mean 2-D axisymmetric Mexican Hat.
% The kernel is normalized to support scale comparisons across s = 1, 2, 4.

    Ny = frameSize(1);
    Nx = frameSize(2);
    y = (-floor(Ny / 2)):(ceil(Ny / 2) - 1);
    x = (-floor(Nx / 2)):(ceil(Nx / 2) - 1);
    [xGrid, yGrid] = meshgrid(x, y);
    radialSquared = (xGrid .^ 2 + yGrid .^ 2) / max(scale ^ 2, eps);
    kernel = (2 - radialSquared) .* exp(-0.5 * radialSquared);
    kernel = kernel - mean(kernel(:));
    kernel = kernel / max(norm(kernel(:)), eps);
end
