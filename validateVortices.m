function validation = validateVortices(featureFrames, tracks, inputData, config)
%VALIDATEVORTICES Validate candidate structures against lambda_2 and beta^2.
% Scars use bwskel-based geodesic centers to avoid centroid drift for highly
% elongated structures. Metrics are aggregated using area S(t), because area
% is more robust than count for fragmented scar signatures.

    arguments
        featureFrames (:) struct
        tracks (:) struct
        inputData struct
        config struct
    end

    Nt = numel(featureFrames);
    [Ny, Nx, ~] = size(inputData.etaStack);
    lambda2Stack = [];
    if isfield(inputData, 'lambda2Stack')
        lambda2Stack = inputData.lambda2Stack;
    elseif isfield(inputData, 'velocityGradientStack')
        lambda2Stack = computeLambda2FromGradients(inputData.velocityGradientStack);
    end

    scarArea = zeros(Nt, 1);
    dimpleArea = zeros(Nt, 1);
    scarCount = zeros(Nt, 1);
    dimpleCount = zeros(Nt, 1);
    validatedScarArea = zeros(Nt, 1);
    validatedDimpleArea = zeros(Nt, 1);

    for frameIndex = 1:Nt
        features = featureFrames(frameIndex).features;
        for featureIndex = 1:numel(features)
            feature = features(featureIndex);
            if ~feature.isPersistent
                continue;
            end

            mask = false(Ny, Nx);
            mask(feature.pixelIdxList) = true;

            geodesicCenter = [];
            if strcmp(feature.classLabel, 'scar')
                geodesicCenter = findScarGeodesicCenter(mask, feature.centroid);
                scarArea(frameIndex) = scarArea(frameIndex) + feature.area;
                scarCount(frameIndex) = scarCount(frameIndex) + 1;
            else
                dimpleArea(frameIndex) = dimpleArea(frameIndex) + feature.area;
                dimpleCount(frameIndex) = dimpleCount(frameIndex) + 1;
            end

            isLambda2Validated = true;
            if ~isempty(lambda2Stack)
                lambda2Values = lambda2Stack(:, :, frameIndex);
                isLambda2Validated = mean(lambda2Values(feature.pixelIdxList), 'omitnan') < config.lambda2Threshold;
            end

            featureFrames(frameIndex).features(featureIndex).geodesicCenter = geodesicCenter;
            featureFrames(frameIndex).features(featureIndex).lambda2Validated = isLambda2Validated;

            if isLambda2Validated
                if strcmp(feature.classLabel, 'scar')
                    validatedScarArea(frameIndex) = validatedScarArea(frameIndex) + feature.area;
                else
                    validatedDimpleArea(frameIndex) = validatedDimpleArea(frameIndex) + feature.area;
                end
            end
        end
    end

    referenceLagFrames = round(0.8 * config.Tinf / config.dt);
    [beta2ScarSignal, beta2DimpleSignal] = buildReferenceSignals(inputData, config, referenceLagFrames, Nt);

    validation = struct();
    validation.featureFrames = featureFrames;
    validation.lambda2Threshold = config.lambda2Threshold;
    validation.Lnu = config.Lnu;
    validation.metrics = struct();
    validation.metrics.referenceLagFrames = referenceLagFrames;
    validation.metrics.scarArea = scarArea;
    validation.metrics.dimpleArea = dimpleArea;
    validation.metrics.scarCount = scarCount;
    validation.metrics.dimpleCount = dimpleCount;
    validation.metrics.validatedScarArea = validatedScarArea;
    validation.metrics.validatedDimpleArea = validatedDimpleArea;
    validation.metrics.beta2ScarSignal = beta2ScarSignal;
    validation.metrics.beta2DimpleSignal = beta2DimpleSignal;
    validation.metrics.scarAreaNCC = computeNormalizedCrossCorrelation(validatedScarArea, beta2ScarSignal);
    validation.metrics.dimpleAreaNCC = computeNormalizedCrossCorrelation(validatedDimpleArea, beta2DimpleSignal);
end

function lambda2Stack = computeLambda2FromGradients(gradientStack)
%COMPUTELAMBDA2FROMGRADIENTS Compute lambda_2 for a 2x2 in-plane gradient tensor.
% The accepted input shapes are [Ny x Nx x 4 x Nt] using the order
% [dudx, dudy, dvdx, dvdy] or [Ny x Nx x 2 x 2 x Nt].

    dims = ndims(gradientStack);
    if dims == 4 && size(gradientStack, 3) == 4
        dudx = gradientStack(:, :, 1, :);
        dudy = gradientStack(:, :, 2, :);
        dvdx = gradientStack(:, :, 3, :);
        dvdy = gradientStack(:, :, 4, :);
    elseif dims == 5 && all(size(gradientStack, 3:4) == [2 2])
        dudx = gradientStack(:, :, 1, 1, :);
        dudy = gradientStack(:, :, 1, 2, :);
        dvdx = gradientStack(:, :, 2, 1, :);
        dvdy = gradientStack(:, :, 2, 2, :);
    else
        error('velocityGradientStack must be [Ny x Nx x 4 x Nt] or [Ny x Nx x 2 x 2 x Nt].');
    end

    S11 = dudx;
    S22 = dvdy;
    S12 = 0.5 * (dudy + dvdx);
    O12 = 0.5 * (dudy - dvdx);

    M11 = S11 .^ 2 + S12 .^ 2 - O12 .^ 2;
    M22 = S22 .^ 2 + S12 .^ 2 - O12 .^ 2;
    M12 = S12 .* (S11 + S22);

    traceM = M11 + M22;
    discriminant = sqrt(max(traceM .^ 2 - 4 * (M11 .* M22 - M12 .^ 2), 0));
    lambda2Stack = 0.5 * (traceM - discriminant);
    lambda2Stack = squeeze(lambda2Stack);
end

function geodesicCenter = findScarGeodesicCenter(mask, fallbackCentroid)
%FINDSCARGEODESICCENTER Use bwskel and geodesic distance to find a scar center.
% The selected pixel maximizes the minimum geodesic distance to the detected
% skeleton endpoints, which is a practical approximation to the geodesic
% center required for elongated scars.

    skeleton = bwskel(mask);
    if ~any(skeleton(:))
        geodesicCenter = fallbackCentroid;
        return;
    end

    endpoints = bwmorph(skeleton, 'endpoints');
    [endY, endX] = find(endpoints);
    [skelY, skelX] = find(skeleton);

    if numel(endX) < 2
        centerIndex = round(numel(skelX) / 2);
        geodesicCenter = [skelX(centerIndex), skelY(centerIndex)];
        return;
    end

    minDistance = inf(numel(skelX), 1);
    for endpointIndex = 1:numel(endX)
        distMap = bwdistgeodesic(skeleton, endX(endpointIndex), endY(endpointIndex), 'quasi-euclidean');
        endpointDistances = distMap(sub2ind(size(mask), skelY, skelX));
        minDistance = min(minDistance, endpointDistances);
    end

    [~, centerIndex] = max(minDistance);
    geodesicCenter = [skelX(centerIndex), skelY(centerIndex)];
end

function [beta2ScarSignal, beta2DimpleSignal] = buildReferenceSignals(inputData, config, referenceLagFrames, Nt)
%BUILDREFERENCESIGNALS Build lagged beta^2 references for scars and dimples.
% Scars sample beta^2 at z = L_nu, whereas dimples use a Gaussian depth
% weighting to emulate the surface-attached decay of vertical vortices.

    if isfield(inputData, 'beta2ByDepth') && isfield(inputData, 'zLevels')
        beta2ByDepth = squeeze(inputData.beta2ByDepth);
        if isvector(beta2ByDepth)
            beta2ByDepth = reshape(beta2ByDepth, [], Nt);
        end
        if size(beta2ByDepth, 2) ~= Nt
            beta2ByDepth = reshape(beta2ByDepth, [], Nt);
        end

        zLevels = inputData.zLevels(:);
        [~, scarDepthIndex] = min(abs(zLevels - config.Lnu));
        beta2ScarSignal = beta2ByDepth(scarDepthIndex, :).';

        dimpleWeights = exp(-0.5 * (zLevels / config.gaussianDecaySigma) .^ 2);
        dimpleWeights = dimpleWeights / sum(dimpleWeights);
        beta2DimpleSignal = (dimpleWeights.' * beta2ByDepth).';
    elseif isfield(inputData, 'beta2Stack')
        beta2Surface = squeeze(mean(mean(inputData.beta2Stack, 1, 'omitnan'), 2, 'omitnan'));
        beta2ScarSignal = beta2Surface(:);
        beta2DimpleSignal = beta2Surface(:);
    else
        beta2ScarSignal = zeros(Nt, 1);
        beta2DimpleSignal = zeros(Nt, 1);
    end

    beta2ScarSignal = circshift(beta2ScarSignal(:), referenceLagFrames);
    beta2DimpleSignal = circshift(beta2DimpleSignal(:), referenceLagFrames);
end

function nccValue = computeNormalizedCrossCorrelation(signalA, signalB)
%COMPUTENORMALIZEDCROSSCORRELATION Return the zero-lag NCC between signals.

    signalA = signalA(:);
    signalB = signalB(:);
    validMask = isfinite(signalA) & isfinite(signalB);
    signalA = signalA(validMask);
    signalB = signalB(validMask);

    if isempty(signalA) || std(signalA) == 0 || std(signalB) == 0
        nccValue = NaN;
        return;
    end

    signalA = signalA - mean(signalA);
    signalB = signalB - mean(signalB);
    nccValue = sum(signalA .* signalB) / sqrt(sum(signalA .^ 2) * sum(signalB .^ 2));
end
