function [tracks, featureFrames, diagnostics] = trackFeatures(featureFrames, config)
%TRACKFEATURES Link detections through time using overlapping connected masks.
% Features are linked frame-to-frame when they share enough pixels and match
% class labels. Persistence is enforced through the minimum lifetime
% tau_min = 0.166 * T_inf.

    arguments
        featureFrames (:) struct
        config struct
    end

    Nt = numel(featureFrames);
    minLifetimeFrames = max(1, ceil(config.tauMin / config.dt));

    tracks = struct( ...
        'id', {}, ...
        'classLabel', {}, ...
        'frames', {}, ...
        'featureIndices', {}, ...
        'areas', {}, ...
        'centroids', {}, ...
        'lifetimeFrames', {}, ...
        'isPersistent', {});

    activeTrackIds = [];
    nextTrackId = 1;

    for frameIndex = 1:Nt
        currentFeatures = featureFrames(frameIndex).features;
        if isempty(currentFeatures)
            activeTrackIds = [];
            continue;
        end

        assignedTrackIds = nan(numel(currentFeatures), 1);
        availableTrackIds = activeTrackIds;

        for featureIndex = 1:numel(currentFeatures)
            currentFeature = currentFeatures(featureIndex);
            bestTrackId = NaN;
            bestOverlap = -Inf;

            for candidateIdx = 1:numel(availableTrackIds)
                candidateTrackId = availableTrackIds(candidateIdx);
                candidateTrack = tracks(candidateTrackId);
                if ~strcmp(candidateTrack.classLabel, currentFeature.classLabel)
                    continue;
                end

                previousFrame = candidateTrack.frames(end);
                if previousFrame ~= frameIndex - 1
                    continue;
                end

                previousFeatureIndex = candidateTrack.featureIndices(end);
                previousFeature = featureFrames(previousFrame).features(previousFeatureIndex);
                overlapScore = computeOverlapScore(previousFeature.pixelIdxList, currentFeature.pixelIdxList);

                if overlapScore > bestOverlap
                    bestOverlap = overlapScore;
                    bestTrackId = candidateTrackId;
                end
            end

            if ~isnan(bestTrackId) && bestOverlap >= config.overlapThreshold
                assignedTrackIds(featureIndex) = bestTrackId;
                track = tracks(bestTrackId);
                track.frames(end + 1) = frameIndex;
                track.featureIndices(end + 1) = featureIndex;
                track.areas(end + 1) = currentFeature.area;
                track.centroids(end + 1, :) = currentFeature.centroid;
                tracks(bestTrackId) = track;
                availableTrackIds(availableTrackIds == bestTrackId) = [];
            else
                newTrack = struct( ...
                    'id', nextTrackId, ...
                    'classLabel', currentFeature.classLabel, ...
                    'frames', frameIndex, ...
                    'featureIndices', featureIndex, ...
                    'areas', currentFeature.area, ...
                    'centroids', currentFeature.centroid, ...
                    'lifetimeFrames', 1, ...
                    'isPersistent', false);
                tracks(nextTrackId) = newTrack; %#ok<AGROW>
                assignedTrackIds(featureIndex) = nextTrackId;
                nextTrackId = nextTrackId + 1;
            end
        end

        for featureIndex = 1:numel(currentFeatures)
            featureFrames(frameIndex).features(featureIndex).trackId = assignedTrackIds(featureIndex);
        end

        activeTrackIds = assignedTrackIds(~isnan(assignedTrackIds)).';
    end

    for trackIndex = 1:numel(tracks)
        tracks(trackIndex).lifetimeFrames = numel(tracks(trackIndex).frames);
        tracks(trackIndex).isPersistent = tracks(trackIndex).lifetimeFrames >= minLifetimeFrames;
        for sampleIndex = 1:numel(tracks(trackIndex).frames)
            frameIndex = tracks(trackIndex).frames(sampleIndex);
            featureIndex = tracks(trackIndex).featureIndices(sampleIndex);
            featureFrames(frameIndex).features(featureIndex).isPersistent = tracks(trackIndex).isPersistent;
        end
    end

    diagnostics = struct();
    diagnostics.minLifetimeFrames = minLifetimeFrames;
    diagnostics.totalTracks = numel(tracks);
    diagnostics.persistentTracks = nnz([tracks.isPersistent]);
    diagnostics.meanTrackLifetimeFrames = mean([tracks.lifetimeFrames], 'omitnan');
end

function overlapScore = computeOverlapScore(previousPixelIdx, currentPixelIdx)
%COMPUTEOVERLAPSCORE Normalized intersection-over-min-area overlap measure.
% This overlap criterion is robust for fragmented scars because the
% denominator is the smaller of the two connected-component areas.

    if isempty(previousPixelIdx) || isempty(currentPixelIdx)
        overlapScore = 0;
        return;
    end

    overlapCount = numel(intersect(previousPixelIdx, currentPixelIdx));
    overlapScore = overlapCount / max(min(numel(previousPixelIdx), numel(currentPixelIdx)), 1);
end
