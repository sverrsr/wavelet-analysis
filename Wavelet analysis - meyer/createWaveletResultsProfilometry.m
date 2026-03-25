% SCRIPT TO CREATE RESULTS FOR COMPUTER VISION TECHNIQUE ON PROFILOMETRY
% DATA. FILES ARE SAVED AND CAN BE ACCESSED LATER!

% Written by Herman Martens Meyer

% Applied and edited by Sverre Steinnes Romøren
    % Applied to DNS instead of profilometry
        % allowedDistance = 10;
    % proximity-based instead of convection-based

%%
clearvars; clc;

% I Read Video and Prepare Timestamp
%load("D:\sverrsr\Documents\SYNC\flow-statistics\surfElev_RE2500_WEINF.mat");

%%
shadow_path = "E:\SYNC\path-trace-for-free-surface-flow\h5\WEINF\caustics_out_4p00pi.h5";
surfElev2 = h5read(shadow_path, "/caustic_frames");


%%
%
x = 1:256;
y = 1:256;

[X, Y] = meshgrid(x, y);

%%
% FULL DATASET. 
% The data should be on form [x_dim, y_dim, ~]
% surfElev is not, shadows are
%surfElev2   = permute(surfElev, [2, 3, 1]); 
surfElev2   = surfElev2(:,:, 1:2000); % Test size
data1       = surfElev2(:, :, :); % Cropped


%%
% eta_meansub = data1.surfData;
eta_meansub = data1;
numFrames = size(eta_meansub, 3);

%% II. Wavelet Analysis and Filtering
% Parameters for wavelet filtering
% scales = 1:8;
% selected_scale = 8;
eccentricity_threshold  = 0.85; % 1 line 0 circle
solidity_threshold      = 0.6;
[x_dim, y_dim, ~]       = size(eta_meansub);
stamps                  = zeros(1, numFrames);

%% Loop through different W_thr

selected_scale  = 1; % Dette er på en måte blur
W_thr_list = -1e-4;
%W_thr_list = [-5e-4 -2e-4 -1e-4 -5e-5 0 5e-5 1e-4]; %-5e-4

fprintf('Going though different W_thr')

for W_thr = W_thr_list

    
    tracks = struct('id', {}, 'centroids', {}, 'frames', {}, 'active', {}, 'timestamp', {});
    nextTrackId = 1;
    baseYShift = 250*(1/45); % Base offset per time unit. ACCORDING TO MEAN FLOW
    radiusFactor = 1; 
    maxSearchTime = 5; %cannot search for a structure if it hasn't been present for n timesteps
    
    detectionCount = zeros(1, numFrames);
    detectionTime  = zeros(1, numFrames);
   
    binaryDimples = false(x_dim, y_dim, numFrames);
    binaryAll = false(x_dim, y_dim, numFrames);

    for t_index = 1:numFrames
        currentTime = (t_index);
        %currentStamp = stamps(t_index);
        currentStamp = (t_index+35)/45 - 1/45; % accounting the correct startTime
        %disp(t_index)
        % store it
        stamps(t_index)       =  currentStamp;
        %if mod(currentTime,100)==0, disp(currentTime), end
        snapshot = eta_meansub(:, :, t_index);
        % Wavelet transform and filtering (same as before)
        cwt_result = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', selected_scale);
        wavelet_coefficients = cwt_result.cfs;
        % Crop to original image size
        wavelet_coefficients = wavelet_coefficients(1:x_dim, 1:y_dim);
        % Save the raw wavelet coefficients for this frame.
        %wavelet_coeff_all(:, :, t_index) = wavelet_coefficients;
    
        mask = wavelet_coefficients < W_thr;
        %filtered_coefficients = wavelet_coefficients .* mask;
        connected_components = bwconncomp(mask);
        region_props = regionprops(connected_components, 'Area', 'Eccentricity', 'Solidity', 'Centroid');
        validIdx = find([region_props.Eccentricity] <= eccentricity_threshold & ...
        [region_props.Solidity] > solidity_threshold);

        labelImg = labelmatrix(connected_components);
        dimpleMask = ismember(labelImg, validIdx);   % only accepted regions
        binaryDimples(:, :, t_index) = dimpleMask;

        binaryAll(:, :, t_index) = mask;

        %eccentric_regions = ismember(labelmatrix(connected_components), validIdx);
        %filtered_by_eccentricity = wavelet_coefficients .* eccentric_regions;
        %filtered_all_structures(:, :, t_index) = filtered_coefficients;
        %filtered_dimples(:, :, t_index) = filtered_by_eccentricity;
        % Extract centroids of valid regions
        if isempty(validIdx)
            centroids = [];
        else
            %centroids = cat(1, region_props(validIdx).Centroid); % Each row: [x y]
            centroids_px = cat(1, region_props(validIdx).Centroid);
            x_coords = interp2(1:size(X,2), 1:size(X,1), X, centroids_px(:,1), centroids_px(:,2));
            y_coords = interp2(1:size(Y,2), 1:size(Y,1), Y, centroids_px(:,1), centroids_px(:,2));
            centroids = [x_coords, y_coords];
        end
        numDetections = size(centroids, 1);
    
        % STORE COUNT & TIME
        detectionCount(t_index) = numDetections;
        detectionTime (t_index) = currentStamp;
            
        % Identify/Match Structures (Tracking)
        numTracks = length(tracks);
        costMatrix = Inf(numDetections, numTracks);
        for i = 1:numDetections
            for j = 1:numTracks
                if ~tracks(j).active
                    continue; % Skip dead tracks
                end
                lastTime = tracks(j).frames(end);
                dt = currentTime - lastTime;
                % --- Skip if the time gap is more than maxSearchTime frames ---
                if dt > maxSearchTime
                    continue; % costMatrix stays Inf, effectively ignoring this track
                end
                % Predicted position: same x, y shifted by dt*baseYShift
                predicted = tracks(j).centroids(end,:);
                allowedDistance = 50;   % choose a fixed search radius in pixels
                % Compute distance from detection to predicted position
                d = norm(centroids(i,:) - predicted);
                if d <= allowedDistance
                    costMatrix(i,j) = d;
                end
            end
        end
        % Greedy assignment of detections to tracks
        detectionTrackIDs = zeros(numDetections, 1);
        assignments = [];
        if ~isempty(costMatrix)
        while true
        [minVal, idx] = min(costMatrix(:));
        if isinf(minVal), break; end
            [detIdx, trackIdx] = ind2sub(size(costMatrix), idx);
            assignments = [assignments; detIdx, trackIdx, minVal]; %#ok<AGROW>
            detectionTrackIDs(detIdx) = tracks(trackIdx).id;
            costMatrix(detIdx, :) = Inf;
            costMatrix(:, trackIdx) = Inf;
        end
        end
        % Update assigned tracks with new detections
        if ~isempty(assignments)
        for k = 1:size(assignments, 1)
            detIdx = assignments(k, 1);
            trackIdx = assignments(k, 2);
            tracks(trackIdx).centroids(end+1, :) = centroids(detIdx, :);
            tracks(trackIdx).frames(end+1) = currentTime;
            % **Append the current stamp to the new 'timestamp' field**
            tracks(trackIdx).timestamp(end+1) = currentStamp;
        end
        end
        % ***** FIX: Mark tracks not updated in the current frame as inactive *****
        for j = 1:length(tracks)
            if tracks(j).active && tracks(j).frames(end) < currentTime
                tracks(j).active = false;
            end
        end
        % Start new tracks for unassigned detections
        for i = 1:numDetections
            if detectionTrackIDs(i) == 0
                tracks(nextTrackId).id = nextTrackId;
                tracks(nextTrackId).centroids = centroids(i, :);
                tracks(nextTrackId).frames = currentTime;
                % **Initialize the 'timestamp' field with the current timestamp**
                tracks(nextTrackId).timestamp = currentStamp;
                tracks(nextTrackId).active = true;
                detectionTrackIDs(i) = nextTrackId;
                nextTrackId = nextTrackId + 1;
            end
        end
        % Declare lost tracks as dead if not updated in current frame (and if they had multiple updates)
        for j = 1:length(tracks)
            if tracks(j).active && tracks(j).frames(end) < currentTime && numel(tracks(j).frames) >= 2
                tracks(j).active = false;
            end
        end

        % === Optional Visualization === %
        figure(1); clf;
        set(gcf, 'Position', [200, 200, 1000, 500]);
        imagesc(wavelet_coefficients); colormap gray; hold on; axis image
        colorbar;
        if ~isempty(centroids)
        plot(centroids(:,1), centroids(:,2), 'bo', 'MarkerSize', 22);
        for i = 1:numDetections
            % Label the detection with its associated track ID
            %text(centroids(i,1)+22, centroids(i,2)+22, num2str(detectionTrackIDs(i)), ...
            %'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold');
        end
        end
        title(['Time = ' num2str(currentTime) ', Stamp = ' num2str(currentStamp)]);
        set(gca, 'YDir', 'normal');
        set(gca, 'XDir', 'reverse');
        % Capture frame and write to video
    
         % --- free memory from this frame before next iteration ---
        clear snapshot cwt_result wavelet_coefficients mask ...
              connected_components region_props validIdx eccentric_regions ...
              filtered_coefficients filtered_by_eccentricity
    end
    
    % τmin = 0.166* 5.41 = 0,89806 - Babiker 23
    % dt = 0.06 - time / frame
    % frames = τmin / dt = 15 frames
    lifetimeThreshold = 15; % Adjust as needed (in number of frames) 
    numTracks = length(tracks);
    trackInfo = struct('id', {}, 'lifetime', {}, 'coordinates', {});
    for i = 1:numTracks
        if ~isempty(tracks(i).frames)
        % Lifetime defined as number of frames the structure was tracked.
        lifetime = numel(tracks(i).frames);
        trackInfo(end+1) = struct('id', tracks(i).id, 'lifetime', lifetime, 'coordinates', {tracks(i).centroids});
        end
    end
    %disp('Track Information:');
    for i = 1:length(trackInfo)
        if trackInfo(i).lifetime > lifetimeThreshold
        %fprintf('Track %d: Lifetime = %d frames\n', trackInfo(i).id, trackInfo(i).lifetime);
        %disp(trackInfo(i).coordinates);
        end
    end

    % Save for this threshold
    outDir = fullfile(pwd, 'createFigures', 'data');
    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end
    
    fname = fullfile(outDir, sprintf('profilometryTrackingResults_1-10800_Wthr%.3f_s%d_Run2.mat', W_thr, selected_scale));
    save(fname, 'tracks', 'trackInfo', 'stamps', 'baseYShift', ...
           'detectionTime', 'detectionCount', 'binaryDimples');
    fprintf(' W_thr=%.3f  →  saved "%s"\n', W_thr, fname)
    
    % Save vids
    name1 = sprintf('binaryAll_Wthr_%g', W_thr);
    name2 = sprintf('binaryDimples_Wthr_%g', W_thr);

    mat2video(binaryAll, name1);
    mat2video(binaryDimples, name2);

    % Save vids
    name1 = sprintf('binaryAll_Circled_%g', W_thr);
    name2 = sprintf('binaryDimples_Circled_%g', W_thr);
    circle_dimples(binaryAll, name1);
    circle_dimples(binaryDimples, name2);

end



