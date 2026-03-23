clear;clc;

%% I. Read Video and Prepare Timestamps
%data = load ('..\data\SZ_VFD10p5Hz_TimeResolved_Run1_1920x1080_camCal_undistorted_test1.mat');
%data = load ('..\data\SZ_VFD10p5Hz_TimeResolved_Run1_60fps_test1.mat');
%data = load('..\data\SZ_VFD10p5Hz_TimeResolved_Run1_1920x1080.mat');
%data = load ('..\data\SZ_VFD10p5Hz_TimeResolved_Run1_30fps.mat');
%data = load('..\data\SZ_VFD10p5Hz_TimeResolved_Run1_30fps_calibrated.mat');

data = load('..\data\SZ_VFD10p5Hz_TimeResolved_Run1_30fps_25limit_May4th.mat'); % this is new as of May 4th. Green channel limit 25


%%
% video = data.filteredFramesGray;
% times = data.filteredTimeindeces;
% stamps = data.filteredTimestamps;

video = data.calibratedVideo;
times = data.timeIndeces;
stamps = data.timeStamps;
%% Don't ask
disp(stamps(38)-stamps(37)-1/30)
disp(stamps(212)-stamps(37)-1/30)
% disp(14/15-1/15)
% disp(299/15 - 1/15)
% 
% disp(41/45 - 1/45)
% disp(895/45 - 1/45)

%% Short snippet to get data on correct form
[height, width, numFrames] = size(video);

% height=1080;
% width=1920;
% numFrames=2241;

%% Preallocate a 3D matrix for the frames
% eta = zeros(height, width, numFrames, 'uint8'); % Use 'uint8' for grayscale images
% 
% %Populate the 3D matrix
% for t = 1:numFrames
%     eta(:, :, t) = video{t};
% end

eta = video;
disp('Data read and converted to correct form.');

%% MEAN SUBTRACTION TO REMOVE THE BLACK CEILING PANELS
% Convert eta to double, then subtract the mean frame.
mean_frame = mean(eta, 3); % 1080x1920 (double)
eta_meansub = double(eta) - mean_frame;

%% Wavelet Analysis and Filtering
% Parameters for wavelet filtering
%scales = 1:8;
%selected_scale = 8;
eccentricity_threshold = 0.85;
solidity_threshold = 0.6;
[x_dim, y_dim, ~] = size(eta);


%% Loop through different W_thr
selected_scale = 8;
W_thr_list = -16:2:60; % RANGE TO INVESTIGATE

for W_thr = W_thr_list

    tracks = struct('id', {}, 'centroids', {}, 'frames', {}, 'active', {}, 'timestamp', {});
    nextTrackId = 1;
    baseYShift = 16.022; % Base offset per time unit
    radiusFactor = 1; % NOT USED!
    maxSearchTime = 5; %cannot search for a structure if it hasn't been present for n timesteps
    
    detectionCount = zeros(1, numFrames);
    detectionTime  = zeros(1, numFrames);
    
    for t_index = 38:numFrames
        currentTime = times(t_index);
        currentStamp = stamps(t_index);
        currentStamp = currentStamp - stamps(37) - 1/30;
        %disp(currentTime)
        %if mod(currentTime,100)==0, disp(currentTime), end
        snapshot = eta_meansub(:, :, t_index);
        % Wavelet transform and filtering (same as before)
        cwt_result = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', selected_scale);
        wavelet_coefficients = cwt_result.cfs;
        % Save the raw wavelet coefficients for this frame.
        %wavelet_coeff_all(:, :, t_index) = wavelet_coefficients;
    
        mask = wavelet_coefficients > W_thr;
        %filtered_coefficients = wavelet_coefficients .* mask;
        connected_components = bwconncomp(mask);
        region_props = regionprops(connected_components, 'Area', 'Eccentricity', 'Solidity', 'Centroid');
        validIdx = find([region_props.Eccentricity] <= eccentricity_threshold & ...
        [region_props.Solidity] > solidity_threshold);
        %eccentric_regions = ismember(labelmatrix(connected_components), validIdx);
        %filtered_by_eccentricity = wavelet_coefficients .* eccentric_regions;
        %filtered_all_structures(:, :, t_index) = filtered_coefficients;
        %filtered_dimples(:, :, t_index) = filtered_by_eccentricity;
        % Extract centroids of valid regions
        if isempty(validIdx)
            centroids = [];
        else
            centroids = cat(1, region_props(validIdx).Centroid); % Each row: [x y]
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
                predicted = [tracks(j).centroids(end,1), ...
                tracks(j).centroids(end,2) + dt*baseYShift];
                allowedDistance = dt * baseYShift * radiusFactor;
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
        % figure(1); clf;
        % set(gcf, 'Position', [200, 200, 1000, 500]);
        % imagesc(wavelet_coefficients); colormap gray; hold on;
        % colorbar;
        % if ~isempty(centroids)
        % plot(centroids(:,1), centroids(:,2), 'bo', 'MarkerSize', 22);
        % for i = 1:numDetections
        %     % Label the detection with its associated track ID
        %     %text(centroids(i,1)+22, centroids(i,2)+22, num2str(detectionTrackIDs(i)), ...
        %     %'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold');
        % end
        % end
        % title(['Time = ' num2str(currentTime) ', Stamp = ' num2str(currentStamp)]);
        %set(gca, 'YDir', 'normal');
        %set(gca, 'XDir', 'reverse');
        % Capture frame and write to video
    
         % --- free memory from this frame before next iteration ---
        clear snapshot cwt_result wavelet_coefficients mask ...
              connected_components region_props validIdx eccentric_regions ...
              filtered_coefficients filtered_by_eccentricity
    end
    
    lifetimeThreshold = 0; % Adjust as needed (in number of frames)
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
    fname = sprintf('../createFigures/data/ceilingTrackingResults_38-2346_Wthr%d_s%d_Run3.mat', W_thr, selected_scale);
    save(fname, 'tracks', 'trackInfo', 'stamps', 'times', 'baseYShift', 'detectionTime', 'detectionCount');
    fprintf(' W_thr=%d  →  saved "%s"\n', W_thr, fname);

end

% %% for scale = 5
% 
% selected_scale = 5;
% W_thr_list = 10:2:50;
% 
% for W_thr = W_thr_list
% 
%     tracks = struct('id', {}, 'centroids', {}, 'frames', {}, 'active', {}, 'timestamp', {});
%     nextTrackId = 1;
%     baseYShift = 16.022; % Base offset per time unit, should be updated
%     radiusFactor = 1; % Adjusts the factor for allowed distance from predicted position will be baseYShift*radiusFactor
%     maxSearchTime = 5; %cannot search for a structure if it hasn't been present for n timesteps
% 
%     detectionCount = zeros(1, numFrames);
%     detectionTime  = zeros(1, numFrames);
% 
%     for t_index = 38:numFrames
%         currentTime = times(t_index);
%         currentStamp = stamps(t_index);
%         currentStamp = currentStamp - stamps(37) - 1/30;
%         %disp(currentTime)
%         snapshot = eta_meansub(:, :, t_index);
%         % Wavelet transform and filtering (same as before)
%         cwt_result = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', selected_scale);
%         wavelet_coefficients = cwt_result.cfs;
%         % Save the raw wavelet coefficients for this frame.
%         %wavelet_coeff_all(:, :, t_index) = wavelet_coefficients;
% 
%         mask = wavelet_coefficients > W_thr;
%         %filtered_coefficients = wavelet_coefficients .* mask;
%         connected_components = bwconncomp(mask);
%         region_props = regionprops(connected_components, 'Area', 'Eccentricity', 'Solidity', 'Centroid');
%         validIdx = find([region_props.Eccentricity] <= eccentricity_threshold & ...
%         [region_props.Solidity] > solidity_threshold);
%         %eccentric_regions = ismember(labelmatrix(connected_components), validIdx);
%         %filtered_by_eccentricity = wavelet_coefficients .* eccentric_regions;
%         %filtered_all_structures(:, :, t_index) = filtered_coefficients;
%         %filtered_dimples(:, :, t_index) = filtered_by_eccentricity;
%         % Extract centroids of valid regions
%         if isempty(validIdx)
%             centroids = [];
%         else
%             centroids = cat(1, region_props(validIdx).Centroid); % Each row: [x y]
%         end
%         numDetections = size(centroids, 1);
% 
%         % 2) STORE COUNT & TIME
%         detectionCount(t_index) = numDetections;
%         detectionTime (t_index) = currentStamp;
% 
%         % IV. Identify/Match Structures (Tracking)
%         numTracks = length(tracks);
%         costMatrix = Inf(numDetections, numTracks);
%         for i = 1:numDetections
%             for j = 1:numTracks
%                 if ~tracks(j).active
%                     continue; % Skip dead tracks
%                 end
%                 lastTime = tracks(j).frames(end);
%                 dt = currentTime - lastTime;
%                 % --- Skip if the time gap is more than maxSearchTime frames ---
%                 if dt > maxSearchTime
%                     continue; % costMatrix stays Inf, effectively ignoring this track
%                 end
%                 % Predicted position: same x, y shifted by dt*baseYShift
%                 predicted = [tracks(j).centroids(end,1), ...
%                 tracks(j).centroids(end,2) + dt*baseYShift];
%                 allowedDistance = dt * baseYShift * radiusFactor;
%                 % Compute distance from detection to predicted position
%                 d = norm(centroids(i,:) - predicted);
%                 if d <= allowedDistance
%                     costMatrix(i,j) = d;
%                 end
%             end
%         end
%         % Greedy assignment of detections to tracks
%         detectionTrackIDs = zeros(numDetections, 1);
%         assignments = [];
%         if ~isempty(costMatrix)
%         while true
%         [minVal, idx] = min(costMatrix(:));
%         if isinf(minVal), break; end
%             [detIdx, trackIdx] = ind2sub(size(costMatrix), idx);
%             assignments = [assignments; detIdx, trackIdx, minVal]; %#ok<AGROW>
%             detectionTrackIDs(detIdx) = tracks(trackIdx).id;
%             costMatrix(detIdx, :) = Inf;
%             costMatrix(:, trackIdx) = Inf;
%         end
%         end
%         % Update assigned tracks with new detections
%         if ~isempty(assignments)
%         for k = 1:size(assignments, 1)
%             detIdx = assignments(k, 1);
%             trackIdx = assignments(k, 2);
%             tracks(trackIdx).centroids(end+1, :) = centroids(detIdx, :);
%             tracks(trackIdx).frames(end+1) = currentTime;
%             % **Append the current stamp to the new 'timestamp' field**
%             tracks(trackIdx).timestamp(end+1) = currentStamp;
%         end
%         end
%         % ***** FIX: Mark tracks not updated in the current frame as inactive *****
%         for j = 1:length(tracks)
%             if tracks(j).active && tracks(j).frames(end) < currentTime
%                 tracks(j).active = false;
%             end
%         end
%         % Start new tracks for unassigned detections
%         for i = 1:numDetections
%             if detectionTrackIDs(i) == 0
%                 tracks(nextTrackId).id = nextTrackId;
%                 tracks(nextTrackId).centroids = centroids(i, :);
%                 tracks(nextTrackId).frames = currentTime;
%                 % **Initialize the 'timestamp' field with the current timestamp**
%                 tracks(nextTrackId).timestamp = currentStamp;
%                 tracks(nextTrackId).active = true;
%                 detectionTrackIDs(i) = nextTrackId;
%                 nextTrackId = nextTrackId + 1;
%             end
%         end
%         % Declare lost tracks as dead if not updated in current frame (and if they had multiple updates)
%         for j = 1:length(tracks)
%             if tracks(j).active && tracks(j).frames(end) < currentTime && numel(tracks(j).frames) >= 2
%                 tracks(j).active = false;
%             end
%         end
%         % === Optional Visualization === %
%         % figure(1); clf;
%         % set(gcf, 'Position', [200, 200, 1000, 500]);
%         % imagesc(wavelet_coefficients); colormap gray; hold on;
%         % colorbar;
%         % if ~isempty(centroids)
%         % plot(centroids(:,1), centroids(:,2), 'bo', 'MarkerSize', 22);
%         % for i = 1:numDetections
%         %     % Label the detection with its associated track ID
%         %     %text(centroids(i,1)+22, centroids(i,2)+22, num2str(detectionTrackIDs(i)), ...
%         %     %'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold');
%         % end
%         % end
%         % title(['Time = ' num2str(currentTime) ', Stamp = ' num2str(currentStamp)]);
%         %set(gca, 'YDir', 'normal');
%         %set(gca, 'XDir', 'reverse');
%         % Capture frame and write to video
% 
%          % --- free memory from this frame before next iteration ---
%         clear snapshot cwt_result wavelet_coefficients mask ...
%               connected_components region_props validIdx eccentric_regions ...
%               filtered_coefficients filtered_by_eccentricity
%     end
% 
%     lifetimeThreshold = 0; % Adjust as needed (in number of frames)
%     numTracks = length(tracks);
%     trackInfo = struct('id', {}, 'lifetime', {}, 'coordinates', {});
%     for i = 1:numTracks
%         if ~isempty(tracks(i).frames)
%         % Lifetime defined as number of frames the structure was tracked.
%         lifetime = numel(tracks(i).frames);
%         trackInfo(end+1) = struct('id', tracks(i).id, 'lifetime', lifetime, 'coordinates', {tracks(i).centroids});
%         end
%     end
%     %disp('Track Information:');
%     for i = 1:length(trackInfo)
%         if trackInfo(i).lifetime > lifetimeThreshold
%         %fprintf('Track %d: Lifetime = %d frames\n', trackInfo(i).id, trackInfo(i).lifetime);
%         %disp(trackInfo(i).coordinates);
%         end
%     end
% 
%     % Save for this threshold
%     fname = sprintf('../createFigures/data/ceilingTrackingResults_38-2346_Wthr%d_s%d_Run2.mat', W_thr, selected_scale);
%     save(fname, 'tracks', 'trackInfo', 'stamps', 'times', 'baseYShift', 'detectionTime', 'detectionCount');
%     fprintf(' W_thr=%d  →  saved "%s"\n', W_thr, fname);
% 
% end
% 
% %% for scale = 10
% 
% selected_scale = 10;
% W_thr_list = 40:3:90;
% 
% for W_thr = W_thr_list
% 
%     tracks = struct('id', {}, 'centroids', {}, 'frames', {}, 'active', {}, 'timestamp', {});
%     nextTrackId = 1;
%     baseYShift = 16.022; % Base offset per time unit, should be updated
%     radiusFactor = 1; % Adjusts the factor for allowed distance from predicted position will be baseYShift*radiusFactor
%     maxSearchTime = 5; %cannot search for a structure if it hasn't been present for n timesteps
% 
%     detectionCount = zeros(1, numFrames);
%     detectionTime  = zeros(1, numFrames);
% 
%     for t_index = 38:numFrames
%         currentTime = times(t_index);
%         currentStamp = stamps(t_index);
%         currentStamp = currentStamp - stamps(37) - 1/30;
%         %disp(currentTime)
%         snapshot = eta_meansub(:, :, t_index);
%         % Wavelet transform and filtering (same as before)
%         cwt_result = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', selected_scale);
%         wavelet_coefficients = cwt_result.cfs;
%         % Save the raw wavelet coefficients for this frame.
%         %wavelet_coeff_all(:, :, t_index) = wavelet_coefficients;
% 
%         mask = wavelet_coefficients > W_thr;
%         %filtered_coefficients = wavelet_coefficients .* mask;
%         connected_components = bwconncomp(mask);
%         region_props = regionprops(connected_components, 'Area', 'Eccentricity', 'Solidity', 'Centroid');
%         validIdx = find([region_props.Eccentricity] <= eccentricity_threshold & ...
%         [region_props.Solidity] > solidity_threshold);
%         %eccentric_regions = ismember(labelmatrix(connected_components), validIdx);
%         %filtered_by_eccentricity = wavelet_coefficients .* eccentric_regions;
%         %filtered_all_structures(:, :, t_index) = filtered_coefficients;
%         %filtered_dimples(:, :, t_index) = filtered_by_eccentricity;
%         % Extract centroids of valid regions
%         if isempty(validIdx)
%             centroids = [];
%         else
%             centroids = cat(1, region_props(validIdx).Centroid); % Each row: [x y]
%         end
%         numDetections = size(centroids, 1);
% 
%         % 2) STORE COUNT & TIME
%         detectionCount(t_index) = numDetections;
%         detectionTime (t_index) = currentStamp;
% 
%         % IV. Identify/Match Structures (Tracking)
%         numTracks = length(tracks);
%         costMatrix = Inf(numDetections, numTracks);
%         for i = 1:numDetections
%             for j = 1:numTracks
%                 if ~tracks(j).active
%                     continue; % Skip dead tracks
%                 end
%                 lastTime = tracks(j).frames(end);
%                 dt = currentTime - lastTime;
%                 % --- Skip if the time gap is more than maxSearchTime frames ---
%                 if dt > maxSearchTime
%                     continue; % costMatrix stays Inf, effectively ignoring this track
%                 end
%                 % Predicted position: same x, y shifted by dt*baseYShift
%                 predicted = [tracks(j).centroids(end,1), ...
%                 tracks(j).centroids(end,2) + dt*baseYShift];
%                 allowedDistance = dt * baseYShift * radiusFactor;
%                 % Compute distance from detection to predicted position
%                 d = norm(centroids(i,:) - predicted);
%                 if d <= allowedDistance
%                     costMatrix(i,j) = d;
%                 end
%             end
%         end
%         % Greedy assignment of detections to tracks
%         detectionTrackIDs = zeros(numDetections, 1);
%         assignments = [];
%         if ~isempty(costMatrix)
%         while true
%         [minVal, idx] = min(costMatrix(:));
%         if isinf(minVal), break; end
%             [detIdx, trackIdx] = ind2sub(size(costMatrix), idx);
%             assignments = [assignments; detIdx, trackIdx, minVal]; %#ok<AGROW>
%             detectionTrackIDs(detIdx) = tracks(trackIdx).id;
%             costMatrix(detIdx, :) = Inf;
%             costMatrix(:, trackIdx) = Inf;
%         end
%         end
%         % Update assigned tracks with new detections
%         if ~isempty(assignments)
%         for k = 1:size(assignments, 1)
%             detIdx = assignments(k, 1);
%             trackIdx = assignments(k, 2);
%             tracks(trackIdx).centroids(end+1, :) = centroids(detIdx, :);
%             tracks(trackIdx).frames(end+1) = currentTime;
%             % **Append the current stamp to the new 'timestamp' field**
%             tracks(trackIdx).timestamp(end+1) = currentStamp;
%         end
%         end
%         % ***** FIX: Mark tracks not updated in the current frame as inactive *****
%         for j = 1:length(tracks)
%             if tracks(j).active && tracks(j).frames(end) < currentTime
%                 tracks(j).active = false;
%             end
%         end
%         % Start new tracks for unassigned detections
%         for i = 1:numDetections
%             if detectionTrackIDs(i) == 0
%                 tracks(nextTrackId).id = nextTrackId;
%                 tracks(nextTrackId).centroids = centroids(i, :);
%                 tracks(nextTrackId).frames = currentTime;
%                 % **Initialize the 'timestamp' field with the current timestamp**
%                 tracks(nextTrackId).timestamp = currentStamp;
%                 tracks(nextTrackId).active = true;
%                 detectionTrackIDs(i) = nextTrackId;
%                 nextTrackId = nextTrackId + 1;
%             end
%         end
%         % Declare lost tracks as dead if not updated in current frame (and if they had multiple updates)
%         for j = 1:length(tracks)
%             if tracks(j).active && tracks(j).frames(end) < currentTime && numel(tracks(j).frames) >= 2
%                 tracks(j).active = false;
%             end
%         end
%         % === Optional Visualization === %
%         % figure(1); clf;
%         % set(gcf, 'Position', [200, 200, 1000, 500]);
%         % imagesc(wavelet_coefficients); colormap gray; hold on;
%         % colorbar;
%         % if ~isempty(centroids)
%         % plot(centroids(:,1), centroids(:,2), 'bo', 'MarkerSize', 22);
%         % for i = 1:numDetections
%         %     % Label the detection with its associated track ID
%         %     %text(centroids(i,1)+22, centroids(i,2)+22, num2str(detectionTrackIDs(i)), ...
%         %     %'Color', 'y', 'FontSize', 12, 'FontWeight', 'bold');
%         % end
%         % end
%         % title(['Time = ' num2str(currentTime) ', Stamp = ' num2str(currentStamp)]);
%         %set(gca, 'YDir', 'normal');
%         %set(gca, 'XDir', 'reverse');
%         % Capture frame and write to video
% 
%          % --- free memory from this frame before next iteration ---
%         clear snapshot cwt_result wavelet_coefficients mask ...
%               connected_components region_props validIdx eccentric_regions ...
%               filtered_coefficients filtered_by_eccentricity
%     end
% 
%     lifetimeThreshold = 0; % Adjust as needed (in number of frames)
%     numTracks = length(tracks);
%     trackInfo = struct('id', {}, 'lifetime', {}, 'coordinates', {});
%     for i = 1:numTracks
%         if ~isempty(tracks(i).frames)
%         % Lifetime defined as number of frames the structure was tracked.
%         lifetime = numel(tracks(i).frames);
%         trackInfo(end+1) = struct('id', tracks(i).id, 'lifetime', lifetime, 'coordinates', {tracks(i).centroids});
%         end
%     end
%     %disp('Track Information:');
%     for i = 1:length(trackInfo)
%         if trackInfo(i).lifetime > lifetimeThreshold
%         %fprintf('Track %d: Lifetime = %d frames\n', trackInfo(i).id, trackInfo(i).lifetime);
%         %disp(trackInfo(i).coordinates);
%         end
%     end
% 
%     % Save for this threshold
%     fname = sprintf('../createFigures/data/ceilingTrackingResults_38-2346_Wthr%d_s%d_Run2.mat', W_thr, selected_scale);
%     save(fname, 'tracks', 'trackInfo', 'stamps', 'times', 'baseYShift', 'detectionTime', 'detectionCount');
%     fprintf(' W_thr=%d  →  saved "%s"\n', W_thr, fname);
% 
% end