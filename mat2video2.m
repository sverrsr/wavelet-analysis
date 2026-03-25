function mat2video(data, videoName)
% MAT2VIDEO converts a 3D array (height, width, time) into a video.
% Each slice data(:,:,k) becomes one grayscale frame.

arguments
    data (:,:,:) {mustBeNumeric}
    videoName (1,:) char
end

workingDir = pwd;
[h,w,nFrames] = size(data);

% Check finite data
if ~all(isfinite(data(:)))
    error('Input contains NaN or Inf values.');
end

% Normalize safely to [0,1]
dmin = min(data(:));
dmax = max(data(:));

if dmax == dmin
    % Constant array -> write black frames
    data = zeros(size(data), 'double');
else
    data = (double(data) - dmin) ./ (dmax - dmin);
end

% Convert to uint8 explicitly. This is usually more robust.
data = uint8(255 * data);

% Safer codec choice if MPEG-4 gives trouble:
% v = VideoWriter(fullfile(workingDir, videoName + ".avi"), 'Motion JPEG AVI');

v = VideoWriter(fullfile(workingDir, videoName + ".mp4"), 'MPEG-4');
v.FrameRate = 30;
v.Quality = 95;

open(v);

for k = 1:nFrames
    frame = data(:,:,k);   % 2D uint8 grayscale frame
    writeVideo(v, frame);

    if mod(k,250)==0 || k==1 || k==nFrames
        fprintf('Processed %d/%d\n', k, nFrames);
    end
end

close(v);
disp('Done.');
end