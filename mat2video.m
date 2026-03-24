function mat2video(matFile, videoName)
% MAT2VIDEO converts a 3D array (x, matFile, time) into an MP4 video.
% Each slice matFile(:,:,k) is written directly as a grayscale frame.
% No normalization, clamping, or preprocessing is applied.
% The video is saved in the current working directory.

arguments
    matFile
    videoName
end


workingDir = pwd;

[~,~,nFrames] = size(matFile);

% Rescale
gmin = min(matFile(:));
gmax = max(matFile(:));
matFile = (matFile - gmin) ./ (gmax - gmin);

v = VideoWriter(fullfile(workingDir, videoName), 'MPEG-4');
v.FrameRate = 30;
v.Quality = 100;

open(v)

for k = 1:nFrames

    frame = matFile(:,:,k);   % nothing done to the frame
    writeVideo(v, frame)


    if mod(k,250)==0 || k==1 || k==nFrames
        fprintf("Processed %d/%d\n", k, nFrames);
    end
end

close(v)

disp("Done.")

end

