function circle_dimples(binaryMask, videoName)
%CIRCLE_DIMPLES Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    binaryMask 
    videoName
end


[nx, ny, nt] = size(binaryMask);
[X, Y] = meshgrid(1:ny, 1:nx);

circleMask3D = false(nx, ny, nt);

se = strel('disk', 2);   % dilation size

for t = 1:nt
    BW = imdilate(binaryMask(:,:,t), se);              % dilate one frame
    CC = bwconncomp(BW);
    S = regionprops(CC, 'Centroid', 'Area');

    out = false(nx, ny);

    for k = 1:numel(S)
        c = S(k).Centroid;
        % If you want every detection to have the same circle radius,
        % replace r with a constant size
        r = sqrt(S(k).Area/pi);                    % area-matched circle
        out = out | ((X-c(1)).^2 + (Y-c(2)).^2 <= r^2);
        out = bwperim(out);
    end

    circleMask3D(:,:,t) = out;
end

mat2video(circleMask3D, videoName);

% save('circleMask3D.mat', 'circleMask3D')
end