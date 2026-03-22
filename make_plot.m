function make_plot(eta, detTracks, trueMask, tList, outVideoFile)
%==========================================================================
% Make plots in the same style as the Babiker et al. supplemental video:
%   (1) Surface elevation only
%   (2) Surface elevation + tracked vortices (circles) + track trails (dotted)
%   (3) Tracks (circles + dotted trails) + true vortices (blue areas) on white bg
%
% INPUTS
%   eta        : Ny x Nx x Nt (double/single) surface elevation
%   detTracks  : struct array from SOC code (must have fields):
%                  .t   (frame indices)
%                  .pix (cell array of pixel-idx per frame)
%   trueMask   : Ny x Nx x Nt logical mask for "true" vortices (lambda2<thr)
%                (can be [] if you want to skip blue overlay)
%   tList      : vector of frame indices to render (e.g., 1:Nt or 38:Nt)
%   outVideoFile : '' to skip video, or a filename (e.g. 'soc_movie.mp4')
%
% NOTES
% - Circles radius is set from region area: r = sqrt(A/pi)
% - Track trails show last Ntrail points (dotted).
% - Axes are image coordinates (x=columns, y=rows), like your centroid output.
%==========================================================================

if nargin < 5, outVideoFile = ''; end
if nargin < 4 || isempty(tList), tList = 1:size(eta,3); end

[Ny,Nx,Nt] = size(eta);
hasTrue = ~isempty(trueMask);

% ---- style knobs (tweak to match your exact look) ----
Ntrail   = 20;     % how many previous positions to show in dotted trail
alphaTrue = 0.35;  % blue overlay alpha
circleMinR = 2;    % minimum circle radius (pixels) so tiny detections still show

% Figure layout similar to the video (black background with 3 panels)
fig = figure('Color','k','Units','pixels','Position',[100 100 1200 420]);
tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');

% Optional video writer
doVideo = ~isempty(outVideoFile);
if doVideo
    vw = VideoWriter(outVideoFile,'MPEG-4');
    vw.FrameRate = 25;  % adjust to your data
    open(vw);
end

for tt = 1:numel(tList)
    t = tList(tt);
    if t < 1 || t > Nt, continue; end

    % Collect detections at this frame: centroids, radii, and trails
    [C, R, trails] = detections_at_t(detTracks, t, Ny, Nx, Ntrail);

    % ========== Panel 1: Surface elevation only ==========
    nexttile(1); cla;
    imagesc(eta(:,:,t)); axis image off;
    title('Surface elevation only','Color','w','FontSize',11,'FontWeight','normal');

    % ========== Panel 2: Elevation + tracked vortices ==========
    nexttile(2); cla;
    imagesc(eta(:,:,t)); axis image off; hold on;
    title('Surface elevation with tracked vortices','Color','w','FontSize',11,'FontWeight','normal');

    % dotted trails (white)
    for k = 1:numel(trails)
        xy = trails{k};
        plot(xy(:,1), xy(:,2), ':', 'Color','w', 'LineWidth', 1.0);
    end
    % circles (white)
    if ~isempty(C)
        viscircles(C, R, 'Color','w', 'LineWidth', 0.8);
    end
    hold off;

    % ========== Panel 3: Tracks + true vortices ==========
    nexttile(3); cla;
    % white background
    image(ones(Ny,Nx,3)); axis image off; hold on;
    title('Tracked vortices (circles) and true vortices (blue areas)','Color','w','FontSize',11,'FontWeight','normal');

    % true vortices as blue overlay on white
    if hasTrue
        T = trueMask(:,:,t);
        if ~islogical(T), T = T ~= 0; end
        blue = cat(3, zeros(Ny,Nx), zeros(Ny,Nx), ones(Ny,Nx)); % pure blue
        h = image(blue);
        set(h, 'AlphaData', alphaTrue * double(T));
    end

    % dotted trails (gray/black-ish)
    trailColor = [0.35 0.35 0.35];
    for k = 1:numel(trails)
        xy = trails{k};
        plot(xy(:,1), xy(:,2), ':', 'Color',trailColor, 'LineWidth', 1.0);
    end

    % circles (gray)
    if ~isempty(C)
        viscircles(C, R, 'Color',trailColor, 'LineWidth', 0.8);
    end
    hold off;

    drawnow;

    if doVideo
        fr = getframe(fig);
        writeVideo(vw, fr);
    end
end

if doVideo
    close(vw);
    fprintf('Saved video: %s\n', outVideoFile);
end

end

%==========================================================================
% Helper: extract detections (centroids/radii) and trails for a given frame t
%==========================================================================
function [C, R, trails] = detections_at_t(tracks, t, Ny, Nx, Ntrail)

C = [];
R = [];
trails = {};

for i = 1:numel(tracks)
    % tracks(i).t is a vector of frame indices
    idx = find(tracks(i).t == t, 1);
    if isempty(idx), continue; end

    pix = tracks(i).pix{idx};
    if isempty(pix), continue; end

    % centroid from pixels (image coords)
    [yy, xx] = ind2sub([Ny Nx], pix);
    cx = mean(double(xx));
    cy = mean(double(yy));

    % radius from area
    A = numel(pix);
    r = sqrt(A/pi);
    r = max(r, 2); % minimum radius for visibility

    C(end+1,:) = [cx cy]; %#ok<AGROW>
    R(end+1,1) = r;       %#ok<AGROW>

    % trail: last Ntrail centroids up to this time
    j0 = max(1, idx - Ntrail + 1);
    xy = zeros(idx-j0+1, 2);
    c = 0;
    for j = j0:idx
        pj = tracks(i).pix{j};
        if isempty(pj), continue; end
        [yyj, xxj] = ind2sub([Ny Nx], pj);
        c = c + 1;
        xy(c,:) = [mean(double(xxj)), mean(double(yyj))];
    end
    xy = xy(1:c,:);
    if size(xy,1) >= 2
        trails{end+1} = xy; %#ok<AGROW>
    end
end

end