function results = wavelet(inputs)
%==========================================================================
% SOC wavelet vortex detection + lambda2 "true vortex" verification
% Implements the steps/thresholds you listed, using Wavelet Toolbox.
%
% REQUIRED INPUTS (in struct "inputs"):
%   inputs.eta      : Ny x Nx x Nt surface elevation field eta(x,y,t)
%   inputs.dtau     : time step between frames (Delta tau), e.g. 0.060
%   inputs.Tinf     : integral time scale, e.g. 5.41
%
% FOR lambda2 VERIFICATION (also required if you want full Step 5):
%   Either provide velocity gradients at the surface:
%     inputs.grad : Ny x Nx x 3 x 3 x Nt, where grad(:,:,i,j,t) = d u_i / d x_j
%       ordering: u1=u, u2=v, u3=w ; x1=x, x2=y, x3=z
%   OR provide velocities + spacings + optional z-derivatives (less exact):
%     inputs.U, inputs.V, inputs.W : Ny x Nx x Nt surface velocity components
%     inputs.dx, inputs.dy         : grid spacing in x,y
%     (optional) inputs.dUdz, inputs.dVdz, inputs.dWdz : Ny x Nx x Nt
%       If you do NOT provide z-derivatives, they are set to zero (not exact).
%
% OUTPUT (struct "results"):
%   results.detTracksAll        : all wavelet tracks (pre-filter)
%   results.detTracksFiltered   : wavelet tracks passing tau_min + ecc filters
%   results.trueTracksFiltered  : lambda2 true tracks passing tau_min + area>=2
%   results.metrics             : precision/recall etc (if lambda2 available)
%   results.coverage            : coverage per frame for W > Wthr
%
% Thresholds / Parameters (from your spec):
%   Wavelet scale: s=1
%   Wavelet threshold: Wthr = 1.5e-4
%   Coverage target: ~2.5% (reported)
%   Tracking: overlap of connected regions between consecutive frames
%   tau_min = 0.166*Tinf, with Tinf=5.41 (DNS value you gave)
%   Eccentricity: eps_max=0.85
%   Ecc averaging window: preceding tau_min/3 time steps
%   lambda2 true: lambda2 < -2, area>=2 grid points, persist >= tau_min
%   Match: wavelet detection is TP if overlaps true vortex >= 40% of its lifetime
%
% Dependencies:
%   - Wavelet Toolbox (cwtft2)
%   - Image Processing Toolbox (bwconncomp, labelmatrix)
%==========================================================================

%% -------------------- unpack + constants --------------------------------
eta  = inputs.eta;
dtau = inputs.dtau;
Tinf = inputs.Tinf;

[Ny,Nx,Nt] = size(eta);

% SOC parameters (exact values you listed)
s        = 1;
Wthr     = 1.5e-4;
epsMax   = 0.85;
tauMin   = 0.166 * Tinf;    % with Tinf=5.41 -> tauMin ~ 0.898
lambda2thr = -2;

% Eccentricity averaging over preceding tauMin/3
nAvg = max(1, round((tauMin/3)/dtau));

% True vortex pragmatic limits
trueMinAreaPix = 2;

% Connectivity
conn = 8;

% Sign convention:
% If your eta has indentations negative and you want those to map to positive W,
% you might need a sign flip depending on mexh normalization.
% Keep this ON if you expect indentations -> positive W.
doSignFlip = true;

%% -------------------- Step 1: mean removal (optional) -------------------
% If your eta has a static bias, you can remove a temporal mean:
% eta = eta - mean(eta,3);
% (Commented out by default because DNS eta often already centered)

%% -------------------- Step 2-4: wavelet + threshold + CC ----------------
W   = zeros(Ny,Nx,Nt,'like',double(eta));
BW  = false(Ny,Nx,Nt);
CC  = cell(Nt,1);
coverage = zeros(Nt,1);

for t = 1:Nt
    snapshot = double(eta(:,:,t));

    % Wavelet Toolbox 2D CWT at one scale
    cwtRes = cwtft2(snapshot, 'Wavelet', 'mexh', 'Scales', s);
    cfs = cwtRes.cfs;
    if ndims(cfs) == 3
        cfs = cfs(:,:,1);
    end

    if doSignFlip
        cfs = -cfs;
    end

    W(:,:,t) = cfs;

    BW(:,:,t) = (cfs > Wthr);
    coverage(t) = nnz(BW(:,:,t)) / (Ny*Nx);

    CC{t} = bwconncomp(BW(:,:,t), conn);
end

%% -------------------- Step 5: tracking by overlap -----------------------
% Track wavelet detections using overlap of connected regions between frames
detTracksAll = track_by_overlap(CC, W, Ny, Nx, dtau);

%% -------------------- Step 6: physical filters (tau + ecc) --------------
% Apply:
%  - lifetime >= tauMin
%  - max rolling-mean(ecc over preceding nAvg frames) <= epsMax
detTracksFiltered = filter_tracks_soc(detTracksAll, dtau, tauMin, nAvg, epsMax);

%% -------------------- lambda2 "true vortices" (Step 5 verification) -----
hasLambda2 = false;
if isfield(inputs,'grad') || (isfield(inputs,'U') && isfield(inputs,'V') && isfield(inputs,'W') ...
        && isfield(inputs,'dx') && isfield(inputs,'dy'))
    hasLambda2 = true;
end

trueTracksFiltered = struct([]);
lambda2 = [];
trueMaskFiltered = [];

metrics = struct();
if hasLambda2
    % Compute lambda2 field per frame
    lambda2 = compute_lambda2_surface(inputs);

    % Threshold true cores and segment
    CCtrue = cell(Nt,1);
    BWtrue = false(Ny,Nx,Nt);
    for t = 1:Nt
        BWt = (lambda2(:,:,t) < lambda2thr);
        CCt = bwconncomp(BWt, conn);

        % enforce min area >= 2 grid points at the component stage
        keep = cellfun(@numel, CCt.PixelIdxList) >= trueMinAreaPix;
        CCt.PixelIdxList = CCt.PixelIdxList(keep);
        CCt.NumObjects = numel(CCt.PixelIdxList);

        CCtrue{t} = CCt;

        % rebuild BWtrue from kept components
        BWt2 = false(Ny,Nx);
        for k = 1:CCt.NumObjects
            BWt2(CCt.PixelIdxList{k}) = true;
        end
        BWtrue(:,:,t) = BWt2;
    end

    % Track true vortices by overlap
    trueTracksAll = track_by_overlap(CCtrue, [], Ny, Nx, dtau);

    % Filter true tracks by persistence >= tauMin
    trueTracksFiltered = filter_tracks_true(trueTracksAll, dtau, tauMin);

    % Build a filtered true-mask per frame containing only kept true tracks
    trueMaskFiltered = build_mask_from_tracks(trueTracksFiltered, Ny, Nx, Nt);

    % Evaluate wavelet detections vs true vortices
    metrics = evaluate_wavelet_vs_true(detTracksFiltered, trueMaskFiltered, Ny, Nx);
end

%% -------------------- package results -----------------------------------
results = struct();
results.params = struct('s',s,'Wthr',Wthr,'epsMax',epsMax,'tauMin',tauMin,'nAvg',nAvg, ...
                        'lambda2thr',lambda2thr,'dtau',dtau,'Tinf',Tinf,'doSignFlip',doSignFlip);
results.coverage = coverage;                       % fraction per frame
results.coverageMean = mean(coverage);             % should be ~0.025 in their DNS tuning
results.W = W;                                     % wavelet spectrum (can be big)
results.BW = BW;                                   % threshold mask
results.detTracksAll = detTracksAll;
results.detTracksFiltered = detTracksFiltered;

results.lambda2 = lambda2;
results.trueMaskFiltered = trueMaskFiltered;
results.trueTracksFiltered = trueTracksFiltered;

results.metrics = metrics;

% quick console print
fprintf('Wavelet coverage mean = %.4f (%.2f%%)\n', results.coverageMean, 100*results.coverageMean);
fprintf('Wavelet tracks: all=%d, filtered=%d\n', numel(detTracksAll), numel(detTracksFiltered));
if hasLambda2
    fprintf('True(lambda2) tracks filtered=%d\n', numel(trueTracksFiltered));
    fprintf('TP=%d  FP=%d  FN=%d  Precision=%.3f  Recall=%.3f  F1=%.3f\n', ...
        metrics.TP, metrics.FP, metrics.FN, metrics.precision, metrics.recall, metrics.F1);
else
    fprintf('lambda2 inputs not provided -> verification skipped.\n');
end

end

%% ========================================================================
%                               HELPERS
% ========================================================================

function tracks = track_by_overlap(CC, W, Ny, Nx, dtau)
% Tracks connected components frame-to-frame by spatial overlap.
% If W is provided (Ny x Nx x Nt), compute intensity-weighted eccentricity per component.

Nt = numel(CC);
tracks = struct('id',{},'t',{},'pix',{},'area',{},'ecc',{},'alive',{});
nextID = 1;

for t = 1:Nt
    comps = CC{t};
    objPix = comps.PixelIdxList;
    numObj = comps.NumObjects;

    % per-object features
    area = zeros(numObj,1);
    ecc  = zeros(numObj,1);
    for i = 1:numObj
        p = objPix{i};
        area(i) = numel(p);
        if ~isempty(W)
            ecc(i)  = ecc_intensity_cov(p, W(:,:,t), Ny, Nx);
        else
            ecc(i) = NaN;
        end
    end

    if t == 1
        for i = 1:numObj
            tracks(end+1) = newTrack(nextID, t, objPix{i}, area(i), ecc(i)); %#ok<AGROW>
            nextID = nextID + 1;
        end
        continue
    end

    aliveIdx = find([tracks.alive]);
    usedTrack = false(1, numel(aliveIdx));

    assignTo = zeros(numObj,1); % track index in tracks (0 => new)
    assignScore = zeros(numObj,1);

    % find best overlap match for each current object (greedy one-to-one)
    for i = 1:numObj
        cur = objPix{i};
        bestK = 0;
        bestIoU = 0;

        for jj = 1:numel(aliveIdx)
            if usedTrack(jj), continue; end
            k = aliveIdx(jj);
            lastPix = tracks(k).pix{end};

            inter = numel(intersect(cur, lastPix));
            if inter == 0, continue; end
            uni = numel(union(cur, lastPix));
            iou = inter / uni;

            if iou > bestIoU
                bestIoU = iou;
                bestK = k;
            end
        end

        if bestK ~= 0
            assignTo(i) = bestK;
            assignScore(i) = bestIoU;
        end
    end

    % resolve conflicts: if multiple objects picked same track, keep highest IoU
    if any(assignTo > 0)
        uniqueTracks = unique(assignTo(assignTo>0));
        for ut = uniqueTracks(:)'
            idx = find(assignTo == ut);
            if numel(idx) > 1
                [~,m] = max(assignScore(idx));
                keep = idx(m);
                drop = setdiff(idx, keep);
                assignTo(drop) = 0;
            end
        end
    end

    % update matched tracks
    updated = false(1, numel(tracks));
    for i = 1:numObj
        k = assignTo(i);
        if k > 0
            tracks(k) = extendTrack(tracks(k), t, objPix{i}, area(i), ecc(i));
            updated(k) = true;
            % mark that this alive track was used
            jj = find(aliveIdx == k, 1);
            if ~isempty(jj), usedTrack(jj) = true; end
        end
    end

    % mark unmatched alive tracks as dead
    for jj = 1:numel(aliveIdx)
        k = aliveIdx(jj);
        if ~updated(k)
            tracks(k).alive = false;
        end
    end

    % create new tracks for unmatched objects
    for i = 1:numObj
        if assignTo(i) == 0
            tracks(end+1) = newTrack(nextID, t, objPix{i}, area(i), ecc(i)); %#ok<AGROW>
            nextID = nextID + 1;
        end
    end
end

% finalize
for k = 1:numel(tracks)
    tracks(k).alive = false;
end

% store lifetime in tau units for convenience
for k = 1:numel(tracks)
    n = numel(tracks(k).t);
    tracks(k).tau = max(0, (n-1)*dtau);
end

end

function out = filter_tracks_soc(tracks, dtau, tauMin, nAvg, epsMax)
% SOC filter:
%   - tau >= tauMin
%   - rolling average eccentricity over preceding nAvg frames <= epsMax (max over track)

keep = false(1, numel(tracks));
for k = 1:numel(tracks)
    n = numel(tracks(k).t);
    tau = max(0,(n-1)*dtau);
    if tau < tauMin, continue; end

    epsk = tracks(k).ecc;
    if all(isnan(epsk))
        % if no ecc available, fail closed (or set keep=true if you want)
        continue
    end

    % rolling mean over preceding nAvg frames (inclusive)
    roll = zeros(1,n);
    for j = 1:n
        a = max(1, j-nAvg+1);
        roll(j) = mean(epsk(a:j));
    end

    if max(roll) <= epsMax
        keep(k) = true;
    end
end
out = tracks(keep);
end

function out = filter_tracks_true(tracks, dtau, tauMin)
% True-vortex filter: persistence >= tauMin (area>=2 enforced earlier)
keep = false(1, numel(tracks));
for k = 1:numel(tracks)
    n = numel(tracks(k).t);
    tau = max(0,(n-1)*dtau);
    if tau >= tauMin
        keep(k) = true;
    end
end
out = tracks(keep);
end

function mask = build_mask_from_tracks(tracks, Ny, Nx, Nt)
mask = false(Ny,Nx,Nt);
for k = 1:numel(tracks)
    for j = 1:numel(tracks(k).t)
        t = tracks(k).t(j);
        mask(tracks(k).pix{j} + (t-1)*Ny*Nx) = true; % linear into 3D
    end
end
end

function metrics = evaluate_wavelet_vs_true(detTracks, trueMask, Ny, Nx)
% A wavelet detection track is TP if it overlaps trueMask in >=40% of its lifetime frames.

TP = 0; FP = 0;
for k = 1:numel(detTracks)
    n = numel(detTracks(k).t);
    hits = 0;
    for j = 1:n
        t = detTracks(k).t(j);
        pix = detTracks(k).pix{j};
        if any(trueMask(pix + (t-1)*Ny*Nx))
            hits = hits + 1;
        end
    end
    frac = hits / n;
    if frac >= 0.40
        TP = TP + 1;
    else
        FP = FP + 1;
    end
end

% For a simple FN estimate: count true "events" (connected components in trueMask)
% that never overlap any wavelet track at the same time at least once.
% (If you want FN per true track with its own 40% rule, say so and I'll adjust.)
Nt = size(trueMask,3);
everHit = false(Nt,1);
for t = 1:Nt
    anyDetOverlap = false;
    for k = 1:numel(detTracks)
        idx = find(detTracks(k).t == t, 1);
        if isempty(idx), continue; end
        if any(trueMask(detTracks(k).pix{idx} + (t-1)*Ny*Nx))
            anyDetOverlap = true;
            break
        end
    end
    everHit(t) = anyDetOverlap;
end
% FN estimate: frames that have true vortices but no overlaps (frame-based FN)
trueFrames = squeeze(any(reshape(trueMask,Ny*Nx,Nt),1));
FN = nnz(trueFrames & ~everHit);

precision = TP / max(1,(TP+FP));
recall    = TP / max(1,(TP+FN));   % not perfect, but gives a signal
F1        = 2*precision*recall / max(1e-12,(precision+recall));

metrics = struct('TP',TP,'FP',FP,'FN',FN,'precision',precision,'recall',recall,'F1',F1);
end

function e = ecc_intensity_cov(pix, Wframe, Ny, Nx)
% Intensity covariance matrix sigma_ij for a connected region Ac.
% Uses weights = Wframe(pix) (positive wavelet "intensity" within region).

w = double(Wframe(pix));
w(w < 0) = 0;              % just in case
sw = sum(w);
if sw <= 0
    e = 1; % degenerate -> treat as very eccentric
    return
end

[yy, xx] = ind2sub([Ny, Nx], pix);
xx = double(xx); yy = double(yy);

mx = sum(w.*xx)/sw;
my = sum(w.*yy)/sw;

dx = xx - mx; dy = yy - my;

Cxx = sum(w.*dx.^2)/sw;
Cyy = sum(w.*dy.^2)/sw;
Cxy = sum(w.*dx.*dy)/sw;

C = [Cxx Cxy; Cxy Cyy];

lam = sort(eig(C), 'descend');
L1 = lam(1); L2 = lam(2);
if L1 <= 0
    e = 1;
else
    r = max(0, min(1, L2/L1));
    e = sqrt(1 - r);
end
end

function tr = newTrack(id, t, pix, area, ecc)
tr = struct();
tr.id = id;
tr.t = t;
tr.pix = {pix};
tr.area = area;
tr.ecc = ecc;
tr.alive = true;
tr.tau = 0;
end

function tr = extendTrack(tr, t, pix, area, ecc)
tr.t(end+1) = t;
tr.pix{end+1} = pix;
tr.area(end+1) = area;
tr.ecc(end+1) = ecc;
tr.alive = true;
end

function lambda2 = compute_lambda2_surface(inputs)
% Computes Jeong & Hussain lambda2 from surface velocity gradients.
% Requires either:
%   inputs.grad: Ny x Nx x 3 x 3 x Nt
% or velocities + dx,dy (+ optional z-derivatives):
%   inputs.U,V,W, dx,dy, optional dUdz,dVdz,dWdz

if isfield(inputs,'grad')
    G = inputs.grad; % Ny x Nx x 3 x 3 x Nt
    [Ny,Nx,~,~,Nt] = size(G);
else
    U = inputs.U; V = inputs.V; W = inputs.W;
    dx = inputs.dx; dy = inputs.dy;
    [Ny,Nx,Nt] = size(U);

    % x,y derivatives from surface fields
    dudx = zeros(Ny,Nx,Nt); dudy = zeros(Ny,Nx,Nt);
    dvdx = zeros(Ny,Nx,Nt); dvdy = zeros(Ny,Nx,Nt);
    dwdx = zeros(Ny,Nx,Nt); dwdy = zeros(Ny,Nx,Nt);

    for t = 1:Nt
        [dUdy, dUdx] = gradient(double(U(:,:,t)), dy, dx);
        [dVdy, dVdx] = gradient(double(V(:,:,t)), dy, dx);
        [dWdy, dWdx] = gradient(double(W(:,:,t)), dy, dx);
        dudx(:,:,t) = dUdx; dudy(:,:,t) = dUdy;
        dvdx(:,:,t) = dVdx; dvdy(:,:,t) = dVdy;
        dwdx(:,:,t) = dWdx; dwdy(:,:,t) = dWdy;
    end

    % z-derivatives: MUST come from near-surface 3D data for exact lambda2
    if isfield(inputs,'dUdz'), dudz = inputs.dUdz; else, dudz = zeros(Ny,Nx,Nt); end
    if isfield(inputs,'dVdz'), dvdz = inputs.dVdz; else, dvdz = zeros(Ny,Nx,Nt); end
    if isfield(inputs,'dWdz'), dwdz = inputs.dWdz; else, dwdz = zeros(Ny,Nx,Nt); end

    % pack into G
    G = zeros(Ny,Nx,3,3,Nt);
    % row 1: grad of u
    G(:,:,1,1,:) = dudx;  G(:,:,1,2,:) = dudy;  G(:,:,1,3,:) = dudz;
    % row 2: grad of v
    G(:,:,2,1,:) = dvdx;  G(:,:,2,2,:) = dvdy;  G(:,:,2,3,:) = dvdz;
    % row 3: grad of w
    G(:,:,3,1,:) = dwdx;  G(:,:,3,2,:) = dwdy;  G(:,:,3,3,:) = dwdz;
end

% lambda2 per point: median eigenvalue of (S^2 + Omega^2)
lambda2 = zeros(Ny,Nx,Nt);

for t = 1:Nt
    for iy = 1:Ny
        for ix = 1:Nx
            J = squeeze(G(iy,ix,:,:,t));     % 3x3
            S = 0.5*(J + J.');
            O = 0.5*(J - J.');
            M = S*S + O*O;
            ev = sort(eig(M), 'ascend');      % ev(2) is median
            lambda2(iy,ix,t) = ev(2);
        end
    end
end

end