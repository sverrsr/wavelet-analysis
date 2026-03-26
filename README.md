# wavelet-analysis

Wavelet-based MATLAB tooling for detecting periodic free-surface vortex signatures (dimples and scars) from a `256 x 256 x 12500` surface-elevation stack.

## Files

- `Main_Vortex_Analysis.m` orchestrates loading data, setting physical defaults, executing detection/tracking/validation, and saving `wavelet_vortex_analysis_results.mat`.
- `detectFeatures.m` denoises each frame with Wavelet Toolbox (`wavedec2`, `wthresh`, `waverec2`) and applies a periodic 2-D Mexican Hat convolution to classify dimples vs scars.
- `trackFeatures.m` links connected components through time using overlap-based persistence filtering.
- `validateVortices.m` applies `\lambda_2`-based validation, estimates scar geodesic centers with `bwskel`, and computes area-based NCC metrics against `\beta^2`.
- `detectDimplesOnly.m` is the new dimple-only entry point: it accepts a `surfElev` stack, skips `\lambda_2` checks, and saves the persistent dimple binary map, total dimple area per frame, and dimple count per frame.

## Expected input

Place `vortex_input.mat` in the repository root, `input/`, or `data/` with:

- required: `etaStack` shaped `[Ny x Nx x Nt]` (production target: `256 x 256 x 12500`), periodic over both spatial directions;
- optional: `beta2Stack`, `beta2ByDepth`, `zLevels`, `lambda2Stack`, or `velocityGradientStack`.

If no input file is present, the main script generates a compact synthetic periodic dataset so the workflow remains runnable.

## MATLAB requirements

- MATLAB R2022b or later.
- Wavelet Toolbox.
- Image Processing Toolbox.

## Usage

```matlab
Main_Vortex_Analysis
```

The script saves a `wavelet_vortex_analysis_results.mat` file containing the configuration, frame-level detections, track summaries, and NCC metrics.

## Dimple-only workflow

Use `detectDimplesOnly` when you only want dimples and do not want the `\lambda_2` validation stage. The function accepts `surfElev` directly, so the common call is:

```matlab
results = detectDimplesOnly(surfElev, 'OutputFile', 'dimple_detection_results.mat');
```

Saved variables in the MAT-file:

- `dimpleBinaryMap`: logical `Ny x Nx x Nt` map of persistent dimples.
- `dimpleAreaByFrame`: `1 x Nt` total dimple area in pixels for each frame.
- `dimpleCountByFrame`: `1 x Nt` number of persistent dimples in each frame.
- `componentTable`: per-dimple metadata (frame, track id, area, eccentricity, centroid, scale).
