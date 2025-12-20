# Flagella analysis pipeline

MATLAB and Fiji (ImageJ) scripts associated with 'Axonemal dynein contributions to flagellar beat types and waveforms' by S. Fochler, M. H. Doran, T. Beneke, J. Smith, C. Fort, B. J. Walker, A. Brown, E. Gluenz, R. J. Wheeler.

## Installation

Download the repository as a .zip and extract. This should take less than 1 second on most machines.

## Usage


### Waveform tracing from image data

To trace flagellar waveforms from `.tif` files containing single cells, run `Flagellum_Phase_TraceRawCoords.ijm` from the `ImageJ_scripts` directory. When prompted, select the folder containing `.tif` imaging data.

The script will generate a subdirectory for each processed file, outputting pixel coordinates of the flagellum for each frame in `rawcoordinates.txt` and saving the identified flagellar length in each frame in `traceProfile.png`.

Custom preprocessing parameters can be input on a per-directory basis by providing `threshold.txt` in the directory. An example is provided in the `example_data` directory. It is recommended to experiment with parameters before bulk analysing a new dataset.

### Waveform analysis

To perform automated waveform analysis, run `batchProcess.m` in MATLAB. You will be prompted to select a directory containing traced images (i.e. the same directory that was processed in the previous step). Images will be processed in parallel (first-time runs may take longer than expected due to spinning up of the parallel pool).

A range of plots (`.png`) and summary statistics (`summary.txt` and `output.mat`) will be created in the identified subdirectories. Statistics can be aggregated across files by running `aggregateData.m` in MATLAB and selecting a processed directory.

Descriptions of all computed quantities are saved in `output.mat` and viewable in `MATLAB_scripts/assignDescriptionsAndUnitsToTable.m`.

## Demo

1. Run `ImageJ_scripts/Flagellum_Phase_TraceRawCoords.ijm` on the `example_data` directory in FIJI (ImageJ). This should take a few seconds on most machines.

1. Run `batchProcess.m` using the (now-processed) `example_data` directory. This should take a few seconds on most machines. First-time runs may take longer than expected due to spinning up of the (here unused) parallel pool.

1. This will generate a large number of plot `.png` and a summary file in `example.tif_analysis`. To verify correctness, compare the generated `summary.txt` with the top-level `example_summary.txt`. All entries should be exact matches, excluding the filepath at the top of the file.

## System requirements

MATLAB scripts require MATLAB (R2022a or above) with Curve Fitting Toolbox. Tested on MATLAB R2025b on macOS 15.6.1.

Fiji scripts require [Fiji](https://imagej.net/software/fiji/) (v2.0.0 or above). Tested on Fijji 2.16.0 on macOS 15.6.1.

## Pseudocode explanation

### `Flagellum_Phase_TraceRawCoords.ijm`:

1. **Initialisation**: Set default filtering and thresholding parameters, or overwrite them by reading `threshold.txt` if available.

2. **File Iteration**: Loop through all `.tif` files in the directory.


3. **Setup**:
* Create specific output directories; skip analysis if `rawcoordinates.txt` already exists.
* Load metadata to define the flagellum base coordinates and video crop regions.


4. **Frame Processing (Loop)**:
* **Pre-process**: Duplicate frame, apply Gaussian blur, Unsharp Mask, and background subtraction.
* **Threshold**: Apply auto-thresholding and refined binary masking.
* **Clean**: Remove bright artefacts and keep only the largest binary particle (the cell body/flagellum).


5. **Skeletonisation & Tracing**:
* Generate a skeleton of the binary shape and prune short branches.
* Identify the skeleton terminus closest to the defined base point.
* Trace the skeleton pixels pixel-by-pixel until a branch or end is reached.

6. **Measurement**: Calculate the width of the flagellum at every traced point using a Euclidean distance map.


7. **Output**: Save the traced coordinates and widths to `rawcoordinates.txt` and save a summary plot of traced points per frame.

### `processFile.m`:

1. **Initialise**: Set default parameters (filepath, sampling frequency, spatial scale, thresholds).
2. **Load Data**: Import raw coordinates and width data from the text file.
3. **Pre-process (Per Frame)**:
* Separate flagellum from body using width thresholds.
* Translate flagellum base to origin .
* Rotate coordinates to align the cell orientation vector horizontally.
* Compute arclengths and smooth spatial traces.


4. **Filter & Resample**:
* Identify "bad frames" (failed traces/short lengths).
* Interpolate spatially to valid uniform arclength points.
* Interpolate temporally to fill gaps in the "best range" of frames.


5. **Geometric Analysis**: Compute tangent angles and curvature fields across time and space.
6. **Fourier Analysis**:
* Perform FFT on tangent angles.
* Extract dominant frequency, amplitude profiles, and phase linearity.


7. **Reconstruction**: Synthesise a theoretical beat from the dominant frequency and compare with actual data to assess quality of fit.
8. **Beat Metrics**: Calculate beat period (via autocorrelation), max/min curvatures, shear velocities, and wave amplitudes.
9. **Output**: Save computed statistics to `summary.txt`, generate diagnostic plots (kymographs, spectra, waveforms), and save workspace to `output.mat`.