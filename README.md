# Flagella analysis pipeline

MATLAB and Fiji (ImageJ) scripts associated with 'Axonemal dynein contributions to flagellar beat types and waveforms' by S. Fochler, M. H. Doran, T. Beneke, J. Smith, C. Fort, B. J. Walker, A. Brown, E. Gluenz, R. J. Wheeler.

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

Fiji scripts require [Fiji](https://imagej.net/software/fiji/). Tested on Fijji 2.16.0 on macOS 15.6.1.