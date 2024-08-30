function batchProcess(sourceDir)

%% Batch script for processing a directory of output.
% Frame rate, in Hz.
frameRate = 200;
% Spatial scale, in um per pixel.
spatialScale = 1;
% The threshold below which we reject frames, based on what proportion of the
% overall flagellum length is traced from the frame. 0.9 will reject frames
% where less than 90% of the flagellum has been traced. NOTE: we'll truncate
% all of the captured flagella at this proportion!
thresholdForFrameRejection = 0.9;
% The number of interpolated sample points along the analysed flagellum.
numArclengthSamples = 200;


oldpath = addpath(['.',filesep,'MATLAB_scripts']);

if nargin < 1
    % Request the input directory to scour for folders ending in .tif_analysis.
    sourceDir = uigetdir('.','Select the directory containing the .tif_analysis directories');
end

if sourceDir == 0
    disp('No directory provided. Aborting.')
    return
else
    disp(['Source directory: ',sourceDir])
end

% This will get all the directories ending in *.tif_analysis in dataForMATLAB, recursively.
[parentFolders, folders] = traverseDirs(sourceDir, '*.tif_analysis*');

% Set the filename of the coordinates files in every subfolder. MUST BE AN
% EXACT MATCH.
filename = 'rawcoordinates.txt';
tic
parfor i = 1 : length(folders)
    disp([num2str(i), ' / ',num2str(length(folders)),'. ',folders{i},'.'])
    try
        processFile([parentFolders{i}, filesep, folders{i}], filename, frameRate, spatialScale, thresholdForFrameRejection, numArclengthSamples)
    catch exception
        disp('Error. Analysis aborted.')
        disp(exception)
    end
end
toc

end