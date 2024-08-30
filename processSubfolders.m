% Request the input directory to scour for folders, each of which will contain
% files ending in .tif_analysis.
sourceDir = uigetdir('.','Select the directory containing folders that contain the .tif_analysis directories');
if sourceDir == 0
    disp('No directory provided. Aborting.')
    return
else
    disp(['Source directory: ',sourceDir])
end

% This will get all the subdirectories.
directories = dir(sourceDir);
% remove all files (isdir property is 0)
directories = directories([directories(:).isdir]);
% remove '.' and '..' 
directories = directories(~ismember({directories(:).name},{'.','..'}));

aggregate = table();
for i = 1 : length(directories)
    dirName = [directories(i).folder,filesep,directories(i).name];
    disp(dirName)
    try
        batchProcess(dirName);
        aggregate = [aggregate; aggregateData(dirName)];
    catch
        disp(['Error on ',directories(i).name,'.'])
    end
end

save([sourceDir,filesep,'aggregate.mat'],'aggregate')