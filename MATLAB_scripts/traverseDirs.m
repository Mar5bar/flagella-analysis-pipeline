function [folders, filenames] = traverseDirs(sourceDir, matchString)
    dataDirs = dir([sourceDir,filesep,matchString]);
    filenames = [{dataDirs.name}];
    folders = [{dataDirs.folder}];

    subDirs = dir(sourceDir);
    subDirs = subDirs([subDirs.isdir]);
    for i = 1 : length(subDirs)
        if strcmp(subDirs(i).name, '.') || strcmp(subDirs(i).name, '..')
            continue
        else
            [newFolders, newFilenames] = traverseDirs([subDirs(i).folder, filesep, subDirs(i).name], matchString);
            folders = [folders, newFolders];
            filenames = [filenames, newFilenames];
        end
    end
end