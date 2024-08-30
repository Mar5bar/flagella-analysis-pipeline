function [x,y,widths,numFrames] = loadRawCoordinates(filename,spatialScale)

    % Load each line in a cell array of strings.
    data = readlines(filename);

    % There is an empty final row.
    numFrames = length(data)-1;

    % Initialise x,y,widths as cell arrays.
    x = cell(numFrames,1);
    y = cell(numFrames,1);
    widths = cell(numFrames,1);

    for frame = 1 : numFrames
        nums = str2num(data{frame});
        x{frame} = nums(1:3:end)*spatialScale;
        y{frame} = nums(2:3:end)*spatialScale;
        widths{frame} = nums(3:3:end);
    end
end
