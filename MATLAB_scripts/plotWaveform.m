function h = plotWaveform(spatialX, spatialY, numFrames)
    if nargin < 3
        numFrames = size(spatialX,1);
    elseif strcmp(numFrames,'all')
        numFrames = size(spatialX,1);
    end
    h = plot(spatialX(1:numFrames,:)', spatialY(1:numFrames,:)');
    axis equal
end