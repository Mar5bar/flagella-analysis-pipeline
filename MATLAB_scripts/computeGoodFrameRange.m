function [frameRange, badFlag] = computeGoodFrameRange(badFrameMask)
% Find a range of good frames.
    
    goodFrameMask = ~badFrameMask;
    % We'll join up any small gaps in goodFrameMask. If a frame is within
    % windowWidth frames of a good frame on BOTH sides, we allow it.
    windowWidth = 2;
    allowable = false(size(badFrameMask));
    for frame = 1 : length(allowable)
        allowable(frame) = any(goodFrameMask(max(1,frame-windowWidth):frame)) &...
                            any(goodFrameMask(frame:min(length(goodFrameMask),frame+windowWidth)));
    end
    addedFrames = xor(goodFrameMask, allowable);
    goodFrameMask = goodFrameMask | allowable;

    % goodFrameMask = goodFrameMask | logical(filloutliers(double(goodFrameMask),'nearest','movmedian',round(length(goodFrameMask)/10)));

    % Find the largest connected block in goodFrameMask.
    zeroPositions = find(~[0;goodFrameMask;0]); % Locations of zeros, with buffer marking start and end.
    [~, blockLoc] = max(diff(zeroPositions)); % Location of largest gap between zeros ie the block we want.
    goodFrameMask(:) = false;
    goodFrameMask(zeroPositions(blockLoc):zeroPositions(blockLoc+1)-2) = true;
    frameRange = find(goodFrameMask,1,'first') : find(goodFrameMask,1,'last');
       
    % If more than 20% of frames in the range are to be interpolated, this
    % is a bad file.
    badFlag = sum(addedFrames(frameRange)) > 0.2 * length(frameRange);

end