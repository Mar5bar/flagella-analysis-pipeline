function [flagMask, bodyMask, flagOnRight] = computeFlagMask(widths)
% Compute a flagMask from the width profile data that identifies the flagellum.
    
    thresh = min(mode(widths(1:round(end/2))),mode(widths(round(end/2):end)));
    flagMask = widths <= thresh+1;
    bodyMask = widths > thresh+1;

    % We'll join up any small gaps in flagMask. We use a moving median filter
    % with a window length of 20% of signal. OR prevents removal of good points.
    flagMask = flagMask | logical(filloutliers(double(flagMask),'nearest','movmedian', max(round(length(flagMask)/5),1)));

    % The flagellum should now be the largest connected block in the flagMask, which we
    % find.
    zeroPositions = find(~[0,flagMask,0]); % Locations of zeros, with buffer marking start and end.
    [~, blockLoc] = max(diff(zeroPositions)); % Location of largest gap between zeros ie the block we want.
    flagMask(:) = false;
    flagMask(zeroPositions(blockLoc):zeroPositions(blockLoc+1)-2) = true;

    % Do the same thing to find the body without outliers.
    zeroPositions = find(~[0,bodyMask,0]); % Locations of zeros, with buffer marking start and end.
    [~, blockLoc] = max(diff(zeroPositions)); % Location of largest gap between zeros ie the block we want.
    bodyMask(:) = false;
    bodyMask(zeroPositions(blockLoc):zeroPositions(blockLoc+1)-2) = true;

    % Having identified the flagellum in flagMask, ~flagMask will contain the body
    % (+outliers). We'll check if the body is on the left or the right of the
    % flagellum, returning flagOnRight = true if the flagellum is on the right
    % of the body.
    flagLoc = median(find(flagMask));
    bodyLoc = median(find(bodyMask));
    flagOnRight = bodyLoc < flagLoc;

    % Now, make sure that the masks extendsall the way to the relevant
    % ends.
    if flagOnRight
        flagMask(find(flagMask,1,'first') : end) = true;
        bodyMask(1 : find(bodyMask,1,'first')) = true;
    else
        flagMask(1 : find(flagMask,1,'first')) = true;
        bodyMask(find(bodyMask,1,'last') : end) = true;
    end
    
    % We appear to be able to robustly identify the body, whilst the
    % flagellum is sometimes truncated due to wider sections in the
    % flagellum. In an attempt to overcome this, we'll simply invert the
    % body mask to find the flagellum.
    flagMask = ~bodyMask;

    if ~any(bodyMask) || ~any(flagMask)
        warning('Unable to distinguish between the body and flagellum. Body is likely very large.')
    end

end