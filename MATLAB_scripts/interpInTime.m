function [x, y, angles] = interpInTime(x, y, arclengths, badFrameMask)
    % Throw away data from bad frames.
    x(badFrameMask,:) = NaN;
    y(badFrameMask,:) = NaN;
    x = fillmissing(x,'pchip',1);
    y = fillmissing(y,'pchip',1);
    angles = atan2(gradient(y,arclengths,1), gradient(x,arclengths,1));
    angles = unwrap(angles,[],2);
    
end