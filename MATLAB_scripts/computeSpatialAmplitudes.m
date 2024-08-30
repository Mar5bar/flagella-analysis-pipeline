function amplitudes = computeSpatialAmplitudes(x,y) 
    % Compute the amplitudes of material points with coordinates (x,y) over
    % the given data. We define this to be half the largest distance between
    % (x,y) coords for each arclength over the data.

    % Go through each material point, one by one.
    amplitudes = zeros(1,size(x,2));
    for i = 1 : size(x,2)
        X = x(:,i); Y = y(:,i);
        d = ((X - X').^2 + (Y - Y').^2).^(0.5);
        amplitudes(i) = max(d(:))/2;
    end

end