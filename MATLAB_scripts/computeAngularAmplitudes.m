function amplitudes = computeAngularAmplitudes(angles) 
    % Compute the angular amplitudes of material points. We define this to be
    % half the range of theta for each arclength over the data.

    amplitudes = (max(angles,[],1) - min(angles,[],1)) / 2;

end