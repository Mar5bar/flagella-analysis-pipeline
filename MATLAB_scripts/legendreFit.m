function [total, components, coeffs] = legendreFit(data)
%% Compute and return the Legendre polynomial fit to data, assuming it is on
% [0,1] and equally spaced.
    
    coeffs = legendreCoeffs(data);
    components = coeffs .* legendrePolys(linspace(0,1,numel(data)));
    total = sum(components, 2);
end