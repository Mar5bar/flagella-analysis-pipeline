function [spatialX, spatialY] = anglesToSpatial(angles,arclengths)
    if nargin < 2
        arclengths = linspace(0,1,size(angles,2));
    end
    spatialX = cumtrapz(arclengths, cos(angles),2);
    spatialY = cumtrapz(arclengths, sin(angles),2);
end
