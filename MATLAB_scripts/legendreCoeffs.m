function g = fitLegendre(data)
%% Fit shifted Legendre polynomials to the data (1xn or nx1), assumed to be
% evaluated at equally space points on [0,1] (inc. endpoints). The returned
% coefficients will be in ascending order, starting with the mean.

    data = data(:);

    % Define the scalings needed to evaluate the coefficient.
    scales = [1,3,5];

    xs = linspace(0,1,numel(data))';

    g = scales .* trapz(xs, data .* legendrePolys(xs));

end