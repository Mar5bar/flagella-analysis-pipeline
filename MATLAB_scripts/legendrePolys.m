function l = legendrePolys(s)
%% Evaluate shifted Legendre polynomials on s.

    s = s(:);

    l = [1 + 0*s,...
         2*s - 1,...
         6*s.^2 - 6*s + 1];
end