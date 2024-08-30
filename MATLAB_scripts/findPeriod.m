function [period, startFrame, lag] = findPeriod(angles, timepoints)
    % We will take in all of the angle data and compute the temporal
    % autocorrelation at each arclength and sum them to give a single
    % autocorrelation curve. Timepoints must be uniformly spaced.
    if nargin < 2
        timepoints = linspace(0,1,size(angles,1));
    end

    % Sum up the autocorrelations of all the material points.
    autoCorr = zeros(2*size(angles,1)-1,1);
    for i = 1 : size(angles,2)
        [a, l] = xcorr(angles(:,i));
        autoCorr = autoCorr + a;
    end
    mask = l > 0;
    autoCorr = autoCorr(mask); l = l(mask);
    [~, lag] = findpeaks(autoCorr,l,'NPeaks',1,'SortStr','descend');
    if isempty(lag)
        lag = NaN;
        period = NaN;
        startFrame = NaN;
        return
    end
    
    % We compute a precise period by fitting a spline to the autocorrelation.
    ls = linspace(min(l),max(l),10*numel(l));
    f = fit(l(:),double(autoCorr),'spline');
    [~,highresLag] = findpeaks(f(ls),ls,'NPeaks',1,'SortStr','descend');
    period = interp1(l,timepoints(2:end),highresLag,'linear');

    % We find the starting point so as to minimise the jump in the tangent
    % angle over the period.
    d = angles(1:end-lag,:) - angles(lag+1:end,:);
    [~, startFrame] = min(max(abs(d),[],2));
end