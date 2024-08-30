function [overallScore, xFitOptimal, yFitOptimal, optimalRotation] = computeQualityOfFitXY(x,y,xFit,yFit,arclengths,timestamps,maxIterations)
%% Compute the quality of fit of xFit,yFit to x,y, allowing for
% rotation about (0,0).
    if nargin < 7
        maxIterations= 10;
    end
    assert(all(size(x) == size(xFit)),'Fitted data must be the same shape as the original.')
    xFitOptimal = zeros(size(xFit));
    yFitOptimal = zeros(size(yFit));
    optimalRotation = zeros(size(xFit,1),1);
    scores = zeros(size(xFit,1),1) + Inf;

    opts = optimoptions('fminunc','MaxIterations',maxIterations,'Display','none');
    for frame = 1 : size(xFit,1)
        f = @(theta) computeScore(x(frame,:),y(frame,:),xFit(frame,:),yFit(frame,:),arclengths,theta);
        [optimalRotation(frame), scores(frame)] = fminunc(f, 0, opts);
        [~, xFitOptimal(frame,:), yFitOptimal(frame,:)] = f(optimalRotation(frame));
    end
    overallScore = sqrt(trapz(timestamps,trapz(arclengths,(x - xFitOptimal).^2 + (y - yFitOptimal).^2,2),1) / trapz(timestamps,trapz(arclengths,x.^2 + y.^2,2),1));
end

function [score, xTrial, yTrial] = computeScore(x,y,xFit,yFit,arclengths,theta) 
    [xTrial,yTrial] = rotateCoords(xFit,yFit,theta);
    score = trapz(arclengths,((xTrial - x).^2 + (yTrial - y).^2).^0.5) / sum((x.^2 + y.^2).^0.5);
end