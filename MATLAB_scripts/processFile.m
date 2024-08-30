function processFile(filepath, filename, samplingFrequency, spatialScale, arclengthThreshold, numArclengthSamples)
    warning('off','imageio:tiffmexutils:libtiffWarning')

    % Defaults.
    if nargin < 1
        filepath = [pwd,filesep,'dataForMATLAB',filesep,'testFolder1.tif_analysis'];
    end
    if nargin < 2
        filename = 'rawcoordinates.txt';
    end
    if nargin < 3
        samplingFrequency = 100;
    end
    if nargin < 4
        spatialScale = 1;
    end
    if nargin < 5
        arclengthThreshold = 0.9;
    end
    if nargin < 6
        numArclengthSamples = 200;
    end

    message = "";

    %% Load in the coordinate and width data from the file.
    [x,y,widths,numFrames] = loadRawCoordinates([filepath,filesep,filename], spatialScale);

    timestamps = 0:1/samplingFrequency:(numFrames-1)/samplingFrequency;

    % For each frame, find the body and record where it is. Then, remove the
    % body and flip the flagellum data if necessary. Translate and rotate the
    % remaining data.
    overallMaxArclengths = zeros(numFrames,1);
    arclengths = cell(numFrames,1);
    maxArclengths = zeros(numFrames,1);
    for frame = 1 : numFrames
        % Compute the overall arclength of the skeleton, trying to identify
        % frames that just contain the body.
        overallMaxArclengths(frame) = sum((sum(diff([x{frame}(:),y{frame}(:)],1,1).^2,2)).^0.5);
        [flagMask, bodyMask, flagOnRight] = computeFlagMask(widths{frame});
        % Compute the body COM.
        bodyCOM = [mean(x{frame}(bodyMask)), mean(y{frame}(bodyMask))];
        % Retain only the flagellar x,y data.
        x{frame} = x{frame}(flagMask);
        y{frame} = y{frame}(flagMask);
        % Ensure that the flagellum base is the first entry (rather than
        % last).
        if ~flagOnRight
            x{frame} = fliplr(x{frame});
            y{frame} = fliplr(y{frame});
        end
        % Compute the original arclengths.
        arclengths{frame} = [0,cumsum((sum(diff([x{frame}(:),y{frame}(:)],1,1).^2,2)).^0.5)'];
        % Smooth the data in space with a smoothing spline.
        if length(x{frame}) > 1
            [~, x{frame}] = spaps(arclengths{frame},x{frame},0.2*max(arclengths{frame}));
            [~, y{frame}] = spaps(arclengths{frame},y{frame},0.2*max(arclengths{frame}));
        end
        flagBase = [x{frame}(1), y{frame}(1)];
        % Translate the flagellum base to the origin.
        x{frame} = x{frame} - flagBase(1);
        y{frame} = y{frame} - flagBase(2);
        % Compute the cell orientation.
        cellOrientationVector = flagBase - bodyCOM;
        cellOrientationAngle = atan2(cellOrientationVector(2),cellOrientationVector(1));
        % Rotate the coordinates by minus this angle to align with horizontal.
        [x{frame},y{frame}] = rotateCoords(x{frame},y{frame},-cellOrientationAngle);
        % Compute a vector of arclengths for this frame.
        arclengths{frame} = [0,cumsum((sum(diff([x{frame}(:),y{frame}(:)],1,1).^2,2)).^0.5)'];
        maxArclengths(frame) = arclengths{frame}(end);
    end

    % Define a bad frame to be a frame where the flagellum length is less
    % than some proportion of the maximum observed flagellum length. Here,
    % we'll opt for 90% and truncate to that amount.
    overallLength = max(movmean(overallMaxArclengths,round(numFrames/10),'omitnan'));
    flagLength = max(movmean(maxArclengths,round(numFrames/10),'omitnan'));
    badFrameMask = maxArclengths < arclengthThreshold * flagLength | overallMaxArclengths < arclengthThreshold * overallLength;
    numBadFramesOrig = sum(badFrameMask);
    
    % Create a plot that shows how the tracing has done wrt the flagellum,
    % including the estimated flagellum length and the truncation length.
    close all
    f1 = figure('visible','off');

    clf
    plot(1:numFrames,maxArclengths)
    xlabel('Frame','Interpreter','latex')
    ylabel('Flag. length (um)','Interpreter','latex')
    ylim([0,1.1*max(ylim)])
    yline(flagLength,'LineStyle','--')
    yline(arclengthThreshold * flagLength)
    grid on
    set(gca,'FontSize',24)
    text(0.1*max(xlim),0.93*max(ylim), 'Estimated length','FontSize',24)
    text(0.1*max(xlim),0.9*arclengthThreshold * flagLength, 'Truncated length','FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'flagLengthOverTime.png'])

    % Generate the same plot, but showing the length of the whole trace (inc. body).
    clf
    plot(1:numFrames,overallMaxArclengths)
    xlabel('Frame','Interpreter','latex')
    ylabel('Flag. length (um)','Interpreter','latex')
    ylim([0,1.1*max(ylim)])
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'wholeTraceLengthOverTime.png'])

    % Generate uniform arclengths up to arclengthThreshold * flagLength.
    uniformArclengths = linspace(0,arclengthThreshold * flagLength, numArclengthSamples);

    % Resample each frame at the new uniform arclengths, inserting NaNs where
    % no data is available for that frame.
    for frame = 1 : numFrames
        if sum(~isnan(x{frame})) > 1 && sum(~isnan(y{frame})) > 1
            x{frame} = interp1(arclengths{frame},x{frame},uniformArclengths,'linear',NaN);
            y{frame} = interp1(arclengths{frame},y{frame},uniformArclengths,'linear',NaN);
        else
            x{frame} = NaN * uniformArclengths;
            y{frame} = NaN * uniformArclengths;
        end
    end
    % Now that each frame has the same number of samples, assemble in to
    % matrices. Each frame is now a row.
    x = cell2mat(x);
    y = cell2mat(y);

    % Overwrite original stored arclengths with new uniform arclengths.
    arclengths = uniformArclengths;
    
    % Compute the angles again.
    angles = atan2(gradient(y,arclengths,1), gradient(x,arclengths,1));
    % Remove any jumps by 2pi.
    angles = unwrap(angles,[],2);

    % There are potentially lots of NaN values in the data from failed traces,
    % which will ruin a Fourier transform. We'll try to detect the largest
    % block of frames with less than, say, 20% of frames missing, and
    % interpolate through them. We'll interpolate via the angle
    % parameterisation. If such a block doesn't exist, we'll just take all of
    % the frames and raise a flag.

    % Determine best frame range.
    [bestFrameRange, badFileFlag] = computeGoodFrameRange(badFrameMask);
    if badFileFlag
        message = [message, "More than 20% of frames used interpolated data. "];
    end
    numFramesOrig = numFrames;
    numFrames = length(bestFrameRange);
    % If we have fewer than 30 frames, deem this a bad file.
    if numFrames < 30
        message = [message, 'Fewer than 30 good frames. '];
        badFileFlag = true;
    end
    
    xOrig = x;
    yOrig = y;

    %% Only keep data inside the bestFrameRange.
    x = x(bestFrameRange,:);
    y = y(bestFrameRange,:);
    angles = angles(bestFrameRange,:);
    timestamps = timestamps(bestFrameRange);
    badFrameMask = badFrameMask(bestFrameRange);

    xBeforeInterp = x;
    yBeforeInterp = y;
    anglesBeforeInterp = angles;

    [~, ~, angles] = interpInTime(x, y, arclengths, badFrameMask);

    % Now, smooth angles and reconstruct x and y.
    for frame = 1 : numFrames
        [~, angles(frame,:)] = spaps(arclengths, angles(frame,:), 0.01*(arclengths(2)-arclengths(1)));
    end
    [x, y] = anglesToSpatial(angles, arclengths);

    % Compute the curvatures.
    [curvatures, ~] = gradient(angles,arclengths,1);

    % Maximum curvature

    % d(curvature)/dt.
    [~,curvatureVelocity] = gradient(curvatures,1,timestamps);

    % Compute the Fourier spectra.
    [complexSpectrum, powerSpectrum, frequencyDomain] = fourierSpectra(angles, samplingFrequency);
    frequencyBinWidth = frequencyDomain(2) - frequencyDomain(1); % The size of the frequency bins.

    % Create a windowed spectrum, averaging the amplitudes in a centred window of width 2*windowWidth
    windowHalfWidthHz = 1; % The approximate halfwidth of the window, in Hz.
    windowHalfWidthBins = floor(windowHalfWidthHz / frequencyBinWidth);
    windowedPowerSpectra = movmean(powerSpectrum,2*windowHalfWidthBins+1,1) * (2*windowHalfWidthBins+1);

    % Aggregate the Fourier amplitude spectrum over arclengths.
    summedPowerSpectrum = trapz(arclengths,powerSpectrum,2);
    windowedSummedPowerSpectrum = trapz(arclengths,windowedPowerSpectra,2);
    proximalMask = arclengths < max(arclengths)/3;
    distalMask = arclengths >= 2*max(arclengths)/3;
    centralMask = ~(proximalMask | distalMask);
    windowedSummedPowerSpectrum = trapz(arclengths,windowedPowerSpectra,2);
    windowedSummedPowerSpectrumProximal = trapz(arclengths(proximalMask),windowedPowerSpectra(:,proximalMask),2);
    windowedSummedPowerSpectrumCentral = trapz(arclengths(centralMask),windowedPowerSpectra(:,centralMask),2);
    windowedSummedPowerSpectrumDistal = trapz(arclengths(distalMask),windowedPowerSpectra(:,distalMask),2);

    % Find the dominant frequency, ignoring the static mode and the lowest bins aggregated with it.
    [~, domFreqInd] = max(windowedSummedPowerSpectrum(windowHalfWidthBins+2:end)); domFreqInd = domFreqInd+windowHalfWidthBins+1;
    domFreq = frequencyDomain(domFreqInd);

    [~, domFreqIndProximal] = max(windowedSummedPowerSpectrumProximal(windowHalfWidthBins+2:end)); domFreqIndProximal = domFreqIndProximal+windowHalfWidthBins+1;
    domFreqProximal = frequencyDomain(domFreqIndProximal);
    [~, domFreqIndCentral] = max(windowedSummedPowerSpectrumCentral(windowHalfWidthBins+2:end)); domFreqIndCentral = domFreqIndCentral+windowHalfWidthBins+1;
    domFreqCentral = frequencyDomain(domFreqIndCentral);
    [~, domFreqIndDistal] = max(windowedSummedPowerSpectrumDistal(windowHalfWidthBins+2:end)); domFreqIndDistal = domFreqIndDistal+windowHalfWidthBins+1;
    domFreqDistal = frequencyDomain(domFreqIndDistal);

    % Aggregated amplitude of the dominant frequency.
    domFreqAmpSummed = windowedSummedPowerSpectrum(domFreqInd);

    % Get the phase of the dominant frequency along the arclength.
    domFreqPhase = unwrap(angle(complexSpectrum(domFreqInd,:)));

    % Compute a goodness-of-linear-fit metric for domFreqPhase.
    dummyX = [ones(numel(arclengths),1), arclengths(:)];
    dummyY = domFreqPhase(:);
    coeffs = dummyX \ dummyY;
    pred = dummyX * coeffs;
    phaseLinearityRsq = 1 - sum((dummyY - pred).^2)/sum((dummyY - mean(dummyY)).^2);
    if all(dummyY == mean(dummyY))
        phaseLinearityRsq = 1;
    end

    % Windowed amplitude of dominant frequency as a function of arclength.
    domFreqAmp = windowedPowerSpectra(domFreqInd,:);

    %% Fit Legendre polynomials to domFreqPhase and domFreqAmp.
    [legendreFitDomFreqPhase, legendreComponentsDomFreqPhase, legendreCoeffsDomFreqPhase] = legendreFit(domFreqPhase);
    [legendreFitDomFreqAmp, legendreComponentsDomFreqAmp, legendreCoeffsDomFreqAmp] = legendreFit(domFreqAmp);

    %% Continue with the Fourier analysis.

    % Windowed amplitude of static component as a function of arclength
    % (note that this is NOT static curvature, merely the average tangent
    % angle at each arclength).
    staticAngle = windowedPowerSpectra(1,:);
    staticAngleSummed = trapz(arclengths, staticAngle);

    % Frequency dominance of dominant frequency, measured as a proportion of
    % total power in spectrum without static part.
    domFreqDominance = domFreqAmpSummed / (sum(summedPowerSpectrum) - staticAngleSummed);

    % Reconstruct a beat using dominant frequency band. We'll do this by
    % removing all other non-static frequencies from the complex spectrum and
    % inverting. This retains the phase and amplitude information at each
    % arclength, but only for this single frequency +- windowHalfWidthHz.
    reducedComplexSpectra = 0*complexSpectrum;
    mask = [1,...
    domFreqInd + (-windowHalfWidthBins:windowHalfWidthBins),...
    2*length(frequencyDomain)-domFreqInd + (-windowHalfWidthBins:windowHalfWidthBins),...
    2:(windowHalfWidthBins+1),...
    2*length(frequencyDomain)-2-(0:windowHalfWidthBins-1)...
    ];
    reducedComplexSpectra(mask,:) = complexSpectrum(mask,:);
    anglesReconstructed = ifft(reducedComplexSpectra);
    % If we needed to pad by one, remove the final faux timepoint.
    anglesReconstructed = anglesReconstructed(1:numFrames,:);
    [xReconstructed,yReconstructed] = anglesToSpatial(anglesReconstructed,arclengths);

    % Quality of single-frequency fit, judged by angle. Equivalent to a normalised double integral of difference^2 over arclength and time.
    qualityOfFitAngle = sqrt(trapz(timestamps,trapz(arclengths,(angles - anglesReconstructed).^2,2),1) / trapz(timestamps,trapz(arclengths,angles.^2,2),1));

    % Quality of single-frequency fit, judged by xy coordinates, optimally rotated at the base.
    [qualityOfFitXYOptimal, xReconstructedOptimal, yReconstructedOptimal, optimalRotations] = computeQualityOfFitXY(x,y,xReconstructed,yReconstructed,arclengths,timestamps);

    % Number of phase cycles along the flagellum, crudely approximated as min-max.
    numPhaseCycles = (max(domFreqPhase) - min(domFreqPhase)) / (2*pi);

    % Crude approximation to beat period, 1/dominant frequency.
    periodFourier = 1/domFreq;

    % Try to find the period using the autocorrelation method, to compare.
    periodAutoCorr = NaN;
    periodStart = NaN;
    periodFramesAutoCorr = NaN;
    try
        [periodAutoCorr, periodStart, lag] = findPeriod(angles, timestamps);
        periodFramesAutoCorr = lag;
        findPeriodSuccess = true;
    catch
        findPeriodSuccess = false;
    end
    
    % NOTE I recommend using periodStart and periodFramesAutoCorr for period
    % calculation, but I've implemented using the Fourier analysis result.

    % Compute the number of bad frames.
    numBadFrames = sum(badFrameMask);

    % Approximate number of  frames per beat.
    framesPerBeatFourier = find(timestamps - timestamps(1) < periodFourier,1,'last');
    shortFileFlag = 5*framesPerBeatFourier > numFrames;
    singleAndFiveFrameAnalysisDone = ~(shortFileFlag || isempty(framesPerBeatFourier));
    if ~singleAndFiveFrameAnalysisDone
        if ~isempty(framesPerBeatFourier)
            message = [message, 'Warning: based on the estimated period, fewer than 5 good consecutive periods were captured. '];
        else
            message = [message, 'Warning: beat period less than one frame difference. Query period estimation. '];
        end
        disp(strjoin(message))
        startFrameSinglePeriod = 1;
        framesPerBeatFourier = 1;
        singlePeriodFrameRange = 1:numFrames;
        fivePeriodFrameRange = 1:numFrames;

        principalCurv = NaN;
        principalShearVel = NaN;
        reverseCurv = NaN;
        reverseShearVel = NaN;
        amplitudesSpatialFivePeriods = {NaN};
        amplitudeSpatialFivePeriods = NaN;
        amplitudesAngularFivePeriods = {NaN};
        amplitudeAngularFivePeriods = NaN;
        maxCurvSinglePeriodSmoothed = NaN;
        minCurvSinglePeriodSmoothed = NaN;
        smoothedAngles = {NaN};
        curvaturesWithinSingleBeat = {NaN};
        maxShearVelSinglePeriodSmoothed = NaN;
        minShearVelSinglePeriodSmoothed = NaN;

    else
        % Find a run of non-bad frames to use for visualisation/curvature, if one exists.
        startFrameSinglePeriod = find(conv(badFrameMask,ones(1,framesPerBeatFourier),'same')==0, 1, 'first');
        if isempty(startFrameSinglePeriod)
            [~, startFrameSinglePeriod] = min(conv(badFrameMask,ones(1,framesPerBeatFourier),'same'));
        end
        if (startFrameSinglePeriod + framesPerBeatFourier > numFrames)
            startFrameSinglePeriod = round(numFrames/2);
        end
        startFrameSinglePeriod = startFrameSinglePeriod - 1;
        singlePeriodFrameRange = startFrameSinglePeriod + (1:framesPerBeatFourier);
        fivePeriodFrameRange = startFrameSinglePeriod + (1:(5*framesPerBeatFourier));
        % Ensure that the fivePeriodFrameRange is within the limits.
        if fivePeriodFrameRange(end) > numFrames
            fivePeriodFrameRange = fivePeriodFrameRange - (fivePeriodFrameRange(end) - size(x,1));
        end

        % Find max and min curvaturesWithinSingleBeat. First, fit a smoothing spline to eliminate
        % noisy oscillations. Note that this fit is supposed to be crude, caring
        % only about the linear parts of the curves and not about the peaks in
        % angle. We'll do this only for the run of non-bad frames.
        smoothedAngles = zeros(framesPerBeatFourier,size(angles,2));
        curvaturesWithinSingleBeat = zeros(framesPerBeatFourier,size(angles,2));
        for i = singlePeriodFrameRange
            % We can only fit to data without NaNs.
            mask = ~isnan(angles(i,:));
            try
                f = fit(arclengths(mask)', double(angles(i,mask)'), 'smoothingspline','SmoothingParam',0.999);
                smoothedAngles(i-startFrameSinglePeriod,:) = f(arclengths)';
                curvaturesWithinSingleBeat(i-startFrameSinglePeriod,:) = gradient(f(arclengths),arclengths)';
            catch
                % There is not enough data to do the fit. Likely a bad frame.
            end
            curvaturesWithinSingleBeat(i-startFrameSinglePeriod,~mask) = NaN;
        end
        [maxCurvSinglePeriodSmoothed, maxInd] = max(curvaturesWithinSingleBeat(:));
        [maxFrameInd, maxArcInd] = ind2sub(size(curvaturesWithinSingleBeat),maxInd);
        [minCurvSinglePeriodSmoothed, minInd] = min(curvaturesWithinSingleBeat(:));
        [minFrameInd, minArcInd] = ind2sub(size(curvaturesWithinSingleBeat),minInd);
        
        % Compute the shear velocity (assuming constant curvature model).
        maxShearVelSinglePeriodSmoothed = gradient(smoothedAngles(:,maxArcInd), timestamps(singlePeriodFrameRange)); maxShearVelSinglePeriodSmoothed = abs(maxShearVelSinglePeriodSmoothed(maxFrameInd));
        minShearVelSinglePeriodSmoothed = gradient(smoothedAngles(:,minArcInd), timestamps(singlePeriodFrameRange)); minShearVelSinglePeriodSmoothed = abs(minShearVelSinglePeriodSmoothed(minFrameInd));
        
        % Distinguish between principal and reverse bends (principal has
        % highest abs curvature.
        if abs(maxCurvSinglePeriodSmoothed) >= abs(minCurvSinglePeriodSmoothed)
            principalCurv = abs(maxCurvSinglePeriodSmoothed);
            principalShearVel = maxShearVelSinglePeriodSmoothed;
            reverseCurv = abs(minCurvSinglePeriodSmoothed);
            reverseShearVel = minShearVelSinglePeriodSmoothed;
        else
            principalCurv = abs(minCurvSinglePeriodSmoothed);
            principalShearVel = minShearVelSinglePeriodSmoothed;
            reverseCurv = abs(maxCurvSinglePeriodSmoothed);
            reverseShearVel = maxShearVelSinglePeriodSmoothed;
        end

        amplitudesSpatialFivePeriods = computeSpatialAmplitudes(x(fivePeriodFrameRange,:), y(fivePeriodFrameRange,:));
        amplitudeSpatialFivePeriods = max(amplitudesSpatialFivePeriods);
        amplitudesAngularFivePeriods = computeAngularAmplitudes(angles(fivePeriodFrameRange,:));
        amplitudeAngularFivePeriods = max(amplitudesAngularFivePeriods);
    end

    % Compute the amplitudes as a function of arclength. This is from a single period.
    amplitudesSpatialSinglePeriod = computeSpatialAmplitudes(x(singlePeriodFrameRange,:), y(singlePeriodFrameRange,:));
    amplitudeSpatialSinglePeriod = max(amplitudesSpatialSinglePeriod);
    amplitudesAngularSinglePeriod = computeAngularAmplitudes(angles(singlePeriodFrameRange,:));
    amplitudeAngularSinglePeriod = max(amplitudesAngularSinglePeriod);

    amplitudesSpatialAllData = computeSpatialAmplitudes(x, y);
    amplitudeSpatialAllData = max(amplitudesSpatialAllData);
    amplitudesAngularAllData = computeAngularAmplitudes(angles);
    amplitudeAngularAllData = max(amplitudesAngularAllData);

    amplitudesAngularAllDataProximal = computeAngularAmplitudes(angles(:,proximalMask));
    amplitudeAngularAllDataProximal = max(amplitudesAngularAllDataProximal);
    amplitudesAngularAllDataCentral = computeAngularAmplitudes(angles(:,centralMask));
    amplitudeAngularAllDataCentral = max(amplitudesAngularAllDataCentral);
    amplitudesAngularAllDataDistal = computeAngularAmplitudes(angles(:,distalMask));
    amplitudeAngularAllDataDistal = max(amplitudesAngularAllDataDistal);
    
    message = strjoin(message);
    %% Save output summary.
    fh = fopen([filepath,filesep,'summary.txt'],'w');
    fprintf(fh,'%s/%s\n',filepath,filename);
    fprintf(fh,'Bad file (< 30 good frames or lots of interpolation):  %s.\n',mat2str(badFileFlag));
    fprintf(fh,'Short file (fewer than 5 beats captured):  %s.\n',mat2str(shortFileFlag));
    fprintf(fh,'Messages: %s\n',message);
    fprintf(fh,'Number of frames (total): %f\n',numFramesOrig);
    fprintf(fh,'Number of bad frames in file: %f\n',numBadFramesOrig);
    fprintf(fh,'Proportion of bad frames in file: %f\n',numBadFramesOrig / numFramesOrig);
    fprintf(fh,'Number of frames (analysed): %f\n',numFrames);
    fprintf(fh,'Number of bad frames (analysed): %f\n',numBadFrames);
    fprintf(fh,'Proportion of bad frames (analysed): %f\n',numBadFrames / numFrames);
    fprintf(fh,'Start frame: %f\n',bestFrameRange(1));
    fprintf(fh,'End frame: %f\n',bestFrameRange(end));
    fprintf(fh,'Single period start frame: %f\n',singlePeriodFrameRange(1));
    fprintf(fh,'Single period end frame: %f\n',singlePeriodFrameRange(end));
    fprintf(fh,'Flagellum length: %f\n',flagLength);
    fprintf(fh,'Analysed flagellum length: %f\n', arclengthThreshold*flagLength);
    fprintf(fh,'Sampling frequency (Hz): %f\n',samplingFrequency);
    fprintf(fh,'Dominant frequency (Hz): %f\n',domFreq);
    fprintf(fh,'Dominant frequency (proximal) (Hz): %f\n',domFreqProximal);
    fprintf(fh,'Dominant frequency (central) (Hz): %f\n',domFreqCentral);
    fprintf(fh,'Dominant frequency (distal) (Hz): %f\n',domFreqDistal);
    fprintf(fh,['Legendre coeffs of amplitude of dominant frequency: ' num2str(legendreCoeffsDomFreqAmp),'\n']);
    fprintf(fh,['Legendre coeffs of phase of dominant frequency: ' num2str(legendreCoeffsDomFreqPhase),'\n']);
    fprintf(fh,'Linearity of phase as a function of arclength: %f\n',phaseLinearityRsq);
    fprintf(fh,'Approximate period (1/frequency) (s): %f\n',periodFourier);
    fprintf(fh,'Approximate frames per beat: %f\n',framesPerBeatFourier);
    fprintf(fh,'Overall aggregated angular amplitude of dominant frequency: %f\n',domFreqAmpSummed);
    fprintf(fh,'Amp. dominant frequency / all (non-static) amp: %f\n',domFreqDominance);
    fprintf(fh,'Halfwidth of frequency bin (Hz): %f\n',windowHalfWidthHz);
    fprintf(fh,'Overall aggregated static angular amplitude: %f\n',staticAngleSummed);
    fprintf(fh,'Relative error in reconstructed angle (0 is perfect): %f\n',qualityOfFitAngle);
    fprintf(fh,'Overall error in reconstructed beat (0 is perfect): %f\n',qualityOfFitXYOptimal);
    fprintf(fh,'Number of phase cycles along flagellum: %f\n',numPhaseCycles);
    if singleAndFiveFrameAnalysisDone
        fprintf(fh,'Max signed curvature: %f\n',maxCurvSinglePeriodSmoothed);
        fprintf(fh,'Min signed curvature: %f\n',minCurvSinglePeriodSmoothed);
        fprintf(fh,'kappa_P: %f\n',principalCurv);
        fprintf(fh,'kappa_R: %f\n',reverseCurv);
        fprintf(fh,'V_P: %f\n',principalShearVel);
        fprintf(fh,'V_R: %f\n',reverseShearVel);
    end
    if findPeriodSuccess
        fprintf(fh,'Frequency from autocorrelation (Hz): %f\n',1/periodAutoCorr);
        fprintf(fh,'Period from autocorrelation (s): %f\n',periodAutoCorr);
    else
        fprintf(fh,'Not enough good data for autocorrelation method (< 3 frames)');
    end
    fprintf(fh,'Spatial amplitude, single period: %f\n',amplitudeSpatialSinglePeriod);
    if ~shortFileFlag
        fprintf(fh,'Spatial amplitude, five periods: %f\n',amplitudeSpatialFivePeriods);
    else
        fprintf(fh,'Less than five periods long');
    end
    fprintf(fh,'Spatial amplitude, all data: %f\n',amplitudeSpatialAllData);
    fprintf(fh,'Fourier amplitude, single period (rad): %f\n',amplitudeAngularSinglePeriod);
    if ~shortFileFlag
        fprintf(fh,'Fourier amplitude, five periods (rad): %f\n',amplitudeAngularFivePeriods);
    else
        fprintf(fh,'Less than five periods long');
    end
    fprintf(fh,'Fourier amplitude, all data (rad): %f\n',amplitudeAngularAllData);
    fclose(fh); 

    %% Make and save plots.

    close all
    f1 = figure('visible','off');

    clf
    plotWaveform(x,y);
    xlims = xlim; ylims = ylim;
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'beat.png'])%,'ContentType','vector')

    clf
    plotWaveform(xOrig,yOrig);
    xlim(xlims), ylim(ylims)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'beatIncBadData.png'])%,'ContentType','vector')

    clf
    plotWaveform(xReconstructed,yReconstructed);
    xlim(xlims), ylim(ylims)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'reconstructedFromDomFreqBeat.png'])%,'ContentType','vector')

    clf
    tiledlayout(1,2,'TileSpacing','compact')
    nexttile()
    plotWaveform(x,y);
    xlim(xlims), ylim(ylims)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    title('Full')
    nexttile()
    plotWaveform(xReconstructed,yReconstructed);
    xlim(xlims), ylim(ylims)
    xlabel('$x$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24,'YTickLabel',[])
    title('Reconstructed')
    exportgraphics(gcf,[filepath,filesep,'beatVsReconstructedFromDomFreqBeat.png'])%,'ContentType','vector')

    clf
    tiledlayout(1,2,'TileSpacing','compact')
    nexttile()
    plotWaveform(x,y);
    xlim(xlims), ylim(ylims)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    title('Full')
    nexttile()
    plotWaveform(xReconstructedOptimal,yReconstructedOptimal);
    xlim(xlims), ylim(ylims)
    xlabel('$x$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24,'YTickLabel',[])
    th = title({'Reconstructed (optimally rotated).';['Score: ',num2str(qualityOfFitXYOptimal)]});
    th.FontSize = 16;
    exportgraphics(gcf,[filepath,filesep,'beatVsReconstructedFromDomFreqBeatOptimallyRotated.png'])%,'ContentType','vector')

    clf
    plotWaveform(x(startFrameSinglePeriod+(1:framesPerBeatFourier),:),y(startFrameSinglePeriod+(1:framesPerBeatFourier),:));
    xlim(xlims), ylim(ylims)
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'singlePeriodBeat.png'])%,'ContentType','vector')    

    clf
    plot(angles(startFrameSinglePeriod+(1:framesPerBeatFourier),:)')
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Shear angle (rad)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'singlePeriodShearAngle.png'])%,'ContentType','vector')

    clf
    plot(arclengths, domFreqPhase,'Color','black')
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Phase (rad)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'fourierPhaseAlongArclength.png'])%,'ContentType','vector')

    clf
    plot(arclengths, domFreqAmp,'Color','black')
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Amplitude (rad, windowed)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'fourierAmplitudeAlongArclength.png'])%,'ContentType','vector')

    clf
    plot(frequencyDomain, summedPowerSpectrum,'Color','black')
    xline(domFreq)
    xlabel('Frequency (Hz)','Interpreter','latex')
    ylabel('Amplitude (rad)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'summedPowerSpectrum.png'])%,'ContentType','vector')

    clf
    [meshF, meshA] = meshgrid(frequencyDomain, arclengths);
    surf(meshF, meshA, powerSpectrum','LineStyle','none')
    xlabel('Frequency (Hz)','Interpreter','latex')
    ylabel('Arclength ($\mu$m)','Interpreter','latex')
    view(0,90)
    shading interp
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    try
        colormap(viridis)
    end
    axis tight
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'fourierAmplitudeSpectrum.png'])%,'ContentType','vector')

    clf
    plot(arclengths, amplitudesSpatialSinglePeriod)
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Amplitude ($\mu$m)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'spatialAmplitudesSinglePeriod.png'])%,'ContentType','vector')

    clf
    plot(arclengths, amplitudesSpatialAllData)
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Amplitude ($\mu$m)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'spatialAmplitudesAllData.png'])%,'ContentType','vector')

    clf
    plot(arclengths, amplitudesAngularSinglePeriod)
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Amplitude (rad)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'fourierAmplitudesSinglePeriod.png'])%,'ContentType','vector')

    clf
    plot(arclengths, amplitudesAngularAllData)
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Amplitude (rad)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'fourierAmplitudesAllData.png'])%,'ContentType','vector')

    clf
    surf(arclengths,timestamps,curvatures)
    xlabel('Arclength ($\mu$m)','Interpreter','latex')
    ylabel('Timestamp','Interpreter','latex')
    view(0,90)
    shading interp
    c = colorbar;
    c.TickLabelInterpreter = 'latex';
    c.Label.String = 'Curvature';
    c.Label.Interpreter = 'latex';
    caxis(0.5*[-1,1])
    try
        colormap(viridis)
    end
    axis tight
    set(gca,'FontSize',24,'YDir','reverse')
    exportgraphics(gcf,[filepath,filesep,'curvatureKymograph.png'])%,'ContentType','vector')

    clf
    hold on
    tempSes = linspace(0,1,numel(domFreqAmp));
    plot(tempSes, domFreqAmp,'o','LineWidth',1,'Color','black')
    plot(tempSes, legendreFitDomFreqAmp,'LineWidth',2,'Color','black')
    plot(tempSes, legendreComponentsDomFreqAmp)
    xlabel('Arclength / Flag Length','Interpreter','latex')
    ylabel('Amplitude (rad)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'legendreDomFreqAmp.png'])%,'ContentType','vector')

    clf
    hold on
    tempSes = linspace(0,1,numel(domFreqPhase));
    plot(tempSes, domFreqPhase,'o','LineWidth',1,'Color','black')
    plot(tempSes, legendreFitDomFreqPhase,'LineWidth',2,'Color','black')
    plot(tempSes, legendreComponentsDomFreqPhase)
    xlabel('Arclength / Flag Length','Interpreter','latex')
    ylabel('Phase (rad)','Interpreter','latex')
    grid on
    set(gca,'FontSize',24)
    exportgraphics(gcf,[filepath,filesep,'legendreDomFreqPhase.png'])%,'ContentType','vector')

    close(f1)

    clear ans bodyCOM bodyMask c cellOrientationAngle cellOrientationVector domFreqInd f f1 fh flagBase flagMask flagOnRight frame i lag mask maxArcInd maxFrameInd maxInd meshA meshF minArcInd minFrameInd minInd tempSes th uniformArclengths windowHalfWidthBins windowHalfWidthHz xlims ylims amplitudesAngularAllDataCentral amplitudesAngularAllDataDistal amplitudesAngularAllDataProximal dummyX dummyY coeffs pred
    save([filepath,filesep,'output.mat'])

end


% TODO

% Speed of phase (phase of dominant frequency) travelling along flagellum