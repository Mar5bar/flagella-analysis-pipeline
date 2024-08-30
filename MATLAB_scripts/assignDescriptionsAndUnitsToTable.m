function T = assignDescriptionsAndUnitsToTable(T)
%% Assign units and descriptions to variables in the table, which contains
% data generated by processFile.

%% Units.
T.Properties.VariableUnits{'amplitudeAngularAllData'} = 'radians';
T.Properties.VariableUnits{'amplitudeAngularFivePeriods'} = 'radians';
T.Properties.VariableUnits{'amplitudeAngularSinglePeriod'} = 'radians';
T.Properties.VariableUnits{'amplitudesAngularAllData'} = 'radians';
T.Properties.VariableUnits{'amplitudesAngularFivePeriods'} = 'radians';
T.Properties.VariableUnits{'amplitudesAngularSinglePeriod'} = 'radians';
T.Properties.VariableUnits{'amplitudeSpatialAllData'} = 'um';
T.Properties.VariableUnits{'amplitudeSpatialFivePeriods'} = 'um';
T.Properties.VariableUnits{'amplitudeSpatialSinglePeriod'} = 'um';
T.Properties.VariableUnits{'amplitudesSpatialAllData'} = 'um';
T.Properties.VariableUnits{'amplitudesSpatialFivePeriods'} = 'um';
T.Properties.VariableUnits{'amplitudesSpatialSinglePeriod'} = 'um';
T.Properties.VariableUnits{'angles'} = 'radians';
T.Properties.VariableUnits{'anglesBeforeInterp'} = 'radians';
T.Properties.VariableUnits{'anglesReconstructed'} = 'radians';
T.Properties.VariableUnits{'arclengths'} = 'um';
T.Properties.VariableUnits{'bestFrameRange'} = 'frames';
T.Properties.VariableUnits{'curvatures'} = 'radians / um';
T.Properties.VariableUnits{'curvaturesWithinSingleBeat'} = 'radians / um';
T.Properties.VariableUnits{'curvatureVelocity'} = 'radians / um / s';
T.Properties.VariableUnits{'domFreq'} = 'Hz';
T.Properties.VariableUnits{'domFreqCentral'} = 'Hz';
T.Properties.VariableUnits{'domFreqDistal'} = 'Hz';
T.Properties.VariableUnits{'domFreqProximal'} = 'Hz';
T.Properties.VariableUnits{'domFreqAmp'} = 'radians';
T.Properties.VariableUnits{'domFreqAmpSummed'} = 'radians um';
T.Properties.VariableUnits{'domFreqDominance'} = '1';
T.Properties.VariableUnits{'flagLength'} = 'um';
T.Properties.VariableUnits{'framesPerBeatFourier'} = 'frames';
T.Properties.VariableUnits{'frequencyBinWidth'} = 'Hz';
T.Properties.VariableUnits{'frequencyDomain'} = 'Hz';
T.Properties.VariableUnits{'legendreCoeffsDomFreqAmp'} = 'various';
T.Properties.VariableUnits{'legendreCoeffsDomFreqPhase'} = 'various';
T.Properties.VariableUnits{'legendreComponentsDomFreqAmp'} = 'radians';
T.Properties.VariableUnits{'legendreComponentsDomFreqPhase'} = 'radians';
T.Properties.VariableUnits{'legendreFitDomFreqAmp'} = 'radians';
T.Properties.VariableUnits{'legendreFitDomFreqPhase'} = 'radians';
T.Properties.VariableUnits{'maxArclengths'} = 'um';
T.Properties.VariableUnits{'maxCurvSinglePeriodSmoothed'} = 'radians / um';
T.Properties.VariableUnits{'maxShearVelSinglePeriodSmoothed'} = 'radians / second';
T.Properties.VariableUnits{'minCurvSinglePeriodSmoothed'} = 'radians / um';
T.Properties.VariableUnits{'minShearVelSinglePeriodSmoothed'} = 'radians / second';
T.Properties.VariableUnits{'numArclengthSamples'} = 'points';
T.Properties.VariableUnits{'numBadFrames'} = 'frames';
T.Properties.VariableUnits{'numBadFramesOrig'} = 'frames';
T.Properties.VariableUnits{'numFrames'} = 'frames';
T.Properties.VariableUnits{'numFramesOrig'} = 'frames';
T.Properties.VariableUnits{'optimalRotations'} = 'radians';
T.Properties.VariableUnits{'overallLength'} = 'um';
T.Properties.VariableUnits{'overallMaxArclengths'} = 'um';
T.Properties.VariableUnits{'periodAutoCorr'} = 'seconds';
T.Properties.VariableUnits{'periodStart'} = 'frames';
T.Properties.VariableUnits{'periodFramesAutoCorr'} = 'frames';
T.Properties.VariableUnits{'periodFourier'} = 'seconds';
T.Properties.VariableUnits{'phaseLinearityRsq'} = 'none';
T.Properties.VariableUnits{'principalCurv'} = 'radians / um';
T.Properties.VariableUnits{'principalShearVel'} = 'radians / second';
T.Properties.VariableUnits{'reverseCurv'} = 'radians / um';
T.Properties.VariableUnits{'reverseShearVel'} = 'radians / second';
T.Properties.VariableUnits{'samplingFrequency'} = 'Hz';
T.Properties.VariableUnits{'singlePeriodFrameRange'} = 'frames';
T.Properties.VariableUnits{'smoothedAngles'} = 'radians';
T.Properties.VariableUnits{'spatialScale'} = 'um per pixel';
T.Properties.VariableUnits{'startFrameSinglePeriod'} = 'frames';
T.Properties.VariableUnits{'staticAngle'} = 'radians';
T.Properties.VariableUnits{'staticAngleSummed'} = 'radians um';
T.Properties.VariableUnits{'timestamps'} = 'seconds';
T.Properties.VariableUnits{'widths'} = 'pixels';
T.Properties.VariableUnits{'x'} = 'um';
T.Properties.VariableUnits{'xBeforeInterp'} = 'um';
T.Properties.VariableUnits{'xOrig'} = 'um';
T.Properties.VariableUnits{'xReconstructed'} = 'um';
T.Properties.VariableUnits{'xReconstructedOptimal'} = 'um';
T.Properties.VariableUnits{'y'} = 'um';
T.Properties.VariableUnits{'yBeforeInterp'} = 'um';
T.Properties.VariableUnits{'yOrig'} = 'um';
T.Properties.VariableUnits{'yReconstructed'} = 'um';
T.Properties.VariableUnits{'yReconstructedOptimal'} = 'um';


%% Descriptions.
T.Properties.VariableDescriptions{'amplitudeAngularAllData'} = 'Maximum amplitude of tangent angle oscillations, computed from all the data.';
T.Properties.VariableDescriptions{'amplitudeAngularFivePeriods'} = 'Maximum amplitude of tangent angle oscillations, computed from five periods.';
T.Properties.VariableDescriptions{'amplitudeAngularSinglePeriod'} = 'Maximum amplitude of tangent angle oscillations, computed from a single period.';
T.Properties.VariableDescriptions{'amplitudesAngularAllData'} = 'Amplitude of tangent angle oscillations at material points, computed from all the data.';
T.Properties.VariableDescriptions{'amplitudesAngularFivePeriods'} = 'Amplitude of tangent angle oscillations at material points, computed from five periods.';
T.Properties.VariableDescriptions{'amplitudesAngularSinglePeriod'} = 'Amplitude of tangent angle oscillations at material points, computed from a single period.';
T.Properties.VariableDescriptions{'amplitudeSpatialAllData'} = 'Maximum amplitude of spatial oscillations, computed from all the data.';
T.Properties.VariableDescriptions{'amplitudeSpatialFivePeriods'} = 'Maximum amplitude of spatial oscillations, computed from five periods.';
T.Properties.VariableDescriptions{'amplitudeSpatialSinglePeriod'} = 'Maximum amplitude of spatial oscillations, computed from a single period.';
T.Properties.VariableDescriptions{'amplitudesSpatialAllData'} = 'Amplitude of spatial oscillations of material points, computed from all the data.';
T.Properties.VariableDescriptions{'amplitudesSpatialFivePeriods'} = 'Amplitude of spatial oscillations of material points, computed from five periods.';
T.Properties.VariableDescriptions{'amplitudesSpatialSinglePeriod'} = 'Amplitude of spatial oscillations of material points, computed from a single period.';
T.Properties.VariableDescriptions{'angles'} = 'Tangent angle parameterisation of the flagellum, smoothed and interpolated.';
T.Properties.VariableDescriptions{'anglesBeforeInterp'} = 'Tangent angle parameterisation of the flagellum, before interpolation in time.';
T.Properties.VariableDescriptions{'anglesReconstructed'} = 'Tangent angle parameterisation of the flagellum, reconstructed from the dominant frequency.';
T.Properties.VariableDescriptions{'arclengths'} = 'Arclengths of material points along the flagellum, uniformly spaced.';
T.Properties.VariableDescriptions{'badFileFlag'} = 'Is the file too short / contains too many bad frames?';
T.Properties.VariableDescriptions{'badFrameMask'} = 'Logical mask identifying bad frames from the original data, where a bad frame is one capturing less than the threshold proportion of both the flagellum length and the total trace length.';
T.Properties.VariableDescriptions{'bestFrameRange'} = 'Longest range of frames identified as having minimal bad frames, which the data is truncated to.';
T.Properties.VariableDescriptions{'complexSpectrum'} = 'Complex Fourier spectrum of the time-truncated angle parameterisation.';
T.Properties.VariableDescriptions{'curvatures'} = 'Curvature of the flagellum.';
T.Properties.VariableDescriptions{'curvaturesWithinSingleBeat'} = 'Curvature of the flagellum over a single period.';
T.Properties.VariableDescriptions{'curvatureVelocity'} = 'Signed speed of curvature moving along the flagellum (base to tip).';
T.Properties.VariableDescriptions{'domFreq'} = 'Dominant frequency identified by Fourier analysis of the tangent-angle parameterisation.';
T.Properties.VariableDescriptions{'domFreqCentral'} = 'Dominant frequency identified by Fourier analysis of the tangent-angle parameterisation on the central third of the flagellum.';
T.Properties.VariableDescriptions{'domFreqDistal'} = 'Dominant frequency identified by Fourier analysis of the tangent-angle parameterisation on the distal third of the flagellum.';
T.Properties.VariableDescriptions{'domFreqProximal'} = 'Dominant frequency identified by Fourier analysis of the tangent-angle parameterisation on the proximal third of the flagellum.';
T.Properties.VariableDescriptions{'domFreqAmp'} = 'Amplitude at each arclength of the dominant frequency identified by Fourier analysis of the tangent-angle parameterisation.';
T.Properties.VariableDescriptions{'domFreqAmpSummed'} = 'Amplitude of the dominant frequency identified by Fourier analysis of the tangent-angle parameterisation, integrated over the flagellum.';
T.Properties.VariableDescriptions{'domFreqDominance'} = 'Measure of the dominance of the dominant frequency identified by Fourier analysis of the tangent-angle parameterisation, computed as the ratio between the power of this frequency band and the total non-static power of the spectrum.';
T.Properties.VariableDescriptions{'domFreqPhase'} = 'Phase of the dominant frequency identified by Fourier analysis of the tangent-angle parameterisation.';
T.Properties.VariableDescriptions{'filename'} = 'File name from which coordinate and width data was loaded.';
T.Properties.VariableDescriptions{'filepath'} = 'Directory name from which data has been loaded, processed, and saved.';
T.Properties.VariableDescriptions{'fivePeriodFrameRange'} = 'Range of good frames covering approximately five periods of the flagellar beat.';
T.Properties.VariableDescriptions{'flagLength'} = 'Estimated length of the flagellum, before truncation.';
T.Properties.VariableDescriptions{'framesPerBeatFourier'} = 'Estimated number of frames per beat period, computed from the Fourier estimate of the beat frequency.';
T.Properties.VariableDescriptions{'frequencyBinWidth'} = 'Width of the bins used to group nearby frequencies when doing Fourier analysis.';
T.Properties.VariableDescriptions{'frequencyDomain'} = 'Frequency domain used for Fourier analysis.';
T.Properties.VariableDescriptions{'legendreCoeffsDomFreqAmp'} = 'Best-fit Legendre coefficients computed from domFreqAmp.';
T.Properties.VariableDescriptions{'legendreCoeffsDomFreqPhase'} = 'Best-fit Legendre coefficients computed from domFreqPhase.';
T.Properties.VariableDescriptions{'legendreComponentsDomFreqAmp'} = 'Best-fit Legendre components (constant, linear, quadratic) computed from domFreqAmp.';
T.Properties.VariableDescriptions{'legendreComponentsDomFreqPhase'} = 'Best-fit Legendre components (constant, linear, quadratic) computed from domFreqPhase.';
T.Properties.VariableDescriptions{'legendreFitDomFreqAmp'} = 'Best-fit curve to domFreqAmp, made up of Legendre components (constant, linear, quadratic).';
T.Properties.VariableDescriptions{'legendreFitDomFreqPhase'} = 'Best-fit curve to domFreqPhase, made up of Legendre components (constant, linear, quadratic).';
T.Properties.VariableDescriptions{'maxArclengths'} = 'Maximum arclength recorded in each frame of raw data.';
T.Properties.VariableDescriptions{'maxCurvSinglePeriodSmoothed'} = 'Maximum absolute curvature recorded in a single period from smoothed angular data.';
T.Properties.VariableDescriptions{'maxShearVelSinglePeriodSmoothed'} = 'Maximum shear velocity recorded in a single period from smoothed angular data.';
T.Properties.VariableDescriptions{'minCurvSinglePeriodSmoothed'} = 'Minimum absolute curvature recorded in a single period from smoothed angular data.';
T.Properties.VariableDescriptions{'minShearVelSinglePeriodSmoothed'} = 'Minimum shear velocity recorded in a single period from smoothed angular data.';
T.Properties.VariableDescriptions{'numArclengthSamples'} = 'Number of points sampled along the flagellum for analysis.';
T.Properties.VariableDescriptions{'numBadFrames'} = 'Number of rejected frames in the truncated frame range.';
T.Properties.VariableDescriptions{'numBadFramesOrig'} = 'Number of rejected frames in the original frame range.';
T.Properties.VariableDescriptions{'numFrames'} = 'Number of frames in the truncated frame range.';
T.Properties.VariableDescriptions{'numFramesOrig'} = 'Number of frames in the original frame range.';
T.Properties.VariableDescriptions{'numPhaseCycles'} = 'Number of full cycles of the phase of domFreq in the Fourier analysis between the base and the tip of the flagellum.';
T.Properties.VariableDescriptions{'optimalRotations'} = 'Optimal angles by which to rotate the reconstructed spatial data to give the best fit to the full smoothed spatial data in (x,y).';
T.Properties.VariableDescriptions{'overallLength'} = 'Estimated length of the whole cell, including the body and flagellum.';
T.Properties.VariableDescriptions{'overallMaxArclengths'} = 'Estimated length of the whole cell in each frame, including the body and flagellum.';
try
    T.Properties.VariableDescriptions{'periodAutoCorr'} = 'Estimated beating period, computed using the autocorrelation of the angle parameterisation.';
    T.Properties.VariableDescriptions{'periodFramesAutoCorr'} = 'Estimated number of frames in a beating period, computed using the autocorrelation of the angle parameterisation.';
    T.Properties.VariableDescriptions{'periodStart'} = 'Estimated starting frame of the most periodic period of beating, computed using the autocorrelation of the angle parameterisation.';
end
T.Properties.VariableDescriptions{'periodFourier'} = 'Estimated beating period, computed using the frequency identified by Fourier analysis.';
T.Properties.VariableDescriptions{'phaseLinearityRsq'} = 'Measure of linearity (R^2) of phase as a function of arclength.';
T.Properties.VariableDescriptions{'powerSpectrum'} = 'Fourier power spectrum of the time-truncated angle parameterisation.';
T.Properties.VariableDescriptions{'principalCurv'} = 'Absolute value of the maximum curvature over a single beating period.';
T.Properties.VariableDescriptions{'principalShearVel'} = 'Maximum shear velocity recorded in a single period from smoothed angular data.';
T.Properties.VariableDescriptions{'qualityOfFitXYOptimal'} = 'Measure of how poor the single-frequency reconstruction fits the full data (x,y), allowing for rotation of each frame, measured using spatial coordinates. Zero is a perfect fit.';
T.Properties.VariableDescriptions{'qualityOfFitAngle'} = 'Measure of how poor the single-frequency reconstruction fits the full data (x,y), measured using the angle parameterisation. Zero is a perfect fit.';
T.Properties.VariableDescriptions{'reverseCurv'} = 'Absolute value of the minimum curvature over a single beating period.';
T.Properties.VariableDescriptions{'reverseShearVel'} = 'Minimum shear velocity recorded in a single period from smoothed angular data.';
T.Properties.VariableDescriptions{'samplingFrequency'} = 'Frame rate of the imaging data.';
T.Properties.VariableDescriptions{'shortFileFlag'} = 'Does the file have fewer than five periods of good data?';
T.Properties.VariableDescriptions{'singleAndFiveFrameAnalysisDone'} = 'Have we been able to identify and analyse both a single period and a five-period range?';
T.Properties.VariableDescriptions{'singlePeriodFrameRange'} = 'Range of frames corresponding to a single period.';
T.Properties.VariableDescriptions{'smoothedAngles'} = 'Tangent angle parameterisation of the flagellum over a single period, smoothed again for curvature computation.';
T.Properties.VariableDescriptions{'spatialScale'} = 'Number of micrometres per pixel in the raw data.';
T.Properties.VariableDescriptions{'startFrameSinglePeriod'} = 'Starting frame for a single period.';
T.Properties.VariableDescriptions{'staticAngle'} = 'Zero-frequency component of the Fourier spectrum, computed from the angle parameterisation.';
T.Properties.VariableDescriptions{'staticAngleSummed'} = 'Zero-frequency component of the Fourier spectrum, computed from the angle parameterisation and integrated over the flagellum.';
T.Properties.VariableDescriptions{'timestamps'} = 'Timestamps of the truncated frame range.';
T.Properties.VariableDescriptions{'widths'} = 'Widths of the original imported trace, including the body and the flagellum.';
T.Properties.VariableDescriptions{'x'} = 'x coordinates of the beat over the truncated frame range, uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'xBeforeInterp'} = 'x coordinates of the beat over the truncated frame range, before interpolation in time, uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'xOrig'} = 'x coordinates of the beat over the entire frame range, uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'xReconstructed'} = 'x coordinates of the beat over the truncated frame range, reconstructed from the dominant frequency and uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'xReconstructedOptimal'} = 'x coordinates of the beat over the truncated frame range, reconstructed from the dominant frequency, optimally reoriented to match x, and uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'y'} = 'y coordinates of the beat over the truncated frame range, uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'yBeforeInterp'} = 'y coordinates of the beat over the truncated frame range, before interpolation in time, uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'yOrig'} = 'y coordinates of the beat over the entire frame range, uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'yReconstructed'} = 'y coordinates of the beat over the truncated frame range, reconstructed from the dominant frequency and uniformly sampled in arclength.';
T.Properties.VariableDescriptions{'yReconstructedOptimal'} = 'y coordinates of the beat over the truncated frame range, reconstructed from the dominant frequency, optimally reoriented to match y, and uniformly sampled in arclength.';

end