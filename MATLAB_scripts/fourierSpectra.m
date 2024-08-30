function [complexSpectrum, amplitudeSpectrumOneSided, frequencyDomain] = fourierSpectra(data, sampleRate);
    if nargin < 2
        sampleRate = 1;
    end

    L = size(data,1);
    if mod(L,2) == 1
        L = L+1;
    end
    frequencyDomain = sampleRate*(0:(L/2))/L;

    % Compute the Fourier transform wrt time of the data.
    complexSpectrum = fft(data,L,1);

    % Compute the 2-sided spectrum and convert to 1-sided power spectrum.
    amplitudeSpectrumTwoSided = abs(complexSpectrum/L);
    amplitudeSpectrumOneSided = amplitudeSpectrumTwoSided(1:L/2 + 1,:); 
    amplitudeSpectrumOneSided(2:end-1,:) = 2 * amplitudeSpectrumOneSided(2:end-1,:);

end