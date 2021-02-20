numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center

% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;     % side-lobe attenuation, dB

bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM

snrdB = 15;              % SNR in Db
symbolsIn = qamMapper(inpData(:));

% Process all sub-bands together
offset = subbandOffset;
symbolsInOFDM = [zeros(offset, 1); symbolsIn; ...
                 zeros(numFFT-offset-subbandSize*numSubbands, 1)];
ifftOut = sqrt(numFFT).*ifft(ifftshift(symbolsInOFDM));

% Plot power spectral density (PSD) over all subcarriers
[psd,f] = periodogram(ifftOut, rectwin(length(ifftOut)), numFFT*2, ...
                      1, 'centered');
hFig1 = figure;
plot(f,10*log10(psd));
grid on
axis([-0.5 0.5 -100 20]);
xlabel('Normalized frequency');
ylabel('PSD (dBW/Hz)')
title(['OFDM, ' num2str(numSubbands*subbandSize) ' Subcarriers'])
set(hFig1, 'Position', figposition([46 50 25 30]));

% Compute peak-to-average-power ratio (PAPR)
PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
[~,~,paprOFDM] = PAPR2(ifftOut);
disp(['Peak-to-Average-Power-Ratio (PAPR) for OFDM = ' num2str(paprOFDM) ' dB']);
