clc
clear all
s = rng(211);       % Set RNG state for repeatability
numFFT = 512;        % number of FFT points
subbandSize = 20;    % must be > 1
numSubbands = 10;    % numSubbands*subbandSize <= numFFT
subbandOffset = 156; % numFFT/2-subbandSize*numSubbands/2 for band center
 
% Dolph-Chebyshev window design parameters
filterLen = 43;      % similar to cyclic prefix length
slobeAtten = 40;       % side-lobe attenuation, dB
 
bitsPerSubCarrier = 4;   % 2: 4QAM, 4: 16QAM, 6: 64QAM, 8: 256QAM
snrdB = 2;              % SNR in dB
% Design window with specified attenuation
prototypeFilter = chebwin(filterLen, slobeAtten);
% QAM Symbol mapper
qamMapper = comm.RectangularQAMModulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');
% Transmit-end processing
%  Initialize arrays
inpData = zeros(bitsPerSubCarrier*subbandSize, numSubbands);
txSig = complex(zeros(numFFT+filterLen-1, 1));
%  Loop over each subband
for bandIdx = 1:numSubbands
 
    bitsIn = randi([0 1], bitsPerSubCarrier*subbandSize, 1);
    symbolsIn = qamMapper(bitsIn);
    inpData(:,bandIdx) = bitsIn; % log bits for comparison
 
    % Pack subband data into an OFDM symbol
    offset = subbandOffset+(bandIdx-1)*subbandSize;
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT-offset-subbandSize, 1)];
    ifftOut = ifft(ifftshift(symbolsInOFDM));
 
    % Filter for each subband is shifted in frequency
    bandFilter = prototypeFilter.*exp( 1i*2*pi*(0:filterLen-1)'/numFFT* ...
                 ((bandIdx-1/2)*subbandSize+0.5+subbandOffset+numFFT/2) );
    filterOut = conv(bandFilter,ifftOut);
 
    % Plot power spectral density (PSD) per subband
    [psd,f] = periodogram(filterOut, rectwin(length(filterOut)), ...
                          numFFT*2, 1, 'centered');
    plot(f,10*log10(psd));
 
    % Sum the filtered subband responses to form the aggregate transmit
    % signal
    txSig = txSig + filterOut;
end
 
% Add WGN
rxSig = awgn(txSig, snrdB, 'measured');
% Pad receive vector to twice the FFT Length (note use of txSig as input)
%   No windowing or additional filtering adopted
yRxPadded = [rxSig; zeros(2*numFFT-numel(txSig),1)];
 
% Perform FFT and downsample by 2
RxSymbols2x = fftshift(fft(yRxPadded));
RxSymbols = RxSymbols2x(1:2:end);
 
% Select data subcarriers
dataRxSymbols = RxSymbols(subbandOffset+(1:numSubbands*subbandSize));
% Use zero-forcing equalizer after OFDM demodulation
rxf = [prototypeFilter.*exp(1i*2*pi*0.5*(0:filterLen-1)'/numFFT); ...
       zeros(numFFT-filterLen,1)];
prototypeFilterFreq = fftshift(fft(rxf));
prototypeFilterInv = 1./prototypeFilterFreq(numFFT/2-subbandSize/2+(1:subbandSize));
% Equalize per subband - undo the filter distortion
dataRxSymbolsMat = reshape(dataRxSymbols,subbandSize,numSubbands);
EqualizedRxSymbolsMat = bsxfun(@times,dataRxSymbolsMat,prototypeFilterInv);
EqualizedRxSymbols = EqualizedRxSymbolsMat(:);
% Demapping and BER computation
qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
    2^bitsPerSubCarrier, 'BitOutput', true, ...
    'NormalizationMethod', 'Average power');
BER = comm.ErrorRate;
% Perform hard decision and measure errors
rxBits = qamDemod(EqualizedRxSymbols);
ber = BER(inpData(:), rxBits);
 
disp(['UFMC Reception, BER = ' num2str(ber(1)) ' at SNR = ' ...
    num2str(snrdB) ' dB']);
 
% Restore RNG state
rng(s);
 
%% Plotting BER vs SNR
semilogy(ber,snrdB,'-ok');
grid;
title('ufmc Bit Error Rate .VS. Signal To Noise Ratio');
ylabel('BER');
xlabel('SNR [dB]');
