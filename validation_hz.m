%% FIGURE 10 - STATIONARY 10 Hz VALIDATION
% MSc Dissertation
% Validation of Filter-Hilbert, STFT and Multitaper methods

clear;
clc;
close all;

%% ---------------------------------------------------------
% 1. CREATE SYNTHETIC 10 Hz SIGNAL
% ---------------------------------------------------------

fs = 500;                   % Sampling frequency (Hz)
duration = 5;               % Signal duration (seconds)

t = 0:1/fs:duration-1/fs;

f0 = 10;                    % Known signal frequency (Hz)

% Pure stationary 10 Hz sinusoid
x = sin(2*pi*f0*t);


%% ---------------------------------------------------------
% 2. FILTER-HILBERT ANALYSIS
% ---------------------------------------------------------

% Band-pass range surrounding 10 Hz
lowerFreq = 4;
upperFreq = 15;

% FIR filter
filterOrder = 192;

b = fir1(filterOrder, ...
    [lowerFreq upperFreq]/(fs/2), ...
    'bandpass');

filteredSignal = filtfilt(b,1,x);

% Hilbert transform
analyticSignal = hilbert(filteredSignal);

instantaneousAmplitude = abs(analyticSignal);
instantaneousPhase = unwrap(angle(analyticSignal));

% Instantaneous frequency
instantaneousFrequency = ...
    [NaN diff(instantaneousPhase)*fs/(2*pi)];

% Avoid filter-edge regions when calculating mean frequency
validRegion = t >= 1 & t <= duration-1;

hilbertPeak = mean( ...
    instantaneousFrequency(validRegion), ...
    'omitnan');


%% ---------------------------------------------------------
% 3. STFT ANALYSIS
% ---------------------------------------------------------

windowLength = round(0.5*fs);        % 500 ms
overlap = round(0.50*windowLength);  % 50% overlap
nfft = 256;

window = hann(windowLength);

[S,F,T] = spectrogram( ...
    x, ...
    window, ...
    overlap, ...
    nfft, ...
    fs);

STFTPower = abs(S).^2;

% Restrict analysis to 1-40 Hz
freqMask = F >= 1 & F <= 40;

Fstft = F(freqMask);
Pstft = STFTPower(freqMask,:);

% Mean power across time
meanSTFTPower = mean(Pstft,2);

[~,idx] = max(meanSTFTPower);

stftPeak = Fstft(idx);


%% ---------------------------------------------------------
% 4. MULTITAPER ANALYSIS - 1 SECOND WINDOW, FIVE TAPERS
% ---------------------------------------------------------

% Multitaper parameters
mtWindowLength = round(1.0*fs);      % 1 second = 500 samples
NW = 3;                              % Time-bandwidth product
K = 5;                               % Number of DPSS tapers
nfftMT = 512;

% Select a 1-second segment from the centre of the signal
startSample = round(2*fs) + 1;
endSample = startSample + mtWindowLength - 1;

xSegment = x(startSample:endSample);

% Generate five DPSS tapers
[tapers,~] = dpss(mtWindowLength,NW,K);

% Store individual spectra
mtSpectra = zeros(nfftMT/2+1,K);

for k = 1:K

    % Apply kth DPSS taper
    taperedSignal = xSegment .* tapers(:,k)';

    % Fourier transform
    X = fft(taperedSignal,nfftMT);

    % Power spectrum
    P = abs(X).^2;

    % Keep positive frequencies
    mtSpectra(:,k) = P(1:nfftMT/2+1);

end

% Average the five spectral estimates
Pmt = mean(mtSpectra,2);

% Frequency vector
Fmt = (0:nfftMT/2)' * fs/nfftMT;

% Restrict analysis to 1-40 Hz
mtMask = Fmt >= 1 & Fmt <= 40;

FmtPlot = Fmt(mtMask);
PmtPlot = Pmt(mtMask);

% Find maximum spectral value
[~,idxMT] = max(PmtPlot);

mtPeak = FmtPlot(idxMT);

%% ---------------------------------------------------------
% 5. DISPLAY FREQUENCY RESULTS
% ---------------------------------------------------------

fprintf('\n--------------------------------------------\n');
fprintf('10 Hz SYNTHETIC SIGNAL VALIDATION\n');
fprintf('--------------------------------------------\n');
fprintf('True frequency        : %.2f Hz\n',f0);
fprintf('Filter-Hilbert        : %.2f Hz\n',hilbertPeak);
fprintf('STFT                  : %.2f Hz\n',stftPeak);
fprintf('Multitaper            : %.2f Hz\n',mtPeak);
fprintf('--------------------------------------------\n');


%% ---------------------------------------------------------
% 6. CREATE FIGURE 10
% ---------------------------------------------------------

figure('Color','w', ...
       'Position',[100 100 1100 750]);

tiledlayout(2,2, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% (a) Synthetic 10 Hz signal
nexttile;

plot(t,x,'LineWidth',1);

xlabel('Time (s)');
ylabel('Amplitude');

title('(a) Stationary 10 Hz Synthetic Signal');

xlim([0 2]);
grid on;


%% (b) Filter-Hilbert result
nexttile;

yyaxis left

plot(t,filteredSignal,'LineWidth',1);

ylabel('Filtered Amplitude');

yyaxis right

plot(t,instantaneousFrequency,'LineWidth',1);

ylabel('Instantaneous Frequency (Hz)');

xlabel('Time (s)');

title(sprintf( ...
    '(b) Filter-Hilbert: Estimated Frequency = %.2f Hz', ...
    hilbertPeak));

xlim([0 duration]);
ylim([0 20]);

grid on;


%% (c) STFT result
nexttile;

imagesc(T,Fstft,10*log10(Pstft + eps));

axis xy;

xlabel('Time (s)');
ylabel('Frequency (Hz)');

title(sprintf( ...
    '(c) STFT: Peak Frequency = %.2f Hz', ...
    stftPeak));

ylim([1 30]);

colorbar;


%% (d) Multitaper result
nexttile;

plot(FmtPlot, ...
     10*log10(PmtPlot + eps), ...
     'LineWidth',1.3);

xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');

title(sprintf( ...
    '(d) Multitaper: Peak Frequency = %.2f Hz', ...
    mtPeak));

xlim([1 30]);

grid on;


%% ---------------------------------------------------------
% 7. EXPORT HIGH-RESOLUTION FIGURE
% ---------------------------------------------------------

exportgraphics( ...
    gcf, ...
    'Fig10_10Hz_Validation.png', ...
    'Resolution',300);

exportgraphics( ...
    gcf, ...
    'Fig10_10Hz_Validation.pdf', ...
    'ContentType','vector');

disp('Fig. 10 has been exported successfully.');