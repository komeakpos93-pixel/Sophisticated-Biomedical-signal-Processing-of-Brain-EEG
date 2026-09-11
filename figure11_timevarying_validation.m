%% FIGURE 11 - TIME-VARYING SYNTHETIC SIGNAL VALIDATION
% 6 Hz -> 10 Hz -> 20 Hz
% STFT and Multitaper comparison

clear;
clc;
close all;

%% ---------------------------------------------------------
% 1. CREATE TIME-VARYING SYNTHETIC SIGNAL
% ---------------------------------------------------------

fs = 500;                 % Sampling frequency (Hz)
duration = 6;             % Total duration (s)
t = 0:1/fs:duration-1/fs;

rng(2);                   % Reproducible noise
noiseLevel = 0.4;

x = zeros(size(t));

% 0-2 s: 6 Hz
idx1 = t < 2;
x(idx1) = sin(2*pi*6*t(idx1));

% 2-4 s: 10 Hz
idx2 = t >= 2 & t < 4;
x(idx2) = sin(2*pi*10*t(idx2));

% 4-6 s: 20 Hz
idx3 = t >= 4;
x(idx3) = sin(2*pi*20*t(idx3));

% Add noise
x = x + noiseLevel*randn(size(x));


%% ---------------------------------------------------------
% 2. STFT ANALYSIS
% ---------------------------------------------------------

stftWindowLength = round(0.5*fs);       % 500 ms
stftOverlap = round(0.50*stftWindowLength);
stftNFFT = 512;

stftWindow = hann(stftWindowLength);

[S,Fstft,Tstft] = spectrogram( ...
    x, ...
    stftWindow, ...
    stftOverlap, ...
    stftNFFT, ...
    fs);

STFTPower = abs(S).^2;

% Limit to 1-30 Hz
stftMask = Fstft >= 1 & Fstft <= 30;

FstftPlot = Fstft(stftMask);
PstftPlot = STFTPower(stftMask,:);


%% ---------------------------------------------------------
% 3. MULTITAPER TIME-FREQUENCY ANALYSIS
% ---------------------------------------------------------

mtWindowLength = round(0.5*fs);      % 500 ms
mtOverlap = round(0.50*mtWindowLength);
mtStep = mtWindowLength - mtOverlap;

NW = 2;
K = 3;
mtNFFT = 512;

[tapers,~] = dpss(mtWindowLength,NW,K);

% Number of windows
numWindows = floor((length(x)-mtWindowLength)/mtStep) + 1;

% Frequency vector
Fmt = (0:mtNFFT/2)' * fs/mtNFFT;

% Restrict to 1-30 Hz
mtMask = Fmt >= 1 & Fmt <= 30;
FmtPlot = Fmt(mtMask);

MTpower = zeros(sum(mtMask),numWindows);
Tmt = zeros(1,numWindows);

for w = 1:numWindows

    startIdx = (w-1)*mtStep + 1;
    endIdx = startIdx + mtWindowLength - 1;

    segment = x(startIdx:endIdx);

    taperSpectra = zeros(mtNFFT/2+1,K);

    for k = 1:K

        taperedSegment = segment .* tapers(:,k)';

        X = fft(taperedSegment,mtNFFT);

        P = abs(X).^2;

        taperSpectra(:,k) = P(1:mtNFFT/2+1);

    end

    meanSpectrum = mean(taperSpectra,2);

    MTpower(:,w) = meanSpectrum(mtMask);

    Tmt(w) = ((startIdx-1) + mtWindowLength/2)/fs;

end


%% ---------------------------------------------------------
% 4. CREATE FIGURE 11
% ---------------------------------------------------------

figure('Color','w', ...
       'Position',[100 100 1100 850]);

tiledlayout(3,1, ...
    'TileSpacing','compact', ...
    'Padding','compact');


%% (a) Time-domain signal
nexttile;

plot(t,x,'LineWidth',1);

xlabel('Time (s)');
ylabel('Amplitude');

title('(a) Time-Varying Synthetic Signal');

xlim([0 6]);

xline(2,'--');
xline(4,'--');

text(0.9,1.8,'6 Hz');
text(2.9,1.8,'10 Hz');
text(4.9,1.8,'20 Hz');

grid on;


%% (b) STFT spectrogram
nexttile;

imagesc( ...
    Tstft, ...
    FstftPlot, ...
    10*log10(PstftPlot + eps));

axis xy;

xlabel('Time (s)');
ylabel('Frequency (Hz)');

title('(b) STFT Time-Frequency Representation');

ylim([1 30]);

xline(2,'--w');
xline(4,'--w');

colorbar;


%% (c) Multitaper spectrogram
nexttile;

imagesc( ...
    Tmt, ...
    FmtPlot, ...
    10*log10(MTpower + eps));

axis xy;

xlabel('Time (s)');
ylabel('Frequency (Hz)');

title('(c) Multitaper Time-Frequency Representation');

ylim([1 30]);

xline(2,'--w');
xline(4,'--w');

colorbar;


%% ---------------------------------------------------------
% 5. EXPORT
% ---------------------------------------------------------

exportgraphics( ...
    gcf, ...
    'Fig11_TimeVarying_Validation.png', ...
    'Resolution',300);

exportgraphics( ...
    gcf, ...
    'Fig11_TimeVarying_Validation.pdf', ...
    'ContentType','vector');

disp('Fig. 11 has been exported successfully.');