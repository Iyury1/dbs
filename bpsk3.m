% MATLAB Script: Compare BPSK Demod (Direct RF Sampling vs Zero IF)

clear; close all; clc;

%% Parameters
Fs_base = 1e6;       % Baseband sampling rate (for zero-IF path)
Fc = 10e6;           % RF carrier (for direct sampling simulation)
oversampling = 8;    % Over-sampling factor for the RF path
Fs_rf = Fc * 2;      % Example: Nyquist rate ~ 2*Fc (just for demonstration)
% (You could do Fc * oversampling if you want a certain oversampling ratio at RF)

symbolRate = 100e3;  % Symbol rate for BPSK
M = 2;               % BPSK => M = 2
numSymbols = 1e3;    % Number of symbols

% Derived parameters
Tsym = 1/symbolRate;
T_base = 1/Fs_base;
t_base = (0 : numSymbols*Fs_base/symbolRate-1)*T_base;

%% 1) Generate Random Bits and BPSK Baseband Signal
dataIn = randi([0 M-1], numSymbols, 1);
% Map bits to +1/-1
bpskSymbols = 2*dataIn - 1;  % 0 -> -1, 1 -> +1 (or vice versa)

% Upsample BPSK to baseband sampling rate
sps_base = Fs_base / symbolRate;  % samples per symbol in baseband
txBaseband = upsample(bpskSymbols, sps_base);

% Apply a simple pulse shaping filter (rectangular for demonstration)
txBasebandFilt = filter(ones(sps_base,1), 1, txBaseband);

% Truncate transient
txBasebandFilt = txBasebandFilt( (length(ones(sps_base,1))-1)/2 + 1 : end );

% Now we have a baseband signal at Fs_base
t_base = (0 : length(txBasebandFilt)-1)*T_base;

%% 2) Zero-IF Simulation (Direct Baseband)
% Add some AWGN
SNR_dB = 10;
rxBaseband = awgn(txBasebandFilt, SNR_dB, 'measured');

% Demod for zero-IF
% For BPSK in an ideal scenario, we can just look at the sign of the real part
% Downsample and detect
rxBaseband_ds = rxBaseband(1:sps_base:end);  % downsample back to symbol rate
dataDetected_baseband = rxBaseband_ds > 0;

% Convert back to bits  0-> +1  => dataIn=1 is +1? We used 2*bits-1 => +1 => bits=1
% So sign>0 => 1, sign<0 => 0
numErrorsBase = sum(dataDetected_baseband ~= dataIn);
BER_base = numErrorsBase / numSymbols;
fprintf('Zero-IF Approach: BER = %g\n', BER_base);

%% 3) Direct RF Sampling Simulation
% We have the same baseband signal. Now upconvert it to RF:
t_rf = (0 : (length(txBasebandFilt)-1))' / Fs_base;  % time vector at baseband sample times
% But we want to produce a signal sampled at Fs_rf (which might be different from Fs_base).

% Let's re-sample the baseband signal to a new time vector that fits Fs_rf
t_rf_new = (0 : 1 : (ceil(length(txBasebandFilt)*Fs_rf/Fs_base)-1))'/Fs_rf;
% Interpolate baseband signal to match Fs_rf sampling rate
txBasebandInterp = interp1(t_rf, txBasebandFilt, t_rf_new, 'linear', 0);

% Now upconvert to carrier Fc
txRF = real(txBasebandInterp .* exp(1j*2*pi*Fc*t_rf_new)); 
% BPSK is real, but let’s do a full complex multiplication and keep real part.

% Add noise at RF
rxRF = awgn(txRF, SNR_dB, 'measured');

% ---- Digital Downconversion ----
% Multiply by e^{-j2pi Fc t}, then low-pass filter
rxRF_complex = rxRF .* exp(-1j*2*pi*Fc*t_rf_new);

% Design or define a low-pass filter ~ symbolRate/2 (just for demonstration)
lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 0.1, ...
    'StopbandFrequency', 0.2, 'PassbandRipple', 1, ...
    'StopbandAttenuation', 60, 'SampleRate', Fs_rf);

rxBaseband_dig = filter(lpFilt, rxRF_complex);

% Now downsample to baseband rate
decimFactor = round(Fs_rf / Fs_base);
rxBaseband_dig_ds = downsample(rxBaseband_dig, decimFactor);

% Time-align or remove filter delay if necessary
grpDelay = round(mean(grpdelay(lpFilt)));
rxBaseband_dig_ds = rxBaseband_dig_ds(grpDelay+1:end);

% Now we have a complex baseband. For ideal BPSK, the imaginary part should be near zero.
rxBaseband_dig_ds = real(rxBaseband_dig_ds);

% Downsample again to symbol rate if needed
rxBaseband_symbols = downsample(rxBaseband_dig_ds, sps_base);

% Make a sign decision
dataDetected_rf = rxBaseband_symbols > 0;
numErrorsRF = sum(dataDetected_rf(1:numSymbols) ~= dataIn); % match lengths carefully
BER_rf = numErrorsRF / numSymbols;
fprintf('Direct RF Sampling Approach: BER = %g\n', BER_rf);

%% Plot Comparison
figure;
subplot(2,1,1);
plot(real(rxBaseband(1:200))); grid on;
title('Zero-IF Baseband (I component, first 200 samples)');
xlabel('Sample index'); ylabel('Amplitude');

subplot(2,1,2);
plot(real(rxBaseband_dig(1:200))); grid on;
title('Direct RF Sampling, Digitally Downconverted (I component, first 200 samples)');
xlabel('Sample index'); ylabel('Amplitude');