clear; clc; close all;

%% Parameters
fc        = 24e6;       % "Carrier frequency" (use 24 GHz if you really want, but be cautious of huge data)
Rb        = 1e6;        % Bit rate
oversamp  = 10;         % Oversampling factor (samples per bit)
Fs        = Rb * oversamp; % Sampling frequency

numBits   = 1000;       % Number of bits
t_bit     = 1/Rb;       % Bit duration
t_total   = numBits * t_bit;

% Time vector for the entire signal
t = linspace(0, t_total, numBits*oversamp);

%% Generate random bits and map to phases
bits       = randi([0 1], numBits, 1);
symbols    = 2*bits - 1;     % Map 0->-1, 1->+1

% Upsample the symbols to match oversamp
symbols_upsampled = repmat(symbols.', oversamp, 1); 
symbols_upsampled = symbols_upsampled(:).';  % row vector

%% Generate carrier
% BPSK: phase is either 0 or pi, so multiply carrier by +1 or -1
carrier = cos(2*pi*fc*t);

% Create passband BPSK signal
bpsk_passband = symbols_upsampled .* carrier;

%% (Optional) Add noise
% Let's define an SNR
SNR_dB = 10;
rx_noisy = awgn(bpsk_passband, SNR_dB, 'measured');

%% Downconversion via Mixer
% In a hardware receiver, you'd have a local oscillator (LO) at ~fc
LO = cos(2*pi*fc*t);

% Multiply received signal by LO (assume perfect phase lock for simplicity)
mixed = rx_noisy .* LO;

% Low-pass filter to get rid of the high-frequency component
% Design a simple LPF with passband ~ 1.2 * (Rb/2)
lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 0.6*(Rb/Fs), ...
    'StopbandFrequency', 0.8*(Rb/Fs), 'PassbandRipple', 1, ...
    'StopbandAttenuation', 60, 'SampleRate', Fs);

baseband = filter(lpFilt, mixed);

%% BPSK Demodulation
% We assume the baseband signal is roughly +/- amplitude for 1/0
% We'll sample at symbol boundaries.

% Reshape baseband to figure out the symbol decisions
baseband_reshaped = reshape(baseband, oversamp, numBits);

% Integrate (or just look at the sign at midpoint)
decisions = mean(baseband_reshaped, 1);  % integrate each symbol
demod_bits = decisions > 0;             % 1 if positive, 0 if negative

%% Compute Bit Error Rate
numErrors = sum(demod_bits(:) ~= bits);
BER = numErrors / numBits;

fprintf('BER = %g at SNR = %g dB\n', BER, SNR_dB);

% Plot a small snippet of waveforms (for visualization)
figure;
subplot(3,1,1);
plot(t(1:500), bpsk_passband(1:500));
title('Transmitted Passband BPSK (first 500 samples)');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,2);
plot(t(1:500), rx_noisy(1:500));
title('Received + Noise');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(3,1,3);
plot(t(1:500), baseband(1:500));
title('Downconverted Baseband (filtered)');
xlabel('Time (s)'); ylabel('Amplitude');