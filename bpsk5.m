clear; clc; close all;

%% Parameters
Fc = 10e6;       % Carrier frequency (for passband sim)
Fs = 80e6;       % Sampling rate (>= 2*Fc, oversample some)
N  = 1e5;        % Number of samples to simulate
SNR_dB = 10;     % SNR for the channel

% BPSK symbol rate, etc.
symbolRate = 1e6;
Ts = 1  / Fs;


%% 1) Generate BPSK signal with Communications Toolbox
% 1. Create random binary data
dataBits = randi([0 1], 1000, 1);

% 2. Create a BPSK modulator System object
bpskMod = comm.BPSKModulator;

% 3. Modulate the bits
bpskSignal = bpskMod(dataBits);
% Optional: Add noise (requires Communications Toolbox)
bpskSignalNoisy = awgn(bpskSignal, 20, 'measured');

% Example: view constellation
scatterplot(bpskSignalNoisy)
title('BPSK Constellation with AWGN')

%% 1) Generate a passband BPSK-like signal at Fs
t = (0:N-1)'/Fs;

% For simplicity, let's create random +/- 1 baseband samples
% at a symbol rate of 1e6, upsampled to Fs=80e6 => 80 samples/symbol
sps = Fs / symbolRate;
dataBits = randi([0 1], N/sps, 1);
sym = 2*dataBits - 1;  % Map to +/-1


% Upsample
symUpsampled = upsample(sym, sps);

% Mix up to carrier
txRF = real( txBB .* exp(1j*2*pi*Fc*t) );

%% 2) Add noise
rxRF_ideal = awgn(txRF, SNR_dB, 'measured'); % "Infinite-bit" approach, so no quant. yet

% Also create a 1-bit version
rxRF_1bit = sign(rxRF_ideal);  % Simply take sign


%% 3) Sample the RF signals

step = Ts / Tct;

sample_ideal = rxRF_ideal(1:step:end);

sample_1bit = rxRF_1bit(1:step:end);


%A simple way to detect the moment of flip

% Lock onto the wave’s fundamental period (if it is stable). For each cycle (or half-cycle),
% sample the sign (HIGH/LOW). If suddenly in the next cycle that same portion of the wave is
% the opposite sign, that implies a 180° phase flip occurred in between. You can time‐stamp
% that flip by noticing when the output logic changes out of sync with the usual transitions.
% Example using a D flip-flop or XOR
% XOR method: If you have a reference clock synchronized to the same frequency, an XOR gate
% with the input wave can reveal phase shifts. A sudden 180° shift will make the XOR output
% change state abruptly.
% Delay & compare: Another approach is to delay (by one period) the digital wave and compare it
% with the current wave. If they differ consistently where they used to match, a phase flip occurred.


% Go through each symbol, sample at a convenient location on the symbol
% boundary and determine the sign.

