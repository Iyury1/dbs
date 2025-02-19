clear; clc; close all;

%% Parameters
Fc = 10e6;       % Carrier frequency (for passband sim)
Fs = 80e6;       % Sampling rate (>= 2*Fc, oversample some)
N  = 1e5;        % Number of samples to simulate
SNR_dB = 10;     % SNR for the channel

% BPSK symbol rate, etc. (Optional for more advanced sim)
symbolRate = 1e6;

%% 1) Generate a passband BPSK-like signal at Fs
t = (0:N-1)'/Fs;

% For simplicity, let's create random +/- 1 baseband samples
% at a symbol rate of 1e6, upsampled to Fs=80e6 => 80 samples/symbol
sps = Fs / symbolRate;
dataBits = randi([0 1], N/sps, 1);
sym = 2*dataBits - 1;  % Map to +/-1

% Upsample
symUpsampled = upsample(sym, sps);

% Pulse shape (rectangular for simplicity)
txBB = filter(ones(sps,1), 1, symUpsampled);

% Mix up to carrier
txRF = real( txBB .* exp(1j*2*pi*Fc*t) );

%% 2) Add noise
rxRF_ideal = awgn(txRF, SNR_dB, 'measured'); % "Infinite-bit" approach, so no quant. yet

% Also create a 1-bit version
rxRF_1bit = sign(rxRF_ideal);  % Simply take sign

%% 3) Digital Downconversion

% ----- (A) Infinite-bit path -----
% Multiply by e^{-j2*pi*Fc*t}, lowpass filter to baseband
rxMixed_ideal = rxRF_ideal .* exp(-1j*2*pi*Fc*t);

Fpass = 370;
Fstop = 430;
Ap = 1;
Ast = 30;
Fs = 2000;
d = designfilt("lowpassfir", ...
    PassbandFrequency=Fpass, ...
    StopbandFrequency=Fstop, ...
    PassbandRipple=Ap, ...
    StopbandAttenuation=Ast, ...
    SampleRate=Fs);

% Example: simple lowpass FIR
lpFilt = designfilt('lowpassfir', 'PassbandFrequency', 0.1, ...
    'StopbandFrequency', 0.2, 'PassbandRipple', 1, ...
    'StopbandAttenuation', 60, 'SampleRate', Fs);
rxBB_ideal = filter(lpFilt, rxMixed_ideal);

%%  ---- (B) 1-bit path -----
% For a "true" 1-bit approach, we need to decide how to do mixing:
%   sign(a*b) = sign(a)*sign(b)  (for real signals)
% But here, we have complex exponent. Let's do a simpler approach:
%   multiply in float, then sign again.
% A more "pure" 1-bit sim would do sign-multiplies, etc.

loComplex = exp(-1j*2*pi*Fc*t);
rxMixed_temp = rxRF_1bit .* loComplex;  % This is effectively float x sign
rxMixed_1bit = sign(real(rxMixed_temp)) + 1j*sign(imag(rxMixed_temp));

% Lowpass filter in 1-bit domain is tricky; 
% let's approximate by letting the filter be float, then re-quantize:
rxBB_1bit_float = filter(lpFilt, rxMixed_1bit);
rxBB_1bit = rxBB_1bit_float;  % or re-quantize: rxBB_1bit = sign(rxBB_1bit_float);

%% 4) Downsample to symbol rate and detect

% For the ideal path:
decimFactor = sps;  % we used 'sps' upsample
rxBB_ideal_ds = downsample(rxBB_ideal, decimFactor);
% We only look at real part for BPSK:
rxSym_ideal = real(rxBB_ideal_ds);

detBits_ideal = rxSym_ideal > 0;
dataBitsShort = dataBits(1:length(detBits_ideal)); % align lengths
BER_ideal = mean(detBits_ideal ~= dataBitsShort);

% For the 1-bit path:
rxBB_1bit_ds = downsample(rxBB_1bit, decimFactor);
rxSym_1bit = real(rxBB_1bit_ds);

detBits_1bit = rxSym_1bit > 0;
BER_1bit = mean(detBits_1bit ~= dataBitsShort);

fprintf('BER (Infinite-bit) = %g\n', BER_ideal);
fprintf('BER (1-bit)        = %g\n', BER_1bit);