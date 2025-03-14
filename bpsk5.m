clear; clc; close all;

%% Parameters
Fc = 24e9;       % Carrier frequency (for passband sim)
Fs = 1e9;       % Sampling rate (>= 2*Fc, oversample some)
upsample = 4;
Fct = Fs * upsample;      % Continuous time sampling rate

N  = 1e3;        % Number of samples to simulate
SNR_dB = 10;     % SNR for the channel

% BPSK symbol rate, etc.
Fsym = 1e4;
Ts = 1  / Fs;
Tct = 1 / Fct;


%% 1) Generate a passband BPSK-like signal at Fs
t1 = (0:N-1)' * Ts; % N samples Ts seconds apart
t2 = linspace(0, t1(end), N * upsample)'; 
%t2 = (0:(N*upsample)-1)' * Tct; % N samples Ts seconds apart

disp('Length of t1:'); disp(length(t1));
disp('Length of t2:'); disp(length(t2));
disp('Max time in t1:'); disp(t1(end));
disp('Max time in t2:'); disp(t2(end));

% For simplicity, let's create random +/- 1 baseband samples
% at a symbol rate of 1e6, upsampled to Fs=80e6 => 80 samples/symbol
sps_ct = Fct / Fsym; % CT samples per symbol
M = (N*upsample)/sps_ct; % num symbols
dataBits = randi([0 1], M, 1);
sym = 2*dataBits - 1;  % Map to +/-1

% Upsample
symUpsampled = repelem(sym, sps_ct);

disp('Length of symUpsampled (should match length of t2):');
disp(length(symUpsampled));

% Mix up to carrier
txRF = real( symUpsampled .* exp(1j*2*pi*Fc*t2) );

%% 2) Add noise
rxRF_ideal = awgn(txRF, SNR_dB, 'measured'); % "Infinite-bit" approach, so no quant. yet
%rxRF_ideal = txRF;
% Also create a 1-bit version
rxRF_1bit = sign(rxRF_ideal);  % Simply take sign
figure;

subplot(2,2,1); plot(rxRF_ideal); title('Plot 1');
subplot(2,2,2); plot(rxRF_1bit); title('Plot 1');
disp(length(rxRF_ideal));
%% 3) Sample the RF signals

step = Ts / Tct;
disp(['Step ratio (Ts/Tct): ' num2str(step)]);

sample_ideal = rxRF_ideal(1:step:end);

sample_1bit = rxRF_1bit(1:step:end);
subplot(2,2,1); plot(sample_ideal); title('Infinite Bit ADC');
subplot(2,2,2); plot(sample_1bit); title('1 Bit ADC');

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

offset = floor(sps_ct/2);  % e.g., 2000 if sps_ct=4000

% Indices at which we sample each of the M symbols
symbolIndices = offset + (0:M-1)*sps_ct;  % 0-based indexing

% Make sure we don't exceed array bounds:
if symbolIndices(end) > length(rxRF_ideal)
    error('Symbol sampling index exceeds signal length!');
end

% Extract symbol values from both signals
% For "ideal" we take the sign to make a hard decision
rxSym_ideal = sign(rxRF_ideal(symbolIndices));

% For "1-bit" data, it's already ±1, so just pick
rxSym_1bit = rxRF_1bit(symbolIndices);

firstBit = dataBits(1)
firstSign_ideal = sign(rxSym_ideal(1))
firstSign_1bit = sign(rxSym_1bit(1))

secondBit = dataBits(2)
secondSign_ideal = sign(rxSym_ideal(2))
secondSign_1bit = sign(rxSym_1bit(2))


% Convert received symbols back to bits:
%   +1 => bit=1
%   -1 => bit=0
rxBits_ideal = (rxSym_ideal + 1)/2;
rxBits_1bit  = (rxSym_1bit  + 1)/2;

% Compare with the original bits
numErrors_ideal = sum(rxBits_ideal ~= dataBits);
numErrors_1bit  = sum(rxBits_1bit  ~= dataBits);


BER_ideal = numErrors_ideal / M;
BER_1bit  = numErrors_1bit  / M;

fprintf('\n=== BER RESULTS ===\n');
fprintf('BER (ideal) : %g\n', BER_ideal);
fprintf('BER (1-bit) : %g\n', BER_1bit);
