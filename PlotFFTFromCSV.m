%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Using READMATRIX (Recommended for MATLAB R2019b+)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

filename = 'C:\Users\ianyu\OneDrive\MScEng\Chip Testing\iladata.csv';

% (A) Detect import options and set data lines to start from row 3
opts = detectImportOptions(filename);
% The 3rd row in typical (human) counting is row #3.
% DataLines = [3 Inf] means start reading from line 3 to the end of file.
opts.DataLines = [3 Inf];  

% (B) Read the CSV as a matrix according to these options
dataAll = readmatrix(filename, opts);

% (C) Extract only the 4th column
data = dataAll(:,4);

% (D) Define sampling frequency (adjust to your actual Fs)
Fs = 1e9;

% (E) Compute the FFT and power spectrum
N = length(data);
Y = fft(data);
P2 = abs(Y/N);
P1 = P2(1:N/2+1);
powerSpectrum = P1.^2;

% (F) Frequency axis
f = Fs*(0:(N/2))/N;

% (G) Plot the power spectrum
figure;
plot(f, powerSpectrum, 'LineWidth', 1.2);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum');
grid on;


% (9) (Optional) Plot power in decibels
% powdB = 10*log10(powerSpectrum);
% figure;
% plot(f, powdB, 'r', 'LineWidth', 1.5);
% grid on;
% title('Power Spectrum (dB) of Imported Signal');
% xlabel('Frequency (Hz)');
% ylabel('Power (dB)');