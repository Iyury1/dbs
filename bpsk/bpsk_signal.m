function bpsk_signal = generateBPSKBaseband(M, Bd)
    % generateBPSKBaseband generates a BPSK baseband signal.
    %
    % Inputs:
    %   M    - Number of symbols
    %   Fsym - Symbol rate (in Hz)
    %
    % Output:
    %   bpsk_signal - The generated BPSK baseband signal (vector)
    
    % Oversampling factor (number of samples per symbol)
    oversampleFactor = 10;  
    % Sampling frequency based on the symbol rate and oversampling factor
    Fct = oversampleFactor * Bd;
    
    % Generate M random bits (0 or 1)
    bits = randi([0,1], M, 1);
    
    % Map bits to BPSK symbols: 0 -> -1, 1 -> +1
    symbols = 2 * bits - 1;

    % Create the baseband signal by repeating each symbol
    % 'kron' replicates each symbol oversampleFactor times.
    baseband = kron(symbols, ones(oversampleFactor, 1));

    % Create a time vector for the entire signal (in seconds)
    t = (0:length(baseband)-1) / fs;
    
    % Generate the carrier signal: a cosine wave
    carrier = cos(2*pi*fc*t)';
    
    % Modulate the carrier with the BPSK baseband signal
    m = baseband .* carrier;

    Fs = 1e9;
    Ts = 1 / Fs;

    s = m(0:Ts:Bd )