% If we are doing iterative runs, this can be helpful to save preprocessing time
clearvars -except chirpZeroCrossings DBSSampleTimes DBSSampleTimesWithEDError DBSSampeTimePhases
% Define system contants
frameRate = 60;
chirpTime = 1/frameRate;
bandwidth = 1e9;
chirpRate = bandwidth/chirpTime;
global sampleFrequency
sampleFrequency = 300e6;
c = 3e8;
f0 = 76e9;
f1 = f0 + bandwidth;
freq = (f1-f0)/2+f0;
centerFreq = (f1-f0)/2+f0;
quarterCenterFreqTime = (1/centerFreq)/4;
lambda = c/freq;
k = (2*pi)/lambda;
pitch = lambda/2;
noisePower = -35; %dBm
%% Preprocessing calculations to determine the DBS Latch Sample Times
% First, determine the rising edge times
if (~exist('chirpZeroCrossings','var'))
fprintf("Starting Chirp Zero Crossing Calculations\n");
chirpZeroCrossings = zeros(1, round(chirpTime/(1/freq)));
for kthCrossing = 2:2:2*length(chirpZeroCrossings)
chirpZeroCrossings(kthCrossing/2) = calculateChirpZeroCrossing(f0, f1, chirpTime, kthCrossing);
end
end
% Then calculate the DBS Latch Sample Times
if (~exist('DBSSampleTimes', 'var'))
fprintf("Starting DBS Latch Sample Time Calculations\n");
DBSSampleTimes = calculateDBSSampleTimes(chirpZeroCrossings, sampleFrequency);
end
clear chirpZeroCrossings; % We can clear the chirp zero crossings to free memory as we don't need them anymore
%% Setup Antennas and Reflectors
global arrayRows;
global arrayCols;
arrayRows = 64;
arrayCols = 64;
%sr1 = SuperReflector(50, 50, 0, 0, 30, 1);
superReflectors = [];
% To create a list of super reflectors: superReflectors = [sr1];, but we start with it empty
% Reflecting objects to simulate
r1 = Reflector(25, 25, 0, 1);
r2 = Reflector(100, -100, 25, 1);
r3 = Reflector(45, 18, 80, 1);
global reflectors
reflectors = [r1 r2 r3];
for superReflector = superReflectors
    reflectors = [reflectors, superReflector.ReflectorArray];
end
reflectorCount = 0;
for reflector = reflectors
    reflectorCount = reflectorCount + 1;
    fprintf("Reflector %d: Theta = %0.2f, Phi = %0.2f, Distance = %0.2f\n", reflectorCount, reflector.theta, reflector.phi, reflector.r);
end
fprintf("Number of super reflectors: %d\n", length(superReflectors))
fprintf("Number of reflectors: %d\n", reflectorCount)
%%
fprintf("Starting Antenna Simulation\n");
close all;
fractionOfChirp = 0.25;
global Nsamples;
Nsamples = round(size(DBSSampleTimes, 2)*fractionOfChirp); % To reduce runtime and memory useage, only compute a fraction of the chirp
randPhases = rand(arrayRows, arrayCols)*360;
antennasReal = zeros(arrayRows, arrayCols, Nsamples, 'single');
antennasImag = zeros(arrayRows, arrayCols, Nsamples, 'single');
antennasRealNoisy = zeros(arrayRows, arrayCols, Nsamples, 'single');
antennasImagNoisy = zeros(arrayRows, arrayCols, Nsamples, 'single');
Niterations = arrayRows * arrayCols * size(reflectors,2);
curIteration = 0;
% Iterate over all antenna elements, and for each sample time, determine the contribution of each reflector
for row = 1 : arrayRows
    for col = 1 : arrayCols
        for reflector = 1 : size(reflectors,2)
            % Calculate the inphase and quadrature phase of the return signal. The phase depends on the location of the antenna in the grid
            realPhase = -(k*pitch*((row-1)*sind(reflectors(reflector).phi) + (col-1)*sind(reflectors(reflector).theta)))*180/pi;
            imagPhase = realPhase + 90;
            % Determine the index of time that the signal will return to the antenna grid, based on its round trip time
            reflectorDBSSampleTimeIdx = find(DBSSampleTimes > reflectors(reflector).roundTripTime, 1);
            % Create the time vector for the return chirp. We need this because the sample times are not purely periodic, but change with the instantaneous phase of the Tx signal
            t = DBSSampleTimes(reflectorDBSSampleTimeIdx:Nsamples)'-reflectors(reflector).roundTripTime;
            % Create the return signal for inphase and quadrature samplers
            reflectorChirpReal = chirp(t, f0, chirpTime, f1, 'linear', realPhase);
            reflectorChirpImag = chirp(t, f0, chirpTime, f1, 'linear', imagPhase);
            % Add return chirp to antenna. Multiple reflecting objects will have their contributions sum at the antenna
            antennasReal(row, col , reflectorDBSSampleTimeIdx:end) = squeeze(antennasReal(row, col ,reflectorDBSSampleTimeIdx:end))' + reflectorChirpReal';
            antennasImag(row, col , reflectorDBSSampleTimeIdx:end) = squeeze(antennasImag(row, col ,reflectorDBSSampleTimeIdx:end))' + reflectorChirpImag';
            % To simulate errors in sampling, pulled from schematic simulations of 65nm CML/StrongARM latch
            noiseReal = wgn(length(t), 1, noisePower);
            noiseImag = wgn(length(t), 1, noisePower);
            antennasRealNoisy(row, col, reflectorDBSSampleTimeIdx:end) = squeeze(antennasRealNoisy(row, col, reflectorDBSSampleTimeIdx:end))' + reflectorChirpReal' + noiseReal';
            antennasImagNoisy(row, col, reflectorDBSSampleTimeIdx:end) = squeeze(antennasImagNoisy(row, col, reflectorDBSSampleTimeIdx:end))' + reflectorChirpImag' + noiseImag';
            curIteration = curIteration + 1;
            fprintf("Progress: %0.2f%s\n", (curIteration/Niterations)*100, "%");
        end
    end
end
fprintf("Done!\n")
% Digitize the analog signals at the antennas
DBSReal = zeros(arrayRows, arrayCols, Nsamples,'single');
DBSImag = zeros(arrayRows, arrayCols, Nsamples,'single');
DBSReal(antennasReal >= 0) = 1;
DBSReal(antennasReal < 0) = -1;
DBSImag(antennasImag >= 0) = 1;
DBSImag(antennasImag < 0) = -1;
DBSRealNoisy = zeros(arrayRows, arrayCols, Nsamples,'single');
DBSImagNoisy = zeros(arrayRows, arrayCols, Nsamples,'single');
DBSRealNoisy(antennasRealNoisy >= 0) = 1;
DBSRealNoisy(antennasRealNoisy < 0) = -1;
DBSImagNoisy(antennasImagNoisy >= 0) = 1;
DBSImagNoisy(antennasImagNoisy < 0) = -1;
%% Accumulation
downsamplingRatio = 1220; % This downsampling ratio needs to be calculated based on the desired number of samples per chirp and sampling clock rate
global Naccumulations;
global accumulationProgress;
Naccumulations = floor(size(DBSReal, 3)/downsamplingRatio);
accumulatorInPhase = zeros(Naccumulations, arrayRows, arrayCols, 'single');
accumulatorQuadrature = zeros(Naccumulations, arrayRows, arrayCols, 'single');
accumulatorInPhaseNoisy = zeros(Naccumulations, arrayRows, arrayCols, 'single');
accumulatorQuadratureNoisy = zeros(Naccumulations, arrayRows, arrayCols, 'single');
accumulationSamplesOutCounter = 1;
index = 1;
fprintf("Starting Accumulation...\n");
for accumulationSamplesOutCounter = 1 : Naccumulations
    for accumulationCounter = 1 : downsamplingRatio
        index = downsamplingRatio*(accumulationSamplesOutCounter-1) + accumulationCounter;
        for row = 1 : arrayRows
            for col = 1 : arrayCols
                % Add the digital value from the Direct Binary Sampler to the accumulator register, until there are "downsamplingRatio" samples
                accumulatorInPhase(accumulationSamplesOutCounter, row, col) = accumulatorInPhase(accumulationSamplesOutCounter, row, col) + DBSReal(row, col, index);
                accumulatorQuadrature(accumulationSamplesOutCounter, row, col) = accumulatorQuadrature(accumulationSamplesOutCounter, row, col) + DBSImag(row, col, index);
                accumulatorInPhaseNoisy(accumulationSamplesOutCounter, row, col) = accumulatorInPhaseNoisy(accumulationSamplesOutCounter, row, col) + DBSRealNoisy(row, col, index);
                accumulatorQuadratureNoisy(accumulationSamplesOutCounter, row, col) = accumulatorQuadratureNoisy(accumulationSamplesOutCounter, row, col) + DBSImagNoisy(row, col, index);
            end
        end
    end
    fprintf("Accumulation Progress: %0.2f%s\n", ((accumulationSamplesOutCounter)/Naccumulations)*100,"%");
end
fprintf("Done Accumulation!\n");
% Normalize
for row = 1 : arrayRows
    for col = 1 : arrayCols
        % Shift to center around 0
        accumulatorInPhase(:,row,col) = accumulatorInPhase(:,row,col)-((downsamplingRatio)/2);
        accumulatorQuadrature(:,row,col) = accumulatorQuadrature(:,row,col)-((downsamplingRatio)/2);
        accumulatorInPhaseNoisy(:,row,col) = accumulatorInPhaseNoisy(:,row,col)-((downsamplingRatio)/2);
        accumulatorQuadratureNoisy(:,row,col) = accumulatorQuadratureNoisy(:,row,col)-((downsamplingRatio)/2);
        % Normalize
        accumulatorInPhase(:,row,col) = 0.5*(accumulatorInPhase(:,row,col)/max(accumulatorInPhase(:,row,col)))+0.5;
        accumulatorQuadrature(:,row,col) = 0.5*(accumulatorQuadrature(:,row,col)/max(accumulatorQuadrature(:, row,col)))+0.5;
        accumulatorInPhaseNoisy(:,row,col) = 0.5*(accumulatorInPhaseNoisy(:,row,col)/max(accumulatorInPhaseNoisy(:,row,col)))+0.5;
        accumulatorQuadratureNoisy(:,row,col) = 0.5*(accumulatorQuadratureNoisy(:,row,col)/max(accumulatorQuadratureNoisy(:,row,col)))+0.5;
    end
end
%% Beamforming and Distance Estimation
close all;
% FFT Based beamforming of accumulator output
% 2D FFT (make sure to include inphase and quadrature components). Take instantaneous value of all DBS (accumulator actually) values
beamformerAccumulatorFFT = zeros(arrayRows, arrayCols, size(accumulatorInPhase,1), 'single');
beamformerAccumulatorNoisyFFT = zeros(arrayRows, arrayCols, size(accumulatorInPhase,1), 'single');
for accCycle = 1 : Naccumulations
    beamformerAccumulatorFFT(:,:,accCycle) = fft2((squeeze(accumulatorInPhase(accCycle,:,:))) + (1i*(squeeze(accumulatorQuadrature(accCycle,:,:)))));
    beamformerAccumulatorNoisyFFT(:,:,accCycle) = fft2((squeeze(accumulatorInPhaseNoisy(accCycle,:,:))) + (1i*(squeeze(accumulatorQuadratureNoisy(accCycle,:,:)))));
end
% FFT based distance estimate from beamforming FFT
% For a single 2D index (specifical azimuth + elevation), take the values from all 2D beamforming FFTs and do a time domain FFT to determine object distance in that direction
distanceAccumulatorFFT = zeros(arrayRows, arrayRows, size(beamformerAccumulatorFFT,3)/2 + 1, 'single');
distanceAccumulatorNoisyFFT = zeros(arrayRows, arrayRows, size(beamformerAccumulatorFFT,3)/2 + 1, 'single');
for row = 1 : arrayRows
    for col = 1 : arrayCols
        [f, distanceAccumulatorFFT(row,col,:)] = calculateFFT(squeeze(beamformerAccumulatorFFT(row, col, :)), sampleFrequency/downsamplingRatio);
        [f, distanceAccumulatorNoisyFFT(row,col,:)] = calculateFFT(squeeze(beamformerAccumulatorNoisyFFT(row, col, :)), sampleFrequency/downsamplingRatio);
    end
end
%% DSP to extract distance from important beam angles
% For each distanceAccumulatorFFT, determine the "x" value for the max "y",
% the calculate the range from that.
for row = 1 : arrayRows
    for col = 1 : arrayRows
        [val, idx] = max(distanceAccumulatorFFT(row,col,:));
        maxFrequency = f(idx);
        R = (c*maxFrequency)/(2*chirpRate);
        distanceAccumulatorFFT_vals(row, col) = val;
        distanceAccumulatorFFT_ranges(row, col) = R;
        [val, idx] = max(distanceAccumulatorNoisyFFT(row,col,:));
        maxFrequency = f(idx);
        R = (c*maxFrequency)/(2*chirpRate);
        distanceAccumulatorFFT_Noisy vals(row, col) = val;
        distanceAccumulatorFFT_Noisy ranges(row, col) = R;
    end
end
% Average out the beamforming across all accumulator samples
beamformerAccumulatorFFTAverage = zeros(arrayRows, arrayRows);
beamformerAccumulatorFFTAverage_Noisy = zeros(arrayRows, arrayRows);
for row = 1 : arrayRows
    for col = 1 : arrayRows
        beamformerAccumulatorFFTAverage(row, col) = mean(abs(beamformerAccumulatorFFT(row, col, :)));
        beamformerAccumulatorFFTAverage_Noisy(row, col) = mean(abs(beamformerAccumulatorNoisyFFT(row,col, :)));
    end
end
beamformingFFTMax = max(squeeze(abs(beamformerAccumulatorFFTAverage(:,:))), [], "all");
beamformingFFTMin = min(squeeze(abs(beamformerAccumulatorFFTAverage(:,:))), [], "all");
beamformingFFTMid = (beamformingFFTMax-beamformingFFTMin)/2;
beamformingFFTMax_Noisy = max(squeeze(abs(beamformerAccumulatorFFTAverage_Noisy(:,:))), [], "all");
beamformingFFTMin_Noisy = min(squeeze(abs(beamformerAccumulatorFFTAverage_Noisy(:,:))), [], "all");
beamformingFFTMid_Noisy = (beamformingFFTMax_Noisy-beamformingFFTMin_Noisy)/2;
for row = 1 : arrayRows
for col = 1 : arrayRows
if abs(beamformerAccumulatorFFTAverage(row,col)) < beamformingFFTMid
distanceAccumulatorFFT_ranges(row, col) = 0;
end
if abs(beamformerAccumulatorFFTAverage_Noisy(row,col)) < beamformingFFTMid_Noisy
distanceAccumulatorFFT_Noisy ranges(row, col) = 0;
end
end
end
%% Save Beamforming and Range Data
fprintf("Starting beamforming and range file save...\n");
filename = sprintf("%d_Rows_%d_Cols_%0.2f_MHz_Sampling_%d_Reflectors_%d_NoisePower_BEAMFORM_AND_DISTANCE.mat", arrayRows, arrayCols, sampleFrequency/1e6, length(reflectors), noisePower);
save(filename, ...
    "beamformerAccumulatorFFTAverage",...
    "distanceAccumulatorFFT ranges",...
    "beamformerAccumulatorFFTAverage Noisy",...
    "distanceAccumulatorFFT Noisy ranges",...
    "arrayRows", ...
    "arrayCols", ...
    "reflectors",...
    "sampleFrequency", ...
    "noisePower", ...
    "-v7.3")
fprintf("Done beamforming and range file save!\n");
%% Plotting of Beamforming and Distance Information
close all;
plotBeamformingDistance(arrayRows, arrayCols, beamformerAccumulatorFFTAverage, distanceAccumulatorFFT_ranges, "Ideal");
plotBeamformingDistance(arrayRows, arrayCols, beamformerAccumulatorFFTAverage_Noisy, distanceAccumulatorFFT_Noisy_ranges, "Noisy");
plotBeamformingDistanceThreshold(arrayRows, arrayCols, distanceAccumulatorFFT_ranges, "Ideal");
plotBeamformingDistanceThreshold(arrayRows, arrayCols, distanceAccumulatorFFT_Noisy_ranges, "Noisy");
% Reminder of reflectors
reflectorCount = 0;
for reflector = reflectors
reflectorCount = reflectorCount + 1;
fprintf("Reflector %d: Theta = %0.2f, Phi = %0.2f, Distance = %0.2f\n", reflectorCount, reflector. theta, reflector.phi, reflector.r);
end
%% Functions
function ChirpZeroCrossing = calculateChirpZeroCrossing(f1, f2, T, k)
% Mathematically calculate the zero crossing times for a chirped signal
% f1 : starting chirp frequency [float]
% f2 : ending chirp frequency [float]
% T : period of chirp [float]
% k : Get the k'th zero crossing. We want this to be only even or odd to get only rising or falling edges
% returns:
% - ChirpZeroCrossing : vector of rising or falling edge times
a = (f2-f1)/(2*T);
b = f1;
c = -((k/2)-(1/4));
ChirpZeroCrossing = (-b + sqrt((b)^2 - (4*a*c)))/(2*a);
end
function DBSSampleTimes = calculateDBSSampleTimes(risingEdgeTimes, sampleFrequency)
% Calculates sample times for the sampler. This is needed to align/synchronize the phase of the clock
% to the phase of the transmit signal
% risingEdgeTimes : all transmit chirp zero crossing times [vector of float]
% sampleFrequency : sample clock frequency [float]
% returns:
% - DBSSampleTimes : Synchronized/aligned sample times for the sampler. These will be
% 1/sampleFrequency + time to next transmit rising edge
latchIdxs = [];
curSampleTime = 1/sampleFrequency;
for risingEdgeIdx = 1 : numel(risingEdgeTimes)
risingEdgeTime = risingEdgeTimes(risingEdgeIdx);
if(risingEdgeTime > curSampleTime)
latchIdxs(end+1) = risingEdgeTimes(risingEdgeIdx);
curSampleTime = curSampleTime + (1/sampleFrequency);
end
end
DBSSampleTimes = latchIdxs;
end
function [f, data_FFT] = calculateFFT(data, Fs)
% Applies a single dimension FFT
% data : data to apply FFT to [vector of float]
% Fs : sample frequency of "data" [float]
% returns:
% - f : vector of frequency bins
% - data_FFT : FFT data for each frequency bin
L = length(data);
data_FFT = abs(fft(data));
data_FFT = data_FFT(1:L/2+1);
data_FFT(2:end-1) = 2*data_FFT(2:end-1);
f = Fs*(0:(L/2))/L;
end
function plotBeamformingDistance(arrayRows, arrayCols, beamformingData, distanceData, titleName)
    % Creates and saves a 2-D plot of reflecting objects locations and their distances
    % arrayRows : number of antennas per row [integer]
    % arrayCols : number of antennas per column [integer]
    % beamformingData : preprocessed 2D array of beamforming data [2D x float]
    % distanceData : preprocessed 2D array of distance data [2D x float]
    % titleName : name of plot, also used in filename when saving plot
    % X-axis = Azimuth
    % Y-axis = Elevation
    % Depth = magnitude/strength of detection
    % Text overlaid on elegible directions for distance
    N = arrayRows;
    angles = asin((2 * (0:(N-1)) / N) -1);
    angles_deg = rad2deg(angles);
    finalFigure = figure('Position', get(0, 'Screensize'));
    BFD1 = flip(fftshift(abs(beamformingData), 1), 1);
    BFD2 = fftshift(abs(BFD1), 2);
    imagesc(BFD2);
    title(titleName)
    colorbar;
    colormap(parula);
    DD1 = flip(fftshift(distanceData', 1), 2);
    DD2 = fftshift(DD1, 2);
    for row = 1 : arrayRows
        for col = 1 : arrayRows
            if DD2(row, col) > 0
                text(row, col, num2str(DD2(row, col)), 'Color', 'black', 'HorizontalAlignment','center');
            end
        end
    end
    X = compose('%0.2f', angles_deg);
    Y = compose('%0.2f', flip(angles_deg));
    xticks(1:1:arrayRows);
    yticks(1:1:arrayRows);
    xticklabels(X);
    yticklabels(Y);
    xlabel("Azimuth (\theta)");
    ylabel("Elevation (\phi)");
    global arrayRows
    global arrayCols
    global sampleFrequency
    global reflectors
    filename = sprintf('%dx%d_Antennas_%0.2fMHz_SampFreq_%d_Reflectors_%s', arrayRows, arrayCols, sampleFrequency/1e6, size(reflectors, 2), titleName);
    savefig(finalFigure, filename + ".fig");
    saveas(finalFigure, filename + ".png");
end
function plotBeamformingDistanceThreshold(arrayRows, arrayCols, distanceData, titleName)
    % Creates and saves a 2-D plot of reflecting objects locations and their distances
    % The distance data is pre-processed to threshold shown data. Custom colormaps can be useful to show white noise background
    % arrayRows : number of antennas per row [integer]
    % arrayCols : number of antennas per column [integer]
    % distanceData : preprocessed 2D array of distance data [2D x float]
    % titleName : name of plot, also used in filename when saving plot
    % X-axis = Azimuth
    % Y-axis = Elevation
    % Depth = Object distance
    N = arrayRows;
    angles = asin((2 * (0:(N-1)) / N) -1);
    angles_deg = rad2deg(angles);
    finalFigure = figure('Position', get(0, 'Screensize'));
    DD1 = flip(fftshift(distanceData, 1), 1);
    DD2 = fftshift(DD1, 2);
    imagesc(DD2);
    title(titleName)
    colorbar;
    colormap(parula);
    X = compose('%0.2f', angles_deg);
    Y = compose('%0.2f', flip(angles_deg));
    xticks(1:1:arrayRows);
    yticks(1:1:arrayRows);
    xticklabels(X);
    yticklabels(Y);
    xlabel('Azimuth (\theta)');
    ylabel('Elevation (\phi)');
    global arrayRows
    global arrayCols
    global sampleFrequency
    global reflectors
    filename = sprintf('%dx%d_Antennas_%0.2fMHz_SampFreq_%d_Reflectors_%s_Threshold', arrayRows, arrayCols, sampleFrequency/1e6, size(reflectors, 2), titleName);
    savefig(finalFigure, filename + ".fig");
    saveas(finalFigure, filename + ".png");
end
