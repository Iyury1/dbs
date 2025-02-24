%% Callback function: updatePlots
function updatePlots(~,~)
    % Get current slider values
    fs = round(get(hSlider_fs, 'Value'));  % Ensure integer number of samples
    fc = get(hSlider_fc, 'Value');
    
    % Define number of points for the two time vectors
    N_sampled = fs;                % Number of samples for the "sampled" signal
    N_over = oversampleFactor * fs;  % Oversampled points
    
    % Use a time grid that spans one period: 0 to 2*pi.
    % Using increments so that the number of points exactly equals N_sampled and N_over.
    dt_sampled = 2*pi/N_sampled;
    dt_over = 2*pi/N_over;
    t_sampled = 0:dt_sampled:(2*pi - dt_sampled);
    t_over = 0:dt_over:(2*pi - dt_over);
    
    %% Define the modulated carrier signal.
    % Here we use a simple amplitude modulation:
    %   s(t) = [1 + 0.5*sin(t)] .* sin(fc*t)
    % The modulator (sin(t)) has a period of 2pi.
    oversampledSignal = (1 + 0.5*sin(t_over)) .* sin(fc*t_over);
    % Extract the sampled version by taking every oversampleFactor-th point.
    sampledSignal = oversampledSignal(1:oversampleFactor:end);
    
    %% Plot time domain signals
    % Oversampled time domain
    plot(ax1, t_over, oversampledSignal, 'b');
    ax1.XLim = [0 2*pi];
    title(ax1, sprintf('Oversampled Signal (Time Domain) [fs = %d, fc = %.2f]', fs, fc));
    
    % Sampled time domain
    plot(ax2, t_sampled, sampledSignal, 'r.-','MarkerSize',20);
    ax2.XLim = [0 2*pi];
    title(ax2, sprintf('Sampled Signal (Time Domain) [fs = %d, fc = %.2f]', fs, fc));
    
    %% Frequency domain calculations
    % For the oversampled signal:
    Y_over = fftshift(fft(oversampledSignal))/N_over;
    % Define frequency axis: we display from -((oversampleFactor*fs)/2) to +((oversampleFactor*fs)/2)
    f_over = linspace(-oversampleFactor*fs/2, oversampleFactor*fs/2, N_over);
    
    % For the sampled signal:
    Y_sampled = fftshift(fft(sampledSignal))/N_sampled;
    f_sampled = linspace(-fs/2, fs/2, N_sampled);
    
    %% Plot frequency domain signals
    % Oversampled frequency domain
    plot(ax3, f_over, abs(Y_over), 'b');
    ax3.XLim = [-oversampleFactor*fs/2, oversampleFactor*fs/2];
    title(ax3, sprintf('Oversampled Signal (Frequency Domain) [fs_{over} = %d]', oversampleFactor*fs));
    
    % Sampled frequency domain
    plot(ax4, f_sampled, abs(Y_sampled), 'r');
    ax4.XLim = [-fs/2, fs/2];
    title(ax4, sprintf('Sampled Signal (Frequency Domain) [fs = %d]', fs));
end