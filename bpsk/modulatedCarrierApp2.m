function modulatedCarrierApp
    %% Default parameters
    fs_default = 1;       % Default sample frequency (samples per period)
    fc_default = 24;        % Default carrier frequency
    bd_default = 1;
    Tsym = 1 / bd_default;
    oversampleFactor = 10; % Oversampling factor

    %% Create the main uifigure
    fig = uifigure('Name','Modulated Carrier Simulation',...
                   'Position',[100 100 1200 800]);
               
    %% Create main grid layout: two rows
    % Row 1 will contain the axes; Row 2 will contain the slider controls.
    mainLayout = uigridlayout(fig, [2, 1]);
    mainLayout.RowHeight = {'4x','1x'};  % Adjust as needed
    mainLayout.ColumnWidth = {'1x'};
    
    %% Create panel and grid layout for axes (plots)
    axesPanel = uipanel(mainLayout);
    axesPanel.Layout.Row = 1;
    axesPanel.Layout.Column = 1;
    
    % Create a grid layout for 6 axes (3 rows, 2 columns)
    axesGrid = uigridlayout(axesPanel, [3,2]);
    axesGrid.RowHeight = {'1x','1x','1x'};
    axesGrid.ColumnWidth = {'1x','1x'};
    
    % Create six uiaxes
    ax_baseband_t = uiaxes(axesGrid);
    ax_baseband_t.Layout.Row = 1; ax_baseband_t.Layout.Column = 1;
    title(ax_baseband_t, 'Baseband Signal (Time Domain)');
    xlabel(ax_baseband_t, 'Time (s)');
    ylabel(ax_baseband_t, 'Amplitude');
    
    ax_baseband_f = uiaxes(axesGrid);
    ax_baseband_f.Layout.Row = 1; ax_baseband_f.Layout.Column = 2;
    title(ax_baseband_f, 'Baseband Signal (Frequency Domain)');
    xlabel(ax_baseband_f, 'Frequency (normalized)');
    ylabel(ax_baseband_f, 'Magnitude');


    ax_carrier_t = uiaxes(axesGrid);
    ax_carrier_t.Layout.Row = 2; ax_carrier_t.Layout.Column = 1;
    ax_carrier_t.XLim = [0, 2/fc_default];
    title(ax_carrier_t, 'Carrier Signal (Time Domain)');
    xlabel(ax_carrier_t, 'Time (s)');
    ylabel(ax_carrier_t, 'Amplitude');

    ax_carrier_f = uiaxes(axesGrid);
    ax_carrier_f.Layout.Row = 2; ax_carrier_f.Layout.Column = 2;
    title(ax_carrier_f, 'Carrier Signal (Frequency Domain)');
    xlabel(ax_carrier_f, 'Frequency (normalized)');
    ylabel(ax_carrier_f, 'Magnitude');
    
    ax_sampled_t = uiaxes(axesGrid);
    ax_sampled_t.Layout.Row = 3; ax_sampled_t.Layout.Column = 1;
    ax_sampled_t.XLim = [0, 2*Tsym];

    title(ax_sampled_t, 'Sampled Signal (Time Domain)');
    xlabel(ax_sampled_t, 'Time (s)');
    ylabel(ax_sampled_t, 'Amplitude');

    ax_sampled_f = uiaxes(axesGrid);
    ax_sampled_f.Layout.Row = 3; ax_sampled_f.Layout.Column = 2;
    title(ax_sampled_f, 'Sampled Signal (Frequency Domain)');
    xlabel(ax_sampled_f, 'Frequency (normalized)');
    ylabel(ax_sampled_f, 'Magnitude');

    
    %% Create panel and grid layout for slider controls
    controlsPanel = uipanel(mainLayout);
    controlsPanel.Layout.Row = 2;
    controlsPanel.Layout.Column = 1;
    
    % Create a grid layout for the controls: 3 rows (one per slider), 2 columns (slider and label)
    controlsGrid = uigridlayout(controlsPanel, [3,2]);
    controlsGrid.RowHeight = {'1x','1x','1x'};
    controlsGrid.ColumnWidth = {'3x','1x'};
    
    %% Create Sliders

    % Slider for Sample Frequency (fs)
    hSlider_fs = uislider(controlsGrid);
    hSlider_fs.Limits = [0.1 1];
    hSlider_fs.MajorTicks = 0.1:0.1:1;             % Display ticks at every 0.1 increment
    hSlider_fs.MinorTicks = [];                    % (Optional) Clear minor ticks if not needed
    hSlider_fs.MinorTicksMode = 'manual';

    hSlider_fs.Value = fs_default;
    hSlider_fs.Layout.Row = 1; 
    hSlider_fs.Layout.Column = 1;
    

    hSlider_fs.ValueChangedFcn = @(sld, event) hSlider_fsCallback(sld);
    function hSlider_fsCallback(sld)
        set(sld, 'Value', round((sld.Value - 0.1)/0.1)*0.1 + 0.1);
        % subtracts the lower limit (0.1), divides by the step size (0.1),
        % rounds to the nearest integer, and then reconstructs the snapped
        % value by multiplying back by 0.1 and adding the lower limit back.
        % Call updatePlots
        updatePlots();
    end

    label_fs = uilabel(controlsGrid);
    label_fs.Layout.Row = 1; 
    label_fs.Layout.Column = 2;
    label_fs.Text = sprintf('Fs: %f (GHz)', hSlider_fs.Value);
    


    % Slider for Carrier Frequency (fc)
    hSlider_fc = uislider(controlsGrid);
    hSlider_fc.Limits = [24 25];
    hSlider_fc.MajorTicks = 24:0.1:25;
    hSlider_fc.MinorTicks = 24:0.005:25;
    hSlider_fc.MinorTicksMode = 'manual';
    hSlider_fc.Value = fc_default;
    hSlider_fc.Layout.Row = 2; 
    hSlider_fc.Layout.Column = 1;
    
    hSlider_fc.ValueChangedFcn = @(sld, event) hSlider_fcCallback(sld);
    function hSlider_fcCallback(sld)
        % Snap the slider value
        set(sld, 'Value', round((sld.Value - 24)/0.005)*0.005 + 24);
        % Call updatePlots
        updatePlots();
    end

    label_fc = uilabel(controlsGrid);
    label_fc.Layout.Row = 2; 
    label_fc.Layout.Column = 2;
    label_fc.Text = sprintf('Fc: %f (GHz)', hSlider_fc.Value);
    
    % Slider for Baud Rate (Bd)
    hSlider_fsym = uislider(controlsGrid);
    hSlider_fsym.Limits = [1 100];
    hSlider_fsym.MajorTicks = 0:10:100;
    hSlider_fsym.MinorTicks = 1:1:100;
    hSlider_fsym.Value = bd_default;
    hSlider_fsym.Layout.Row = 3; 
    hSlider_fsym.Layout.Column = 1;
    
    hSlider_fsym.ValueChangedFcn = @(sld, event) hSlider_fsymCallback(sld);
        function hSlider_fsymCallback(sld)
        % Snap the slider value
        set(sld, 'Value', round((sld.Value - 1)/1)*1 + 1);
        % Call updatePlots
        updatePlots();
    end


    label_fsym = uilabel(controlsGrid);
    label_fsym.Layout.Row = 3; 
    label_fsym.Layout.Column = 2;
    label_fsym.Text = sprintf('Bd: %f (MHz)', hSlider_fsym.Value);
    
    %% Set slider callbacks to update plots when values change
    %hSlider_fs.ValueChangedFcn = @(sld,event) updatePlots();
    %hSlider_fc.ValueChangedFcn = @(sld,event) updatePlots();
    %hSlider_fsym.ValueChangedFcn = @(sld,event) updatePlots();

    %% Define updatePlots function
    function updatePlots()
        % Update slider labels with current values
        label_fs.Text = sprintf('Fs: %f (GHz)', hSlider_fs.Value);
        label_fc.Text = sprintf('Fc: %f (GHz)', hSlider_fc.Value);
        label_fsym.Text = sprintf('Bd: %f (MHz)', hSlider_fsym.Value);
        
        % Retrieve current slider values
        fs = hSlider_fs.Value * 1e9;
        fc = hSlider_fc.Value * 1e9;
        bd = hSlider_fsym.Value * 1e6;
        % Fsym is not used in this example update, but you can incorporate it as needed.
        

        % Oversampling factor (number of samples per symbol)
        % Sampling frequency based on the symbol rate and oversampling factor
        f_over = 8*fc;

        Tbd = 1 / bd;
        
        M = 10;
        SNR_dB = 100;

        % Generate M random bits (0 or 1)
        bits = randi([0,1], M, 1);
        
        % Map bits to BPSK symbols: 0 -> -1, 1 -> +1
        symbols = 2 * bits - 1;
    
        oversampleFactor = round(f_over/bd);

        % Create the baseband signal by repeating each symbol
        % 'kron' replicates each symbol oversampleFactor times.
        baseband = kron(symbols, ones(oversampleFactor, 1));
    
        % Create a time vector for the entire signal (in seconds)
        t_over = (0:length(baseband)-1) / f_over;

        % Generate the carrier signal: a cosine wave
        carrier = cos(2*pi*fc*t_over);
        
        % Modulate the carrier with the BPSK baseband signal
        m = baseband .* carrier';
        
        %% 2) Add noise
        rxRF_ideal = awgn(m, SNR_dB, 'measured'); % "Infinite-bit" approach, so no quant. yet
        
        % Also create a 1-bit version
        rxRF_1bit = sign(rxRF_ideal);  % Simply take sign

        %% 3) Sample the RF signals
        
        step = round(f_over / fs);

        
        sample_ideal = rxRF_ideal(1:step:end);
        
        sample_1bit = rxRF_1bit(1:step:end);

        t_sampled = (0:length(sample_1bit)-1) / fs;

        % Update time domain plot


        ax_baseband_t.XLim = [0, 5*Tbd];
        ax_baseband_t.YLim = [-2, 2];

        plot(ax_baseband_t, t_over, baseband, 'b');
        title(ax_baseband_t, sprintf('Baseband Signal (Time Domain) [Bd = %f MHz]', bd / 1e6));
        xlabel(ax_baseband_t, 'Time (s)'); ylabel(ax_baseband_t, 'Amplitude');

        %t_over(1:bd_idx), carrier(1:bd_idx)

        ax_carrier_t.XLim = [0, 5/f_over];
        plot(ax_carrier_t, t_over, m, 'r');
        title(ax_carrier_t, sprintf('Carrier Signal (Time Domain) [Fc = %.2f GHz]', fc / 1e9));
        xlabel(ax_carrier_t, 'Time (s)'); ylabel(ax_carrier_t, 'Amplitude');
        
        ax_sampled_t.XLim = [0, 5*Tbd];
        stem(ax_sampled_t, t_sampled, sample_1bit, 'm');
        title(ax_sampled_t, sprintf('Sampled Signal (Time Domain) [Fs = %.2f GHz]', fs / 1e9));
        xlabel(ax_sampled_t, 'Time (s)'); ylabel(ax_sampled_t, 'Amplitude');


        % Update frequency domain plots for the sampled signal

        Y_baseband = fftshift(fft(baseband, 1024))/1024;
        f_baseband = linspace(-f_over/2, f_over/2, 1024);

        disp(length(Y_baseband))
        disp(length(f_baseband))
        plot(ax_baseband_f, f_baseband, abs(Y_baseband), 'r');
        title(ax_baseband_f, sprintf('Sampled Signal (Frequency Domain) [fs = %f]', fs));
        xlabel(ax_baseband_f, 'Frequency (normalized)'); ylabel(ax_carrier_f, 'Magnitude');
        % Update frequency domain plots for the oversampled signal

        Y_m = fftshift(fft(m, 1024))/1024;
        f_m = linspace(-f_over/2, f_over/2, 1024);
        plot(ax_carrier_f, f_m, abs(Y_m), 'b');
        title(ax_carrier_f, sprintf('Carrier Signal (Frequency Domain) [fs_{over} = %f]', f_over));
        xlabel(ax_carrier_f, 'Frequency (normalized)'); ylabel(ax_baseband_f, 'Magnitude');

        % Update frequency domain plots for the sampled signal
        Y_sampled = fftshift(fft(sample_1bit, 1024))/1024;
        f_sampled = linspace(-fs/2, fs/2, 1024);
        plot(ax_sampled_f, f_sampled, abs(Y_sampled), 'r');
        title(ax_sampled_f, sprintf('Sampled Signal (Frequency Domain) [fs = %f]', fs));
        xlabel(ax_sampled_f, 'Frequency (normalized)'); ylabel(ax_carrier_f, 'Magnitude');



    end

    %% Initial update of plots
    updatePlots();
end
