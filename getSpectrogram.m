function [time, f, spectrogram] = getSpectrogram(time, x, sr, varargin)
%getSpectrogram calculates spectrogram
%
%   getSpectrogram(time, x, sr, [t, [iol, [n, [wol]]]]) gets spectrogram of
%   signal x, sampled at sr Hz, using intervals of length t with iol %
%   overlap, and within each interval, using Welch's method to estimate the
%   Power Spectral Density (PSD) using n Hann windows with wol % overlap.
%   The default settings are 10-second intervals with 50% overlap and 2
%   Hann windows with 50% overlap. The datetime values are passed to the
%   function in the time variable to generate the x-axis labels for the
%   plot
%

% x = time series
% sr = sample rate (Hz)
% T_interval = duration of interval to use (sec)
% OL_interval = overlap between consecutive intervals (%)
% N_window = number of Hanning windows to use within each interval
% OL_window = overlap between consectuive Hanning windows (%)

% Nomenclature
% T = duration (sec)
% t = array of time values (sec)
% L = length of vector (samples)

% Charalambous et al., 2020
    
    T_interval_default = 10;
    OL_interval_default = 50;
    N_window_default = 2;
    OL_window_default = 50;
    
    if nargin == 3
        T_interval  = T_interval_default;
        OL_interval = OL_interval_default;
        N_window    = N_window_default;
        OL_window   = OL_window_default;
    elseif nargin == 4
        T_interval  = varargin{1};
        OL_interval = OL_interval_default;
        N_window    = N_window_default;
        OL_window   = OL_window_default;
    elseif nargin == 5
        T_interval  = varargin{1};
        OL_interval = varargin{2};
        N_window    = N_window_default;
        OL_window   = 50;
    elseif nargin == 6
        T_interval  = varargin{1};
        OL_interval = varargin{2};
        N_window    = varargin{3};
        OL_window   = OL_window_default;
    elseif nargin == 7
        T_interval  = varargin{1};
        OL_interval = varargin{2};
        N_window    = varargin{3};
        OL_window   = varargin{4};
    end

    disp(join(['[SPspectrogram]: Calculating spectrogram for ', inputname(2), ...
               ', sample rate = ', num2str(sr), ' Hz']));

    disp(join(['[SPspectrogram]: Using ', ...
               num2str(T_interval), ' sec intervals with ', ...
               num2str(OL_interval), '% overlap and ' , ...
               num2str(N_window), ' Hann windows with ' , ...
               num2str(OL_window), '% overlap']));
        
    % Calculate length of each time series interval
    L_interval = T_interval*sr;

    % Calculate length of Hanning window within each time series interval
    L_window = floor(L_interval/N_window);

    % Generate Hanning window
    w = hann(L_window); 

    % Generate the time of each sample in each time series interval
    t_interval = [0:1/sr:T_interval]';

    % Calculate the spacing between intervals (default is 50%)
    OL_interval = 1 - OL_interval/100;
    T_interval_spacing = OL_interval*T_interval; 

    % Calculate the overlaps between Hanning windows (default is 50%)
    OL_window = OL_window/100;
    L_window_overlap = OL_window*L_window;

    % Initialise data arrays
    tsearch = [];
    spectrogram = [];
    n_interval = 1;

    % Calculate end time of time series
    T_max = length(x)/sr;
    
    % Progress bar
    maxt = floor(T_max-T_interval); 

    % Calculate noise in each interval
    for t_start = 1:T_interval_spacing:floor(T_max-T_interval)

        % Display progress
        % disp(join(['Progress: ', num2str(100*t_start/maxt), '%']));
        % Calculate the stop time for this interval
        t_stop = t_start + T_interval;
        
        % Extract the signal samples for this interval
        x_start = round(sr*t_start + 1);
        x_stop  = round(sr*t_stop + 1); 
        
        if x_stop > length(x)
            % If x_stop exceeds the last index of x, then set limit the
            % value of x_stop to the length of x
            x_stop = length(x);
        end
        
        x_interval = x(x_start:x_stop);
        
        if x_stop == length(x)
            % If x_stop exceeds the last index of x, and x_interval has
            % therefore been shortened, need to also shorten t_interval so
            % that polyfit will work
            t_interval = t_interval(1:length(x_interval));
        end

        % Detrend the signal samples for this interval
        % Set opol = 0 for no trend removal
        opol = 2;
        [p, ~, mu] = polyfit(t_interval, x_interval, opol);
        x_fit = polyval(p, t_interval, [], mu);
        x_interval = x_interval - x_fit;
                                        
%         % Calculate the number of discrete Fourier transform (DFT) points to use in
%         % the PSD estimate. The default nfft is the greater of 256 or the next
%         % power of 2 greater than the length of the windows (pow2aboveL_window)
%         pow2aboveL_window = ceil(log(L_window)/log(2));
%         nfft = max(256, 2^pow2aboveL_window);
% 
%         % Calculate the PSD for this interval
%         [PSD_interval, f_interval] = pwelch(x_interval, ...
%                                             w, ...
%                                             L_window_overlap, ...
%                                             nfft, ...
%                                             sr);
                                        
        % Calculate the PSD for this interval
        % (Swapped in original way of calculating
        % see WTP email 13/04/19)
        SampleNumber = L_interval;
        [PSD_interval, f_interval] = pwelch(x_interval, ...
                                            w, ...
                                            floor(SampleNumber/N_window/2), ...
                                            floor(SampleNumber/N_window), ...
                                            sr);
                                        
        % Calculate ASD for this interval
        ASD_interval = PSD_interval;%sqrt(PSD_interval);

        % Append the interval start time to a vector recording the interval
        % start times
        tsearch = [tsearch; t_start];

        % Append the spectrogram
        spectrogram = horzcat(spectrogram, ASD_interval);

        n_interval = n_interval + 1;
    end

    t = tsearch + T_interval/2;
    f = f_interval;

         
    time = time(1) + seconds(tsearch + T_interval/2);


end