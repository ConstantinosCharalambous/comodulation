function x_out = removeLow(x, sr, corner_freq)
%SPREMOVELOW Remove low-frequency components from a signal
%
%   y = removeLow(x, sr, fc) removes frequencies that pass a low-pass
%   filter with a corner frequency of fc from the signal x sampled at a
%   rate of sr Hz, returning the modified signal y.

    type = 'low';
    
    disp(join(['[SPlowpassfilt]: Filtering ', inputname(1), ...
        ', sample rate = ', num2str(sr), ' Hz', ...
        ', using a ', type, '-pass Chebyshev Type I filter', ...
        ', with a corner frequency of ', num2str(corner_freq)]));

    [b_h, a_h] = cheby1(3, 0.1, corner_freq/(sr*.5), type);
    x_filtered = filtfilt(b_h, a_h, x);
    x_out = x - x_filtered;
    
end