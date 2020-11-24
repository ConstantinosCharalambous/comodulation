%% ----------------------------------------------------------------------- %%
%                                                                          %
% Supplementary code for "A Comodulation Analysis of Atmospheric Energy    %
% Injection into the Ground Motion at InSight, Mars".                      %
%                                                                          %
% This code produces the relationships shown in Figure 8. It uses the      %
% envelope of the VBB Zground acceleration to obtain predictions of the    %
% wind speed over various frequency bands. It then plots the predictions   %
% against the measured wind speed. Finaly it shows how the mean and        %
%, i.e. the first two moments used to matched the time series, vary across %
% the frequency range.                                                     %
%                                                                          %
% Developed for InSight mission to Mars. No warranty is implied.           %
%                                                                          %
%%----------------------------------------------------------------------- %%

%---------------------------------------%
% Authors: C. Charalambous et al., 2020 %
%---------------------------------------%  

% Figure 8 
% The procedures below are shown specificaly for Figura 8a i.e. wind speed predicted
% from the VBB Z power. However, the same procedures can be applied to obtain 
% Figure 8b, by replacing VBB Z with the pressure time series countsPRE.
% Note that pressure is at 10 sps contrary to 20 sps of VBB Z, so the
% frequency range of this analysis should reach up to 4 Hz.

clear all
load sols_237_239_vbb_wind_pressure.mat

ftSz = 16.5; % font size

%% ------------------------------------------------------------------- %%
% This section extracts wind speed predictions from the sqrt of 
% the VBB Z power of ground acceleration over various frequency bands
% through each of the iterations in the loop. It also plots the wind
% speed predictions against measured wind speeds. 
% -------------------------------------------------------------------- %%

meanPower = [];
stdPower = [];

% Initialize paramaters
tinterval = 100; % time window in seconds for which we extract power
nAvg = 1.2;
sampleNumber = tinterval*srVBB; % samples of the signal per window
w  = hann(floor(sampleNumber)); % hanning window

freqStart = 0.1; % lower frequency to start with
freqBand = 1; % frequency range for each window applied
numBands = 8; % number of frequency bands to iterate

CM = flipud(parula(numBands)); % apply your preferred colormap

figure(8)
clf

% Iterate through the frequency bands requested
for i = 1 : numBands
    
    disp(['Running: ', num2str(i), '/', num2str(numBands), ' frequency band'])
    
    fLow = freqStart+(i-1)*freqBand;
    fHigh = freqStart+freqBand+(i-1)*freqBand;
    
    %-------------------- VBB Z ---------------------------%
    [bandZ,d3] = bandpass(Zvbb,[fLow fHigh], srVBB, 'ImpulseResponse', 'iir', 'Steepness', 0.95);
    
    [~, fVBB, tVBB, pVBB] = (spectrogram(bandZ, w, ...
        floor(sampleNumber/nAvg), ...
        floor(sampleNumber/nAvg), ...
        srVBB,'yaxis'));
    
    timeSpectro = timeVBB(1)+seconds(tVBB);
    
    
    % Obtain moments (mean and variance) from Wind speed and the sqrt of the VBB Z power
    meanWS = mean(speedWS);
    stdWS = std(speedWS);
    df = fVBB(2) - fVBB(1);
    meanPrVBB = mean(sqrt(sqrt(sum(pVBB.*df))));
    stdPrVBB = std(sqrt(sqrt(sum(pVBB.*df))));
    meanPower(i) = meanPrVBB;
    stdPower(i) = stdPrVBB;
    pE_MoM =( sqrt(sqrt(sum(pVBB.*df))) - meanPrVBB )/stdPrVBB * stdWS + meanWS;
    
    
    % Plot measured wind speed time series 
    figure(8)
    ax354a = subplot(1,7,[1 3])
    hold on
    if i == 1
        scatter(timeWS, medfilt1(speedWS,5), 3, [204/255 204/255 204/255], 'filled')
    end
    
    
    % Plot the predicted wind speed from seismic data 
    scatter(timeSpectro, pE_MoM, 3, CM(i,:), 'MarkerFaceAlpha', 0.8, 'MarkerEdgeAlpha', 0.3)
    xlim([timeWS(1) timeWS(end)])
    box on
    
    
    % Create timetable and sync in order to get a scatter plot of predicted
    % vs measured wind speed
    peTT = timetable(timeSpectro', pE_MoM');
    windTT = timetable(timeWS, speedWS);
    
    vbbModeTempTT = [];
    vbbModeTempTT = synchronize(windTT, peTT, timeSpectro,'linear');
    vbbModeTempTT.Properties.VariableNames = {'windSpeed' 'powerVBB'};
    
    ylabel('\bf{Wind Speed (ms$^{-1}$)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    xlabel('\bf{Time}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    set(gca,'linewidth',1.1)
    set(gca,'TickLabelInterpreter','latex')
    set(gca,'FontSize',ftSz)
    ylim([0 20])
    xlim
    
    subplot(1,7,[5 6])
    hold on
    scatter(vbbModeTempTT.windSpeed, vbbModeTempTT.powerVBB, 3, CM(i,:), 'MarkerFaceAlpha', 0.8,'MarkerEdgeAlpha', 0.3)

end

box on
ylabel('\bf{Predicted Wind Speed (ms$^{-1}$)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
xlabel('\bf{Wind Speed (ms$^{-1}$)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
set(gca,'linewidth',1.1)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',ftSz)
xlim([0 20])
ylim([0 20])


h = colorbar
colormap(CM)
set(h, 'YTick',[0 0.25 0.5 0.75 1],'YTickLabel', {'0','2','4','6','8'});
h.Label.String = '\bf{Frequency Band (Hz)}';
h.Label.FontWeight = 'bold';
h.Label.Interpreter = 'latex';
h.TickLabelInterpreter = 'latex';


%% ------------------------------------------------------------------- %%
% Plots the mean and variance of the sqrt of the VBB Z power that is  
% applied in the method of moments for various frequency bands. 
% This section is a repeat of the previous section but resolves moments
% through smaller frequency bands to observe the evolution and variability
% of finer changes e.g. lander modes at 4 hz, or tick noise at 1 Hz
% -------------------------------------------------------------------- %%


meanPower = [];
stdPower = [];

% Initialize paramaters
freqStart = 0.1; % lower frequency to start with
freqBand = 0.1; % frequency range for each band applied
numBands = 80; % number of frequency bands to iterate

CM = flipud(parula(numBands)); % apply your preferred colormap

% Iterate through the frequency bands requested
for i = 1 : numBands
   
    disp(['Running: ', num2str(i), '/', num2str(numBands), ' frequency band'])

    fLow = freqStart+(i-1)*freqBand;
    fHigh = freqStart+(i-1)*freqBand + 0.1;
    
    %-------------------- VBB Z  ---------------------------%
    [bandZ, ~] = bandpass(Zvbb,[fLow fHigh],srVBB,'ImpulseResponse','iir','Steepness',0.95);
    
    [~, fVBB, tVBB, pVBB] = (spectrogram(bandZ, w, ...
        floor(sampleNumber/nAvg), ...
        floor(sampleNumber/nAvg), ...
        srVBB,'yaxis'));
    
    df = fVBB(2) - fVBB(1);
     
    
    % Obtain moments (mean and variance) from Wind speed and the sqrt of the VBB Z power
    meanWS = mean(speedWS);
    stdWS = std(speedWS);
    meanPrVBB = mean(sqrt(sqrt(sum(pVBB.*df))));
    stdPrVBB = std(sqrt(sqrt(sum(pVBB.*df))));
    meanPower(i) = meanPrVBB;
    stdPower(i) = stdPrVBB;
    
end


figure(9)
hold on
yyaxis left
lins = 1:1:numBands;
plot(freqStart + (lins - 1) * freqBand, meanPower,'-', 'LineWidth', 2, 'Color', [204/255 204/255 204/255])
ylabel('\bf{Mean}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
yyaxis right
ylabel('\bf{Variance}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
plot(freqStart + (lins-1) * freqBand, stdPower, '-', 'LineWidth', 2, 'Color', [80/255 80/255 80/255])
ax = gca
ax.YAxis(1).Color = [130/255 130/255 130/255];
ax.YAxis(2).Color = [40/255 40/255 40/255];
box on

xlabel('\bf{Frequency (Hz)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
set(gca,'linewidth',1.1)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',ftSz)