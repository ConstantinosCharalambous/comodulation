%% ---------------------------------------------------------------------- %%
%                                                                         %
% Supplementary code for "A Comodulation Analysis of Atmospheric Energy   %
% Injection into the Ground Motion at InSight, Mars".                     %
%                                                                         %
% This code produces the relationships shown in Figure 3c. By             %
% demonstrating here how these relationships are obtained, the user will  %
% also be able to adapt these procedures to recreate the relationships in %
% figures 2c,d 4 & 5                                                      %
%                                                                         %
% Developed for InSight mission to Mars. No warranty is implied.          %
%                                                                         %
%%----------------------------------------------------------------------- %%

%---------------------------------------%
% Authors: C. Charalambous et al., 2020 %
%---------------------------------------%  

clear all
load sols_237_239_vbb_wind_pressure.mat

%% Figure 2 Example - Relationships between ground acceleration, pressure and wind

figure(2)
clf

numBands = 4; % number of frequency bands
tInterval = 100; % time window in seconds
nAvg = 1.2;

%----- Set Frequency bands ------%

fLow(1) = 0.01;
fHigh(1) = 0.1;

fLow(2) = 0.1;
fHigh(2) = 0.5;

fLow(3) = 0.5;
fHigh(3) = 0.99;

fLow(4) = 1.01;
fHigh(4) = 4;


%%  1. Calculate VBB Z power, pressure energy and then plot relationships

for band = 1 : numBands
         
    sampleNumber = tInterval*srVBB;
    w  = hann(floor(sampleNumber)); % hanning window

    [bandZvbb, ~] = bandpass(Zvbb, [fLow(band) fHigh(band)], srVBB, 'ImpulseResponse', 'iir', 'Steepness',0.95);
    
    % Get power values from VBB Z
    [~,fVBB,tVBB,pVBB] = (spectrogram(bandZvbb, w, ...
        floor(sampleNumber/nAvg), ...
        floor(sampleNumber/nAvg), ...
        srVBB,'yaxis'));
    
    df = fVBB(2) - fVBB(1);
    powerVBB = sqrt(sum(pVBB).*df);
    powerVBBTime = timeVBB(1)+seconds(tVBB);
    
    if fHigh(band) < srPRE/2.5
        
        [bandPRE, ~] = bandpass(countsPRE,[fLow(band) fHigh(band)], srPRE, 'ImpulseResponse', 'iir', 'Steepness', 0.95);
    else
        [bandPRE, ~] = bandpass(countsPRE,[fLow(band) srPRE/2.5], srPRE, 'ImpulseResponse', 'iir','Steepness', 0.95);
    end
    
    sampleNumber = tInterval*srPRE;
    w  = hann(floor(sampleNumber)); % hanning window
    
    % Get power values from Pressure
    [~, fPre, tPre, pPre ] = (spectrogram(bandPRE, w, ...
        floor(sampleNumber/nAvg), ...
        floor(sampleNumber/nAvg), ...
        srPRE,'yaxis'));
    
    df = fPre(2)-fPre(1);
    meanPowerPRE= sqrt(sum(pPre).*df);
    timePREspect = timePRE(1)+seconds(tPre);
    
    
    % Create timetables
    windTT = [];
    vbbTT = [];
    preTT = [];
    
    windTT = timetable(timeWS,speedWS);
    vbbTT = timetable(powerVBBTime',((powerVBB)'));
    preTT = timetable(timePREspect',((meanPowerPRE)'));
    
    % Sync timetables
    windVBB_TT = [];
    windPRE_TT = [];
    
    windVBB_TT = synchronize(windTT, vbbTT, powerVBBTime', 'linear');
    windPRE_TT = synchronize(windTT, preTT, powerVBBTime', 'linear');
    
    windVBB_TT.Properties.VariableNames = {'windSpeed' 'powerVBB'};
    windPRE_TT.Properties.VariableNames = {'windSpeed' 'powerPRE'};

    
% 2. PLOT RELATIONSHIPS between wind speed, accel energy and pressure energy
    
    ftSz = 18;
    ti = [];
    cm = [];
    TI = hours(minutes(colorVec)); % colorVec contains the color vector for the corresponding Local Mean Solar Time (LMST)
    cn = ceil(max(TI));                                             
    lmstColVec = fix(TI);
    lmstColVec = round(lmstColVec+1);
    
    
    %------------- WIND/VBB ------------------%
    
    ax = subplot(3,numBands,1 + band-1)
    
    scatter(windVBB_TT.windSpeed, windVBB_TT.powerVBB, 2, lmstColVec, 'filled')
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    set(gca,'FontName', 'Times New Roman','FontSize',ftSz)
    ylabel('\bf{Acceleration (ms$^{-2}$)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    xlabel('\bf{Wind Speed (ms$^{-1}$)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    box on
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.02 .02] , ...
        'YGrid'       , 'off'      , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'YMinorTick'  , 'on', ...
        'XMinorTick'  , 'on', ...
        'LineWidth'   , 1         );
    
    set(gca,'linewidth',1.1)
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    ylim([.6E-10 .3E-6])
    yticks([10^-10 10^-8])
    xlim([1 20])
    xticks([1 10])

    set(gca,'TickLabelInterpreter','latex')

    
    %------------- PRESSURE/WIND ------------------%
    
    ax2 = subplot(3,numBands,5 + band-1)
    
    scatter(windPRE_TT.windSpeed, windPRE_TT.powerPRE, 2, lmstColVec, 'filled')
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    set(gca,'FontName', 'Times New Roman','FontSize',ftSz)
    
    xlabel('\bf{Wind Speed (ms$^{-1}$)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    ylabel('\bf{Pressure (Pa)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    
    box on
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.02 .02] , ...
        'YGrid'       , 'off'      , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'YMinorTick'  , 'on', ...
        'XMinorTick'  , 'on', ...
        'LineWidth'   , 1         );
    
    set(gca,'linewidth',1.1)
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    ylim([0.001 0.2])
    yticks([10^-3 10^-2 10^-1])
    xlim([1 20])
    xticks([1 10])
    set(gca,'TickLabelInterpreter','latex')

    colormap(map) % map is the colormap for 
    
    %------------- VBB/PRE ------------------%
    
    ax3 = subplot(3,numBands,9+band-1)
    

    scatter(windPRE_TT.powerPRE, windVBB_TT.powerVBB, 2,lmstColVec,'filled')
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    set(gca,'FontName', 'Times New Roman','FontSize',ftSz)
    ylabel('\bf{Acceleration (ms$^{-2}$)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    xlabel('\bf{Pressure (Pa)}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
    
    box on
    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'in'     , ...
        'TickLength'  , [.02 .02] , ...
        'YGrid'       , 'off'      , ...
        'XColor'      , [.0 .0 .0], ...
        'YColor'      , [.0 .0 .0], ...
        'YMinorTick'  , 'on', ...
        'XMinorTick'  , 'on', ...
        'LineWidth'   , 1         );
    
    set(gca,'linewidth',1.1)
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    xlim([.001 0.2])
    set(gca,'TickLabelInterpreter','latex')
    
    ylim([.6E-10 .3E-6])
    xticks([10^-3 10^-2 10^-1])
    yticks([10^-10 10^-8])

    pos = ax3.Position;

    colormap(map)
    
end

set(gcf, 'Position', get(0, 'Screensize'));