%% ---------------------------------------------------------------------- %%
%                                                                         %
% Supplementary code for "A Comodulation Analysis of Atmospheric Energy   %
% Injection into the Ground Motion at InSight, Mars".                     %
%                                                                         %
% This code produces the matched-moment time series of the envelopes and  %
% wind speed time series as shown in Figure 3b. However, it may easily    %
% adapted to reproduce  other figures in the paper  which incorportate the%
% the approach of the method of moments to match time series, i.e.        %
% Figure 1c, 2b & 14.                                                     %
%                                                                         %
% Developed for InSight mission to Mars. No warranty is implied.          %
%                                                                         %
%%----------------------------------------------------------------------- %%

%---------------------------------------%
% Authors: C. Charalambous et al., 2020 %
%---------------------------------------%  

%% 1. Plot the moment-matched time series of the pressure power, vertical ground acceleration power and wind speed

clear all
load sols_237_239_vbb_wind_pressure.mat

ftSz = 16.5;

% Calculate vbb energy, pressure energy

tinterval = 50;
nAvg = 1.2;

sampleNumber = tinterval*srVBB;
w  = hann(floor(sampleNumber)); % hanning window

fLow = 0.1;
fHigh = srVBB/2.5;

%-------------------- VBB Z ---------------------------%

[Zvbb_band, ~] = bandpass(Zvbb,[fLow fHigh],srVBB,'ImpulseResponse','iir','Steepness',0.95);

[~,fVBB,tVBB,pVBB] = (spectrogram(Zvbb_band, w, ...
    floor(sampleNumber/nAvg), ...
    floor(sampleNumber/nAvg), ...
    srVBB,'yaxis'));

df = fVBB(2)-fVBB(1);
powerVBB = (sqrt(sum(pVBB)*df));
powerVBBTime = timeVBB(1)+seconds(tVBB);

%-------------------- Pressure ---------------------------%

if fHigh < srPRE/2.5
    
    [PRE_band,d3] = bandpass(countsPRE, [fLow fHigh], srPRE,'ImpulseResponse','iir','Steepness',0.95);
else
    
    [PRE_band,d3] = bandpass(countsPRE, [fLow srPRE/2.5], srPRE,'ImpulseResponse','iir','Steepness',0.95);
end

sampleNumber = tinterval*srPRE;
w  = hann(floor(sampleNumber)); % hanning window

[~, fPre, tPre, pPre] = (spectrogram(PRE_band, w, ...
    floor(sampleNumber/nAvg), ...
    floor(sampleNumber/nAvg), ...
    srPRE,'yaxis'));

df = fPre(2)-fPre(1);
powerPRE = (sqrt(sum(pPre)*df));
powerPRETime = timePRE(1)+seconds(tPre);

% Create timetables for the derived time series of power 

windTT = [];
vbbTT =[];
preTT = [];

windTT = timetable(timeWS,speedWS);
vbbTT = timetable(powerVBBTime',((powerVBB)'));
preTT = timetable(powerPRETime',((powerPRE)'));

% Sync derived time series of power through timetables

windVBB_TT = [];
windPRE_TT = [];

windVBB_TT = synchronize(windTT,vbbTT,powerVBBTime','linear');
windPRE_TT = synchronize(windTT,preTT,powerVBBTime','linear');

windVBB_TT.Properties.VariableNames = {'windSpeed' 'powerVBB'};
windPRE_TT.Properties.VariableNames = {'windSpeed' 'powerPRE'};

% Get mean and variance for  wind, vbb and pressure

meanWS = mean(windVBB_TT.windSpeed);
stdWS = std(windVBB_TT.windSpeed);

meanVBB = mean(windVBB_TT.powerVBB);
stdVBB = std(windVBB_TT.powerVBB);

meanPRE = mean(windPRE_TT.powerPRE);
stdPRE = std(windPRE_TT.powerPRE);

% Match the first two moments of VBB and PRE to wind's first two moments (mean and variance) 

matchVBB = ((windVBB_TT.powerVBB - meanVBB)/stdVBB) * stdWS + meanWS;
matchPRE = ((windPRE_TT.powerPRE - meanPRE)/stdPRE) * stdWS + meanWS;

% Plot the moment-matched time series (single mean and variance matching) 

figure(1238231)
clf

subplot(1,3,[1 3])
hold on

scatter(windVBB_TT.timeWS, medfilt1(windVBB_TT.windSpeed,5),3,'MarkerFaceColor',[0.6 0.6 0.6],'MarkerEdgeColor','none','MarkerFaceAlpha',.5)
plot(windVBB_TT.timeWS, matchPRE,'Color', 'k')
plot(windVBB_TT.timeWS, matchVBB, 'Color', [77/255 190/255 238/255])

ylim([0 max(medfilt1(windVBB_TT.windSpeed,5))])
xlim([min(windVBB_TT.timeWS) max(windVBB_TT.timeWS)])

ylabel('\textbf{Wind Speed $\bf{(ms^{-1})}$}', 'FontSize', ftSz, 'FontWeight', 'bold','interpreter','latex');
box on
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.01 .01] , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.0 .0 .0], ...
    'YColor'      , [.0 .0 .0], ...
    'YMinorTick'  , 'off', ...
    'XMinorTick'  , 'off', ...
    'LineWidth'   , 1         );

set(gca,'linewidth',1.1)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',ftSz)
leg = legend('Wind Speed', 'VBB Z','Pressure')
set(leg,'Interpreter','latex')

