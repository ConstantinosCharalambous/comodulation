%% -------------------------------------------------------------------- %%
%                                                                       %
% Supplementary code for "A Comodulation Analysis of Atmospheric Energy %
% Injection into the Ground Motion at InSight, Mars".                   %
%                                                                       %
% This code plots the environmental SNR given example datasets for      % 
% Marsquake events. The output is identical to Figures 11 & 12          %
% of the paper.                                                         %
%                                                                       %
% Developed for InSight mission to Mars. No warranty is implied.        %
%                                                                       %
%%--------------------------------------------------------------------- %%

%---------------------------------------%
% Authors: C. Charalambous et al., 2020 %
%---------------------------------------%  

clear all
close all

%-----------------------------------------------------------------------------------------------------------------------------------------------------%
% a. To obtain the operational aspect of the environmental SNR metric implemented in Marsquake catalogue V2 (doi:10.12686/a7): (i) for the LF and BB  %
% the power is calculated for fmin = 0.2Hz and fmax = 0.5Hz, (ii) for the HF, VF and 2.4Hz events the power is taken at the 2.4Hz resonance, for which% 
% we use the bandwidth fmin = 2.2Hz and fmax = 2.6Hz.                                                                                                 %
% b. For significant events that vary beyond this frequency range such as S0173a, S0128a and S0235b shown here, the corresponding frequency range is  %
% adapted from Giardini et al., 2020 and subsequent papers based on Marsquake on newly released catalogues.     
%-----------------------------------------------------------------------------------------------------------------------------------------------------%


% Uncomment to run another of the example events and change parameters in lines 93-105 accordingly
load event_S0173a.mat % S0173a LF Event, therefore fmin = 0.18, fmax = 0.48 (Giardini et al., 2020), KMM = 1000 s (change below in lines 93-105)
%load event_S0128a.mat % S0128a VF Event, therefore fmin = 2.2, fmax = 2.6 (Giardini et al., 2020), KMM = 500 s (change below in lines 93-105)
%load event_S0235b.mat % S0235b BB Event, therefore fmin = 0.15, fmax = 0.91 (Giardini et al., 2020), KMM = 1000 s (change below in lines 93-105)

%%
for i = 1:length(e)
    
sr = max([e(i).spsrmax, e(i).vbbsrmax]); % sample rate of data


T_interval = 50; % spectrogram window length -> Wlen

OL_interval = 90; % spectrogram window overlap -> OLwin
N_window = 1.1; % number of PSD averages -> Nav
OL_window = 90; % PSD overlap -> OLpsd

e(i).vbbzne(1).spect = [];
e(i).vbbzne(2).spect = [];
e(i).vbbzne(3).spect = [];
e(i).spzne(1).spect = [];
e(i).spzne(2).spect = [];
e(i).spzne(3).spect = [];
e(i).spzne(1).freq_spec = [];
e(i).spzne(1).tspec = [];
e(i).vbbzne(1).freq_spec = [];
e(i).vbbzne(1).tspec = [];
e(i).spzne(2).freq_spec = [];
e(i).spzne(2).tspec = [];
e(i).vbbzne(2).freq_spec = [];
e(i).vbbzne(2).tspec = [];
e(i).spzne(3).freq_spec = [];
e(i).spzne(3).tspec = [];
e(i).vbbzne(3).freq_spec = [];
e(i).vbbzne(3).tspec = [];
e(i).pressdata.spect = [];
e(i).pressdata.freq_spec = [];
e(i).pressdata.tspec = [];
e(i).winddata.spect = [];
e(i).winddata.freq_spec = [];
e(i).winddata.tspec = [];

%--------------------------------------------------%
% Time is the global synchronized time vector for pressure, wind and VBB ZNE acceleration data
%--------------------------------------------------%
Time = (e(i).sp(1).t);

%--------------------------------------------------%
% e(i).vbbzne(1).a is VBB Z acceleration data in m/s^2
% e(i).vbbzne(2).a is VBB N acceleration data in m/s^2
% e(i).vbbzne(3).a is VBB E acceleration data in m/s^2
% e(i).pressdata.data is the pressure data in Pascals
% e(i).winddata.data_speed is the wind data in m/s
%--------------------------------------------------%

[e(i).pressdata.tspec, e(i).pressdata.freq_spec, e(i).pressdata.spect] = getSpectrogram(Time, (e(i).pressdata.data), sr, ...
    T_interval, OL_interval, N_window, OL_window);
[e(i).vbbzne(1).tspec, e(i).vbbzne(1).freq_spec, e(i).vbbzne(1).spect] = getSpectrogram(Time, e(i).vbbzne(1).a, sr, ...
    T_interval, OL_interval, N_window, OL_window);
[e(i).vbbzne(2).tspec, e(i).vbbzne(2).freq_spec, e(i).vbbzne(2).spect] = getSpectrogram(Time, e(i).vbbzne(2).a, sr, ...
    T_interval, OL_interval, N_window, OL_window);
[e(i).vbbzne(3).tspec, e(i).vbbzne(3).freq_spec, e(i).vbbzne(3).spect] = getSpectrogram(Time, e(i).vbbzne(3).a, sr, ...
    T_interval, OL_interval, N_window, OL_window);


end

%% This code on calculates the SNR

%-----------------------------------------------------------------------------------------------------------------------------------------------------%
% KMM = 1000 s and LMM = 0 for longer period LF/BB events and KMM = 500 s and LMM = 0 for the HF/VF/2.4Hz events
%-----------------------------------------------------------------------------------------------------------------------------------------------------%

TwinK= 1000; %-> KMM - moving moment matching window length prior. 
TwinL = 0; %-> LMM - moving moment matching window length post
TwinMeanK = 500; %-> KSNR - SNR averaging window length prio
TwinMeanL = 500; %->  LSNR - SNR window length post
Nsigma = 5;

%-----------------------------------------------------------------------------------------------------------------------------------------------------%
% a. For the operational aspect of the environmental SNR metric: (i) for the LF and BB the power is calculated for fmin = 0.2Hz and fmax = 0.5Hz and  %                              
% (ii ) for the HF, VF and 2.4Hz events the power is taken at the 2.4Hz resonance, for which we can use the bandwidth fmin = 2.2Hz and fmax = 2.6Hz.
% b/ For significant events that vary beyond this frequency range, e.g. Giardini et al., 2020, see frequency range of events suggested therein        

LowF1 = 1/5; %-> parameter fmin 
HiF1 = 1/2; %-> parameter fmax 

%-----------------------------------------------------------------------------------------------------------------------------------------------------%

%toggle plot spectrogram
plotspect = 0; 

sr = max([e(i).spsrmax, e(i).vbbsrmax]);
f_1 = e(i).vbbzne(1).freq_spec;
Time = e(i).vbb(1).t; % synced time vector for all variables
tspec = e(i).vbbzne(1).tspec;


IndStart1 = find(f_1>=LowF1);
IndEnd1 = find(f_1>=HiF1);
df = f_1(3)-f_1(2);

Z_env_full = sqrt(sum(e(i).vbbzne(1).spect(IndStart1(1):IndEnd1(1),:))*df);
N_env_full = sqrt(sum(e(i).vbbzne(2).spect(IndStart1(1):IndEnd1(1),:))*df);
E_env_full = sqrt(sum(e(i).vbbzne(3).spect(IndStart1(1):IndEnd1(1),:))*df);

Press_env = sqrt(sum(e(i).pressdata.spect(:,:))*df);
Wind_env = interp1(Time,e(i).winddata.data_speed,tspec,'nearest');


% 
Xlimits = [min(e(i).vbbzne(1).tspec) max(e(i).vbbzne(1).tspec)];
Ylimits = [N_window/T_interval, sr/2.5];
climits = [-10 -7];
fontsize = 14;
fontsizebar = 12;
Ax = 'log';
if plotspect == 1
figure(i)
clf
subplot(4,1,1)
hold on;
colormap jet;
surf(e(i).vbbzne(1).tspec, e(i).vbbzne(1).freq_spec, log10(sqrt(e(i).vbbzne(1).spect)), 'EdgeColor', 'none')
shading(gca, 'flat');
set(gca,'YScale',Ax);
caxis(climits);
set(gca, 'layer', 'top');
colorbar;
h = colorbar;

% x-axis
% xlabel('Time (UTC)');
xlim(Xlimits);

% y-axis
ylabel('Frequency (Hz)','FontSize', fontsize, 'FontWeight', 'bold');
ylim(Ylimits);
% z-axis
ylabel(h, 'Accel (log(m s^{-2} Hz^{-1/2}))','FontSize', fontsizebar, 'FontWeight', 'bold')
box on
title('Z','FontSize', fontsize, 'FontWeight', 'bold')

subplot(4,1,2)
hold on;
colormap jet;
surf(e(i).vbbzne(2).tspec, e(i).vbbzne(2).freq_spec, log10(sqrt(e(i).vbbzne(2).spect)), 'EdgeColor', 'none')
shading(gca, 'flat');
set(gca,'YScale',Ax);
caxis(climits);
set(gca, 'layer', 'top');
colorbar;
h = colorbar;

% x-axis
xlim(Xlimits);

% y-axis
ylabel('Frequency (Hz)','FontSize', fontsize, 'FontWeight', 'bold');
ylim(Ylimits);
% z-axis
ylabel(h, 'Accel (log(m s^{-2} Hz^{-1/2}))','FontSize', fontsizebar, 'FontWeight', 'bold')
box on
title('N','FontSize', fontsize, 'FontWeight', 'bold')

subplot(4,1,3)
hold on;
colormap jet;
surf(e(i).vbbzne(3).tspec, e(i).vbbzne(3).freq_spec, log10(sqrt(e(i).vbbzne(3).spect)), 'EdgeColor', 'none')
shading(gca, 'flat');
set(gca,'YScale',Ax);
caxis(climits);
set(gca, 'layer', 'top');
colorbar;
h = colorbar;

% x-axis
xlim(Xlimits);

% y-axis
ylabel('Frequency (Hz)','FontSize', fontsize, 'FontWeight', 'bold');
ylim(Ylimits);
% z-axis
ylabel(h, 'Accel (log(m s^{-2} Hz^{-1/2}))','FontSize', fontsizebar, 'FontWeight', 'bold')
box on
title('E','FontSize', fontsize, 'FontWeight', 'bold')

subplot(4,1,4)
plot(tspec, Z_env_full)
hold on
plot(tspec, N_env_full)
plot(tspec, E_env_full)
xlim(Xlimits);
h = colorbar;
ylabel('Accel (m s^{-2})','FontSize', fontsizebar, 'FontWeight', 'bold')
xlabel('Time (UTC)', 'FontSize', fontsizebar, 'FontWeight', 'bold')
legend('Z', 'N', 'E')
title(['Acceleration envelopes bandpassed between f: ', num2str(LowF1), ' - ', num2str(HiF1), ' Hz'])
h.Visible = 'off';
end


%------------------- Match the moments and SNR ---------------------------%

wind = Wind_env';
VBBN = N_env_full;
VBBE = E_env_full;
VBBZ = Z_env_full;
press = Press_env;

time = tspec;


duration = min([length(VBBZ)...
    length(VBBN) length(VBBE)]);

wind = wind(1:duration);

if e(i).NumPressDataRates>0
press = press(1:duration);
else
    press = ones(1,length(wind));
end

VBBZ = VBBZ(1:duration);
VBBN = VBBN(1:duration);
VBBE = VBBE(1:duration);
time = time(1:duration);

% e(i).t_event is the time of the detected seismic event as given in doi:10.12686/a7
% e(i).duration is the duration of the detected seismic event as given in doi:10.12686/a7
IdxSNRstart_1 = find(time>=e(i).t_event-minutes(2));
IdxSNRend_1 = find(time>=(e(i).t_event+e(i).duration));

e(i).WindMean = mean(wind(IdxSNRstart_1(1):IdxSNRend_1(1)),'omitnan');

e(i).PressMean = mean(press(IdxSNRstart_1(1):IdxSNRend_1(1)),'omitnan');


wind(wind <= 0) = NaN;

% Produce moment-matched data to selected seis axis

%Select Vertical, north or east component
seis = VBBZ;
% Set seis to a minimum offset

% Match in log
seis = log(seis);
wind = log(wind);
press = log(press);
pressNorm  = press;
seisNorm  = seis;


% Window length for moment matching
Tdelta = seconds(time(2)-time(1))


n=1;

winK = floor(TwinK/Tdelta);% window in samples
winL = floor(TwinL/Tdelta);% window in samples
delay = floor(winK);

% Seis
% --------
% Calculate running means and variances, offset by win 
seisMean  = movmean(seis, [winK-1 winL]); 
seisMean = [ones(1, winK)*seisMean(1) seisMean(1:end-winK)];

seisVar = movvar(seis, [winK-1 winL]);
seisVar = [ones(1, delay)*seisVar(1) seisVar(1:end-delay)];

% Set to press to NaN for N sigma excursions   
seisSD  = exp((seis - seisMean)./sqrt(seisVar));
idx = seisSD > Nsigma;
seisNorm(idx) = NaN;
seisNaNs = sum(isnan(seisNorm))/length(seisNorm)

%  Recalculate running means and variances
seisMean  = movmean(seisNorm, [winK-1 0], 'omitnan');
seisMean = [ones(1, delay)*seisMean(1) seisMean(1:end-delay)];

seisVar = movvar(seisNorm, [winK-1 0], 'omitnan');
seisVar = [ones(1, delay)*seisVar(1) seisVar(1:end-delay)];

% Wind
% --------
windMean  = movmean(wind, [winK-1 winL], 'omitnan');
windMean = [ones(1, delay)*windMean(1) windMean(1:end-delay)];

windVar = movvar(wind, [winK-1 winL], 'omitnan');
windVar = [ones(1, delay)*windVar(1) windVar(1:end-delay)];


% Pressure
% --------
% Calculate running means and variances, offset by win 
pressMean  = movmean(press, [winK-1 winL]);
pressMean = [ones(1, delay)*pressMean(1) pressMean(1:end-delay)];

pressVar = movvar(press, [winK-1 winL]);
pressVar = [ones(1, delay)*pressVar(1) pressVar(1:end-delay)];

% Set to press to NaN for N sigma excursions   
pressSD  = exp((press - pressMean)./sqrt(pressVar));
idx = pressSD > Nsigma;
pressNorm(idx) = NaN;

pressNaNs = sum(isnan(pressNorm))/length(pressNorm)

%  Recalculate pressure running means and variances

pressMean  = movmean(pressNorm, [winK-1 winL], 'omitnan');
pressMean = [ones(1, delay)*pressMean(1) pressMean(1:end-delay)];

pressVar = movvar(pressNorm, [winK-1 winL], 'omitnan');
pressVar = [ones(1, delay)*pressVar(1) pressVar(1:end-delay)];

% Moment matching

windMM = (wind - windMean).*sqrt(seisVar./windVar)+seisMean;
pressMM = (press - pressMean).*sqrt(seisVar./pressVar)+seisMean;

%SNR1 calc
winMean = 1;
SNR1_wind = movmean(exp(2*(seis-windMM)),winMean, 'omitnan');
SNR1_press = movmean(exp(2*(seis-pressMM)),winMean, 'omitnan');


%SNR2 calc
winMeanK = floor(TwinMeanK/Tdelta);
winMeanL = floor(TwinMeanL/Tdelta);
SNR2_wind = movmean(exp(2*(seis-windMM)),[winMeanK-1, winMeanL] , 'omitnan');
SNR2_press = movmean(exp(2*(seis-pressMM)),[winMeanK-1, winMeanL], 'omitnan');


IdxSNRstart = find(time>=e(i).t_event);
IdxSNRend = find(time>=(e(i).t_event+e(i).duration/2));


[MaxVal_SNR1wind, IndMax_SNR1wind] = max(SNR1_wind(IdxSNRstart(1):IdxSNRend(1)));
[MaxVal_SNR1press, IndMax_SNR1press] = max(SNR1_press(IdxSNRstart(1):IdxSNRend(1)));

[MaxVal_SNR2wind, IndMax_SNR2wind] = max(SNR2_wind(IdxSNRstart(1):IdxSNRend(1)));
[MaxVal_SNR2press, IndMax_SNR2press] = max(SNR2_press(IdxSNRstart(1):IdxSNRend(1)));

e(i).SNR1wind = MaxVal_SNR1wind;
e(i).SNR1press = MaxVal_SNR1press;

e(i).SNR2wind = MaxVal_SNR2wind;
e(i).SNR2press = MaxVal_SNR2press;

figure(100)
clf

ax1 = subplot(3,1,1)
plot(time, exp(seis), 'k')
hold on
plot(time, exp(windMM), '.b')
plot(time, exp(pressMM), 'r')
plot(time(IdxSNRstart_1(1):IdxSNRend_1(1)), exp(seis(IdxSNRstart_1(1):IdxSNRend_1(1))), 'k','LineWidth',2)

ylim([0 1.5*max(exp(seis))])
legend('seis', 'wind MM', 'press MM')
title('Moment-matched injection')
set(gca, 'YScale', 'log')
hold off

ax2 = subplot(3,1,2)

hold on
plot(time, SNR1_wind, '.b')
plot(time, SNR1_press, 'r')
plot(time(IndMax_SNR1wind+IdxSNRstart(1)),MaxVal_SNR1wind,'*k','MarkerSize',20)
plot(time(IndMax_SNR1press+IdxSNRstart(1)),MaxVal_SNR1press,'*k','MarkerSize',20)
ylim([0.1 150])
set(gca, 'YScale', 'log')
legend( 'wind SNR1','press SNR1',...
    'SNR1 max val')
title('SNR1')


ax3 = subplot(3,1,3)

semilogy(time, SNR2_wind, 'b.')
hold on
semilogy(time, SNR2_press, 'r')
plot(time(IndMax_SNR2wind+IdxSNRstart(1)),MaxVal_SNR2wind,'*k','MarkerSize',20)
plot(time(IndMax_SNR2press+IdxSNRstart(1)),MaxVal_SNR2press,'*k','MarkerSize',20)
title(['SNR2 over ' num2str(TwinMeanK+TwinMeanL) ' s'])
legend( 'wind SNR2','press SNR2',...
    'SNR2 max val')
ylim([0.1 100])


linkaxes([ax1 ax2 ax3 ], 'x')
xlim([e(i).t_event-hours(1), e(i).t_event+hours(1)])

