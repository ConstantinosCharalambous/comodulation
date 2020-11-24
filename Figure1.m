%% Figure 1 Example

clear all
load sols_237_239_vbb_wind_pressure.mat

%% Calculate VBB and Pressure spectrograms

tinterval = 200; % window length in seconds
nAvg = 1.5;
sampleNumber = tinterval*srVBB;
w  = hann(floor(sampleNumber)); % hanning window

fLow = 0.01; % set frequency for filtering (high-pass)

%-------------------- E/W ---------------------------%

Evbb_low = removeLow(Evbb,srVBB,fLow);

[~,f,t,pE] = (spectrogram(Evbb_low, w, ...
    floor(sampleNumber/nAvg), ...
    floor(sampleNumber/nAvg), ...
    srVBB,'yaxis'));

meanPowerVBBE_time = timeVBB(1)+seconds(t);

%-------------------- N/S ---------------------------%

Nvbb_low = removeLow(Nvbb,srVBB,fLow);

[~,f,t,pN] = (spectrogram(Nvbb_low, w, ...
    floor(sampleNumber/nAvg), ...
    floor(sampleNumber/nAvg), ...
    srVBB,'yaxis'));

meanPowerVBBN_time = timeVBB(1)+seconds(t);

%-------------------- Z ---------------------------%

Zvbb_low = removeLow(Zvbb,srVBB,fLow);

[~,f,t,pZ] = (spectrogram(Zvbb_low, w, ...
    floor(sampleNumber/nAvg), ...
    floor(sampleNumber/nAvg), ...
    srVBB,'yaxis'));

meanPowerVBBZ_time = timeVBB(1)+seconds(t);

%------------ PRE ---------------------------------%

PRE_low = removeLow(countsPRE,srPRE,fLow);

[~,f,t,pPRE] = (spectrogram(PRE_low, w, ...
    floor(sampleNumber/nAvg), ...
    floor(sampleNumber/nAvg), ...
    srVBB,'yaxis'));

meanPowerPRE_time = timePRE(1)+seconds(t);

%% Plot spectrogram 
figure(1)

%-------------------------------- VBB E/W ------------------------------------%

ax1 = subplot(4,1,1)
colormap jet;
surf(meanPowerVBBE_time,f,log10(sqrt(pE)),'edgecolor','none');
h = colorbar
view(0,90);
caxis([-9.5 -7]);
shading(gca,'flat');
ylabel('Frequency (Hz)')
set(gca, 'YScale', 'log');
ylim([0.1 srVBB/2])
ylabel('\textbf{Acceleration (m/$s^2$)}','FontWeight', 'bold','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
ylabel(h, '{Accel log(m/$s^2$/$\sqrt(Hz)$}','FontSize', 10, 'FontWeight', 'bold','interpreter','latex')
title('\textbf{VBB E/W}','interpreter','latex','FontWeight', 'bold') 
box on

%-------------------------------- VBB N/S ------------------------------------%

ax2 = subplot(4,1,2)
colormap jet;
surf(meanPowerVBBN_time,f,log10(sqrt(pN)),'edgecolor','none');
h = colorbar
view(0,90);
caxis([-9.5 -7]);
shading(gca,'flat');
ylabel('Frequency (Hz)')
set(gca, 'YScale', 'log');
ylim([0.1 srVBB/2])
ylabel('\textbf{Acceleration (m/$s^2$)}','FontWeight', 'bold','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
ylabel(h, '{Accel log(m/$s^2$/$\sqrt(Hz)$}','FontSize', 10, 'FontWeight', 'bold','interpreter','latex')
title('\textbf{VBB N/S}','interpreter','latex','FontWeight', 'bold') 
box on

%-------------------------------- VBB Z ------------------------------------%

ax3 = subplot(4,1,3)
colormap jet;
surf(meanPowerVBBZ_time,f,log10(sqrt(pZ)),'edgecolor','none');
h = colorbar
view(0,90);
caxis([-9.5 -7]);
shading(gca,'flat');
ylabel('Frequency (Hz)')
set(gca, 'YScale', 'log');
ylim([0.1 srVBB/2])
ylabel('\textbf{Acceleration (m/$s^2$)}','FontWeight', 'bold','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
ylabel(h, '{Accel log(m/$s^2$/$\sqrt(Hz)$}','FontSize', 10, 'FontWeight', 'bold','interpreter','latex')
title('\textbf{VBB Z}','interpreter','latex','FontWeight', 'bold') 
box on

%-------------------------------- PRESSURE ------------------------------------%
ax4 = subplot(4,1,4)
colormap jet;
surf(meanPowerPRE_time+seconds(t),f,log10(sqrt(pPRE)),'edgecolor','none');
h = colorbar
view(0,90);
caxis([-3 -1.6])
shading(gca,'flat');
ylabel('Frequency (Hz)')
set(gca, 'YScale', 'log');
ylim([0.1 srPRE/2])
ylabel('\textbf{Pressure (Pa)}','FontWeight', 'bold','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
ylabel(h, '{Pressure log(Pa/$\sqrt(Hz)$}','FontSize', 10, 'FontWeight', 'bold','interpreter','latex')
title('\textbf{Pressure}','interpreter','latex','FontWeight', 'bold') 
box on

load turbo
colormap(turbo)

linkaxes([ax1, ax2, ax3, ax4],'x');
xlim([meanPowerVBBZ_time(1) meanPowerVBBZ_time(end)])