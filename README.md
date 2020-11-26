# Comodulation_InSight

Processing software for "A Comodulation Analysis of Atmospheric Energy Injection into the Ground Motion at InSight, Mars" Charalambous, C., Stott, A. E., Pike, T., McClean, J., Warren, T., Spiga, A., ... & Navarro, S. (2020).  Earth Space Sci. Open Arch., 48 (2020), 10.1002/essoar.10503206.1

The code provided here can be run by obtaining the sample datasets as provided in [doi:10.18715/IPGP.2020.KHW87ROO](https://doi.org/10.18715/IPGP.2020.KHW87ROO). These event files have been pre-processed in order to synchronize the timestamps of pressure, wind and ground acceleration and remove data gaps. All raw data and datasets used in the paper are available for download either at the IRIS Data Management Center, or through the [Mars SEIS Data service MSDS](https://www.seis-insight.eu/en/science/science-summary "SEIS InSight Homepage") doi:10.18715/SEIS.INSIGHT.XB_2016

The authoritative source for the archived InSight data, including all instruments is the Planetary Data System. SEIS data is archived at the [Geosciences node](https://pds-geosciences.wustl.edu/missions/insight/index.htm). Wind and pressure data from APSS sensors is archived at the [Atmospheres node](https://atmos.nmsu.edu/data_and_services/atmospheres_data/INSIGHT/insight.html).

THe following MATLAB codes are included for producing the main figures in the paper. The code can easily be adapted to produce other figures given the parameters suggested in the paper, and allows further experimentation by obtaining data from the repositories suggested above.

Figure1.m - Produces the spectrograms of Figure 1 (data set required from [doi:10.18715/IPGP.2020.KHW87ROO](https://doi.org/10.18715/IPGP.2020.KHW87ROO)). This code can be adapted to produce the spectrograms of Figure 3, 4, 9, 10, 11 & 14. This code requires the removeLow.m function to run.

Figure3b.m - Produces moment-matched time-series and recreates Figure 3b (data set required from [doi:10.18715/IPGP.2020.KHW87ROO](https://doi.org/10.18715/IPGP.2020.KHW87ROO)). This code can be adapted to reproduce Figures 1c, 2b & 14 in the paper which incorporate the approach of the method of moments to match time series.

Figure3c.m - Produces comodulation relationships over various frequency bands as shown in Figure 3c (data set required from [doi:10.18715/IPGP.2020.KHW87ROO](https://doi.org/10.18715/IPGP.2020.KHW87ROO)). This code can be adapted to recreate the relationships in Figures 2c,d 4 & 5.

Figure8.m - This code produces the relationships of wind speed predictions, mean and variance as shown in Figure 8.

SNRAnalysiscode.m - This code calculates and plots the environmental SNR given example datasets for Marsquake events to reproduce Figures 11 & 12 of the paper (sample data set required from [doi:10.18715/IPGP.2020.KHW87ROO](https://doi.org/10.18715/IPGP.2020.KHW87ROO)). By obtaining data from doi:10.18715/SEIS.INSIGHT.XB_2016, the code can be adapted to obtain environmental SNRs for all Marsquake events as listed in the Marsquake Catalogue V2 doi:10.12686/a7. This code requires the getSpectrogram.m function to run.

All software in this repository is licensed using GNU General Public License version 3.

Affiliation: Imperial College London, UK. This research is in support of InSight Contribution Number 170.

Colormaps used in MATLAB obtained from: Stephen Cobeldick (2020). MatPlotLib Perceptually Uniform Colormaps (https://www.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps), MATLAB Central File Exchange. Retrieved November 24, 2020.
