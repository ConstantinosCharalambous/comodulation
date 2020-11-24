# Comodulation_InSight

Processing software for  "A Comodulation Analysis of Atmospheric Energy Injection into the Ground Motion at InSight, Mars" Charalambous, C., Stott, A. E., Pike, T., McClean, J., Warren, T., Spiga, A., ... & Navarro López, S. (2020). A Comodulation Analysis of Atmospheric Energy Injection into the Ground Motion at InSight, Mars. Earth
1784 Space Sci. Open Arch., 48 (2020), 10.1002/essoar.10503206.1

The code provided here can be run by obtaining the sample datasets as provided in seis-insight... These event files have been pre-processed in order to synchronize the timestamps of pressure, wind and ground acceleration and remove data gaps. All raw data and datasets used in the paper is available for download either at the IRIS Data Management Center, or through the [Mars SEIS Data service MSDS](https://www.seis-insight.eu/en/science/science-summary "SEIS InSight Homepage") doi:10.18715/SEIS.INSIGHT.XB_2016

The authoritative source for the archived InSight data, including all instruments is the Planetary Data System. SEIS data is archived at the [Geosciences node](https://pds-geosciences.wustl.edu/missions/insight/index.htm). Wind and pressure data from APSS sensors is archived at the [Atmospheres node](https://atmos.nmsu.edu/data_and_services/atmospheres_data/INSIGHT/insight.html).

THe following MATLAB codes are included for producing the main figures in the paper. The code can easily be adapted to produce the rest of the figures given the parameters suggested in the paper:

Figure1a.m - Panel A of figure 2 (Note that the data to run the code are obtained from doi:)

Figure3.m - Panel B of figure 2 (Note that the data to run the code are obtained from doi:)

Figure8.m - Figure 3 (Note that the data to run the code are obtained from doi:)


The specific data used in these codes can also be downloaded at the NASA Open Data Portal archive preserved at this archive. This data should be saved accordingly to the variable names as provided by the data in , in order to use the MATLAB codes as written.

All software in this repository is licensed using GNU General Public License version 3.

Affiliation: Imperial College London, UK. This research is in support of InSight Contribution Number 170.

Colormaps used in MATLAB obtained from: Stephen Cobeldick (2020). MatPlotLib Perceptually Uniform Colormaps (https://www.mathworks.com/matlabcentral/fileexchange/62729-matplotlib-perceptually-uniform-colormaps), MATLAB Central File Exchange. Retrieved November 24, 2020.
