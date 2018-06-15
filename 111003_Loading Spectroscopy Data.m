clear all


Spectroscopy_data=importdata('D:\Spectroscopy data - Green.txt');

wavelength_Spectroscopy_micron=Spectroscopy_data(:,1)/1000;

power_Spectroscopy=Spectroscopy_data(:,2);

power_spectroscopy_interp=interp1(wavelength_Spectroscopy_micron,power_Spectroscopy,wavelength_micron);