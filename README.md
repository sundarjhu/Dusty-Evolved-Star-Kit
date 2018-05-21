<img src="https://3ff009b6523c0c1f8b94-091582592a5f780b4bac3b68414d35fd.ssl.cf5.rackcdn.com/default/_superImage/President-Hayes-Desk.jpg" width="200">

#DESK#

SED-fitting python script for fitting data from evolved stars (photometry or spectra) with DUSTY models. Package contains the fitting script, a simple plotting script, and a script for converting the output from DUSTY to two fits files. 

**Input**: csv with first column wavelength in um and second column flux in Jy (File can have other columns).

**Output**: Results file with best fit model, as well as a results file with specifics for plotting the output. 

**Options**: In the sed_fitting.py you can specify:
	The model grid
	distance (in kpc)
	the wavelength range to fit
	normalizations range to try 
	the number of values in that normalization range

Several grids are in the models directory (change using the model_grid variable), but you can also create your own model grid.

1. Run dusty
2. Put all outputs (spectra files .s* and output files *.out) into a directory of the same name (see example grid directories)
3. Run the dusty_to_grid.py script

This will create to fits files containing all spectra (*directoryname*_models.fits), and all outputs (*directoryname*_outputs.fits)
