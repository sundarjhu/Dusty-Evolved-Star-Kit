
Dusty-Evolved-Star-Kit
======================

<img src="docs/gold.jpg"  width="170" height="100">

SED-fitting python scripts for fitting data from evolved stars (photometry or spectra) with DUSTY models. Package contains scripts for:
1. Running DUSTY in a multiprocessing batch mode
2. Converting the output from DUSTY to two fits files
3. Least square fitting of the models to data
4. Plotting the results

**Input**: csv with first column wavelength in um and second column flux in Jy (File can have other columns).

**Output**: Results file with best fit model, as well as a results file with specifics for plotting the output. 

**Options**: In the sed_fitting.py you can specify:
 * The model grid
 * distance (in kpc)
 * the wavelength range to fit
 * normalizations range to try
 * the number of values in that normalization range

Several grids are in the models directory (change using the model_grid variable), but you can also create your own model grid. To to this: 

1. Run dusty
2. Put all outputs (spectra files .s* and output files *.out) into a directory of the same name (see example grid directories)
3. Run the dusty_to_grid.py script

This will create two fits files containing all spectra (*directoryname*_models.fits), and all outputs (*directoryname*_outputs.fits).

USING DESK
==========

All of the important script files can be found in "/gold-fit/python_scripts"

Just add the csv data files you want to fit to the *put_target_data_here* directory and run sed_fitting.py.

<img src="gold-fit/output_seds.png"  width="400" height="500">

Under development
-----------------

License
-------

This project is Copyright (c) Steven Goldman and licensed under
the terms of the BSD 3-Clause license. This package is based upon
the `Astropy package template <https://github.com/astropy/package-template>`_
which is licensed under the BSD 3-clause licence. See the licenses folder for
more information.
