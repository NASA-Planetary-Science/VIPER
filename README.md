# VIPER : Voyager Ion PLS Experiment Response

This is the IDL code that was used to fit plasma parameters (ions & electrons; density, flow, temperature) to the Voyager PLS data at Uranus and Neptune.  The fitted plasma parameters, and raw data, may be found on the PDS in the "**vg-pls-reprocess-bagenal2021pdart**" bundle.

This archive has the doi [10.5281/zenodo.22802843](https://doi.org/10.5281/zenodo.22802843) (via Zenodo), please cite with this doi if you use this software in your work.


## Instructions to run VIPER for Uranus and Neptune PLS data ##

Opening and running `BATCH_URANUS.pro` and `BATCH_NEPTUNE.pro` in IDL will
recreate the fits that were turned into the PDS format. It will also
calculate 1-sigma uncertainties and produce plots in user-specified directories on
their local machine.

Depending how paths are set up on your IDL, you may have to copy the contents of the `CSV_Files` directory to your working directory first, or edit the batch pro files' `Input_Filename` keyword of each command to include the path to each CSV.

In `Data/{Planet}/Fitting_tables`, you will find all possible spectra in the same
format as the files that `BATCH_URANUS/NEPTUNE.pro` calls on. If there are other
regions of interest to fit, or if you believe that some spectra have been omitted from
our analysis, they can be included by the user.

## VIPER for Jupiter PLS data ##

Jupiter analysis was completed in 2017 and led to three publications, and used an older format of the VIPER
codebase.  A copy of that code is provided in the **Legacy_VIPER_at_Jupiter_Code** directory here, for completeness of a VIPER repository.
However we expect people to only use the newer (i.e. as updated for Uranus and Neptune) code for any further work.

## VIPER for Saturn PLS data ##

The VIPER code will need some updates to work with Saturn data, it was not set up or tested to run on Saturn PLS data since Cassini provided nearly 300 orbits of coverage.

## VIPER analysis on records of your choice ##

To run your own fits (instead of the batch jobs above), you need to make an input CSV file, with one time record to process listed per row.<br>
We suggest looking at the `BATCH_URANUS/NEPTUNE.pro` files to find an existing CSV for the data mode and planet of your interest, saving that to a new name ('my.csv'), then in your copy delete the rows you do not want, and do any other editing (like number of plasma species or plasma species m/q) to that CSV file.  

If you just want to process the whole input CSV, then (after compiling all) you'd run this one line (the 'u' is for Uranus, or use 'n' for Neptune, or 'j' for Jupiter):<br> 
```
; Example Fits to all rows of data in my.csv
VIPER, 'u', Input_Filename = 'my.csv', OUTPUT_FILENAME='myFits.csv', POD='/my_output_files', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE
```
If you only want to run one record from the CSV, e.g. row 50, then:<br>
```
; Example Fit to fit a specific row of data
; CSV_Indices is what is displayed in Excel, with the header row = 1, first row of data = 2, etc.
VIPER, 'u', Input_Filename = 'my.csv', OUTPUT_FILENAME='myFits.csv', POD='/my_output_files', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, CSV_INDICES = [50]
```

## Final Note ##

Although the I in VIPER is for Ion, this version will also work on PLS electrons.
