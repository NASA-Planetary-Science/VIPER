# Legacy VIPER code and documentation for Jupiter Analysis

The documents and zip files contained in the directory are for the legacy VIPER code for PLS's Jupiter analysis as written and used by the original student team.  The results from which lead to these three papers:
- Bagenal, F., L. P. Dougherty, K. M. Bodisch, J. D. Richardson, and J. M. Belcher (2017), Survey of Voyager plasma science ions at Jupiter: 1. Analysis method, J. Geophys. Res. Space Physics, 122, 8241–8256, doi:[10.1002/2016JA023797](https://doi.org/10.1002/2016JA023797).
- Dougherty L. P., K. M. Bodisch, and F. Bagenal (2017), Survey of Voyager plasma science ions at Jupiter: 2. Heavy ions, J. Geophys. Res. Space Physics, 122, 8257–8276, doi:[10.1002/2017JA024053](https://doi.org/10.1002/2017JA024053).
- Bodisch, K. M., L. P. Dougherty, and F. Bagenal (2017), Survey of Voyager plasma science ions at Jupiter: 3. Protons and minor ions, J. Geophys. Res. Space Physics, 122, 8277–8294, doi:[10.1002/2017JA024148](https://doi.org/10.1002/2017JA024148).

This readme is a copy of instructions the students put on a webpage at the time, with minor formatting updates.<BR>
It is in the IDL language, and paths are given assuming the user is on a Mac or Linux machine (i.e. path separator of /, not \\).<br>
*This code does not work for any planet but Jupiter.  For other planets, use the newer code that is one directory back.*

All the files in this directory were extracted from [https://lasp.colorado.edu/mop/missions/voyager/viper/](https://lasp.colorado.edu/mop/missions/voyager/viper/) with only a little tidying up and removal of some unwanted files.


## VIPER – Voyager Ion PLS Experiment Response

### In this directory:
- VIPER_CODE.zip, Zip file of the VIPER code used at the time for Jupiter (only) analysis
- VOYAGER1.zip, Zip file of the PLS data used as input for Voyager 1
- VOYAGER2.zip, Zip file of the PLS data used as input for Voyager 2
- Voyager PLS Analysis Manual (Document, from 2017)
- Voyager PLS Error Analysis (Document, likely from 2017 too)

They are input data files for individual days.<BR>
Code to unpack these raw data files is found in the VIPER code.<BR>
*The data files in VOYAGER1.zip and VOYAGER2.zip are a different format to that used in later work.*

### Instructions for *(this version of)* VIPER:


1. Download and unzip all necessary folders (VIPER_CODE.zip, VOYAGER1.zip, VOYAGER2.zip)

2. Place the contents of all unzipped folders into one folder on your machine.

3. Once the downloaded files have been unzipped and saved to one folder on your local machine, Open IDL.

4. Upon opening IDL, the user must CHANGE CURRENT DIRECTORY. This can be done two ways. If option A does not work, use option B.

    (A) Type into the IDL command line :    CD,  'path/to/folder' <br>
   e.g.:    CD, '/Users/username/Desktop/Voyager_PLS/VIPER'

    (B) Above the IDL Command line window there is a yellow folder button, click this button and you will then be able to look through the folders on the machine and select the appropriate working directory.

5. You should now be in the correct working directory and be able to run the code.

    To run the code type into the command line: @VIPER

6. VIPER will now compile necessary procedures and begin to run the program.

7. A CSV file will  now automatically open in your computer's default CSV editor (Excel or similar spreadsheet programs provide the most clean output).

8. Now edit the CSV file to change the time or manipulate parameters and save it.

    Press "Continue" if, when saving, the CSV a prompt appears about the format containing special features. The IDL code will work regardless.

9. The command prompt will now tell you 'type any key', now hit a button on the keyboard, with the cursor on the command line. Not all keys will work such as shift. Alpha-numeric keys will always work though.

10. The program will now run and plot the output of the response function along with the original data.

### A second set of instructions... ###
*We found this second set of instructions from 2017.  Including for completeness, even thought not all directories exist anymore.*

For Jupiter:
1. Download "VIPER_CODE.zip"
2. Unzip and place in directory of choice
3. From the IDL window run the following command<br>
    VIPER, 'j', /INITIAL <br>
    This compiles all necessary functions and runs a test fit. It should complete within a few minutes. You can stop the code once programs are compiled, if desired.
4. Open "batch_fits.pro"<BR>
  Has examples of different ways to use the code given the CSV files of the data

Notes:
- The code currently works with data from Jupiter only.
- Plots of fits are output into "/Output/Plots"
- "Data/Data_Preparation" contains conversions of SPL files from Digital Numbers (0-255) into currents (fA) for data at planets other than Jupiter.
- "Data/Jupiter/PLS_Data" and "Data/Jupiter/Trajectories" contain the following:<BR>
    "Day__.txt" files that are an older format of the SPL file that existed before VIPER was created. We only have these files for Jupiter.<BR>
    ".ssedr.rob" and "SpiceRotations" that give spacecraft position and pointing for every time measurement.  These are not in the SPL files.
- "/CSV_Files" includes all of the files that we use for fitting. The files can be named anything so long as they are passed into the "VIPER, 'j', INPUT_FILENAME = '_____' " at the command line
- The bulk of the code lies in “VOYAGER_PLS.pro” and “ADDITIONAL_PLS_FUNCTIONS.pro”
- There is older documentation from February 10th, 2017. It is not the most up to date, but describes more of the inner workings of the code as needed.  *That is the **VoyagerPLSAnalysisManual.pdf** (in this GitHub Directory).*
