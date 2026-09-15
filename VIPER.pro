PRO VIPER, PlanetName, MAKE_PLOT=MAKE_PLOT, UNCERTAINTIES=UNCERTAINTIES, WRITE_UNCERTAINTIES_IN_FILE=WRITE_UNCERTAINTIES_IN_FILE, $
           CALC_FITTING_RESIDUALS = CALC_FITTING_RESIDUALS, INPUT_FILENAME=INPUT_FILENAME, CSV_INDICES=CSV_INDICES, $
           FULL_CSV=FULL_CSV, NO_CSV_RETURN=NO_CSV_RETURN, OUTPUT_FILENAME=OUTPUT_FILENAME, CLEAR=CLEAR, FORCE_SAVE_JPG=FORCE_SAVE_JPG, $
           NO_FIT=NO_FIT, FIT=FIT, CASSINI_COMPOSITION=CASSINI_COMPOSITION, COMPUTE_CHI_SPACE=COMPUTE_CHI_SPACE, N_CHI_PTS=N_CHI_PTS, $
           N_CHI_SIG=N_CHI_SIG, INPUT_UNCERTAINTY_FILENAME=INPUT_UNCERTAINTY_FILENAME, PLOT_QUIETLY=PLOT_QUIETLY, $
           DELAMERE_COMPOSITION=DELAMERE_COMPOSITION, NO_PLOT_SAVE=NO_PLOT_SAVE, POD=POD, ANALYTICAL_VPHI=ANALYTICAL_VPHI, $
           CUPA=CUPA, CUPB=CUPB, CUPC=CUPC, CUPD=CUPD, ELECTRON_DATA_PRODUCT=ELECTRON_DATA_PRODUCT, RETURN_ELECTRON_POWERLAW=RETURN_ELECTRON_POWERLAW, $
           REMOVE_BAD_KAPPA=REMOVE_BAD_KAPPA, INITIAL=INITIAL, CENTRIFUGAL_EQUATOR=CENTRIFUGAL_EQUATOR
           
;
;     
;     
; Measure time of program
TIC                                                                                             
; Voyager Ion PLS Experiment Response
;
; DESCRIPTION:
; - Analyzes data from the Voyager PLS instrument
; Currently focuses on data from the Jupiter epoch, but the code is being written so that it can recognize multiple planetary systems or the solar wind
; - Nov 2023: Uranus added as option. To be combined with GX's notes on how to prepare the data
; 
; 
;
; INPUTS
; PlanetName                 - String that determines what planet the spacecraft is around (see the strings below, at 'Check which planet it is' section)
; 
; Keyword Arugments
; MAKE_PLOT                  - Keyword argument that plots the data. Additional argument in the CSV, or keyword below, to save the plot automatically
; UNCERTAINTIES              - Keyword argument that when set will calculate the uncertainties of the parameters that are being fit
; WRITE_UNCERTAINTIES_IN_FILE- Keyword argument that saves the uncertainties in the same file. SEE NOTES BELOW
; CALC_FITTING_RESIDUALS     - Keyword argument that calculates the residuals between fitting and raw data, and saves them on a file
; INPUT_FILENAME             - Keyword arugment for the input file that should be used for running the code
; CSV_INDICES                - Keyword argument that will tell the program what rows of the CSV file to use. Note that CSVs starts on row 2
; FULL_CSV                   - Keyword argument to process the full csv. This is the default option if the above isn't specified
; NO_CSV_RETURN              - Keyword argument to prevent the return of a csv file - Useful for just plotting the data
; OUTPUT_FILENAME            - Keyword argument to specify the output name of the csv file being returned
; CLEAR                      - Keyword argument that will clear any existing graphic windows
; FORCE_SAVE_JPG             - Keyword argument that overrides what is written in the CSV and automatically saves a JPG of the plot. May not plot on screen
; NO_FIT                     - Keyword argument to explicitly turn off fitting procedure - Useful for some batch jobs to manipulate CSV files from console
; FIT                        - Keyword argument to explicitly turn on fitting procedure
; CASSINI_COMPOSITION        - Keyword argument to use the CASSINI UVIS observations for constraints on the composition of ion species. Obtained from analysis by Edward Nerney
; INPUT_UNCERTAINTY_FILENAME - Keyword argument with name of uncertainty file. For use with the chi space computation since it needs the best fit value and the 1 sigma value.
; COMPUTE_CHI_SPACE          - Keyword argument to calculate chi space for the given parameters (this must be run with fitting enabled or have a input_uncertainty_file provided)
; N_CHI_PTS                  - Keyword argument for number of points to use (Creates a square array of that size)
; N_CHI_SIG                  - Keyword argument for number of sigma to deviate from best fit
;                              If no name is entered, then it will assume the autosaved filename for uncertainties: 'filename_Unc.csv'
; PLOT_QUIETLY               - Keyword argument that produces and saves the plots but does so without popping them up on the computer screen
;                              Does this by adding 100 to the Save_Plot structure tag that the code then understands (there are some flaws in this if you do it more than once.)
; DELAMERE_COMPOSITION       - Keyword argument to use the composition from Delamere's 2004 Paper. Roughly 30 points were extrapolated from his plots of mixing ratios and then interpolated to give our tables.
;                              These values are in Data/Jupiter/Chemistry_Ratios: 'Delamere_Ratios.txt'. Columns are in the following order: Rj, nO+/nO+, nO++/nO+, nS+++/nO+, nS++/nO+, nS+/nO+
;                              The chemistry model is valid for 6-9 Rj. Before this point, the calculation has not been done. Past this point, it is assumed that mixing has ;stopped and that the density
;                              ratios are the same as the values at 9 Rj.
; NO_PLOT_SAVE               - Keyword argument to turn off plot saving
; POD                        - Keyword argument to specify Plot Output Directory (POD). [GX: Give the path AFTER the working folder, e.g. '/Uranus/Output_plots' ]
; ANALYTICAL_VPHI            - Keyword argument that specifies an analytical equation to use for the Vphi profile. It defaults to a 1, which is a hyperbolic tangent model
; CUPA (B,C,D)               - Keyword argument which allows you to select individual cup plots to be displayed in output
; Centrifugal Equator        - Not tested and from a old code import - Does not seem to work 2020_08_20 - The CEq file may be returning density*density somehow
;
;
;
; OUTPUTS
; CSV_File                   - Given the date of run in UTC time as the name if OUTPUT_FILENAME keyword argument is not provided. Must include the .csv extention
;
;
; -----------
; EXAMPLE INPUT: Viper, 'Jup', input_filename  = 'V1_JUP_M_Final.csv', csv_indices = 5, /make_plot,/no_fit
;                     -this runs viper using the V1_JUP_M_Final csv and plots a single set of A,B,C and D cup plots without fitting the data
;                Viper, 'Jup', input_filename  = 'V1_JUP_M_Final.csv', csv_indices = 5, /compute_chi_space       
;                     -this computes the chi space for the fit parameters as defined by the V1_JUP_M_Final csv in this process, the parmeters are fit and uncertainties are calculated
;                     note : must provide uncertainty file if you do not want to go through the fitting procedure.
;                VIPER, 'u', /INITIAL, Input_Filename = '/Users/username/Documents/Uranus/PLS/VIPER/VIPER_Code/Data/Uranus/PLS_Data/test_2.csv', POD = '/Data/Uranus/TEST_OUTPUT', /MAKE_PLOT, OUTPUT_FILENAME = 'test_2_OUTPUT.csv'
;
; ** For Windows: as the path dividers are different than Mac, you need to run the command as:
;    VIPER, 'u', /INITIAL, Input_Filename = 'C:\Users\username\Documents\Uranus\PLS\VIPER\VIPER_Code\Data\Uranus\PLS_Data\test_2_mini.csv', POD = '\Data\Uranus\TEST_OUTPUT', /MAKE_PLOT, OUTPUT_FILENAME = 'test_2_OUTPUT.csv', /PLOT_QUIETLY
;    No need to change anything else.
; -----------
; Notes:
;   - Compiling and running viper will compile all necessary files
;   - WRITE_UNCERTAINTIES_IN_FILE works if there are uncertainty columns in the initial csv file (currently they are named "SPECIES_#_N_UNCERT" and "SPECIES_#_T_UNCERT").
;   It only works with /UNCERTAINTIES. It is programmed for Uranus, so for Jupiter you'll need to add columns, or do it the old way: first run VIPER to find the density 
;   and temperature, then give them as inputs and  rerun VIPER with the argument /UNCERTAINTIES
; 
; -------------------------------------
; Written by Kaleb Bodisch and Logan Dougherty
; Modified for Uranus and Neptune by George Xystouris

; ---------------------- GX: ADD STUFF FOR OTHER PLANETS HERE!!!
CONRD
;MODULE_NAMES     = ['CONRD','ADDITIONAL_PLS_FUNCTIONS','VOYAGER_PLS_JUPITER','VOYAGER_PLS_URANUS','VOYAGER_PLS_NEPTUNE','PLS_UNCERTAINTIES','mpfitfun','DATA_PLOTS','Field_Lines_VIPER']
;MAG_MODULE_NAMES = ['dPdtheta','getLegendre','getSchmidtCoefs','kk_2009','latdistVoy','footpoint_logan','latdist_max']
;RESOLVE_ROUTINE, MODULE_NAMES,     /COMPILE_FULL_FILE, /EITHER
;RESOLVE_ROUTINE, MAG_MODULE_NAMES, /COMPILE_FULL_FILE, /EITHER

; RJW shortened this list
MODULE_NAMES     = ['CONRD','ADDITIONAL_PLS_FUNCTIONS','VOYAGER_PLS_JUPITER','VOYAGER_PLS_URANUS','VOYAGER_PLS_NEPTUNE','PLS_UNCERTAINTIES','mpfitfun','DATA_PLOTS',$
'latdistVoy','footpoint_logan' $ ; These two needed ?
]
RESOLVE_ROUTINE, MODULE_NAMES,     /COMPILE_FULL_FILE, /EITHER
; ------------------------ TO HERE

; Set PlanetName to be all upper case letters to determine the planet
PlanetName = STRUPCASE(PlanetName)

; Check which planet it is
JupiterStrings = ['J','JUP','JUPITER']
UranusStrings = ['U','URA','URAN','URANUS']
NeptuneStrings = ['N','NEP','NEPT','NEPTUNE']


; ------------------------------------------------
; STRUCTURE INPUT
; ------------------------------------------------
; Year, DOY, Hour, Minute, and Second of spacecraft to match Data
;
; Spacecraft:         1 for Voyager 1, 2 for Voyager 2
; Planet_Number:      5 for Jupiter, 6 for Saturn, 7 for Uranus, 8 for Neptune, and 10 for Solar Wind
; Response:           0 for CUPINT/DCPINT, 1 for LABCUR/LDCUR
; L_or_M_Mode:        0 for L Mode (16 channels per cup), 1 for M Mode
; Short_or_Long:      1 for long integration time, 0 for short integration time
; Threshold_Current:  OUTPUT ONLY! Input value does not matter
; 
; Save_Plot:          0 doesn't save plot, 1 for .jpg extension, 2 for .png extension, 3 for .eps extension
; Fit:                0 to not fit the data, 1 to fit the data using Chi-Square minimization technique
; Iterations:         Number of maximum iterations to use in the fitting procedure
; Channels_CupX:      A 6 digit number (not string!) to tell what channels to fit. First 3 digits are first channel (001). Last 3 digits
;                     are the final channel (128). Can do any range of channels.
;                       e.g. 1016 means channels from 001 to 016 (leading zeros are omitted)
; Total_number_of_fitted_channels:  OUTPUT ONLY! Input value does not matter
; Chi_Square          OUTPUT ONLY! If set it gives the chi square value of the fitting on the output
; 
; Vary_VX:            Flag to vary the velocity in the fitting procedure. 1 to vary the velocity.
;                     For the Jupiter analysis, cylindrical coordinates are used (V1,V2,V3 = Vr,Vphi,Vz)
; Vary Density:       Flag for allowing the densities to vary. 1 to vary.
; Vary_Temperature:   Flag to allow temperatures of the ion species to vary. 1 to vary.
; VX:                 Initial value of the flowspeed in cylindrical coordinates (Vr,Vphi,Vz)
; Common_Temperature: Uses a common temperature (eV) for all ion species. Supercedes temperature given in Species_X_T,
;                     unless Species_X removes this restriction
;
; Species_X: Flag to process ion species as a parameter or not
;            NOTE: In the code, a parameter means a value that is allowed to vary freely by itself and is
;            not coupled to any other species. A fixed composition of one ion species will not be called a
;            parameter but is rather a fixed constant of the model
;
;            0: This species is not in parameter space and is a constant of the model. It will follow the
;               Common_Temperature value if it is greater than 0 and the species temperature will be neglected.
;            1: This species is not in parameter space and does not follow the Common_Temperature value, but follows
;               its own temperature value.
;            2: This species is in parameter space with its density and temperature as a parameter.
;               If Common_Temperature is greater than 0, then the species temperature will be ignored.
;            3: This species is in parameter space but does not follow the Common_Temperature value, even if
;               it is specified above 0.
;            4: This species is taken indirectly into parameter space. This value means that the density of this species
;               is coupled to other species. Only the first species with a 4 will technically be a parameter. Other
;               species with the value of 4 will have that proportion of density and add to the value of current
;               measured by the PLS sensor, but they will not vary themselves. Follows a common temperature if given.
;            5: Same as above but explicitly has its own temperature as a parameter, although the density is tied.
;
;               The first species with a value of 4 will have its density provided in number/cc, while the second
;               species with a value of 4 will have its density provided as a fraction of the first
; Species_X_A: Mass of ion species (1 for H+, 16 for O+ and O++, etc.)
; Species_X_Z: Charge of ion species (1 for O+, 2 for O++, etc.)
; Species_X_n: Density of ion species, given in number/cc, unless argument of 4 for Species_X, as described above
; Species_X_T: Temperature, in eV, of species. A Common_Temperature will usually supercede this value
; Species_X_n_Unc: Density Uncertainty - if keyword 'WRITE_UNCERTAINTIES_IN_FILE' set [SEE NOTES ABOVE]
; Species_X_T_Unc: Temperature Uncertainty - if keyword 'WRITE_UNCERTAINTIES_IN_FILE' set [SEE NOTES ABOVE]


; ------------------
;     JUPITER
; ------------------
IF TOTAL(STRCMP(PlanetName, JupiterStrings,/FOLD_CASE)) NE 0 THEN BEGIN ; Jupiter case
  Structure = {Year: 1979, DOY: 64, Hour: 10, Minute: 16, Second: 0, Spacecraft: 1, Planet_Number: 5, Response: 0, $
               L_or_M_Mode: 1, Save_Plot: 0,  Fit: 0, Iterations: 25, Channels_CupA: 001128 , Channels_CupB: 001128, $
               Channels_CupC: 001128, Channels_CupD: 001128, Vary_V1: 0, Vary_V2: 0, Vary_V3: 0, Vary_Density: 0, $
               Vary_Temperature: 0, Vary_Kappa: 1, Analytical_Vphi: 0d, V1: 0d, V2: 0d, V3: 0d, Common_Temperature: 0d, Common_Kappa: 0d, Delamere_Composition: 0d, $
               Species_1: 0, Species_1_A: 1d, Species_1_Z: 1d, Species_1_n: 50d, Species_1_T: 3d, Species_1_Kappa: 2d, $
               Species_2: 0, Species_2_A: 1d, Species_2_Z: 1d, Species_2_n: 50d, Species_2_T: 3d, Species_2_Kappa: 2d, $
               Species_3: 0, Species_3_A: 1d, Species_3_Z: 1d, Species_3_n: 50d, Species_3_T: 3d, Species_3_Kappa: 2d, $
               Species_4: 0, Species_4_A: 1d, Species_4_Z: 1d, Species_4_n: 50d, Species_4_T: 3d, Species_4_Kappa: 2d, $
               Species_5: 0, Species_5_A: 1d, Species_5_Z: 1d, Species_5_n: 50d, Species_5_T: 3d, Species_5_Kappa: 2d, $
               Species_6: 0, Species_6_A: 1d, Species_6_Z: 1d, Species_6_n: 50d, Species_6_T: 3d, Species_6_Kappa: 2d, $
               Species_7: 0, Species_7_A: 1d, Species_7_Z: 1d, Species_7_n: 50d, Species_7_T: 3d, Species_7_Kappa: 2d}
  Structure_Names = TAG_NAMES(Structure)
  
  ; Check if there is a specified input file and read data from that to fill the structure and then proceed to analysis.
  ; If no file exists, spawn a CSV file that can be directly edited by the user before continuing.

  IF ISA(Input_Filename) EQ 1 THEN BEGIN
    
    ; Change directory to where CSVs are located
    CD, CURRENT = RESTORE_DIR
    STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
    CD, STR_DIR
    Input_CSV = READ_CSV(Input_Filename, HEADER = Structure_Names)
    CD, RESTORE_DIR
    
    ; If Full CSV keyword is set, use every row of CSV
    IF KEYWORD_SET(FULL_CSV) THEN CSV_Indices = LINDGEN(N_ELEMENTS(Input_CSV.Field01)) ELSE $
    ; If there is not a specific range, default to using every row, otherwise use the selected inputs
    IF ISA(CSV_Indices) EQ 0 THEN CSV_Indices = LINDGEN(N_ELEMENTS(Input_CSV.Field01)) ELSE $
                                  CSV_Indices = CSV_Indices - 2 ; -2 to go from CSV indexing to IDL indexing
    ; Set keyword clear on if quiet plotting is turned on so that it will clear. Can't see or access files anyways, so clear to free RAM 
    IF KEYWORD_SET(PLOT_QUIETLY) THEN CLEAR = 1               
    nloop = 0 ; Counter for row number of CSV file
    
    IF KEYWORD_SET(ELECTRON_DATA_PRODUCT) THEN EDP = 1 ELSE EDP = 0
    IF KEYWORD_SET(RETURN_ELECTRON_POWERLAW) THEN EPL = 1 ELSE EPL = 0
    
    
    FOR i = 0, N_ELEMENTS(CSV_Indices) - 1 DO BEGIN                       ;start loop
      nloop = nloop + 1 ; Augment loop for next row of CSV file
      CSV_Index = CSV_Indices(i)
      
      ; Clear any existing plot
      IF KEYWORD_SET(CLEAR) THEN BEGIN
        wwwwww = GETWINDOWS()
        ; If plots exist, clear them
        IF ISA(wwwwww) EQ 1 THEN FOR wwwwwwi = 0, N_ELEMENTS(wwwwww)-1 DO wwwwww[wwwwwwi].close
      ENDIF

      ; Read in CSV to proper structure format with correct tag names
      TempArr = MAKE_ARRAY(N_TAGS(Input_CSV)) ; Create temporary array to hold structure information
      FOR j = 0, N_TAGS(Input_CSV)-1 DO TempArr[j] = (Input_CSV.(j))[CSV_Index] ; Fill in temporary array
      ; Define Structure in correct format for VIPER code to use
      New_Struct = {col1: Structure_Names, col2: TempArr}
      struct_elements = N_ELEMENTS(New_Struct.col2)
      TempStruct = CREATE_STRUCT(New_Struct.col1[0], DOUBLE(New_Struct.col2[0]))
      FOR k = 1, N_ELEMENTS(New_Struct.col1) - 1 DO BEGIN
        TempStruct = CREATE_STRUCT(TempStruct, New_Struct.col1[k], DOUBLE(New_Struct.col2[k]))
      ENDFOR
      Structure = TempStruct
      
      ; Need to program in a compositional table call
      
      FORWARD_FUNCTION CASSINI_RATIOS
      IF KEYWORD_SET(CASSINI_COMPOSITION) THEN Structure = CASSINI_RATIOS(Structure, ISA(UNCERTAINTIES))
      
      FORWARD_FUNCTION DELAMERE_RATIOS
      ; Want to insert Delamere Ratios as a tag in the structure
      ; So try this:
      IF KEYWORD_SET(DELAMERE_COMPOSITION) THEN Structure.Delamere_Composition = 1
      IF Structure.Delamere_Composition EQ 1 THEN Structure = DELAMERE_RATIOS(Structure, ISA(UNCERTAINTIES))
      ;IF KEYWORD_SET(DELAMERE_COMPOSITION) THEN Structure = DELAMERE_RATIOS(Structure, ISA(UNCERTAINTIES))
      
      IF KEYWORD_SET(ANALYTICAL_VPHI) THEN Structure.Analytical_Vphi = 1d
      IF Structure.Analytical_Vphi GE 1 THEN BEGIN
        FORWARD_FUNCTION VGR_CupVelocity
        FORWARD_FUNCTION VPHI_PROFILE
        Rj_Distance       = VGR_CupVelocity([0d,0d,0d], 'A', Structure, /RETURN_RJ)
        Vphi_Fixed        = VPHI_PROFILE(Rj_Distance, Structure.Analytical_Vphi)
        Structure.V2      = Vphi_Fixed ; Set appropriate value for Vphi
        Structure.Vary_V2 = 0d ; Forces the phi component of velocity to be a constant at that value
      ENDIF
       
      IF KEYWORD_SET(FORCE_SAVE_JPG) THEN BEGIN
        Structure.Save_Plot = 1
        IF ISA(MAKE_PLOT) EQ 0 THEN MAKE_PLOT=1 ; Turn on keyword argument for plotting in this case if it wasn't already on
      ENDIF
      
      IF KEYWORD_SET(POD) THEN BEGIN
        ; Check if the directory exists
        ;CD, CURRENT=RESTORE_DIR
        ;Plot_Out_Dir = STRCOMPRESS(STRING(RESTORE_DIR, POD))
        ;DIR_CHECK = FILE_TEST(Plot_Out_Dir, /DIRECTORY)
        ;IF DIR_CHECK EQ 0 THEN FILE_MKDIR, Plot_Out_Dir
        
        ; Logan 2025_02_27
        CD, CURRENT=RESTORE_DIR
        cutoff_pos = strpos(RESTORE_DIR, 'Documents', /REVERSE_SEARCH)
        IF cutoff_pos NE -1 THEN BEGIN
          Root_User_Path = STRMID(RESTORE_DIR, 0, cutoff_pos + STRLEN('Documents'))
          ; PRINT, 'Root_User_Path' ; Test PRINT command
        ENDIF ELSE BEGIN
          PRINT, 'ERROR in determining Root path'
        ENDELSE
        
        Root_VIPER_Output = STRCOMPRESS(STRING(Root_User_Path, '/VIPER_Output'))
        
        Plot_Out_Dir = STRCOMPRESS(STRING(Root_VIPER_Output, '/Plots', POD))
        ;print, Plot_Out_Dir ; Test Line
        ;Print, POD ; Test Line
        DIR_CHECK = FILE_TEST(Plot_Out_Dir, /DIRECTORY)
        IF DIR_CHECK EQ 0 THEN FILE_MKDIR, Plot_Out_Dir
        ; Logan 2025_02_27
        
      ENDIF ELSE BEGIN
        ; Maintain this ELSE to save to a location in VIPER Git dir if POD not specified
        POD = '/Output/Plots'
        CD, CURRENT=RESTORE_DIR
        Plot_Out_Dir = STRCOMPRESS(STRING(RESTORE_DIR, POD))
      ENDELSE
      
      
      IF KEYWORD_SET(NO_PLOT_SAVE) THEN Structure.Save_Plot = 0
      
      ; Add 100 to Save_Plot argument of structure to add a buffer to the plot and not have it pop up on the screen
      IF KEYWORD_SET(PLOT_QUIETLY) THEN Structure.Save_Plot = Structure.Save_Plot + 100
      
      IF KEYWORD_SET(NO_FIT) THEN Structure.Fit = 0
      IF KEYWORD_SET(FIT) THEN Structure.Fit = 1
      
      cupss = 0
      IF KEYWORD_SET(CUPA) THEN cupss = 1
      IF KEYWORD_SET(CUPB) THEN cupss = 2
      IF KEYWORD_SET(CUPC) THEN cupss = 3
      IF KEYWORD_SET(CUPD) THEN cupss = 4
      
      ; Extrapolate plasma measurements to the Equator
      ; This section of the code does not fit the data, it takes fitted data from where the S/C
      ; made in situ measurements of the plasma and then traces the magnetic field line to determine
      ; the location of the centrifugal equator. Using the code, latdistVoy, the density of the plasma
      ; is then calculated at the location of the centrifugal equator.
      ;
      ; Note: Extensive testing has not been done on this part of the code yet
      IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
        IF KEYWORD_SET(MAKE_PLOT) OR KEYWORD_SET(UNCERTAINTIES) OR KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN
          MESSAGE, 'VIPER code can not calculate centrifugal equator while producing a plot, fitting data, calculating uncertainties, or computing chi-space distributions.'
        ENDIF
      ENDIF
      
      ; Run through model to process the data  
      IF KEYWORD_SET(MAKE_PLOT) AND KEYWORD_SET(UNCERTAINTIES) THEN BEGIN   
        Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
        ;stop
        Return_Structure = VOYAGER_PLS_JUPITER(Structure, Plot_Out_Dir, cupss, EDP, EPL, /MAKE_PLOT, /UNCERTAINTIES)

        ; Run the procedure to calculate uncertainties
      ENDIF ELSE IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN 
        ;stop
        Return_Structure = VOYAGER_PLS_JUPITER(Structure, Plot_Out_Dir, cupss, EDP, EPL, /MAKE_PLOT)

        ; Run the procedure to calculate uncertainties and plot the result
      ENDIF ELSE IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN 
        Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
        ;stop
        Return_Structure = VOYAGER_PLS_JUPITER(Structure, Plot_Out_Dir, cupss, EDP, EPL, /UNCERTAINTIES)

        ; Run the procedure as normal without producing a plot or uncertainties for the parameters - Useful for
        ; running batch fits so that plots don't need to be created and deleted, freeing up RAM,
        ; and it doesn't waste time generating uncertainties in case the batch fit was unsuccessful.
      ENDIF ELSE IF KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN 
        ; Check if there is an input number of points and sigma to use. If not default to 21 points and 3 sigma
        IF ISA(N_CHI_PTS) EQ 0 THEN BEGIN
          N_CHI_PTS = 21
          PRINT, 'No input provided for number of points in Chi space calculation. Defaulting to 21.'
        ENDIF
        IF ISA(N_CHI_SIG) EQ 0 THEN BEGIN
          N_CHI_SIG = 3
          PRINT, 'No input provided for number of sigma in Chi space calculation. Defaulting to 3 sigma.'
        ENDIF
        
        ; Check that there is no input_uncertainty_filename on the first iteration
        IF ((ISA(INPUT_UNCERTAINTY_FILENAME) EQ 0) AND (nloop EQ 1)) THEN UNC_FILE_LOGIC_CHECK = 0
        
        
        IF UNC_FILE_LOGIC_CHECK EQ 0 THEN BEGIN
          ; Do a quick search for the correct name and if a file exists
          strlength = STRLEN(Input_Filename)
          strlength_unc = strlength - 4
          new_string = STRMID(Input_Filename, 0, strlength_unc)
          INPUT_UNCERTAINTY_FILENAME = STRCOMPRESS(STRING(new_string, '_Unc.csv'))
          
          ; Change directory to where CSVs are located
          CD, CURRENT = RESTORE_DIR
          STR_DIR = RESTORE_DIR + PATH_SEP() + 'CSV_Files'
          CD, STR_DIR
          ; Check if the file exists as it should be named
          IF FILE_TEST(INPUT_UNCERTAINTY_FILENAME) EQ 0 THEN MESSAGE, 'Error: No file containing parameter uncertainties was provided or found.'
          CD, RESTORE_DIR
        ENDIF
        
        ; Change directory to where CSVs are located
        CD, CURRENT = RESTORE_DIR
        STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
        CD, STR_DIR
        ; Read in uncertainty filename
        Input_CSV_Unc = READ_CSV(INPUT_UNCERTAINTY_FILENAME, HEADER = Structure_Names)
        CD, RESTORE_DIR
        
        ; Read in CSV to proper structure format with correct tag names
        TempArr = MAKE_ARRAY(N_TAGS(Input_CSV_Unc)) ; Create temporary array to hold structure information
        FOR j = 0, N_TAGS(Input_CSV_Unc)-1 DO TempArr[j] = (Input_CSV_Unc.(j));[CSV_Index] ; Fill in temporary array
        ; Define Structure in correct format for VIPER code to use
        New_Struct = {col1: Structure_Names, col2: TempArr}
        struct_elements = N_ELEMENTS(New_Struct.col2)
        TempStruct = CREATE_STRUCT(New_Struct.col1[0], DOUBLE(New_Struct.col2[0]))
        FOR k = 1, N_ELEMENTS(New_Struct.col1) - 1 DO BEGIN
          TempStruct = CREATE_STRUCT(TempStruct, New_Struct.col1[k], DOUBLE(New_Struct.col2[k]))
        ENDFOR
        Structure_Unc = TempStruct   
        
        ;stop
        Return_Structure = VOYAGER_PLS_JUPITER(Structure, Plot_Out_Dir, Structure_Unc=Structure_Unc, cupss, EDP, EPL, N_CHI_PTS=N_CHI_PTS, N_CHI_SIG=N_CHI_SIG, /COMPUTE_CHI_SPACE)

      ENDIF ELSE IF KEYWORD_SET(REMOVE_BAD_KAPPA) THEN BEGIN
        ; Return just the input (it's only checking for bad values of Kappa, i.e at the lower limit)
        Return_Structure = Structure
      ENDIF ELSE IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
        Return_Structure = VOYAGER_PLS_JUPITER(Structure, Plot_Out_Dir, cupss, /CEq)
      ENDIF ELSE BEGIN
        ;stop
        Return_Structure = VOYAGER_PLS_JUPITER(Structure, Plot_Out_Dir, cupss, EDP, EPL)
      ENDELSE
     
      
      
      ; Save the file - Skip if No Return keyword is set
      IF NOT KEYWORD_SET(NO_CSV_RETURN) THEN BEGIN
        IF nloop EQ 1 THEN BEGIN ; When nloop is 1, it is on the first iteration so it has to handle the saving differently
          OutputStruct = CREATE_STRUCT(Return_Structure)
          Header = TAG_NAMES(Return_Structure)
        ENDIF ELSE BEGIN ; Concactanate the structures into one large array
          RESTORE, 'VIPER_SAVE_FILE.sav' ; Restore save file
          IF KEYWORD_SET(REMOVE_BAD_KAPPA) THEN BEGIN
            ; This is basically a one off feature to parse the data down
            ; This line does not remove the first iteration if it is a 'bad' value. If first value is less than 1.501001 it will require manual deletion (quicker than altering code)
            ; IF Return_Structure.Species_1_Kappa LE 1.501001d THEN OutputStruct = OutputStruct ELSE OutputStruct = [OutputStruct,Return_Structure]
            ; Quick override of above line on 07/09 - Trying to remove bad fits from spectra and this already exists, but doing it on the Maxwellian fits, so I increase kappa from 0 on bad fits
            IF Return_Structure.Species_1_Kappa GE 0.100d THEN OutputStruct = OutputStruct ELSE OutputStruct = [OutputStruct,Return_Structure]
          ENDIF ELSE BEGIN
            OutputStruct = [OutputStruct, Return_Structure]
          ENDELSE
        ENDELSE
                
        
        ; Use the specified filename for output if given. Otherwise, create a name specific to the time
        IF KEYWORD_SET(OUTPUT_FILENAME) THEN BEGIN
          OutputFilename = OUTPUT_FILENAME
        ENDIF ELSE BEGIN
          
          IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
            ; Automatically name uncertainty file if input not given
            strlength = STRLEN(Input_Filename)
            strlength_unc = strlength - 4
            new_string = STRMID(Input_Filename, 0, strlength_unc)
            OutputFilename = STRCOMPRESS(STRING(new_string, '_Unc.csv'))
          ENDIF
          
          IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
            ; Automatically name uncertainty file if input not given
            strlength = STRLEN(Input_Filename)
            strlength_unc = strlength - 4
            new_string = STRMID(Input_Filename, 0, strlength_unc)
            OutputFilename = STRCOMPRESS(STRING(new_string, '_CEq.csv'))
          ENDIF
                
          ; Use IDL systime command to determine the exact UTC time that the program was run so as to avoid creating
          ; duplicate save files in the current IDL directory
          YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
          DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
          IF ULONG(DayOfRun) LT 10 THEN BEGIN
            DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
          ENDIF
          MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
          HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
          MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
          SecondOfRun = STRMID(SYSTIME(/UTC),17,2)

          ; Name the file if not already created
          IF ISA(OutputFilename) EQ 0 THEN OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
                                                                   HourOfRun, MinuteOfRun, '_', SecondOfRun, '.csv', $
                                                                   FORMAT = '(A,A,A,A,A,A,A,A,A,A)')
        ENDELSE
        
        ; Create a save file to free up RAM
        SAVE, OutputStruct, FILENAME = 'VIPER_SAVE_FILE.sav' ; Create a save file
       
       
        ; Change directory to where OUTPUT CSVs are located ; Logan 2025_02_27
        CD, CURRENT=RESTORE_DIR
        cutoff_pos = strpos(RESTORE_DIR, 'Documents', /REVERSE_SEARCH)
        IF cutoff_pos NE -1 THEN BEGIN
          Root_User_Path = STRMID(RESTORE_DIR, 0, cutoff_pos + STRLEN('Documents'))
          ; PRINT, 'Root_User_Path' ; Test PRINT command
        ENDIF ELSE BEGIN
          PRINT, 'ERROR in determining Root path'
        ENDELSE

        CSV_VIPER_Output = STRCOMPRESS(STRING(Root_User_Path, '/VIPER_Output/CSV_Files'))
        CD, CSV_VIPER_Output
        ; Write CSV file
        WRITE_CSV, OutputFilename, OutputStruct, HEADER = Header ; Save structures
        CD, RESTORE_DIR
        OutputStruct = 0 ; Free up RAM
      ENDIF 
    ENDFOR
  
  PRINT, 'Program has finished.'
  ; If no filename is provided, then create a CSV that can be edited and saved
  ENDIF


; ------------------
;     URANUS
; ------------------
ENDIF ELSE IF TOTAL(STRCMP(PlanetName, UranusStrings, /FOLD_CASE)) NE 0 THEN BEGIN 
  
  IF ISA(Input_Filename) EQ 1 THEN BEGIN
    
    ; Change directory to where CSVs are located
    CD, CURRENT = RESTORE_DIR
    STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
    CD, STR_DIR
    Input_CSV = READ_CSV(Input_Filename, HEADER = Structure_Names)
    CD, RESTORE_DIR
    
    ; If Full CSV keyword is set, use every row of CSV
    IF KEYWORD_SET(FULL_CSV) THEN CSV_Indices = LINDGEN(N_ELEMENTS(Input_CSV.Field01)) ELSE $
    ; If there is not a specific range, default to using every row, otherwise use the selected inputs
    IF ISA(CSV_Indices) EQ 0 THEN CSV_Indices = LINDGEN(N_ELEMENTS(Input_CSV.Field01)) ELSE $
                                  CSV_Indices = CSV_Indices - 2 ; -2 to go from CSV indexing to IDL indexing
    ; Set keyword clear on if quiet plotting is turned on so that it will clear. Can't see or access files anyways, so clear to free RAM 
    IF KEYWORD_SET(PLOT_QUIETLY) THEN CLEAR = 1               
    nloop = 0 ; Counter for row number of CSV file
    
    IF KEYWORD_SET(ELECTRON_DATA_PRODUCT) THEN EDP = 1 ELSE EDP = 0
    IF KEYWORD_SET(RETURN_ELECTRON_POWERLAW) THEN EPL = 1 ELSE EPL = 0
    
    FOR i = 0, N_ELEMENTS(CSV_Indices) - 1 DO BEGIN                       ;start loop
      nloop = nloop + 1 ; Augment loop for next row of CSV file
      CSV_Index = CSV_Indices(i)
      
      ; Clear any existing plot
      IF KEYWORD_SET(CLEAR) THEN BEGIN
        wwwwww = GETWINDOWS()
        ; If plots exist, clear them
        IF ISA(wwwwww) EQ 1 THEN FOR wwwwwwi = 0, N_ELEMENTS(wwwwww)-1 DO wwwwww[wwwwwwi].close
      ENDIF

      ; Read in CSV to proper structure format with correct tag names
      TempArr = MAKE_ARRAY(N_TAGS(Input_CSV)) ; Create temporary array to hold structure information
      FOR j = 0, N_TAGS(Input_CSV)-1 DO TempArr[j] = (Input_CSV.(j))[CSV_Index] ; Fill in temporary array
      ; Define Structure in correct format for VIPER code to use
      New_Struct = {col1: Structure_Names, col2: TempArr}
      struct_elements = N_ELEMENTS(New_Struct.col2)
      TempStruct = CREATE_STRUCT(New_Struct.col1[0], DOUBLE(New_Struct.col2[0]))
      FOR k = 1, N_ELEMENTS(New_Struct.col1) - 1 DO BEGIN
        TempStruct = CREATE_STRUCT(TempStruct, New_Struct.col1[k], DOUBLE(New_Struct.col2[k]))
      ENDFOR
      Structure = TempStruct
      
;      STOP ; Logan 2025_04_09 Testing

;                      ;       ----------------------------------------------
;                      ;       LEAVE THESE HERE FOR WHEN WE DEAL WITH SATURN!
;                      ;       ----------------------------------------------
;                      ;      ; Need to program in a compositional table call
;                      ;      
;                            FORWARD_FUNCTION CASSINI_RATIOS
;                            IF KEYWORD_SET(CASSINI_COMPOSITION) THEN Structure = CASSINI_RATIOS(Structure, ISA(UNCERTAINTIES))
;                      ;      
;                            FORWARD_FUNCTION DELAMERE_RATIOS
;                      ;      ; Want to insert Delamere Ratios as a tag in the structure
;                      ;      ; So try this:
;                            IF KEYWORD_SET(DELAMERE_COMPOSITION) THEN Structure.Delamere_Composition = 1
;                            IF Structure.Delamere_Composition EQ 1 THEN Structure = DELAMERE_RATIOS(Structure, ISA(UNCERTAINTIES))
;                            ;IF KEYWORD_SET(DELAMERE_COMPOSITION) THEN Structure = DELAMERE_RATIOS(Structure, ISA(UNCERTAINTIES))
      
      IF KEYWORD_SET(ANALYTICAL_VPHI) THEN Structure.Analytical_Vphi = 1d
      IF Structure.Analytical_Vphi GE 1 THEN BEGIN ; not fixed for Uranus
        FORWARD_FUNCTION VGR_CupVelocity
        FORWARD_FUNCTION VPHI_PROFILE
        Rj_Distance       = VGR_CupVelocity([0d,0d,0d], 'A', Structure, /RETURN_RJ)                                         ; RENAME TO R_PLANET. OR SOMETHING
        Vphi_Fixed        = VPHI_PROFILE(Rj_Distance, Structure.Analytical_Vphi)
        Structure.V2      = Vphi_Fixed ; Set appropriate value for Vphi
        Structure.Vary_V2 = 0d ; Forces the phi component of velocity to be a constant at that value
      ENDIF
       
      IF KEYWORD_SET(FORCE_SAVE_JPG) THEN BEGIN
        Structure.Save_Plot = 1
        IF ISA(MAKE_PLOT) EQ 0 THEN MAKE_PLOT=1 ; Turn on keyword argument for plotting in this case if it wasn't already on
      ENDIF
      
      IF KEYWORD_SET(POD) THEN BEGIN
        ; Check if the directory exists
        ;CD, CURRENT=RESTORE_DIR
        ;Plot_Out_Dir = STRCOMPRESS(STRING(RESTORE_DIR, POD))
        ;DIR_CHECK = FILE_TEST(Plot_Out_Dir, /DIRECTORY)
        ;IF DIR_CHECK EQ 0 THEN FILE_MKDIR, Plot_Out_Dir
        
        ; Logan 2025_02_27
        CD, CURRENT=RESTORE_DIR
        cutoff_pos = strpos(RESTORE_DIR, 'Documents', /REVERSE_SEARCH)
        IF cutoff_pos NE -1 THEN BEGIN
          Root_User_Path = STRMID(RESTORE_DIR, 0, cutoff_pos + STRLEN('Documents'))
          ; PRINT, 'Root_User_Path' ; Test PRINT command
        ENDIF ELSE BEGIN
          PRINT, 'ERROR in determining Root path'
        ENDELSE
        
        Root_VIPER_Output = STRCOMPRESS(STRING(Root_User_Path, '/VIPER_Output'))
        
        Plot_Out_Dir = STRCOMPRESS(STRING(Root_VIPER_Output + '/Plots' + POD)) ; Logan 2025_03_04 Changed from commas to +
        ;print, Plot_Out_Dir ; Test Line
        ;Print, POD ; Test Line
        DIR_CHECK = FILE_TEST(Plot_Out_Dir, /DIRECTORY)
        ;print, plot_out_dir
        IF DIR_CHECK EQ 0 THEN FILE_MKDIR, Plot_Out_Dir
        ; Logan 2025_02_27
        
      ENDIF ELSE BEGIN
        ; Maintain this ELSE to save to a location in VIPER Git dir if POD not specified
        POD = '/Output/Plots'
        CD, CURRENT=RESTORE_DIR
        Plot_Out_Dir = STRCOMPRESS(STRING(RESTORE_DIR, POD))
      ENDELSE
      
      IF KEYWORD_SET(NO_PLOT_SAVE) THEN Structure.Save_Plot = 0
      
      ; Add 100 to Save_Plot argument of structure to add a buffer to the plot and not have it pop up on the screen
      IF KEYWORD_SET(PLOT_QUIETLY) THEN Structure.Save_Plot = Structure.Save_Plot + 100
      
      IF KEYWORD_SET(NO_FIT) THEN Structure.Fit = 0
      IF KEYWORD_SET(FIT) THEN Structure.Fit = 1
      
      cupss = 0
      IF KEYWORD_SET(CUPA) THEN cupss = 1
      IF KEYWORD_SET(CUPB) THEN cupss = 2
      IF KEYWORD_SET(CUPC) THEN cupss = 3
      IF KEYWORD_SET(CUPD) THEN cupss = 4
      
      
      ; Extrapolate plasma measurements to the Equator
      ; This section of the code does not fit the data, it takes fitted data from where the S/C
      ; made in situ measurements of the plasma and then traces the magnetic field line to determine
      ; the location of the centrifugal equator. Using the code, latdistVoy, the density of the plasma
      ; is then calculated at the location of the centrifugal equator.
      ;
      ; Note: Extensive testing has not been done on this part of the code yet
      IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
        IF KEYWORD_SET(MAKE_PLOT) OR KEYWORD_SET(UNCERTAINTIES) OR KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN
          MESSAGE, 'VIPER code can not calculate centrifugal equator while producing a plot, fitting data, calculating uncertainties, or computing chi-space distributions.'
        ENDIF
      ENDIF

      ; Run through model to process the data  
      IF KEYWORD_SET(MAKE_PLOT) AND KEYWORD_SET(UNCERTAINTIES) THEN BEGIN   
        Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
        Return_Structure = VOYAGER_PLS_URANUS(Structure, Plot_Out_Dir, cupss, EDP, EPL, /MAKE_PLOT, /UNCERTAINTIES, WRITE_UNCERTAINTIES_IN_FILE=WRITE_UNCERTAINTIES_IN_FILE, CALC_FITTING_RESIDUALS=CALC_FITTING_RESIDUALS )

        ; Run the procedure to calculate uncertainties
      ENDIF ELSE IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN 
        Return_Structure = VOYAGER_PLS_URANUS(Structure, Plot_Out_Dir, cupss, EDP, EPL, /MAKE_PLOT, CALC_FITTING_RESIDUALS=CALC_FITTING_RESIDUALS)
        ; Run the procedure to calculate uncertainties and plot the result
      ENDIF ELSE IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN 
        Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
        Return_Structure = VOYAGER_PLS_URANUS(Structure, Plot_Out_Dir, cupss, EDP, EPL, /UNCERTAINTIES, WRITE_UNCERTAINTIES_IN_FILE=WRITE_UNCERTAINTIES_IN_FILE, CALC_FITTING_RESIDUALS=CALC_FITTING_RESIDUALS)

        ; Run the procedure as normal without producing a plot or uncertainties for the parameters - Useful for
        ; running batch fits so that plots don't need to be created and deleted, freeing up RAM,
        ; and it doesn't waste time generating uncertainties in case the batch fit was unsuccessful.
      ENDIF ELSE IF KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN 
        ; Check if there is an input number of points and sigma to use. If not default to 21 points and 3 sigma
        IF ISA(N_CHI_PTS) EQ 0 THEN BEGIN
          N_CHI_PTS = 21
          PRINT, 'No input provided for number of points in Chi space calculation. Defaulting to 21.'
        ENDIF
        IF ISA(N_CHI_SIG) EQ 0 THEN BEGIN
          N_CHI_SIG = 3
          PRINT, 'No input provided for number of sigma in Chi space calculation. Defaulting to 3 sigma.'
        ENDIF
        
        ; Check that there is no input_uncertainty_filename on the first iteration
        IF ((ISA(INPUT_UNCERTAINTY_FILENAME) EQ 0) AND (nloop EQ 1)) THEN UNC_FILE_LOGIC_CHECK = 0
        
        
        IF UNC_FILE_LOGIC_CHECK EQ 0 THEN BEGIN
          ; Do a quick search for the correct name and if a file exists
          strlength = STRLEN(Input_Filename)
          strlength_unc = strlength - 4
          new_string = STRMID(Input_Filename, 0, strlength_unc)
          INPUT_UNCERTAINTY_FILENAME = STRCOMPRESS(STRING(new_string, '_Unc.csv'))
          
          ; Change directory to where CSVs are located
          CD, CURRENT = RESTORE_DIR
          STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
          CD, STR_DIR
          ; Check if the file exists as it should be named
          IF FILE_TEST(INPUT_UNCERTAINTY_FILENAME) EQ 0 THEN MESSAGE, 'Error: No file containing parameter uncertainties was provided or found.'
          CD, RESTORE_DIR
        ENDIF
        
        ; Change directory to where CSVs are located
        CD, CURRENT = RESTORE_DIR
        STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
        CD, STR_DIR
        ; Read in uncertainty filename
        Input_CSV_Unc = READ_CSV(INPUT_UNCERTAINTY_FILENAME, HEADER = Structure_Names)
        CD, RESTORE_DIR
        
        ; Read in CSV to proper structure format with correct tag names
        TempArr = MAKE_ARRAY(N_TAGS(Input_CSV_Unc)) ; Create temporary array to hold structure information
        FOR j = 0, N_TAGS(Input_CSV_Unc)-1 DO TempArr[j] = (Input_CSV_Unc.(j));[CSV_Index] ; Fill in temporary array
        ; Define Structure in correct format for VIPER code to use
        New_Struct = {col1: Structure_Names, col2: TempArr}
        struct_elements = N_ELEMENTS(New_Struct.col2)
        TempStruct = CREATE_STRUCT(New_Struct.col1[0], DOUBLE(New_Struct.col2[0]))
        FOR k = 1, N_ELEMENTS(New_Struct.col1) - 1 DO BEGIN
          TempStruct = CREATE_STRUCT(TempStruct, New_Struct.col1[k], DOUBLE(New_Struct.col2[k]))
        ENDFOR
        Structure_Unc = TempStruct   
        

        Return_Structure = VOYAGER_PLS_URANUS(Structure, Plot_Out_Dir, Structure_Unc=Structure_Unc, cupss, EDP, EPL, N_CHI_PTS=N_CHI_PTS, N_CHI_SIG=N_CHI_SIG, /COMPUTE_CHI_SPACE)

      ENDIF ELSE IF KEYWORD_SET(REMOVE_BAD_KAPPA) THEN BEGIN
        ; Return just the input (it's only checking for bad values of Kappa, i.e at the lower limit)
        Return_Structure = Structure
      ENDIF ELSE IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
        Return_Structure = VOYAGER_PLS_URANUS(Structure, Plot_Out_Dir, cupss, /CEq)
      ENDIF ELSE BEGIN
        ;stop
        Return_Structure = VOYAGER_PLS_URANUS(Structure, Plot_Out_Dir, cupss, EDP, EPL)
      ENDELSE
     
      
      
      ; Save the file - Skip if No Return keyword is set
      IF NOT KEYWORD_SET(NO_CSV_RETURN) THEN BEGIN
        IF nloop EQ 1 THEN BEGIN ; When nloop is 1, it is on the first iteration so it has to handle the saving differently
          OutputStruct = CREATE_STRUCT(Return_Structure)
          Header = TAG_NAMES(Return_Structure)
        ENDIF ELSE BEGIN ; Concactanate the structures into one large array
          RESTORE, 'VIPER_SAVE_FILE.sav' ; Restore save file
          IF KEYWORD_SET(REMOVE_BAD_KAPPA) THEN BEGIN
            ; This is basically a one off feature to parse the data down
            ; This line does not remove the first iteration if it is a 'bad' value. If first value is less than 1.501001 it will require manual deletion (quicker than altering code)
            ; IF Return_Structure.Species_1_Kappa LE 1.501001d THEN OutputStruct = OutputStruct ELSE OutputStruct = [OutputStruct,Return_Structure]
            ; Quick override of above line on 07/09 - Trying to remove bad fits from spectra and this already exists, but doing it on the Maxwellian fits, so I increase kappa from 0 on bad fits
            IF Return_Structure.Species_1_Kappa GE 0.100d THEN OutputStruct = OutputStruct ELSE OutputStruct = [OutputStruct,Return_Structure]
          ENDIF ELSE BEGIN
            OutputStruct = [OutputStruct, Return_Structure]
          ENDELSE
        ENDELSE
                
        
        ; Use the specified filename for output if given. Otherwise, create a name specific to the time
        IF KEYWORD_SET(OUTPUT_FILENAME) THEN BEGIN
          OutputFilename = OUTPUT_FILENAME
        ENDIF ELSE BEGIN
          
          IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
            ; Automatically name uncertainty file if input not given
            strlength = STRLEN(Input_Filename)
            strlength_unc = strlength - 4
            new_string = STRMID(Input_Filename, 0, strlength_unc)
            OutputFilename = STRCOMPRESS(STRING(new_string, '_Unc.csv'))
          ENDIF
          
          IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
            ; Automatically name uncertainty file if input not given
            strlength = STRLEN(Input_Filename)
            strlength_unc = strlength - 4
            new_string = STRMID(Input_Filename, 0, strlength_unc)
            OutputFilename = STRCOMPRESS(STRING(new_string, '_CEq.csv'))
          ENDIF
                
          ; Use IDL systime command to determine the exact UTC time that the program was run so as to avoid creating
          ; duplicate save files in the current IDL directory
          YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
          DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
          IF ULONG(DayOfRun) LT 10 THEN BEGIN
            DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
          ENDIF
          MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
          HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
          MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
          SecondOfRun = STRMID(SYSTIME(/UTC),17,2)

          ; Name the file if not already created
          IF ISA(OutputFilename) EQ 0 THEN OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
                                                                   HourOfRun, MinuteOfRun, '_', SecondOfRun, '.csv', $
                                                                   FORMAT = '(A,A,A,A,A,A,A,A,A,A)')
        ENDELSE
        
        ; Create a save file to free up RAM
        SAVE, OutputStruct, FILENAME = 'VIPER_SAVE_FILE.sav' ; Create a save file
       
       
        ; Change directory to where OUTPUT CSVs are located ; Logan 2025_02_27
        CD, CURRENT=RESTORE_DIR
        cutoff_pos = strpos(RESTORE_DIR, 'Documents', /REVERSE_SEARCH)
        IF cutoff_pos NE -1 THEN BEGIN
          Root_User_Path = STRMID(RESTORE_DIR, 0, cutoff_pos + STRLEN('Documents'))
          ; PRINT, 'Root_User_Path' ; Test PRINT command
        ENDIF ELSE BEGIN
          PRINT, 'ERROR in determining Root path'
        ENDELSE

        CSV_VIPER_Output = STRCOMPRESS(STRING(Root_User_Path, '/VIPER_Output/CSV_Files'))
        CD, CSV_VIPER_Output
        ; Write CSV file
        WRITE_CSV, OutputFilename, OutputStruct, HEADER = Header ; Save structures
;stop
        CD, RESTORE_DIR
        OutputStruct = 0 ; Free up RAM
      ENDIF 
    ENDFOR
  
  PRINT, 'Program has finished.'; If no filename is provided, then create a CSV that can be edited and saved
  ENDIF


  ; ------------------
  ;     NEPTUNE
  ; ------------------
ENDIF ELSE IF TOTAL(STRCMP(PlanetName, NeptuneStrings, /FOLD_CASE)) NE 0 THEN BEGIN

  IF ISA(Input_Filename) EQ 1 THEN BEGIN

    ; Change directory to where CSVs are located
    CD, CURRENT = RESTORE_DIR
    STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
    CD, STR_DIR
    Input_CSV = READ_CSV(Input_Filename, HEADER = Structure_Names)
    CD, RESTORE_DIR

    ; If Full CSV keyword is set, use every row of CSV
    IF KEYWORD_SET(FULL_CSV) THEN CSV_Indices = LINDGEN(N_ELEMENTS(Input_CSV.Field01)) ELSE $
      ; If there is not a specific range, default to using every row, otherwise use the selected inputs
      IF ISA(CSV_Indices) EQ 0 THEN CSV_Indices = LINDGEN(N_ELEMENTS(Input_CSV.Field01)) ELSE $
      CSV_Indices = CSV_Indices - 2 ; -2 to go from CSV indexing to IDL indexing
    ; Set keyword clear on if quiet plotting is turned on so that it will clear. Can't see or access files anyways, so clear to free RAM
    IF KEYWORD_SET(PLOT_QUIETLY) THEN CLEAR = 1
    nloop = 0 ; Counter for row number of CSV file

    IF KEYWORD_SET(ELECTRON_DATA_PRODUCT) THEN EDP = 1 ELSE EDP = 0
    IF KEYWORD_SET(RETURN_ELECTRON_POWERLAW) THEN EPL = 1 ELSE EPL = 0


    FOR i = 0, N_ELEMENTS(CSV_Indices) - 1 DO BEGIN                       ;start loop
      nloop = nloop + 1 ; Augment loop for next row of CSV file
      CSV_Index = CSV_Indices(i)

      ; Clear any existing plot
      IF KEYWORD_SET(CLEAR) THEN BEGIN
        wwwwww = GETWINDOWS()
        ; If plots exist, clear them
        IF ISA(wwwwww) EQ 1 THEN FOR wwwwwwi = 0, N_ELEMENTS(wwwwww)-1 DO wwwwww[wwwwwwi].close
      ENDIF

      ; Read in CSV to proper structure format with correct tag names
      TempArr = MAKE_ARRAY(N_TAGS(Input_CSV)) ; Create temporary array to hold structure information
      FOR j = 0, N_TAGS(Input_CSV)-1 DO TempArr[j] = (Input_CSV.(j))[CSV_Index] ; Fill in temporary array
      ; Define Structure in correct format for VIPER code to use
      New_Struct = {col1: Structure_Names, col2: TempArr}
      struct_elements = N_ELEMENTS(New_Struct.col2)
      TempStruct = CREATE_STRUCT(New_Struct.col1[0], DOUBLE(New_Struct.col2[0]))
      FOR k = 1, N_ELEMENTS(New_Struct.col1) - 1 DO BEGIN
        TempStruct = CREATE_STRUCT(TempStruct, New_Struct.col1[k], DOUBLE(New_Struct.col2[k]))
      ENDFOR
      Structure = TempStruct

      IF KEYWORD_SET(ANALYTICAL_VPHI) THEN Structure.Analytical_Vphi = 1d
      IF Structure.Analytical_Vphi GE 1 THEN BEGIN ; not fixed for Neptune
        FORWARD_FUNCTION VGR_CupVelocity
        FORWARD_FUNCTION VPHI_PROFILE
        Rj_Distance       = VGR_CupVelocity([0d,0d,0d], 'A', Structure, /RETURN_RJ)                                         ; RENAME TO R_PLANET. OR SOMETHING
        Vphi_Fixed        = VPHI_PROFILE(Rj_Distance, Structure.Analytical_Vphi)
        Structure.V2      = Vphi_Fixed ; Set appropriate value for Vphi
        Structure.Vary_V2 = 0d ; Forces the phi component of velocity to be a constant at that value
      ENDIF

      IF KEYWORD_SET(FORCE_SAVE_JPG) THEN BEGIN
        Structure.Save_Plot = 1
        IF ISA(MAKE_PLOT) EQ 0 THEN MAKE_PLOT=1 ; Turn on keyword argument for plotting in this case if it wasn't already on
      ENDIF

      IF KEYWORD_SET(POD) THEN BEGIN
        ; Check if the directory exists
        ;CD, CURRENT=RESTORE_DIR
        ;Plot_Out_Dir = STRCOMPRESS(STRING(RESTORE_DIR, POD))
        ;DIR_CHECK = FILE_TEST(Plot_Out_Dir, /DIRECTORY)
        ;IF DIR_CHECK EQ 0 THEN FILE_MKDIR, Plot_Out_Dir
        
        ; Logan 2025_02_27
        CD, CURRENT=RESTORE_DIR
        cutoff_pos = strpos(RESTORE_DIR, 'Documents', /REVERSE_SEARCH)
        IF cutoff_pos NE -1 THEN BEGIN
          Root_User_Path = STRMID(RESTORE_DIR, 0, cutoff_pos + STRLEN('Documents'))
          ; PRINT, 'Root_User_Path' ; Test PRINT command
        ENDIF ELSE BEGIN
          PRINT, 'ERROR in determining Root path'
        ENDELSE
        
        Root_VIPER_Output = STRCOMPRESS(STRING(Root_User_Path, '/VIPER_Output'))
        
        Plot_Out_Dir = STRCOMPRESS(STRING(Root_VIPER_Output, '/Plots', POD))
        ;print, Plot_Out_Dir ; Test Line
        ;Print, POD ; Test Line
        DIR_CHECK = FILE_TEST(Plot_Out_Dir, /DIRECTORY)
        IF DIR_CHECK EQ 0 THEN FILE_MKDIR, Plot_Out_Dir
        ; Logan 2025_02_27
        
      ENDIF ELSE BEGIN
        ; Maintain this ELSE to save to a location in VIPER Git dir if POD not specified
        POD = '/Output/Plots'
        CD, CURRENT=RESTORE_DIR
        Plot_Out_Dir = STRCOMPRESS(STRING(RESTORE_DIR, POD))
      ENDELSE

      IF KEYWORD_SET(NO_PLOT_SAVE) THEN Structure.Save_Plot = 0

      ; Add 100 to Save_Plot argument of structure to add a buffer to the plot and not have it pop up on the screen
      IF KEYWORD_SET(PLOT_QUIETLY) THEN Structure.Save_Plot = Structure.Save_Plot + 100

      IF KEYWORD_SET(NO_FIT) THEN Structure.Fit = 0
      IF KEYWORD_SET(FIT) THEN Structure.Fit = 1

      cupss = 0
      IF KEYWORD_SET(CUPA) THEN cupss = 1
      IF KEYWORD_SET(CUPB) THEN cupss = 2
      IF KEYWORD_SET(CUPC) THEN cupss = 3
      IF KEYWORD_SET(CUPD) THEN cupss = 4

      ; Extrapolate plasma measurements to the Equator
      ; This section of the code does not fit the data, it takes fitted data from where the S/C
      ; made in situ measurements of the plasma and then traces the magnetic field line to determine
      ; the location of the centrifugal equator. Using the code, latdistVoy, the density of the plasma
      ; is then calculated at the location of the centrifugal equator.
      ;
      ; Note: Extensive testing has not been done on this part of the code yet
      IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
        IF KEYWORD_SET(MAKE_PLOT) OR KEYWORD_SET(UNCERTAINTIES) OR KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN
          MESSAGE, 'VIPER code can not calculate centrifugal equator while producing a plot, fitting data, calculating uncertainties, or computing chi-space distributions.'
        ENDIF
      ENDIF

      ;stop
      ; Run through model to process the data
      IF KEYWORD_SET(MAKE_PLOT) AND KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
        Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
        Return_Structure = VOYAGER_PLS_NEPTUNE(Structure, Plot_Out_Dir, cupss, EDP, EPL, /MAKE_PLOT, /UNCERTAINTIES, WRITE_UNCERTAINTIES_IN_FILE=WRITE_UNCERTAINTIES_IN_FILE, CALC_FITTING_RESIDUALS=CALC_FITTING_RESIDUALS )

        ; Run the procedure to calculate uncertainties
      ENDIF ELSE IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN
        Return_Structure = VOYAGER_PLS_NEPTUNE(Structure, Plot_Out_Dir, cupss, EDP, EPL, /MAKE_PLOT, CALC_FITTING_RESIDUALS=CALC_FITTING_RESIDUALS)
        ; Run the procedure to calculate uncertainties and plot the result
      ENDIF ELSE IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
        Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
        Return_Structure = VOYAGER_PLS_NEPTUNE(Structure, Plot_Out_Dir, cupss, EDP, EPL, /UNCERTAINTIES, WRITE_UNCERTAINTIES_IN_FILE=WRITE_UNCERTAINTIES_IN_FILE, CALC_FITTING_RESIDUALS=CALC_FITTING_RESIDUALS)

        ; Run the procedure as normal without producing a plot or uncertainties for the parameters - Useful for
        ; running batch fits so that plots don't need to be created and deleted, freeing up RAM,
        ; and it doesn't waste time generating uncertainties in case the batch fit was unsuccessful.
      ENDIF ELSE IF KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN
        ; Check if there is an input number of points and sigma to use. If not default to 21 points and 3 sigma
        IF ISA(N_CHI_PTS) EQ 0 THEN BEGIN
          N_CHI_PTS = 21
          PRINT, 'No input provided for number of points in Chi space calculation. Defaulting to 21.'
        ENDIF
        IF ISA(N_CHI_SIG) EQ 0 THEN BEGIN
          N_CHI_SIG = 3
          PRINT, 'No input provided for number of sigma in Chi space calculation. Defaulting to 3 sigma.'
        ENDIF

        ; Check that there is no input_uncertainty_filename on the first iteration
        IF ((ISA(INPUT_UNCERTAINTY_FILENAME) EQ 0) AND (nloop EQ 1)) THEN UNC_FILE_LOGIC_CHECK = 0


        IF UNC_FILE_LOGIC_CHECK EQ 0 THEN BEGIN
          ; Do a quick search for the correct name and if a file exists
          strlength = STRLEN(Input_Filename)
          strlength_unc = strlength - 4
          new_string = STRMID(Input_Filename, 0, strlength_unc)
          INPUT_UNCERTAINTY_FILENAME = STRCOMPRESS(STRING(new_string, '_Unc.csv'))

          ; Change directory to where CSVs are located
          CD, CURRENT = RESTORE_DIR
          STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
          CD, STR_DIR
          ; Check if the file exists as it should be named
          IF FILE_TEST(INPUT_UNCERTAINTY_FILENAME) EQ 0 THEN MESSAGE, 'Error: No file containing parameter uncertainties was provided or found.'
          CD, RESTORE_DIR
        ENDIF

        ; Change directory to where CSVs are located
        CD, CURRENT = RESTORE_DIR
        STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
        CD, STR_DIR
        ; Read in uncertainty filename
        Input_CSV_Unc = READ_CSV(INPUT_UNCERTAINTY_FILENAME, HEADER = Structure_Names)
        CD, RESTORE_DIR

        ; Read in CSV to proper structure format with correct tag names
        TempArr = MAKE_ARRAY(N_TAGS(Input_CSV_Unc)) ; Create temporary array to hold structure information
        FOR j = 0, N_TAGS(Input_CSV_Unc)-1 DO TempArr[j] = (Input_CSV_Unc.(j));[CSV_Index] ; Fill in temporary array
        ; Define Structure in correct format for VIPER code to use
        New_Struct = {col1: Structure_Names, col2: TempArr}
        struct_elements = N_ELEMENTS(New_Struct.col2)
        TempStruct = CREATE_STRUCT(New_Struct.col1[0], DOUBLE(New_Struct.col2[0]))
        FOR k = 1, N_ELEMENTS(New_Struct.col1) - 1 DO BEGIN
          TempStruct = CREATE_STRUCT(TempStruct, New_Struct.col1[k], DOUBLE(New_Struct.col2[k]))
        ENDFOR
        Structure_Unc = TempStruct


        Return_Structure = VOYAGER_PLS_NEPTUNE(Structure, Plot_Out_Dir, Structure_Unc=Structure_Unc, cupss, EDP, EPL, N_CHI_PTS=N_CHI_PTS, N_CHI_SIG=N_CHI_SIG, /COMPUTE_CHI_SPACE)

      ENDIF ELSE IF KEYWORD_SET(REMOVE_BAD_KAPPA) THEN BEGIN
        ; Return just the input (it's only checking for bad values of Kappa, i.e at the lower limit)
        Return_Structure = Structure
      ENDIF ELSE IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
        Return_Structure = VOYAGER_PLS_NEPTUNE(Structure, Plot_Out_Dir, cupss, /CEq)
      ENDIF ELSE BEGIN
        ;stop
        Return_Structure = VOYAGER_PLS_NEPTUNE(Structure, Plot_Out_Dir, cupss, EDP, EPL)
      ENDELSE



      ; Save the file - Skip if No Return keyword is set
      IF NOT KEYWORD_SET(NO_CSV_RETURN) THEN BEGIN
        IF nloop EQ 1 THEN BEGIN ; When nloop is 1, it is on the first iteration so it has to handle the saving differently
          OutputStruct = CREATE_STRUCT(Return_Structure)
          Header = TAG_NAMES(Return_Structure)
        ENDIF ELSE BEGIN ; Concactanate the structures into one large array
          RESTORE, 'VIPER_SAVE_FILE.sav' ; Restore save file
          IF KEYWORD_SET(REMOVE_BAD_KAPPA) THEN BEGIN
            ; This is basically a one off feature to parse the data down
            ; This line does not remove the first iteration if it is a 'bad' value. If first value is less than 1.501001 it will require manual deletion (quicker than altering code)
            ; IF Return_Structure.Species_1_Kappa LE 1.501001d THEN OutputStruct = OutputStruct ELSE OutputStruct = [OutputStruct,Return_Structure]
            ; Quick override of above line on 07/09 - Trying to remove bad fits from spectra and this already exists, but doing it on the Maxwellian fits, so I increase kappa from 0 on bad fits
            IF Return_Structure.Species_1_Kappa GE 0.100d THEN OutputStruct = OutputStruct ELSE OutputStruct = [OutputStruct,Return_Structure]
          ENDIF ELSE BEGIN
            OutputStruct = [OutputStruct, Return_Structure]
          ENDELSE
        ENDELSE


        ; Use the specified filename for output if given. Otherwise, create a name specific to the time
        IF KEYWORD_SET(OUTPUT_FILENAME) THEN BEGIN
          OutputFilename = OUTPUT_FILENAME
        ENDIF ELSE BEGIN

          IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
            ; Automatically name uncertainty file if input not given
            strlength = STRLEN(Input_Filename)
            strlength_unc = strlength - 4
            new_string = STRMID(Input_Filename, 0, strlength_unc)
            OutputFilename = STRCOMPRESS(STRING(new_string, '_Unc.csv'))
          ENDIF

          IF KEYWORD_SET(CENTRIFUGAL_EQUATOR) THEN BEGIN
            ; Automatically name uncertainty file if input not given
            strlength = STRLEN(Input_Filename)
            strlength_unc = strlength - 4
            new_string = STRMID(Input_Filename, 0, strlength_unc)
            OutputFilename = STRCOMPRESS(STRING(new_string, '_CEq.csv'))
          ENDIF

          ; Use IDL systime command to determine the exact UTC time that the program was run so as to avoid creating
          ; duplicate save files in the current IDL directory
          YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
          DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
          IF ULONG(DayOfRun) LT 10 THEN BEGIN
            DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
          ENDIF
          MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
          HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
          MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
          SecondOfRun = STRMID(SYSTIME(/UTC),17,2)

          ; Name the file if not already created
          IF ISA(OutputFilename) EQ 0 THEN OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
            HourOfRun, MinuteOfRun, '_', SecondOfRun, '.csv', $
            FORMAT = '(A,A,A,A,A,A,A,A,A,A)')
        ENDELSE

        ; Create a save file to free up RAM
        SAVE, OutputStruct, FILENAME = 'VIPER_SAVE_FILE.sav' ; Create a save file


        ; Change directory to where OUTPUT CSVs are located ; Logan 2025_02_27
        CD, CURRENT=RESTORE_DIR
        cutoff_pos = strpos(RESTORE_DIR, 'Documents', /REVERSE_SEARCH)
        IF cutoff_pos NE -1 THEN BEGIN
          Root_User_Path = STRMID(RESTORE_DIR, 0, cutoff_pos + STRLEN('Documents'))
          ; PRINT, 'Root_User_Path' ; Test PRINT command
        ENDIF ELSE BEGIN
          PRINT, 'ERROR in determining Root path'
        ENDELSE

        CSV_VIPER_Output = STRCOMPRESS(STRING(Root_User_Path, '/VIPER_Output/CSV_Files'))
        CD, CSV_VIPER_Output
        ; Write CSV file
        WRITE_CSV, OutputFilename, OutputStruct, HEADER = Header ; Save structures
        ;stop
        CD, RESTORE_DIR
        OutputStruct = 0 ; Free up RAM
      ENDIF
    ENDFOR

    PRINT, 'Program has finished.'; If no filename is provided, then create a CSV that can be edited and saved
  ENDIF
ENDIF

  
  
;  ; LOGAN 05/14/2020 - The below code was written years ago when the codebase had a different design
;  ;                    It is functional for the provided /INITIAL argument. At some point, the 'change_paramaters.csv'
;  ;                    might be updated, but the above code is more powerful than fitting a single time with manual
;  ;                    CSV entry
;  ENDIF ELSE BEGIN
;    
;    ; Change directories to produce the CSV file
;    CD, CURRENT = RESTORE_DIR
;    STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
;    CD, STR_DIR
;    
;    ; On initial setup it runs through the code once
;    IF KEYWORD_SET(INITIAL) THEN BEGIN
;      
;      ; These three lines set important variables that need to be called
;      IF KEYWORD_SET(RETURN_ELECTRON_POWERLAW) THEN EPL = 1 ELSE EPL = 0
;      IF KEYWORD_SET(RETURN_ELECTRON_POWERLAW) THEN EDP = 1 ELSE EDP = 0
;      cupss = 0
;      INPUT_CSV = READ_CSV('Initial.csv', Header = Structure_Names)
;    ENDIF ELSE BEGIN
;      WRITE_CSV, 'Change_Parameters.csv',  Structure, HEADER = Structure_Names ; Create csv file to edit
;      SPAWN, 'Open Change_Parameters.csv'
;      PRINT, 'Edit Change_Parameters.csv and then save the file. Then press any key to continue.'
;      dum = GET_KBRD(1)               ; Prompt for key input to continue
;    
;      ; Read in CSV to proper structure format with correct tag names
;      Input_CSV = READ_CSV('change_parameters.csv', HEADER = Structure_Names)
;    ENDELSE
;    
;    CD, RESTORE_DIR
;    
;    TempArr = MAKE_ARRAY(N_TAGS(Input_CSV)) ; Create temporary array to hold structure information
;    FOR j = 0, N_TAGS(Input_CSV)-1 DO TempArr[j] = (Input_CSV.(j))[0] ; Fill in temporary array
;    ; Define Structure in correct format for VIPER code to use
;    New_Struct = {col1: Structure_Names, col2: TempArr}
;    struct_elements = N_ELEMENTS(New_Struct.col2)
;    TempStruct = CREATE_STRUCT(New_Struct.col1[0], DOUBLE(New_Struct.col2[0]))
;    FOR k = 1, N_ELEMENTS(New_Struct.col1) - 1 DO BEGIN
;      TempStruct = CREATE_STRUCT(TempStruct, New_Struct.col1[k], DOUBLE(New_Struct.col2[k]))
;    ENDFOR
;    Structure = TempStruct
;    
;    ; Run through model to process the data
;    
;    ; Run the procedure and plot the result
;    IF KEYWORD_SET(MAKE_PLOT) AND KEYWORD_SET(UNCERTAINTIES) THEN BEGIN 
;
;      Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
;      Return_Structure = VOYAGER_PLS(Structure, Plot_Out_Dir,cupss, /MAKE_PLOT, /UNCERTAINTIES)
;
;      ; Run the procedure to calculate uncertainties
;    ENDIF ELSE IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN
;
;      Return_Structure = VOYAGER_PLS(Structure, Plot_Out_Dir, cupss, EDP, EPL, /MAKE_PLOT)
;
;      ; Run the procedure to calculate uncertainties and plot the result
;    ENDIF ELSE IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
;      
;      Structure.Fit = 1 ; Set the Fit tag on so that uncertainties can be calculated
;      Return_Structure = VOYAGER_PLS(Structure, Plot_Out_Dir,cupss, /UNCERTAINTIES)
;
;      ; Run the procedure as normal without producing a plot or uncertainties for the parameters - Useful for
;      ; running batch fits so that plots don't need to be created and deleted and it doesn't waste time generating
;      ; uncertainties in case the batch fit was unsuccessful.
;    ENDIF ELSE BEGIN
;       
;       Return_Structure = VOYAGER_PLS(Structure, Plot_Out_Dir, cupss, EDP, EPL)
;
;    ENDELSE
;    
;  PRINT, 'Program has finished.'
;  ENDELSE
;ENDIF
    

; Measure time of program
TOC
END


