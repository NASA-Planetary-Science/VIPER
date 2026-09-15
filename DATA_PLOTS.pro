FUNCTION DATA_PLOTS
  ; FILLER - Necessary for correct compilation by RESOLVE_ROUTINE
END

PRO DATA_PLOTS_ELECTRONS, input_file, RJ=RJ, CLEAR=CLEAR

  ; Inputs
  ; input_file - The string argument for the CSV being plotted
  ;   Code assumes that the input file is located in the folder "\Master Fits\PDS" as those files are structured appropriately

  ; Optional Keywords
  ; RJ       - Plot versus Rj instead of versus decimal date
  ; CLEAR    - Clear the plots that are on screen at the beginning

  ; Used the command ascii_template('V1_Kappa_PDS_Format.csv') to generate this template for a structure and then copied it in
  VERSION        = FLOAT(1.00000)
  DATASTART      = LONG(1)
  DELIMITER      = BYTE(44)
  MISSINGVALUE   = DOUBLE(-9)
  COMMENTSYMBOL  = ''
  FIELDCOUNT     = LONG(11)
  FIELDTYPES     = LONG([7,3,5,5,5,5,5,5,5,5,5])
  FIELDNAMES     = ['Time', 'Spacecraft', 'RJ', 'X', 'Y', 'Z', 'dn_Cold', 'T_Cold', 'Kappa', $
                    'dn_Hot', 'T_Hot']
  FIELDLOCATIONS = LONG([0,20,22,40,59,77,95,114,133,152,171])
  FIELDGROUPS    = LONG([0,1,2,3,4,5,6,7,8,9,10])
  
  Data_Template = {VERSION:VERSION, DATASTART:DATASTART, DELIMITER:DELIMITER, MISSINGVALUE:MISSINGVALUE, COMMENTSYMBOL:COMMENTSYMBOL, $
                   FIELDCOUNT:FIELDCOUNT, FIELDTYPES:FIELDTYPES, FIELDNAMES:FIELDNAMES, FIELDLOCATIONS:FIELDLOCATIONS, FIELDGROUPS:FIELDGROUPS}

  ; Close existing plots easily
  IF KEYWORD_SET(CLEAR) THEN BEGIN
    wwwwww = GETWINDOWS()
    ; If plots exist, clear them
    IF ISA(wwwwww) EQ 1 THEN FOR wwwwwwi = 0, N_ELEMENTS(wwwwww)-1 DO wwwwww[wwwwwwi].close
  ENDIF
  
  ; Change directories to where CSV file is
  CD, 'D:\Documents\LASP\VIPER Code\Master Fits\PDS'


  ; Read in VIPER_Structure
  Data = READ_ASCII(input_file, template=Data_Template)
   

  ; Define Time array  
  DOY    = DOUBLE(STRMID(Data.Time,6,3))
  Hour   = DOUBLE(STRMID(Data.Time,10,2))
  Minute = DOUBLE(STRMID(Data.Time,13,2))
  Second = DOUBLE(STRMID(Data.Time,16,2))
  xaxis   = DOY + Hour/24d + Minute/(24d*60d) + Second/(24d*60d*60d)
  
  ; Replace X-axis with Rj if requested
  IF KEYWORD_SET(RJ) THEN BEGIN
    xaxis = Data.Rj
  ENDIF ELSE xaxis = Data.Rj ; Was going to create time series originally but only radial profiles were of interest
  
  ; Density and Temperature parameters for plotting
  dn   = [[Data.dn_Cold[*]],[Data.dn_Hot[*]]]
  T    = [[Data.T_Cold[*]],[Data.T_Hot[*]]]
  Kap  = Data.Kappa[*]
  
  ; Split the Data into both spacecraft - NOT YET FUNCTIONAL/DOESN'T DO THIS
  
  ; Define the titles of the plots
  Titles                          = ['$Cold Density vs. R_J$', '$Kappa vs. R_J$', '$Cold Temperature vs. R_J$', '$Hot Density vs. R_J$', '$Hot Temperature vs. R_J$', '$Densities vs. R_J$', '$Temperatures vs. R_J$']
  IF KEYWORD_SET(RJ) THEN Xtitles = '$R_J$' ELSE Xtitles = 'Time (DOY 1979)'
  YTitles                         = ['Density (number/cc)', 'Kappa', 'Temperature (eV)']
  
  XRanges = [5,45]
  ; Plot Naming
  IF Kap[0] GT 0 THEN ModeName = 'Kappa' ELSE ModeName = 'BiMaxwellian' ; If BiMaxwellian Fits, Kappa is 0 for all values so this should auto name
  OutputFilename1 = STRCOMPRESS(STRING(ModeName, '_Densities_', XRanges[0], '-', XRanges[1], 'Rj.jpg'),/REMOVE_ALL)
  OutputFilename2 = STRCOMPRESS(STRING(ModeName, '_Temperatures_', XRanges[0], '-', XRanges[1], 'Rj.jpg'),/REMOVE_ALL)
  
  ; Plots of the Kappa values
  ;p = plot(Xaxis, Kap, symbol = '+', linestyle = 6, YRANGE = [1.4,15], title = Titles[1], xtitle = XTitles, ytitle = YTitles[1], XRange = XRanges)
  
  ; Change directory for plotting and saving
  CD, 'D:\Documents\LASP\VIPER Code\Master Fits'
  p = plot(Xaxis, dn[*,0], symbol = '+', linestyle = 6, title = Titles[5], xtitle = XTitles, ytitle = YTitles[0], yrange = [0,20], XRange = XRanges)
      op1 = plot(Xaxis, dn[*,1], symbol = '+', linestyle = 6, /ylog, yrange = [1d-2,1d3], color = 'red', XRange = XRanges, /OVERPLOT)
  p.SAVE, OutputFilename1, Border = 10, RESOLUTION = 150
  p = plot(Xaxis, T[*,0], symbol = '+', linestyle = 6, /ylog, title = Titles[6], xtitle = XTitles, ytitle = YTitles[2], yrange = [1d1,1d4], XRange = XRanges)
      op1 = plot(Xaxis, T[*,1], symbol = '+', linestyle = 6, /ylog, yrange = [1d1,1d4], color = 'red', XRange = XRanges, /OVERPLOT)
  p.SAVE, OutputFilename2, Border = 10, RESOLUTION = 150
  ; Return to root directory
  CD, 'D:\Documents\LASP\VIPER Code'
END

PRO PDS_DATA_OUTPUT_ELECTRONS, input_file, output_filename
  
  ; Change directories to where CSV file is
  CONRD
  
  CD, 'D:\Documents\LASP\VIPER Code\CSV_Files'
  
  VIPER_Structure = READ_CSV(input_file)

  ; Begin formatting PDS Data
  ; What is needed
  ;   Time
  ;   S/C (1 for V1, 2 for V2)
  ;   Rj (x, y, z)  (Calculated in EDP product so look there for code)
  ;   X_S3
  ;   Y_S3
  ;   Z_S3
  ;   Fit ne_cold
  ;   Fit Te_cold
  ;   Fit Kappa_cold (If fill value then used 2 Maxwellians for fitting procedure)
  ;   Fit ne_hot
  ;   Fit Te_hot
  ;   NO KAPPA HOT (Only ever fit with one kappa distribution)
  ;  
  ; Should be 11 columns if not including Lat and Lon, or 20 columns if I include them
  
  ; Read in the times and format to a proper UTC string
  Year   = LONG(VIPER_Structure.Field01)
  DOY_a  = LONG(VIPER_Structure.Field02) ; Variable name of DOY was calling function in newly added kk_2009.pro called DOY, so changed variable to resolve errors
  Hour   = LONG(VIPER_Structure.Field03)
  Minute = LONG(VIPER_Structure.Field04)
  Second = LONG(VIPER_Structure.Field05)
  
  Times_st = MAKE_ARRAY(N_ELEMENTS(YEAR), /STRING)
  DOY_st   = STRING(DOY_a)
  Hr_st    = STRING(Hour)
  Min_st   = STRING(Minute)
  Sec_st   = STRING(second)
  
  FOR i = 0, N_ELEMENTS(YEAR)-1 DO BEGIN
    ; Add necessary leading zeroes so time format is the same
    IF DOY_a(i)  LT 10  THEN DOY_st(i) = STRCOMPRESS(STRING('00', DOY_st(i)),/REMOVE_ALL) ELSE $
    IF DOY_a(i)  LT 100 THEN DOY_st(i) = STRCOMPRESS(STRING('0', DOY_st(i)), /REMOVE_ALL)
    IF Hour(i)   LT 10  THEN Hr_st(i)  = STRCOMPRESS(STRING('0', Hr_st(i)),  /REMOVE_ALL)
    IF Minute(i) LT 10  THEN Min_st(i) = STRCOMPRESS(STRING('0', Min_st(i)), /REMOVE_ALL)
    IF Second(i) LT 10  THEN Sec_st(i) = STRCOMPRESS(STRING('0', Sec_st(i)), /REMOVE_ALL)
    Times_st(i) = STRCOMPRESS(STRING(Year(i), '-', DOY_st(i), 'T', Hr_st(i), ':', Min_st(i), ':', Sec_st(i), $
                               FORMAT = '(A,A,A,A,A,A,A,A,A)'), /REMOVE_ALL)
  ENDFOR
  
  ; Switch directories to helpwhere other data is
  CD, 'D:\Documents\LASP\VIPER Code'
  
  ; Get Rj, X, Y, and Z
  Rj_arr = MAKE_ARRAY(N_ELEMENTS(Times_st), /DOUBLE)
  X_arr  = Rj_arr
  Y_arr  = Rj_arr
  Z_arr  = Rj_arr

  FOR i = 0, N_ELEMENTS(Rj_arr)-1 DO BEGIN
    IF i mod 10 EQ 0 THEN PRINT, i
    IF VIPER_Structure.Field06[i] EQ 1 THEN BEGIN
      COMMON VGR1Data
      SSEDR_Data = VGR1SSEDR
    ENDIF ELSE IF VIPER_Structure.Field06[i] EQ 2 THEN BEGIN
      COMMON VGR2Data
      SSEDR_Data = VGR2SSEDR
    ENDIF
    
    ; Time in SSEDR file
    SSEDR_Time = SSEDR_Data[1,*] + (SSEDR_Data[2,*] + SSEDR_Data[3,*]/60d + SSEDR_Data[4,*]/3600d)/24d
    
    ; Time of Input
    Time_Requested = DOY_a(i) + (Hour(i) + Minute(i)/60d + Second(i)/3600d)/24d
    ; Time of Data
    ;      IF VIPER_Structure.L_or_M_Mode EQ 0 THEN Data_Time = Dum_Struct.Ltimes_Array ELSE $
    ;      IF VIPER_Structure.L_or_M_Mode EQ 1 THEN Data_Time = Dum_Struct.Mtimes_Array
    
    ; Need to do some large scale thing and make it work for the VIPER_Structure. On Monday stop that Optional I lie because that is just going to make everything, and I mean everything break
    
    ; Recreate structure with appropriate tags that match initial VIPER input in order to extract X, Y, Z, and Rj
    ; Note: Only time and s/c is needed to pull correct data out of file. Included V1, V2, and V3 in case some variance based off of s/c speed is added
    Structure = {Year: VIPER_Structure.Field01[i], DOY: VIPER_Structure.Field02[i], Hour: VIPER_Structure.Field03[i], Minute: VIPER_Structure.Field04[i], $
                 Second: VIPER_Structure.Field05[i], Spacecraft: VIPER_Structure.Field06[i], Planet_Number: VIPER_Structure.Field07[i], Response: VIPER_Structure.Field08[i], $
                 L_or_M_Mode: VIPER_Structure.Field09[i], Save_Plot: 0,  Fit: 0, Iterations: 25, Channels_CupA: 001128 , Channels_CupB: 001128, $
                 Channels_CupC: 001128, Channels_CupD: 001128, Vary_V1: 0, Vary_V2: 0, Vary_V3: 0, Vary_Density: 0, $
                 Vary_Temperature: 0, Vary_Kappa: 1, Analytical_Vphi: 0d, V1: VIPER_Structure.Field24[i], V2: VIPER_Structure.Field25[i], V3: VIPER_Structure.Field26[i], $
                 Common_Temperature: 0d, Common_Kappa: 0d, Delamere_Composition: 0d, $
                 Species_1: 0, Species_1_A: 1d, Species_1_Z: 1d, Species_1_n: 50d, Species_1_T: 3d, Species_1_Kappa: 2d, $
                 Species_2: 0, Species_2_A: 1d, Species_2_Z: 1d, Species_2_n: 50d, Species_2_T: 3d, Species_2_Kappa: 2d, $
                 Species_3: 0, Species_3_A: 1d, Species_3_Z: 1d, Species_3_n: 50d, Species_3_T: 3d, Species_3_Kappa: 2d, $
                 Species_4: 0, Species_4_A: 1d, Species_4_Z: 1d, Species_4_n: 50d, Species_4_T: 3d, Species_4_Kappa: 2d, $
                 Species_5: 0, Species_5_A: 1d, Species_5_Z: 1d, Species_5_n: 50d, Species_5_T: 3d, Species_5_Kappa: 2d, $
                 Species_6: 0, Species_6_A: 1d, Species_6_Z: 1d, Species_6_n: 50d, Species_6_T: 3d, Species_6_Kappa: 2d, $
                 Species_7: 0, Species_7_A: 1d, Species_7_Z: 1d, Species_7_n: 50d, Species_7_T: 3d, Species_7_Kappa: 2d}
    Structure_Names = TAG_NAMES(Structure)
    
    ; This is unfortunately bulky coding and takes a while to run, since it reads in the whole data file and uses one entry each iteration
    ;;; Fortunately, this code is a primer, and the output file includes ALL relevant information for making plots of the data
    Dum_Struct = PLS_DATA_READ(Structure)
    Data_Time = Dum_Struct.times_Array  
    
    ; Create DOY, Hour, Minute, and Second from Time_Data
    Data_DOY  = FIX(Data_Time)
    Data_Hour = FIX((Data_Time - Data_DOY)*24d)
    Data_Min  = FIX(((Data_Time - Data_DOY)*24d - Data_Hour)*60d)
    Data_Sec  = FIX((((Data_Time - Data_DOY)*24d - Data_Hour)*60d - Data_Min)*60d)
    
    ; Find SSEDR index of data to find S/C location
    SSEDR_ind = CLOSEST(TRANSPOSE(SSEDR_Time[0,*]), Data_Time[0])
    ; Get X, Y, and Z locations of spacecraft and compute radial distance (Rj)
    X_arr(i)  = SSEDR_Data[6,SSEDR_ind] ; System III cartesian coordinates
    Y_arr(i)  = SSEDR_Data[7,SSEDR_ind]
    Z_arr(i)  = SSEDR_Data[8,SSEDR_ind]
    Rj_arr(i) = SQRT(X_arr(i)*X_arr(i)+Y_arr(i)*Y_arr(i)+Z_arr(i)*Z_arr(i))/71492d
  ENDFOR
  
  ; Populate parameters
  V1  = DOUBLE(VIPER_Structure.Field24)
  V2  = DOUBLE(VIPER_Structure.Field25)
  V3  = DOUBLE(VIPER_Structure.Field26)
  
  dnc = DOUBLE(VIPER_Structure.Field33)
  Tec = DOUBLE(VIPER_Structure.Field34)
  Kec = DOUBLE(VIPER_Structure.Field35)
  ; Set fill values for Kappa
  FOR i = 0, N_ELEMENTS(Kec)-1 DO BEGIN
    IF Kec[i] LE 1.5 THEN Kec[i] = -9d
  ENDFOR
  
  
  dnh = DOUBLE(VIPER_Structure.Field39)
  Teh = DOUBLE(VIPER_Structure.Field40)
  
  SC = LONG(VIPER_Structure.Field06)

  
  ; Write CSV file for final output
  
  ; The output format of the CSV file is as follows
  ;;; Time       - UTC string
  ;;; Spacecraft - Long integer of 1 or 2 for Voyager 1 and Voyager 2 respectively
  ;;; Rj         - Radius in jovian radii
  ;;; X          - System III X    (km)
  ;;; Y          - System III Y    (km)
  ;;; Z          - System III Z    (km)
  ;;; dn_cold    - Density of cold electrons (number/cc)
  ;;; T_cold     - Temperature of cold electrons (eV)
  ;;; Kappa      - Kappa value if using Kappa fits. -9 if no kappa fit. There are no kappa values for the "hot" because only one kappa was ever used. Kappa > 3/2
  ;;; dn_hot     - Density of hot electrons (number/cc) - The "hot" values are included if data is being fit with 2 Maxwellians
  ;;; T_hot      - Temperature of hot electrons (eV)
  Output_CSV = {Time:Times_st[*], Spacecraft:SC[*], Rj:Rj_arr[*], X:X_arr[*], Y:Y_arr[*], Z:Z_arr[*], dn_cold:dnc[*], T_cold:Tec[*], Kappa:Kec[*], $
                dn_hot:dnh[*], T_hot:Teh[*]}
  Header = TAG_NAMES(Output_CSV)
  
  ; Change to directory for PDS Output
  CD, 'D:\Documents\LASP\VIPER Code\Master Fits\PDS'
  WRITE_CSV, output_filename, Output_CSV, HEADER = Header
  
  ; Return home
  CD, 'D:\Documents\LASP\VIPER Code'

END
