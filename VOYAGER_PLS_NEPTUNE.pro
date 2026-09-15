FUNCTION VOYAGER_PLS_NEPTUNE, Structure, POD, cupss, EDP, EPL, MAKE_PLOT=MAKE_PLOT, UNCERTAINTIES=UNCERTAINTIES, Structure_Unc=Structure_Unc, $
                      N_CHI_PTS=N_CHI_PTS, N_CHI_SIG=N_CHI_SIG, COMPUTE_CHI_SPACE=COMPUTE_CHI_SPACE, CEq = CEq, WRITE_UNCERTAINTIES_IN_FILE=WRITE_UNCERTAINTIES_IN_FILE, $
                      CALC_FITTING_RESIDUALS = CALC_FITTING_RESIDUALS
  ; Inputs
  ; Structure - A structure containing information on how to handle all parameters in the fit
  ; POD               - Directory name for where to save the plots - Creates an output folder if necessary
  ; PLOT              - Keyword argument to produce a plot of the data
  ; UNCERTAINTIES     - Keyword argument to calculate uncertainties for the given parameters
  ; N_CHI_PTS         - Keyword argument for number of points to use (Creates a square array of that size)
  ; N_CHI_SIG         - Keyword argument for number of sigma to deviate from best fit
  ; Structure_Unc     - Keyword argument with uncertainties of structures
  ; WRITE_UNCERTAINTIES_IN_FILE - if it's not set a new file will be created with the uncertainties 
  ;                     in the place of N and T. If it's set it will write the uncertainties in the same file
  ;                     at the N_UNCERTAINTIES and T_UNCERTAINTIES columns. It needs the \UNCERTAINTIES keyword
  ; CALC_FITTING_RESIDUALS - Calculated the residuals between the fitting and the raw data - GX: I created this, but it seems that mpfit can also calculate residuals

  ; Output
  ; Return_Value - Returns a structure that is of the same size as source expression
  ;                Values will have changed to either better fit parameters or uncertainties in the parameters 
  
  ; Unpack Structure for relevant analysis
  ;
  
  ; Notes: V2_NEPTUNE_PLS is the frame closest to the distances of the SEDR file. It was created by Rob Wilson.
  
  ; Written by: Kaleb Bodisch & Logan Dougherty
  ; Adapted for Neptune by: George Xystouris
  
  
  ; Neptune physical characteristics to be used further down the script
  r_Neptune = 24764d     ; Radius of Neptune [km] (acquired using SPICE bodvrd command)
  u_Rotation = 2.683d  ; Rotation speed of Neptune [km/s] (acquired using SPICE deg_rate = cspice_bodvrd('NEPTUNE', 'PM', 3 ) command. deg_rate(2) gives the degrees the planet rotated in 86400s) 
  
  
  Year  = Structure.Year
  DOY   = Structure.DOY
  Hour  = Structure.Hour
  Minute = Structure.Minute
  Second = Structure.Second
  Spacecraft = Structure.Spacecraft
  IF Structure.L_or_M_Mode EQ 1  THEN n_channels = 128
  IF structure.L_or_M_Mode EQ 0  THEN n_channels = 16
  IF structure.L_or_M_Mode EQ -1 THEN n_channels = 28
  
  ; Define Forward Functions so IDL doesn't crash
  FORWARD_FUNCTION PLS_Data_Read
  FORWARD_FUNCTION Generate_PLS_Current_Neptune
  FORWARD_FUNCTION ADDITIONAL_PLS_FUNCTIONS
  FORWARD_FUNCTION CONRD
  FORWARD_FUNCTION Extrapolate_to_Equator  ;it's not working for Neptune at the moment
  
  ; Read in the Data for the timestamp and spacecraft given
  Dum_Struct = PLS_Data_Read(Structure) ; Works similar to Current_Plot_3 I believe (LOGAN 1/5/15)
  ;STOP
  
  ; Combine data into one array, organized as a 4 x Number of Channels array
  PLS_Data_Transp = [Dum_Struct.current_array_A, Dum_Struct.current_array_B, Dum_Struct.current_array_C, Dum_Struct.current_array_D]
  ; Combine data into one array of 1 x (4 x Number of Channels) - It goes A, B, C, and then D cup
  PLS_Data = [TRANSPOSE(PLS_Data_Transp[0,*]),TRANSPOSE(PLS_Data_Transp[1,*]), $
              TRANSPOSE(PLS_Data_Transp[2,*]),TRANSPOSE(PLS_Data_Transp[3,*])]

; ERROR DEFINITION (following Uranus):
;
; Background + Measurement uncertainty: Poisson (square root of measurement)
; Background error is defined as the thermal noise of the electronics. Here we are using the thermal noise 
; as given in Memo 178, during the pre-Neptune strong flyby.
  IF Structure.L_or_M_Mode EQ 1 THEN PLS_Error = SQRT( 40d^2 + PLS_Data ) ELSE BEGIN
    CASE Structure.Short_or_Long OF
      1 : PLS_Error = SQRT( 9.49d^2 + PLS_Data )
      0 : PLS_Error = SQRT( 20d^2 + PLS_Data )
    ENDCASE
  ENDELSE


  ;------------------ GX CHANGE FROM HERE
  ; Extrapolate to the Centrifugal Equator                      ; FOR NOW WE DON'T CARE ABOUT THIS!
  IF KEYWORD_SET(CEq) THEN BEGIN
    ; Get location of S/C in V2_NEPTUNE_PLS coordinates            
    COMMON VGR2Data_NEPTUNE
    SSEDR_Data = VGR2SSEDR_NEPTUNE

    ; Time in SSEDR file                        THEY ARE GETTING THE CLOSEST TIMES BETWEEN THE TWO (S/C & PLS)
    SSEDR_Time = SSEDR_Data[1,*] + (SSEDR_Data[2,*] + SSEDR_Data[3,*]/60d + SSEDR_Data[4,*]/3600d)/24d
    
    ; Time of Measurement according to PLS instrument
    Data_Time = Dum_Struct.times_Array

    ; Find SSEDR index of data to find S/C location
    SSEDR_ind = CLOSEST(TRANSPOSE(SSEDR_Time[0,*]), Data_Time[0])
    ; Get X, Y, and Z locations of spacecraft and compute radial distance (r_RU)
    xSC_V2_NEPTUNE_PLS = SSEDR_Data[6,SSEDR_ind] ; V2_NEPTUNE_PLS cartesian coordinates - see documentation or paper
    ySC_V2_NEPTUNE_PLS = SSEDR_Data[7,SSEDR_ind]
    zSC_V2_NEPTUNE_PLS = SSEDR_Data[8,SSEDR_ind]

    Position_V2_NEPTUNE_PLS = [xSC_V2_NEPTUNE_PLS, ySC_V2_NEPTUNE_PLS, zSC_V2_NEPTUNE_PLS] 
    New_Structure = Extrapolate_To_Equator(Structure, Position_V2_NEPTUNE_PLS, r_Neptune, /ISOTROPIC)
    RETURN, New_Structure 

  ENDIF
  ;------------------------- TO HERE
  
  
  ; Convert Electron data to Reduced Distribution and fit this
  ; (i.e. PHASE SPACE DENSITY FOR 1 DIRECTION - 1 CUP).
  ; Following methodology in Sittler 1983 Technical memorandum (85037)
  ;
  IF Structure.L_or_M_Mode EQ -1 THEN BEGIN
    PLS_electron_data = PLS_Data
    PLS_electron_data_Transp = PLS_Data_Transp ; measured current from raw data
    
    ; Create grid edges potentials
    phij = MAKE_ARRAY(29)
    FOR i = 0, 28 DO BEGIN
      IF i LE 16 THEN BEGIN
        phij(i) = 60.*10.^(double(i)/32.) - 50.
      ENDIF ELSE IF i GE 17 THEN BEGIN
        phij(i) = 60.*10.^(double(i-12)/8.) - 50.
      ENDIF
    ENDFOR
    
    ; Constants
    ECHRGE = 1.6d-19                       ; Electron charge in Coulombs  !!!!!!!!!!!!! CHECK WHETHER YOU NEED TO ADD THE MINUS!!!!!!!!!!
    me = 9.11d-31                           ; Electron mass in kg
    Aeff = 100/1d4                          ; Area of sensor (cm^2) converted to m^2
    TN = 0.56d                              ; normal transmision of sensor (see section 2 of Sittler manual)
    ;stop
    vj = SQRT(2 * ECHRGE * phij / (me))    ; Equivalent electron speed for each grid edge in m/s

    ; Create appropriately sized arrays
    dvj = MAKE_ARRAY(28)
    vjavg = dvj                             ; average electron speed for each grid
    alphaj = dvj                            ;  
    test = dvj
    ; For each step of the grid we calculate:
    FOR loopi = 0, 27 DO BEGIN
      dvj(loopi) = vj(loopi+1) - vj(loopi)          ; electrons speed for this grid
      vjavg(loopi) = 0.5d*(vj(loopi+1) + vj(loopi)) ; average electron speed for each grid
      alphaj(loopi) = ECHRGE * Aeff * TN * !Pi * vjavg(loopi) * vjavg(loopi) *vjavg(loopi) * dvj(loopi)  ; this is the denominator in eq. 23
    ENDFOR
    alphaj = [alphaj,alphaj,alphaj,alphaj]
    ; This is the observed reduced distribution function (eq. 23). The factor 10^-15 is conversion of the measured current from fA to A
    PLS_DATA = (PLS_DATA * 1d-15)/alphaj
    PLS_DATA_Transp[3,*] = (PLS_DATA_Transp[3,*] * 1d-15)/alphaj
    PLS_DATA_Transp_ORDF = PLS_DATA_Transp
    PLS_ERROR_ORDF = (PLS_ERROR * 1d-15)/alphaj
    PLS_ERROR_ORDF = [PLS_ERROR_ORDF[-28:*],PLS_ERROR_ORDF[-28:*],PLS_ERROR_ORDF[-28:*],PLS_ERROR_ORDF[-28:*]] ; MPFIT needs to have the error for all 4 cups, even if A, B, and C have no data
    
    ; FLUX CONVERSION - FILLER
    Delta_E = MAKE_ARRAY(28) ; energy difference for each grid
    Eenergy = MAKE_ARRAY(28) ; average energy of each grid
    FOR kloop = 0, 27 DO BEGIN
      if kloop LE 15 then begin
        Delta_E(kloop) = ((60.*(10.^((double(kloop)+1.)/32.))-50) - (60.*(10.^((double(kloop))/32.))-50))
        Eenergy(kloop) = (((60.*(10.^((double(kloop)+1.)/32.))-50) + (60.*(10.^((double(kloop))/32.))-50)))/2.
      endif else if kloop EQ 16 then begin
        Delta_E(kloop) = (60.*(10.^(5./8.))-50.)-((60.*(10.^(16./32.))-50.))
        Eenergy(kloop) = ((60.*(10.^(5./8.))-50.)+((60.*(10.^(16./32.))-50.)))/2.
      endif else if kloop GE 17 then begin
        Delta_E(kloop) = ((60.*(10.^((double(kloop)+1.-12.)/8.))-50) - (60.*(10.^((double(kloop)-12.)/8.))-50))
        Eenergy(kloop) = (((60.*(10.^((double(kloop)+1.-12.)/8.))-50) + (60.*(10.^((double(kloop)-12.)/8.))-50)))/2.
      endif
    endfor
    ;stop
    dOmega = 2.8                                          ; ??? Solid angle? it corresponds to an half-angle cone of around 56.33 degrees
    FluxConversion = Aeff*dOmega*ECHRGE*TN*Delta_E/1d-15  ; ??? Factor that converts current to flux. Flux = I / this_factor
  ENDIF

 
    ; Calculation of an electron powerlaw (EPL) distribution
  IF EPL EQ 1 THEN BEGIN
    PLS_DATA_FLUX = PLS_DATA_Transp_ORDF[3,*]*alphaj/1d-15/FluxConversion[*] ; flux converted for current in fA
    indices_logfit = lindgen(16)+6
    ind_remove_log = WHERE(PLS_DATA_FLUX EQ 0, nind_remove_log)
    IF nind_remove_log GT 0 THEN BEGIN
      PLS_DATA_FLUX(ind_remove_log) = !VALUES.D_NAN 
      Eenergy(ind_remove_log) = !VALUES.D_NAN 
    ENDIF
    indx = WHERE(ind_remove_log LT 22, nindx)
    IF nindx GT 0 THEN indices_logfit = lindgen(16-nindx)+6
    ylogdata = ALOG(PLS_DATA_FLUX[indices_logfit])
    xlogdata = ALOG(Eenergy[indices_logfit])
    measure_errors = SQRT(ABS(ylogdata))
    result_EPL = linfit(xlogdata, ylogdata, measure_errors=measure_errors)
    alin = result_EPL[0]
    blin = result_EPL[1]
    print, result_epl
    aexp = EXP(alin)
    
    yfitlog = aexp * Eenergy^(blin)
    IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN
      pt = plot(Eenergy, PLS_DATA_FLUX, color = 'black', /ylog, /xlog, yrange = [1d3,1d8], xrange = [10,6000])
      opt = plot(Eenergy, yfitlog, color = 'red', /OVERPLOT)
    ENDIF

;------------------ GX CHANGE FROM HERE
    ; NOW FOR THE RETURN FILE - I CAN MAKE IT LIKE THE EDP PRODUCT
    ; Get constants for saving data
      COMMON VGR2Data_NEPTUNE
      SSEDR_Data = VGR2Data_NEPTUNE

    ; Time in SSEDR file
    SSEDR_Time = SSEDR_Data[1,*] + (SSEDR_Data[2,*] + SSEDR_Data[3,*]/60d + SSEDR_Data[4,*]/3600d)/24d

    ; Time of Input
    Time_Requested = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
    ; Time of Data
    ;      IF Structure.L_or_M_Mode EQ 0 THEN Data_Time = Dum_Struct.Ltimes_Array ELSE $
    ;      IF Structure.L_or_M_Mode EQ 1 THEN Data_Time = Dum_Struct.Mtimes_Array

    Data_Time = Dum_Struct.times_Array

    ; Create DOY, Hour, Minute, and Second from Time_Data
    ; LOGAN 1/7/15 - Check to see if this works
    Data_DOY  = FIX(Data_Time)
    Data_Hour = FIX((Data_Time - Data_DOY)*24d)
    Data_Min  = FIX(((Data_Time - Data_DOY)*24d - Data_Hour)*60d)
    Data_Sec  = FIX((((Data_Time - Data_DOY)*24d - Data_Hour)*60d - Data_Min)*60d)

    ; Find SSEDR index of data to find S/C location
    SSEDR_ind = CLOSEST(TRANSPOSE(SSEDR_Time[0,*]), Data_Time[0])
    ; Get X, Y, and Z locations of spacecraft and compute radial distance (R)
    xSC_V2_NEPTUNE_PLS = SSEDR_Data[6,SSEDR_ind] ; V2_NEPTUNE_PLS cartesian coordinates
    ySC_V2_NEPTUNE_PLS = SSEDR_Data[7,SSEDR_ind]
    zSC_V2_NEPTUNE_PLS = SSEDR_Data[8,SSEDR_ind]
    r_RU = SQRT(xSC_V2_NEPTUNE_PLS*xSC_V2_NEPTUNE_PLS+ySC_V2_NEPTUNE_PLS*ySC_V2_NEPTUNE_PLS+zSC_V2_NEPTUNE_PLS*zSC_V2_NEPTUNE_PLS)/r_Neptune
      
    xSC_RN = xSC_V2_NEPTUNE_PLS/r_Neptune
    ySC_RN = ySC_V2_NEPTUNE_PLS/r_Neptune
    zSC_RN = zSC_V2_NEPTUNE_PLS/r_Neptune
    SC_lat = asin(zSC_RU/r_RU)
    SC_lon = 180./!Pi*ATAN(ySC_RU,xSC_RU)
  
    ; Calculate L-Shell              ; GX: We don't care about this right now for Neptune 
    ; These 2 csv files are from Fran's website                                                       GOOD TO KEEP IN THE FINAL RESULT - BUT WHAT L...?
    V1Mag = READ_CSV('Spacecraft_trajectory_vip4_can.csv')
    V2Mag = READ_CSV('Spacecraft_trajectory_vip4_can_V2.csv')
    MagTimes = [V1Mag.Field01,V2Mag.Field01]
    LShellVals = [V1Mag.Field05,V2Mag.Field05]
    Data_Dec_Time = DOUBLE(Structure.DOY) + DOUBLE(Structure.Hour)/24. + DOUBLE(Structure.Minute)/(24.*60.) + DOUBLE(Structure.Second)/(24.*3600.)
    indMag = closest(MagTimes[*],Data_Dec_Time)
    LShell = LShellVals[indmag]
;------------------ TO HERE
      
    PDF = PLS_DATA_FLUX
    ; Create the output structure
    EPL_Structure = {Year: 1979, DOY: Structure.DOY, Hour: Structure.Hour, Minute: Structure.Minute, Second: Structure.Second, $
                       x_SC:xSC_RU, y_SC:ySC_RU, z_SC:zSC_RU, R: r_RU, Latitude: SC_lat*180d/(!Pi), Longitude: SC_lon*180d/(!Pi), $
                       LShell: LShell, a_coefficient: aexp, b_coefficient: blin, $
                       Channel01:PDF[0], Channel02:PDF[1], Channel03:PDF[2], Channel04:PDF[3], Channel05:PDF[4], Channel06:PDF[5], $
                       Channel07:PDF[6], Channel08:PDF[7], Channel09:PDF[8], Channel10:PDF[9], Channel11:PDF[10], Channel12:PDF[11], $
                       Channel13:PDF[12], Channel14:PDF[13], Channel15:PDF[14], Channel16:PDF[15], Channel17:PDF[16], Channel18:PDF[17], $
                       Channel19:PDF[18], Channel20:PDF[19], Channel21:PDF[20], Channel22:PDF[21], Channel23:PDF[22], Channel24:PDF[23], $
                       Channel25:PDF[24], Channel26:PDF[25], Channel27:PDF[26], Channel28:PDF[27]}
    
    RETURN, EPL_Structure
    
  ENDIF
  
  
  
  ; If fitting procedure is turned on then fit the data
  IF Structure.Fit EQ 1 THEN BEGIN
    ; Remove data points that are not being fit (sets them to NaN, which MPFIT ignores)
    ; Create dummy array called Current_Holder
    Current_Holder = MAKE_ARRAY(n_channels*4,/DOUBLE)
    Current_Holder[*] = !VALUES.D_NAN        ; Initially populate with NaN values and then return true values
    
    IF Structure.L_or_M_Mode GE 0 THEN BEGIN
      ; Use Good_Channels function to find gchn 
      ; LOGAN 1/5/15 - MAKE SURE TO CHANGE MODE TO ACCOUNT FOR 0 AND 1 FOR L AND M
      cup = 0
      gchn = Good_Channels(long(Structure.Channels_CupA), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupA GT 0 THEN Current_Holder[gchn[0]-1 + 0*n_channels:gchn[1]-1 + 0*n_channels] $
                                                          = PLS_Data_Transp[0, gchn[0]-1:gchn[1]-1] 
      cup = 1
      gchn = Good_Channels(long(Structure.Channels_CupB), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupB GT 0 THEN Current_Holder[gchn[0]-1 + 1*n_channels:gchn[1]-1 + 1*n_channels] $
                                                          = PLS_Data_Transp[1, gchn[0]-1:gchn[1]-1]
      cup = 2
      gchn = Good_Channels(long(Structure.Channels_CupC), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupC GT 0 THEN Current_Holder[gchn[0]-1 + 2*n_channels:gchn[1]-1 + 2*n_channels] $
                                                          = PLS_Data_Transp[2, gchn[0]-1:gchn[1]-1]
      cup = 3
      gchn = Good_Channels(long(Structure.Channels_CupD), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupD GT 0 THEN Current_Holder[gchn[0]-1 + 3*n_channels:gchn[1]-1 + 3*n_channels] $
                                                          = PLS_Data_Transp[3, gchn[0]-1:gchn[1]-1]
   ENDIF ELSE BEGIN 
      ; ELECTRONS: using the observed electron distribution function
      cup = 0
      gchn = Good_Channels(long(Structure.Channels_CupA), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupA GT 0 THEN Current_Holder[gchn[0]-1 + 0*n_channels:gchn[1]-1 + 0*n_channels] $
                                                          = PLS_DATA_Transp_ORDF[0, gchn[0]-1:gchn[1]-1] 
      cup = 1
      gchn = Good_Channels(long(Structure.Channels_CupB), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupB GT 0 THEN Current_Holder[gchn[0]-1 + 1*n_channels:gchn[1]-1 + 1*n_channels] $
                                                          = PLS_DATA_Transp_ORDF[1, gchn[0]-1:gchn[1]-1]
      cup = 2
      gchn = Good_Channels(long(Structure.Channels_CupC), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupC GT 0 THEN Current_Holder[gchn[0]-1 + 2*n_channels:gchn[1]-1 + 2*n_channels] $
                                                          = PLS_DATA_Transp_ORDF[2, gchn[0]-1:gchn[1]-1]
      cup = 3
      gchn = Good_Channels(long(Structure.Channels_CupD), cup, Structure.L_or_M_Mode)
      IF structure.Channels_CupD GT 0 THEN Current_Holder[gchn[0]-1 + 3*n_channels:gchn[1]-1 + 3*n_channels] $
                                                          = PLS_DATA_Transp_ORDF[3, gchn[0]-1:gchn[1]-1]
                                                          
                                                          
   ENDELSE
   

; --------------------- CHANGE STUFF FROM HERE (IF NECESSARY)    
    ; THRESHOLD SET AND CHANNELS REMOVAL
    ;   
    ; Set low values to NaNs so that they are not fit in the process
    ; For L and E modes: we set the background depending on whether the specific scan 
    ; has a short or long integration time. For short integration time the background 
    ; is 103.726fA, and for the long it is 3.055fA.
    ; For M mode: there is only one integration time, at 30.549fA.
    ; 
    ; More on the background in Memo 178. While it is for Neptune, the modes and thresholds
    ;  are the same - GX 16/Feb/24
    ;
    Current_Threshold = 3.1d ; set default, changed below if required
    IF TOTAL(STRCMP(TAG_NAMES(Structure),'Short_or_Long', /FOLD_CASE))  EQ 1  THEN BEGIN
      IF Structure.L_or_M_Mode GE 0 THEN BEGIN
        IF Structure.L_or_M_Mode EQ 1 THEN Current_Threshold = 23.421d ELSE BEGIN
          CASE Structure.Short_or_Long OF ; short = 0, long = 1
            0 : Current_Threshold = 103.726d
            1 : Current_Threshold = 3.055d
            ELSE : PRINT, 'Warning: integration time has not been set properly (needs to be 0 for short or 1 for long). Continuing with current threshold (background) set at '+STRTRIM(Current_Threshold,2)+'fA.'
          ENDCASE
        ENDELSE  
      ENDIF ELSE BEGIN
        CASE Structure.Short_or_Long OF ; short = 0, long = 1
          ; for the electrons the threshold error is converted based on the calculated RDF (eq. 23 in Sittler's manual, with I = threshold current)
          0 : Current_Threshold = (103.726d * 1d-15)/alphaj
          1 : Current_Threshold = (3.055d * 1d-15)/alphaj
        ENDCASE
      ENDELSE
    ENDIF ELSE PRINT, 'Warning: integration time has not been set. Continuing with current threshold (background) set at '+STRTRIM(Current_Threshold,2)+'fA.'

    ; Add threshold current to fitting matrix
    ; If we have electrons, we are re-defining the threshold, because above it was defined based on the RDF
    Thresh_StringCompare = STRCMP('THRESHOLD_CURRENT', TAG_NAMES(Structure), /FOLD_CASE)
    Thresh_ind = WHERE(Thresh_StringCompare EQ 1, nind, /NULL)
    IF nind EQ 1 THEN BEGIN
      struct_current = Current_Threshold
      IF Structure.L_or_M_Mode EQ -1 THEN BEGIN
        CASE Structure.Short_or_Long OF
          0 : struct_current = 103.726d
          1 : struct_current = 3.055d
        ENDCASE
      ENDIF
      Structure.(Thresh_ind) = struct_current
    ENDIF ELSE PRINT,'Warning: Threshold_Current column not found.' 
; --------------------- TO HERE
    
    ; Remove all channels under the threshold
    ind_remove = WHERE((Current_Holder LE Current_Threshold) or (FINITE(Current_Holder) EQ 0), /NULL) ; it gives a "Program caused arithmetic error: Floating illegal operand" warning that can be ignored. It's because of the NaN in the Current_Holder
    
    ; Remove channels 1, 16, 17, 18, and 28 from D Cup for Electron mode fitting
    ; 1 and 28 are the beginning and end and are a little discontinuous
    ; 16, 17, and 18 are where E1 and E2 modes cross over and there is a dicontinuity
    ; ** GX: I don't understand why we're removing channel 18 
    IF Structure.L_or_M_Mode EQ -1 THEN BEGIN
      ind_remove = [ind_remove,84,99,100,111]
      ; remove duplicate channels
      sort_ind_remove = SORT(ind_remove)
      sorted_arr = ind_remove[sort_ind_remove]
      unique_ind_remove = UNIQ(sorted_arr)
      ind_remove = sorted_arr[unique_ind_remove]
    ENDIF
    
    
    ; Add total number of fitted channels to fitting matrix
    Chann_num_StringCompare = STRCMP('TOTAL_NUMBER_OF_FITTED_CHANNELS', TAG_NAMES(Structure), /FOLD_CASE)
    Chann_ind = WHERE(Chann_num_StringCompare EQ 1, nind, /NULL)
    ;IF Structure.L_or_M_Mode EQ -1 THEN channels_number = 28 ELSE channels_number = N_ELEMENTS(Current_Holder)
    channels_number = N_ELEMENTS(Current_Holder) ; Logan 2025_03_06 - As we now count the NaN channels, E-Modes are now counted as 4*28 total channels, instead of just D-Cup
    IF nind EQ 1 THEN Structure.(Chann_ind) = channels_number - N_ELEMENTS(ind_remove) ELSE PRINT,'Warning: Total_Number_of_Fitted_Channels column not found.'
    
    ; Remove bad datapoints. These points will also be removed from the calculation for the Chi-Space and the uncertainties
    IF ind_remove NE !NULL THEN Current_Holder[ind_remove] = !VALUES.D_NAN
    
    ; Save to array for MPFIT to use
    PLS_Data_Fit = Current_Holder
    
    ; Create New Structure to house where parameters are located
    P_Locations = Structure
    FOR j = 0, N_TAGS(Structure) - 1 DO P_Locations.(j) = 0d ; Initially set all to zero
    
    ; Determine what velocities are active
    IF Structure.Vary_V1 EQ 1d THEN P_Locations.V1 = 1d
    IF Structure.Vary_V2 EQ 1d THEN P_Locations.V2 = 1d
    IF Structure.Vary_V3 EQ 1d THEN P_Locations.V3 = 1d
    
    ; Determine if both density and flowspeeds are allowed to vary
    IF Structure.Vary_Density EQ 1d THEN Vary_Density = 1d ELSE Vary_Density = 0d
    IF Structure.Vary_Temperature EQ 1d THEN Vary_Temperature = 1d
    IF Structure.Common_Temperature GT 0d AND Structure.Vary_Temperature NE 0 THEN BEGIN 
      Common_Temp = Structure.Common_Temperature
      P_Locations.Common_Temperature = 1d
    ENDIF
    
    ; If using a Common Kappa fit only that
    ; Species Mode of 6 will overwrite that but will default to individual densities and common temperature
    IF Structure.Common_Kappa GT 0d AND Structure.Vary_Kappa NE 0 THEN BEGIN
      Common_Kappa = Structure.Common_Kappa
      P_Locations.Common_Kappa = 1d
    ENDIF
    
    ; Determine how many species there are on the csv file
    taggg_names = TAG_NAMES(Structure)
    StringCompare = STRCMP('Species_', taggg_names, 8, /FOLD_CASE)
    ind = WHERE(StringCompare EQ 1, nind)
    IF nind NE 0 THEN BEGIN ; this if replaced n_species = nind/6  that was initially here for jupiter. We made it so you can have as many columns as you want for variables
      species_columns = TOTAL(STRCMP(taggg_names,'SPECIES_1',9, /FOLD_CASE))
      n_species = FLOAT(nind)/species_columns ; the float for nind is there to convert it from integer. Weird things happen if you divide integer by integer in IDL
      IF (n_species NE ROUND(n_species)) THEN MESSAGE, 'Error: Species do not have the same number of columns.'
      n_species = ROUND(n_species)
    ENDIF ELSE n_species = 0
    Arguments = MAKE_ARRAY(n_species) ; Stores the arguments
    
    
    n_tied = 0 ; Counter for number of tied density parameters
    FOR j = 1, n_species DO BEGIN
      str1 = STRCOMPRESS(STRING('Species_', j), /REMOVE_ALL)
      StringCompare = STRCMP(str1, taggg_names, 10, /FOLD_CASE) ; Extract all variables/parameters for the specific species
      ; If it is species 10 or greater then there is one more letter in the variable name for IDL to check
      IF n_species GE 10 THEN StringCompare = STRCMP(str1, taggg_names, 11, /FOLD_CASE)
      ind = WHERE(StringCompare EQ 1, nind)
      ;stop
      IF nind LT 1 THEN MESSAGE, 'Error: No species were provided for analysis. Quitting program.'
      Arguments(j-1) = Structure.(ind)

      
      ; If 0, then it is not in parameter space, either density or temperature (will follow its own temperature if common
      ; temperature isn't provided)
      IF Structure.(ind) EQ 0 THEN CONTINUE ELSE $
        
      ; If 1, then it is not in parameter space, and explicitly follows its own temperature
      IF Structure.(ind) EQ 1 THEN CONTINUE ELSE $
        
      ; If 2, then it is in parameter space, it will follow common temperature, or its own if common isn't provided
      IF Structure.(ind) EQ 2 THEN BEGIN
        ; Put density into parameter space
        str2 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d
        
        ; Check for common temperature tag - If not using common temperature, then include temperature in parameter space
        IF Structure.Common_Temperature EQ 0d THEN BEGIN
          str3 = STRCOMPRESS(STRING('Species_', j, '_T'), /REMOVE_ALL)
          str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE)
          ind3 = WHERE(str3cmp EQ 1)
          P_Locations.(ind3) = 1d
        ENDIF 
      ENDIF ELSE IF Structure.(ind) EQ 3 THEN BEGIN
        ; Put density into parameter space
        str2 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d

        ; Put temperature into parameter space
        str3 = STRCOMPRESS(STRING('Species_', j, '_T'), /REMOVE_ALL)
        str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind3 = WHERE(str3cmp EQ 1)
        P_Locations.(ind3) = 1d
      ENDIF ELSE IF (Structure.(ind) EQ 4) AND (Vary_Density EQ 1) THEN BEGIN
        ; Only the first tied density is a parameter in parameter space
        IF n_tied GT 0d THEN CONTINUE
        n_tied = n_tied + 1
        str2 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d
      ENDIF ELSE IF Structure.(ind) EQ 5 THEN BEGIN
        ; Put temperature into parameter space
        str2 = STRCOMPRESS(STRING('Species_', j, '_T'), /REMOVE_ALL)
        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
        ind2 = WHERE(str2cmp EQ 1)
        P_Locations.(ind2) = 1d
        
        ; Put density into the parameter space
        IF n_tied GT 0d THEN CONTINUE
        IF Vary_Density EQ 1 THEN BEGIN
          n_tied = n_tied + 1
          str3 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
          str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE)
          ind3 = WHERE(str3cmp EQ 1)
          P_Locations.(ind3) = 1d
        ENDIF
      ENDIF ELSE IF Structure.(ind) EQ 6 THEN BEGIN
        ; Put density into parameter space
        IF Vary_Density NE 0 THEN BEGIN
          str2 = STRCOMPRESS(STRING('Species_', j, '_N'), /REMOVE_ALL)
          str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE)
          ind2 = WHERE(str2cmp EQ 1)
          P_Locations.(ind2) = 1d
        ENDIF

        ; Check for common temperature tag - If not using common temperature, then include temperature in parameter space
        IF Structure.Common_Temperature EQ 0d THEN BEGIN
          str3 = STRCOMPRESS(STRING('Species_', j, '_T'), /REMOVE_ALL)
          str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE)
          ind3 = WHERE(str3cmp EQ 1)
          P_Locations.(ind3) = 1d
        ENDIF
        
        IF P_Locations.Common_Kappa GT 0d THEN BEGIN
          P_Locations.Common_Kappa = 0d ; Overwrite if somehow using a common value for kappa
          Structure.Common_Kappa = 0d
        ENDIF
        ; Note that this is perhaps bulkier than with temperatures, but temperatures were always for multiple species where this could be for just 1 and I don't want weird errors in the code
        str4 = STRCOMPRESS(STRING('Species_', j, '_Kappa'), /REMOVE_ALL)
        str4cmp = STRCMP(str4, TAG_NAMES(Structure), 16, /FOLD_CASE)
        ind4 = WHERE(str4cmp EQ 1)
        P_Locations.(ind4) = 1d
      ENDIF
    ENDFOR
    
    ; Make sure that initial kappa values are at least 1.501
    IF ISA(ind4) EQ 1 THEN BEGIN
      IF Structure.(ind4) LE 1.5 THEN Structure.(ind4) = 1.501
    ENDIF
    ; This is limiting the amount of fits that are provided. Chi-squared space is determined to push Kappa below 1.5, which isn't a physical result
    ; Therefore there are more Maxwellian fits than Kappa fits
    
    
    ; Save information to common blocks so it can be accessed inside of function that MPFIT calls
    Save_Structures, Struct_In_A=Structure, Struct_In_B=P_Locations
    
    ; Determine initial parameters
    dum_arr = MAKE_ARRAY(N_TAGS(Structure))
    FOR k = 0, N_TAGS(Structure) - 1 DO dum_arr(k) = P_locations.(k)
    ind = WHERE(dum_arr EQ 1)
    P = MAKE_ARRAY(N_ELEMENTS(ind), /DOUBLE)
   
    ; Generate parinfo structure for MPFIT to limit densities and temperatures to be above 0
    parinfo = REPLICATE({limited:[0,0],limits:[0.D,0]}, N_ELEMENTS(ind))

    FOR j = 0, N_ELEMENTS(ind) - 1 DO BEGIN
      P[j] = DOUBLE(Structure.(ind[j]))
      ; Check where density and temperature tags are set so that MPFIT can set limits and prevent values
      ; from going negative

      str1_check = '*_N*'
      str2_check = '*_T*'
      str3_check = '*_Kappa*'
      
      str1_check_match = STRMATCH((TAG_NAMES(Structure))[ind[j]], str1_check, /FOLD_CASE)
      str2_check_match = STRMATCH((TAG_NAMES(Structure))[ind[j]], str2_check, /FOLD_CASE)
      str3_check_match = STRMATCH((TAG_NAMES(Structure))[ind[j]], str3_check, /FOLD_CASE)

      ; Check that it is a temperature or it is a density - Constrain it to above 0
      IF str1_check_match EQ 1 OR str2_check_match EQ 1 THEN BEGIN
        parinfo[j].limited[0] = 1
        ; Logan 2025_08_13 - Even lower limit on this for electrons
        parinfo[j].limits[0] = FLOAT(0.0001) ; Lowered it to 0.001 from 0.01 for electron fitting - I think it's hitting a snag on the hot densities and trying to blow up the temperature with them
      ENDIF
      
      ; Kappa values must be above 3/2 to have temperature result represent something physical
      IF str3_check_match EQ 1 THEN BEGIN
        parinfo[j].limited[0] = 1
        parinfo[j].limits[0] = FLOAT(1.501)
      ENDIF
      
      ; Constrain temperatures to be under a certain value
      ;IF str2_check_match EQ 1 THEN BEGIN
      ;  parinfo[j].limited[1] = 1
      ;  parinfo[j].limits[1]  = FLOAT(300) ; 300 eV limit on temperatures in the fitting procedure
      ;ENDIF
    ENDFOR
    
    
 
 
  ; PLASMA PARAMETERS CALCULATIONS
  ; --------------------------------------
     ; Run through MPFIT and calculate the best parameters
  
     ; Fit with MPFIT - Use all float for this part of the calculation - LOGAN 1/6/15 MPFITFUN does not seem to work with DOUBLE
     ; precision. Parameters are not changed how they should be for some reason. The flow velocities seem like they can change, yet
     ; everything else is completely stuck and does not vary. Simply changing to FLOAT precision allows for MPFIT to fit the
     ; parameters as it should, so it is not a limit that is being put on, but rather a limitation of MPFIT. It may be due to the
     ; range of the data for the PLS instrument, but it shouldn't be a coding issue on the end of the VIPER routines
  
     ; change stuff from here
     IF Structure.L_or_M_Mode NE -1 THEN BEGIN     ; IONS
       Best_P = MPFITFUN('Fit_Parameters_Neptune', FLOAT(MAKE_ARRAY(4*n_channels, /DOUBLE)), FLOAT(PLS_Data_Fit), $
       FLOAT(PLS_Error), FLOAT(P), /NOCATCH, /NAN, PARINFO=parinfo, MAXITER=Structure.Iterations, BESTNORM=Chi_Final) ; QUIET Keyword to suppress MPFIT from printing
     ENDIF ELSE BEGIN     ; ELECTRONS: as input it uses the observed reduced distribution function
       Best_P = MPFITFUN('Fit_Parameters_Neptune', FLOAT(MAKE_ARRAY(4*n_channels, /DOUBLE)), FLOAT(PLS_Data_Fit), $
        FLOAT(PLS_ERROR_ORDF), FLOAT(P), /NOCATCH, /NAN, PARINFO=parinfo, MAXITER=Structure.Iterations, BESTNORM=Chi_Final) ; QUIET Keyword to suppress MPFIT from printing
     ENDELSE
     ; Now put the best parameters and chi square back into the structure
     FOR j = 0, N_ELEMENTS(ind) - 1 DO Structure.(ind[j]) = Best_P[j]
     chi_column = STRCMP('Chi_Square', TAG_NAMES(Structure), /FOLD_CASE)
     IF TOTAL(chi_column) EQ 1 THEN Structure.Chi_Square = Chi_Final
 
 
  ; UNCERTAINTIES CALCULATIONS
  ; --------------------------------------
    ; Need to set the parinfo command so that it doesn't set densities and temperatures to zero
    IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
      str_date = STRCOMPRESS(STRING(ULONG(Structure.Year), ':', ULONG(Structure.DOY), ':', ULONG(Structure.Hour), ':', $
                                    ULONG(Structure.Minute), '::', ULONG(Structure.Second), $
                                    FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02)'), /REMOVE_ALL)
      PRINT, STRCOMPRESS(STRING('Calculating parameter uncertainties for ', str_date))
      ; Determine what data points should be ignored in Chi-Square calculation
      index = WHERE(FINITE(PLS_Data_Fit) EQ 0, nindex) ; Indexes of calculation to avoid
      IF nindex EQ 0 THEN index = -9d ; Set as a filler number so that it will account for all data if none is bad
      Delta = 0.01d ; A 1 percent variation if the relative keyword is set
      Model_Name = 'Generate_PLS_Current_Neptune'
      ;
      ; IONS: currents as inputs
      ; ELECTRONS: ERDF as inputs (and errors connected to this)
      CASE Structure.L_or_M_Mode OF
        1 :   PLS_Unc = PLS_UNCERTAINTIES(Structure, P_Locations, PLS_Data_Fit, PLS_Error, Model_Name, index, Delta, Current_Threshold, /RELATIVE) ; 1xn array
        0 :   PLS_Unc = PLS_UNCERTAINTIES(Structure, P_Locations, PLS_Data_Fit, PLS_Error, Model_Name, index, Delta, Current_Threshold, /RELATIVE) ; 1xn array
        ; FOR THE ELECTRONS: use the electron distribution function
        -1 :  PLS_Unc = PLS_UNCERTAINTIES(Structure, P_Locations, PLS_Data_Fit, PLS_ERROR_ORDF, Model_Name, index, Delta, Current_Threshold, /RELATIVE) ; 1xn array
      ENDCASE
      ;
      ; Decide whether to write uncertainties in file or create new
      IF KEYWORD_SET(WRITE_UNCERTAINTIES_IN_FILE) THEN BEGIN ; write uncertainties in the same file        
;        find the length of the PLS_Unc -> divide by 2 -> write variable in species
        Structure_Unc = Structure
        FOR i_species = 1, SIZE(PLS_Unc, /N_ELEMENTS)/2 DO BEGIN  ; write variables for each species
          Unc_StringCompare = STRCMP('Species_' + STRTRIM(i_species,2) + '_N_UNCERT', TAG_NAMES(Structure_Unc), /FOLD_CASE) 
          uncert_ind = WHERE(Unc_StringCompare EQ 1, nind, /NULL) 
          IF nind NE 1 THEN MESSAGE,'Error: column labeled "SPECIES_#_N_UNCERT" not found! Check .csv file'
          Structure_Unc.(uncert_ind) = PLS_Unc(2*i_species - 2) 
          Unc_StringCompare = STRCMP('Species_' + STRTRIM(i_species,2) + '_T_UNCERT', TAG_NAMES(Structure_Unc), /FOLD_CASE)
          uncert_ind = WHERE(Unc_StringCompare EQ 1, nind, /NULL)
          IF nind NE 1 THEN MESSAGE,'Error: column labeled "SPECIES_#_N_UNCERT" not found! Check .csv file'
          Structure_Unc.(uncert_ind) = PLS_Unc(2*i_species - 1)
        ENDFOR
      ENDIF ELSE BEGIN  ; create a new file and write the uncertainties in the place of n and T (that's how it was done for Jupiter)
        Structure_Unc = Structure ; Set a similar structure for the uncertainties
        FOR j = 0, N_ELEMENTS(ind) - 1 DO Structure_Unc.(ind[j]) = PLS_Unc[j] ; Save uncertainties to structure -- ind = columns that uncertainties routine returned values for
      ENDELSE
      ;STOP

    ENDIF ELSE IF KEYWORD_SET(COMPUTE_CHI_SPACE) THEN BEGIN
      str_date = STRCOMPRESS(STRING(ULONG(Structure.Year), ':', ULONG(Structure.DOY), ':', ULONG(Structure.Hour), ':', $
                                    ULONG(Structure.Minute), '::', ULONG(Structure.Second), $
                                    FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02)'), /REMOVE_ALL)
      PRINT, STRCOMPRESS(STRING('Calculating Chi-Square space for ', str_date))
      ; Use underscores for text filename
      str_date = STRCOMPRESS(STRING(ULONG(Structure.Year), '_', ULONG(Structure.DOY), '_', ULONG(Structure.Hour), '_', $
                                    ULONG(Structure.Minute), '_', ULONG(Structure.Second), '_CHI_ARRAYS.txt', $
                                    FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02,A)'), /REMOVE_ALL)
      ; Run through the Chi-Square calculations
      N_Chi_Steps = N_Chi_Pts
      N_Sigma     = N_Chi_Sig
      Model_Name = 'Generate_PLS_Current_Neptune'
      index = WHERE(FINITE(PLS_Data_Fit) EQ 0, nindex) ; Indexes of calculation to avoid
      IF nindex EQ 0 THEN index = -9d ; Set as a filler number so that it will account for all data if none is bad
      
      Structure_Hold = Structure ; Hold a copy of the structure before computing the chi arrays, just in case IDL overwrites it
      FORWARD_FUNCTION CHI_ARRAYS
      Chis = CHI_ARRAYS(Structure, Structure_Unc, P_Locations, PLS_Data_Fit, PLS_Error, Model_Name, index, N_Chi_Steps, N_Sigma, str_date)    
      Structure = Structure_Hold ; Return normal structure after computing chi arrays
    ENDIF
  ENDIF
  
  
  ; Get the fitting residuals in a file, if required
  ; --------------------------------------
  IF KEYWORD_SET(CALC_FITTING_RESIDUALS) THEN BEGIN
      FORWARD_FUNCTION FITTING_RESIDUALS
      PLS_Species_Currents = Generate_PLS_Current_Neptune(Structure)
      PLS_Simulated_Currents = TOTAL(PLS_Species_Currents, 2)
      IF Structure.L_or_M_Mode NE -1 THEN BEGIN
        residuals = FITTING_RESIDUALS(Structure, PLS_Data_Fit, PLS_Simulated_Currents)
      ENDIF ELSE residuals = FITTING_RESIDUALS(Structure, PLS_Data_Fit*alphaj, PLS_Simulated_Currents*alphaj) ; multiply by the factor alphaj to go from REDF to current
  ENDIF
  
  
  ; --------------------------------------
  ; SIMULATE AND PLOT THE CURRENT, BASED ON THE CALCULATED PARAMETERS ABOVE
  ; --------------------------------------
  ; [GX - IN THIS PART: first take into consideration with an if the ions and the elseif is for the electrons]
  ;
  ; IONS:
  IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN
    IF Structure.L_or_M_Mode NE -1 THEN BEGIN ; Plotting procedure for ions in M and L-modes
      ; Generate final result of model
      PLS_Currents = Generate_PLS_Current_Neptune(Structure)
      
      ; Sum up currents to get a total current curve
      IF Structure.L_or_M_Mode EQ 1 THEN n_channels = 128 ELSE n_channels = 16 ; Number of channels per cup
      Total_Current = MAKE_ARRAY(n_channels*4, /DOUBLE)
   
      ; Sum up currents from different ion species
      Modified_Array = TOTAL(PLS_Currents, 2) ; Reduce it to a 4 x Number of Channel Array, summing together the number of species
      
      ; Set up plot arguments
      xtitle = 'Channel Number'
      COMMON GENERATE
      IF (LDIST EQ "FALSE")  THEN ytitle = "Current (fA)"
      IF (LDIST EQ "TRUE") THEN ytitle = "Reduced Distribution Function (FemtoMhos)"

      ; Color titles for the plots
      colors = ['deep pink']
      ; Use Generate PLS Currents to get the plasma properties
      Plasma_Properties = Generate_PLS_Current_Neptune(Structure, /RETURN_PLASMA_PROPERTIES) ; this returns ONLY the plasma properties
      Masses       = Plasma_Properties[0,*]
      Charges      = Plasma_Properties[1,*]
      Densities    = Plasma_Properties[2,*]
      Temperatures = Plasma_Properties[3,*]
      
      ; COLOURS (HARDCODED). If we use two species then we have two colours
      ; (one for each species), otherwise we only have one (given that we are using a single species of H).
      IF N_ELEMENTS(Masses) EQ 2 THEN colors = ['royal blue','purple'] ELSE colors = MAKE_ARRAY(1, N_ELEMENTS(Masses), /STRING, VALUE = 'royal blue')
      
    
      ; Positions on the plot for each faraday cup
      POSITIONS = [[.10,.60,.45,.90],[.60,.60,.95,.90], [.10,.20,.45,.50], [.60,.20,.95,.50]]

      ; Determine a range of currents (logarithmic scaling)
      ;YRANGE = PLOT_RANGES_LOG(PLS_Data_Transp) ; - Automatically
      YRANGE = [1d0,1d5] ; - Manual override
    
      ; Turn off IDL errors for logarithmic plotting of values that may be zero
      !Except = 0
    
      ; Create dual error array so that array elements don't go below 0 - IDL will not plot negative values on a logarithmic scale
      PLS_Error_Array = MAKE_ARRAY(2, N_ELEMENTS(PLS_Error), /DOUBLE)
      
      ind_lt0 = WHERE(PLS_DATA - PLS_Error LE 0, nind_lt0)
      Dum_Arr = PLS_Error
      Dum_Arr[ind_lt0] = PLS_DATA[ind_lt0] - 1d-3 ; Set error values just above zero so that they appear in plot
      PLS_Error_Array(0,*) = Dum_Arr[*]
      PLS_Error_Array(1,*) = PLS_Error[*]
      ; Set errors to be in the proper format of 2 x 4 x number of channels
      Errors = MAKE_ARRAY(2,4,n_channels)
      FOR j = 0, 3 DO BEGIN
        Errors[*,j,*] = PLS_Error_Array[*, j*n_channels:(j+1)*n_channels-1]
      ENDFOR
      Cups = ['A','B','C','D']
    
    
      IF cupss eq 0 then begin
      ; Plot the Data
      ; Begin the loop through each cup
      FOR iii = 0, 3 DO BEGIN
        ; Values greater than 100 are for saving plots without them popping up
        IF Structure.Save_Plot GE 100 THEN BEGIN
          PLOT_BUFFER=1
        ENDIF
        dummy = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720], current=iii,  $
                     xrange = [1,n_channels], POSITION = positions[*,iii], BUFFER=PLOT_BUFFER, yrange = YRANGE, Color = 'Deep Pink', xminor = 1) ; as the x label is every 2 ticks, xminor=1 shows the label in-between 
        dummy_dat = PLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], /OVERPLOT, POSITION = positions[*,iii], thick=2)
        eplot = ERRORPLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], Errors[*,iii,*], $                                                        ;errorbar overplot
            ERRORBAR_COLOR = 'Orange', $
            ERRORBAR_CAPSIZE = 0.0, $
            errorbar_thick = 1, $
            /OVERPLOT)
        ; Calculate the number of species and then overplot them individually
        n_species = N_ELEMENTS(PLS_Currents)/(n_channels*4)
        FOR jjj = 0, n_species-1 DO BEGIN
          pplot= PLOT(FINDGEN(n_channels)+1, PLS_Currents(iii,jjj,*), title = Cups[iii]+' Cup', xtitle = xtitle, color=colors[jjj],$              ; LOGAN 5/5 - Need to make this adhere to the colors rule
                       ytitle = ytitle, yrange = yrange, /ylog, /overplot)
        ENDFOR
        dummy_tot = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720], current=iii,  $
                         xrange = [1,n_channels], POSITION = positions[*,iii], BUFFER=PLOT_BUFFER, yrange = YRANGE, Color = 'Deep Pink', /OVERPLOT, thick=1.5)
      ENDFOR
      ENDIF
    
      IF cupss NE 0 THEN BEGIN
        positionss = [.10,.20,.95,.90]
        IF cupss EQ 1 THEN Cup = 'A'
        IF cupss EQ 2 THEN Cup = 'B'
        IF cupss EQ 3 THEN Cup = 'C'
        IF cupss EQ 4 THEN Cup = 'D'
                        
      

        ; Plot the Data
        ; Begin the loop through each cup
        FOR iii = cupss-1, cupss-1 DO BEGIN
          ; Values greater than 100 are for saving plots without them popping up
          IF Structure.Save_Plot GE 100 THEN BEGIN
            PLOT_BUFFER=1
          ENDIF
          dummy = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720],  $
            xrange = [1,n_channels], BUFFER=PLOT_BUFFER,POSITION = positionss, yrange = YRANGE, Color = 'Deep Pink')
          dummy_dat = PLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], thick=2,/overplot)
          eplot = ERRORPLOT(FINDGEN(n_channels)+1, PLS_Data_Transp[iii,*], Errors[*,iii,*], $                                                        ;errorbar overplot
            ERRORBAR_COLOR = 'Orange', $
            ERRORBAR_CAPSIZE = 0.0, $
            errorbar_thick = 1, $
            /OVERPLOT)
          ; Calculate the number of species and then overplot them individually
          n_species = N_ELEMENTS(PLS_Currents)/(n_channels*4)
          FOR jjj = 0, n_species-1 DO BEGIN
            pplot= PLOT(FINDGEN(n_channels)+1, PLS_Currents(iii,jjj,*), title = Cups[iii]+' Cup', xtitle = xtitle, color=colors[jjj+1],$              ; LOGAN 5/5 - Need to make this adhere to the colors rule
              ytitle = ytitle, yrange = yrange, /ylog, /overplot)
          ENDFOR
          dummy_tot = PLOT(FINDGEN(n_channels)+1, Modified_Array[iii,*], dimensions=[1000,720],  $
            xrange = [1,n_channels], BUFFER=PLOT_BUFFER,POSITION = positionss, yrange = YRANGE, Color = 'Deep Pink', /OVERPLOT, thick=1.5)
        ENDFOR
      ENDIF

      ; Print the relevant information on the graphic
    
      ; Determine Response, Radial Distance, and Time of plot
      IF Structure.Response EQ 0 THEN Response_Name = 'CUPINT/DCPINT' ELSE $
      IF Structure.Response EQ 1 THEN Response_Name = 'LABCUR/LDCUR'
    
      COMMON VGR2Data_NEPTUNE
      SSEDR_Data = VGR2SSEDR_NEPTUNE
    
      ; Time in SSEDR file
      SSEDR_Time = SSEDR_Data[1,*] + (SSEDR_Data[2,*] + SSEDR_Data[3,*]/60d + SSEDR_Data[4,*]/3600d)/24d
    
      ; Time of Input
      Time_Requested = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
      ; Time of Data
;      IF Structure.L_or_M_Mode EQ 0 THEN Data_Time = Dum_Struct.Ltimes_Array ELSE $
;      IF Structure.L_or_M_Mode EQ 1 THEN Data_Time = Dum_Struct.Mtimes_Array
    
      Data_Time = Dum_Struct.times_Array
    
      ; Create DOY, Hour, Minute, and Second from Time_Data
      ; LOGAN 1/7/15 - Check to see if this works
      Data_DOY  = FIX(Data_Time)
      Data_Hour = FIX((Data_Time - Data_DOY)*24d)
      Data_Min  = FIX(((Data_Time - Data_DOY)*24d - Data_Hour)*60d)
      Data_Sec  = FIX((((Data_Time - Data_DOY)*24d - Data_Hour)*60d - Data_Min)*60d)
    
      ; Find SSEDR index of data to find S/C location
      SSEDR_ind = CLOSEST(TRANSPOSE(SSEDR_Time[0,*]), Data_Time[0])
      ; Get X, Y, and Z locations of spacecraft and compute radial distance (R)
      xSC_V2_NEPTUNE_PLS = SSEDR_Data[6,SSEDR_ind] ; V2_NEPTUNE_PLS cartesian coordinates
      ySC_V2_NEPTUNE_PLS = SSEDR_Data[7,SSEDR_ind]
      zSC_V2_NEPTUNE_PLS = SSEDR_Data[8,SSEDR_ind]
      r_RU = SQRT(xSC_V2_NEPTUNE_PLS*xSC_V2_NEPTUNE_PLS+ySC_V2_NEPTUNE_PLS*ySC_V2_NEPTUNE_PLS+zSC_V2_NEPTUNE_PLS*zSC_V2_NEPTUNE_PLS)/r_Neptune
      
      
      IF Structure.SHORT_OR_LONG EQ 0 THEN integr_time = 'SHORT' ELSE $
        IF Structure.SHORT_OR_LONG EQ 1 THEN integr_time = 'LONG'
      
      ; Print information to graphic
;      t = TEXT(0.5,0.98, 'Mode:' + Response_Name + '           Year:' + STRING(ULONG(Structure.Year)) + $
;                         '           Date (DOY, Hour, Min, Second):' + $
;                         STRING(ULONG(Data_DOY)) + STRING(ULONG(Data_Hour)) + STRING(ULONG(Data_Min)) + STRING(ULONG(Data_Sec)) + $
;                         '           r:    ' + STRING(r_RU, format='(F0.3)'),  font_size = 10, color = 'black', $
;                         ALIGNMENT = .5, vertical_ALIGNMENT = .5, /normal)

      t = TEXT(0.03,0.96, 'Year-DOY Hour:Min:Second: ' + STRTRIM(ULONG(Structure.Year),2) + '-'  + STRING(ROUND(Data_DOY), format='(I03)') + '  ' $
                 + STRING(ROUND(Data_Hour), format='(I02)') + ':' + STRING(ROUND(Data_Min), format='(I02)') + ':' + STRING(ROUND(Data_Sec), format='(I02)') + $
        '          r: ' + STRING(r_RU, format='(F0.3)')+' $R_N$' + '          Method: ' + Response_Name + '          Integr. time: ' + integr_time,  $
        font_size = 11, color = 'black', /normal)
;      tt = TEXT(0.05,0.94, ,  font_size = 10, color = 'black', /normal)
      ;stop
      ; Print remaining argument tags
      u = TEXT(.02, .09, 'A (amu), Z (q):', /normal)
      u = TEXT(.02, .06, '$n (cm^{-3}):$',  /normal)
      u = TEXT(.02, .03, 'T (eV):',  /normal)

      ; Flowspeed
      q=text(.02, .12, '$Plasma flowspeed [V_r, V_{\phi}, V_z] (km/s):$')  
      dum_names = [Structure.V1, Structure.V2, Structure.V3] 
;      sep=.06
;      FOR k = 0, 2 DO q=text((1+k)*sep*2d + 0.15,0.12, string(dum_names[k], format='(F0.2)')) 
      q=text(0.3,0.12, '[' + string(dum_names[0], format='(F0.2)') + ',     ' + string(dum_names[1], format='(F0.2)') + ',     ' + string(dum_names[2], format='(F0.2)') + ']')
 
      ; Use Generate PLS Currents to get the plasma properties  
      Plasma_Properties = Generate_PLS_Current_Neptune(Structure, /RETURN_PLASMA_PROPERTIES)
      Masses       = Plasma_Properties[0,*]
      Charges      = Plasma_Properties[1,*]
      Densities    = Plasma_Properties[2,*]
      Temperatures = Plasma_Properties[3,*]
      ; If uncertenties were calculated get the uncertainties here
      IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
        Densities_Unc = TRANSPOSE(PLS_UNC[INDGEN(N_ELEMENTS(Masses)) * 2])
        Temperatures_Unc = TRANSPOSE(PLS_UNC[INDGEN(N_ELEMENTS(Masses)) * 2 + 1])
      ENDIF
      
      ; Determine seperation based on the number of species
      seps = [0.29,0.26,0.23,0.20,0.17,0.15,0.12,0.10,0.08]
      seperation = seps[N_ELEMENTS(Masses)-1]
      
      ; Define numbers format and plot them
      ; Colours: for Jupiter a multi-colour option was used, where it was tying each colour with its species, and it made sense since there 
      ; could be up to 8 different species. For Neptune we are colour-coding each species, to make the difference more prominent. [GX]

      ; Create number format arrays. Format different for H and N. Needs to change if using more than 2 components
      density_format= STRARR(N_ELEMENTS(masses)) & temperature_format = density_format
      density_format[*,0] = '(F5.3)' & temperature_format[*,0] = '(F5.1)' ; Logan 2025_07_31 - Changed from (I 5) on temperature format
      IF N_ELEMENTS(masses) GE 2 THEN BEGIN ; this is specifically for the hot component
         density_format[1,0] = '(F6.4)' & temperature_format[1,0] = '(I 4)'
      ENDIF
      
      
      IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
        FOR k = 0, N_ELEMENTS(Masses)-1 DO BEGIN
          v = TEXT((1+k)*seperation+.047,0.09, STRCOMPRESS(STRING(UINT(Masses(k)))+','+STRING(UINT(Charges(k)))), color=colors[k], /norm)
          v = TEXT((1+k)*seperation+.05,0.06, STRCOMPRESS(STRING(Densities[k],format = density_format[k]) + ' $\pm$ ' + STRING(Densities_Unc[k],format = density_format[k])), color=colors[k], /norm)
          v = TEXT((1+k)*seperation+.05,0.03, STRCOMPRESS(STRING(Temperatures[k], format = temperature_format[k])+ ' $\pm$ ' + STRING(Temperatures_Unc(k),format = temperature_format[k])), color=colors[k], /norm)
        ENDFOR
      ENDIF ELSE BEGIN
        FOR k = 0, N_ELEMENTS(Masses)-1 DO BEGIN
          v = TEXT((1+k)*seperation+.047,0.09, STRCOMPRESS(STRING(UINT(Masses(k)))+','+STRING(UINT(Charges(k)))), color=colors[k], /norm)
          v = TEXT((1+k)*seperation+.05,0.06, STRING(Densities[k],format = density_format[k]), color=colors[k], /norm)
          v = TEXT((1+k)*seperation+.05,0.03, STRING(Temperatures[k], format = temperature_format[k]), color=colors[k], /norm)
        ENDFOR
      ENDELSE
      
;      FOR k = 0, N_ELEMENTS(Masses)-1 DO BEGIN
;        v = TEXT((1+k)*seperation+.047,0.09, STRCOMPRESS(STRING(UINT(Masses(k)))+','+STRING(UINT(Charges(k)))), color=colors[k+1], /norm)
;        v = TEXT((1+k)*seperation+.05,0.06, STRING(Densities(k),format='(F0.4)'), /norm)
;        v = TEXT((1+k)*seperation+.05,0.03, STRING(Temperatures[k], format='(F0.2)'), /norm)
;      ENDFOR
      

      
      ; Save a plot if requested
      IF Structure.Save_Plot GE 1 THEN BEGIN
      
      
; KEEP THIS PART FOR NOW: IT SAVES FILES WITH THE CURRENT DATE
; Rob: this is convenient when you're still in the process of playign with the data, as it doesn't overwrites
; the existent files. Would a flag that gives the option to chose the format be useful?
; 
;        ; Create name for plot based on current UTC time
;        YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
;        DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
;        IF ULONG(DayOfRun) LT 10 THEN BEGIN
;          DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
;        ENDIF
;        MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
;        HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
;        MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
;        SecondOfRun = STRMID(SYSTIME(/UTC),17,2)
;        
;      
;        ; Name the file
;        IF ((Structure.Save_Plot EQ 1) OR (Structure.Save_Plot EQ 101)) THEN $
;        OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
;                                HourOfRun, MinuteOfRun, '_', SecondOfRun, '.jpg', $
;                                FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
;        ELSE IF ((Structure.Save_Plot EQ 2) OR (Structure.Save_Plot EQ 102)) THEN $
;        OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
;                                HourOfRun, MinuteOfRun, '_', SecondOfRun, '.png', $
;                                FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
;        ELSE IF ((Structure.Save_Plot EQ 3) OR (Structure.Save_Plot EQ 103)) THEN $
;        OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
;                                HourOfRun, MinuteOfRun, '_', SecondOfRun, '.eps', $
;                                FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
;        ELSE PRINT, 'Unable to save selected image file type.'


        ; Step 1 - Create name for file that uses "Structure" parameters to indicate measurement
        ; Create as filename the time the meaurement was taken
        ; Mode for the filename
        IF Structure.L_or_M_Mode EQ 1 THEN Mode_for_filename = 'M' ELSE BEGIN
          CASE Structure.Short_or_Long OF
            0: Mode_for_filename = 'L_SHORT'
            1: Mode_for_filename = 'L_LONG'
          ENDCASE
        ENDELSE
        MeasurementFilename = STRING(ULONG(Structure.Year), FORMAT='(I0)') + '-'  + STRING(ULONG(Structure.Doy), FORMAT='(I03)') + $
          '-'  + STRING(ULONG(Structure.Hour), FORMAT='(I02)') + '-'  + STRING(ULONG(Structure.Minute), FORMAT='(I02)') + $
          '-' + STRING(ULONG(Structure.Second), FORMAT='(I02)') + '_' + Mode_for_filename


        ; Step 2 - Create Timestring of current run to append to the file, so that we create unique filenames
        ; Create name for plot based on current UTC time
        YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
        DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
        IF ULONG(DayOfRun) LT 10 THEN BEGIN
          DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
        ENDIF
        MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
        HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
        MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
        SecondOfRun = STRMID(SYSTIME(/UTC),17,2)

        UniqueTimeFilename = STRING('_T-', MonthOfRun, DayOfRun, '-', YearOfRun, '-', $
          HourOfRun, MinuteOfRun, '-', SecondOfRun, FORMAT = '(A,A,A,A,A,A,A,A,A,A)')

        ; Name the file
        IF ((Structure.Save_Plot EQ 1) OR (Structure.Save_Plot EQ 101)) THEN $
          OutputFilename = STRING(MeasurementFilename, UniqueTimeFilename, '.jpg', FORMAT = '(A,A,A)') $
        ELSE IF ((Structure.Save_Plot EQ 2) OR (Structure.Save_Plot EQ 102)) THEN $
          OutputFilename = STRING(MeasurementFilename, UniqueTimeFilename, '.png', FORMAT = '(A,A,A)') $
        ELSE IF ((Structure.Save_Plot EQ 3) OR (Structure.Save_Plot EQ 103)) THEN $
          OutputFilename = STRING(MeasurementFilename, UniqueTimeFilename, '.eps', FORMAT = '(A,A,A)') $
        ELSE PRINT, 'Unable to save selected image file type.'
        
        
        ; Save the file if of the 3 types - Can add different image types if desired at a later date
        IF ISA(OutputFilename) NE 0 THEN BEGIN
          ; Change directories for output of plots
          ; 
          ; LOGAN - Commented out 2025_02_27 - I think this is doing what I did above
          CD, CURRENT = RESTORE_DIR
          CD, POD
;          ; Check if filename exists. If it does, create a number in the end. If a number exists, name the next
;          IF FILE_TEST(OutputFilename) EQ 1 THEN BEGIN
;          file_list = FILE_SEARCH( POD, '*' )
;          count = 0
;          FOR i = 0, N_ELEMENTS(file_list) - 1 DO BEGIN
;            IF STRPOS(file_list[i], OutputFilename) NE -1 THEN count = count + 1
;          ENDFOR
;          ; as the extension changes, get all the characters before the extension
;          parts = STRSPLIT(OutputFilename, '.', /EXTRACT )
;          OutputFilename = STRING(parts[0]) + '(' + STRING(ULONG(count), FORMAT = '(I0)') + ')' + '.' + STRING(parts[1])
;          ENDIF
          dummy.SAVE, OutputFilename, Border = 10, RESOLUTION = 150
          ; Change back directories
          CD, RESTORE_DIR
        ENDIF
      ENDIF

      
    
    ; ELECTRONS
    ; ------------
    ENDIF ELSE IF Structure.L_or_M_Mode EQ -1 THEN BEGIN ; Plot electron data, channels 1-28 - Borrows a lot of above code for ease
      ;STOP

      ;; DATA PREPARATION
      ;---------------------
      ; GATHER DATA FOR CURRENT AND ORDF
      PLS_DATA_ORDF = PLS_DATA_Transp_ORDF[3,*] ; for A (only the 4th column of PLS_DATA_Transp_ORDF has data - the rest are zeros)
      
      ; Fill PLS_DATA_Transp for plotting:
      ; 1st & 2nd column: Observed electron reduced distribution in femtomhos
      ; 3rd column:       electron flux (units?)
      ; 4th column:       electron current in fA
      PLS_DATA_Transp[0,*] = PLS_DATA_ORDF * 1d15
      PLS_DATA_Transp[1,*] = PLS_DATA_ORDF * 1d15
      PLS_DATA_Transp[3,*] = PLS_DATA_Transp[3,*] * alphaj * 1d15
      PLS_DATA_FLUX = PLS_DATA_Transp[3,*]/FluxConversion[*] ; already current in fA
      PLS_DATA_Transp[2,*] = PLS_DATA_Flux
      
      
      ; GENERATE FINAL RESULT OF MODEL
      PLS_Currents = Generate_PLS_Current_Neptune(Structure) ; in RDF for A
      ; Sum up currents to get a total current curve
      n_channels = 28 ; Number of channels in e-
      Total_Current = MAKE_ARRAY(n_channels*4, /DOUBLE)

      ; Fill individual current array
      n_species = N_ELEMENTS(PLS_Currents)/(n_channels*4)
      FOR jloop = 0, n_species-1 DO BEGIN
        ; 1st & 2nd column: Observed electron reduced distribution in femtomhos
        ; 3rd column:       electron flux (units?)
        ; 4th column:       electron current in fA
        PLS_Currents(0,jloop,*) = PLS_Currents(3,jloop,*) * 1d15
        PLS_Currents(1,jloop,*) = PLS_Currents(3,jloop,*) * 1d15
        PLS_Currents(3,jloop,*) = PLS_Currents(3,jloop,*) * alphaj * 1d15 ; in fA
        PLS_Currents(2,jloop,*) = PLS_Currents(3,jloop,*)/FluxConversion[*]
      ENDFOR
      ; Sum up currents from different electron species
      Modified_Array = TOTAL(PLS_Currents, 2) ; Reduce it to a 4 x Number of Channel Array, summing together the number of species


      ; GET PLASMA PROPERTIES USING THE GENERATE CURRENT FUNCTION
      Plasma_Properties = Generate_PLS_Current_Neptune(Structure, /RETURN_PLASMA_PROPERTIES)
      Masses       = Plasma_Properties[0,*]
      Charges      = Plasma_Properties[1,*]
      Densities    = Plasma_Properties[2,*]
      Temperatures = Plasma_Properties[3,*]
      IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
        Densities_Unc = TRANSPOSE(PLS_UNC[INDGEN(N_ELEMENTS(Masses)) * 2])
        Temperatures_Unc = TRANSPOSE(PLS_UNC[INDGEN(N_ELEMENTS(Masses)) * 2 + 1])
      ENDIF
      
      
      ; ERRORS
      ; Only plotting the error for the measured quantities (ORDF, flux, and current), NOT the simulated
      ; A dual error array was created so that array elements don't go below 0 - IDL will not plot negative 
      ; values on a logarithmic scale. If error goes below 0 it changes to a very small number.      
      ;
      ; Error array: 4 first columns for observed electron reduced distribution error (in femtomhos)
      ; 5-6 columns for electron flux error, 7-8 columns for electron current error in fA
      PLS_Error_Array = MAKE_ARRAY(8, N_ELEMENTS(PLS_Error)/4, /DOUBLE)
      ;
      ; ORDF error
      ind_lt0 = WHERE(PLS_DATA_Transp[0,*] - (PLS_ERROR_ORDF[-28:*]*1d15) LE 0, nind_lt0 )  ; factor 1d15 A to fA for the error
      Error_array = PLS_ERROR_ORDF[-28:*]*1d15
      Dum_Arr = Error_array
      Dum_Arr[ind_lt0] = PLS_DATA_Transp[0,ind_lt0] - 1d-20 ; Set error values just above zero so that they appear in plot
      PLS_Error_Array[0,*] = Dum_Arr[*]  &  PLS_Error_Array[2,*] = PLS_Error_Array[0,*]
      PLS_Error_Array[1,*] = Error_array[*]  &  PLS_Error_Array[3,*] = PLS_Error_Array[1,*]
      ; Current error
      ind_lt0 = WHERE(PLS_DATA_Transp[3,*] - PLS_Error[-28:*] LE 0, nind_lt0) ; already in fA
      Error_array = PLS_Error[-28:*]
      Dum_Arr = Error_array
      Dum_Arr[ind_lt0] = PLS_DATA_Transp[3,ind_lt0] - 1d-3
      PLS_Error_Array[6,*] = Dum_Arr[*]
      PLS_Error_Array[7,*] = Error_array
      ; Flux error
      ind_lt0 = WHERE(PLS_DATA_Transp[2,*] - (PLS_Error[-28:*]/FluxConversion) LE 0, nind_lt0) ; already in fA
      Error_array = PLS_Error[-28:*]/FluxConversion
      Dum_Arr = Error_array
      Dum_Arr[ind_lt0] = PLS_DATA_Transp[2,ind_lt0] - 1d-5
      PLS_Error_Array[4,*] = Dum_Arr[*]
      PLS_Error_Array[5,*] = Error_array
      ;
      ; Set errors to be in the proper format of 2 x 4 x number of channels
      Errors = MAKE_ARRAY(2,4,n_channels)
      FOR j = 0, 3 DO BEGIN
        Errors[*,j,*] = PLS_Error_Array[ j*2 : (j*2)+1 ,*]
      ENDFOR
      
      
      
      ; PLOTTING
      ; ------------
      ; Set XRANGES
      xvals = MAKE_ARRAY(4,28)
      xvals(0,*) = SQRT(ECHRGE*2d*Eenergy/me)/1000d
      xvals(1,*) = Eenergy
      xvals(2,*) = Eenergy
      xvals(3,*) = findgen(28) + 1
      
      xranges = MAKE_ARRAY(4,2)
      xranges(0,*) = [2000,45000]
      xranges(1,*) = [10,6000]
      xranges(2,*) = [10,6000]
      xranges(3,*) = [1,28]
      xlogs = [0,1,1,0]
      
      titles  = ["Electron Reduced Distribution Function", "Electron Reduced Distribution Function", "Electron Flux", "Electron Current"]
      xtitles = ["Velocity (km/s)", "Energy (eV)", "Energy (eV)", "Channel Number"]
      ytitles = ["Reduced Distribution Function (FemtoMhos)", "Reduced Distribution Function (FemtoMhos)", "Flux", "Current (fA)"]
      
      ; Set up plot arguments
      xtitle = 'Channel Number'
      COMMON GENERATE
      IF (LDIST EQ "FALSE")  THEN ytitle = "Current (fA)"
      IF (LDIST EQ "TRUE") THEN ytitle = "Reduced Distribution Function (FemtoMhos)"

      ; Color titles for the plots
      colors = ['deep pink','royal blue','purple','red'] ; Logan 2025_08_13 - Added 'red'
      
            
      ; Positions for each faraday cup
      POSITIONS = [[.10,.60,.45,.90],[.60,.60,.95,.90], [.10,.20,.45,.50], [.60,.20,.95,.50]] ; FIX THIS

      ; Automatically Determine a range of currents (logarithmic scaling) that matches well with the data
      ;YRANGE = PLOT_RANGES_LOG(PLS_Data_Transp)

      ; Manually set the y range
      YRANGE = [[1d-10,1d-10,1d0,1d0],[1d1,1d1,1d10,1d5]]
      
      ; Turn off IDL errors for logarithmic plotting of values that may be zero
      !Except = 0


      Cups = ['A','B','C','D']
      
      IF cupss EQ 0 THEN BEGIN
        ; Plot the Data
        ; Begin the loop through each cup
        FOR iii = 0, 3 DO BEGIN ; Just plot the D cup in this instance
          ; Values greater than 100 are for saving plots without them popping up
          IF Structure.Save_Plot GE 100 THEN BEGIN
            PLOT_BUFFER=1
          ENDIF
          ; Model results
          dummy = PLOT(xvals[iii,*], Modified_Array[iii,*], dimensions=[1000,720], current=iii,  $
            xrange = xranges[iii,*], POSITION = positions[*,iii], BUFFER=PLOT_BUFFER, yrange = YRANGE[iii,*], Color = 'Deep Pink', /ylog, xlog = xlogs[iii])
          ; Raw data
          dummy_dat = PLOT(xvals[iii,*], PLS_Data_Transp[iii,*], /OVERPLOT, POSITION = positions[*,iii], thick=2, /ylog, xlog = xlogs[iii])
          eplot = ERRORPLOT(xvals[iii,*], PLS_Data_Transp[iii,*], Errors[*,iii,*], $                                                        ;errorbar overplot
            ERRORBAR_COLOR = 'Orange', $
            ERRORBAR_CAPSIZE = 0.0, $
            errorbar_thick = 1, $
            /OVERPLOT, /ylog, xlog = xlogs[iii])
          ; Overplot number of specied individually
          n_species = N_ELEMENTS(PLS_Currents)/(n_channels*4)
          FOR jjj = 0, n_species-1 DO BEGIN
            pplot= PLOT(xvals[iii,*], PLS_Currents(iii,jjj,*), title = titles[iii], xtitle = xtitles[iii], color=colors[jjj+1],$              ; LOGAN 5/5 - Need to make this adhere to the colors rule
              ytitle = ytitles[iii], yrange = yrange[iii,*], /ylog, xlog = xlogs[iii], /overplot)
;            eplot = ERRORPLOT(xvals[iii,*], PLS_Data_Transp[iii,*], Errors[*,iii,*], $                                                        ;errorbar overplot
;                ERRORBAR_COLOR = 'Orange', $
;                ERRORBAR_CAPSIZE = 0.0, $
;                errorbar_thick = 1, $          xminor = 1) ; as the x label is every 2 ticks, xminor=1 shows the label in-between 
;                /OVERPLOT)
                
          ENDFOR
          dummy_tot = PLOT(xvals[iii,*], Modified_Array[iii,*], dimensions=[1000,720], current=iii,  $
            xrange = xranges[iii,*], POSITION = positions[*,iii], BUFFER=PLOT_BUFFER, yrange = YRANGE[iii,*], Color = 'Deep Pink', /OVERPLOT, thick=1.5)
        ENDFOR
      ENDIF


      IF cupss NE 0 THEN BEGIN
        positionss = [.10,.20,.95,.90]
        IF cupss EQ 1 THEN Cup = 'A'
        IF cupss EQ 2 THEN Cup = 'B'
        IF cupss EQ 3 THEN Cup = 'C'
        IF cupss EQ 4 THEN Cup = 'D'      
      
        ; Plot the Data
        ; Begin the loop through each cup
        FOR iii = cupss-1, cupss-1 DO BEGIN
          ; Values greater than 100 are for saving plots without them popping up
          IF Structure.Save_Plot GE 100 THEN BEGIN
            PLOT_BUFFER=1
          ENDIF
          dummy = PLOT(xvals[iii,*], Modified_Array[iii,*], dimensions=[1000,720],  $
            xrange = xranges[iii,*], BUFFER=PLOT_BUFFER,POSITION = positionss, yrange = YRANGE[iii,*], Color = 'Deep Pink')
          dummy_dat = PLOT(xvals[iii,*], PLS_Data_Transp[iii,*], thick=2,/overplot)
      ;    eplot = ERRORPLOT(xvals[iii,*], PLS_Data_Transp[iii,*], Errors[*,iii,*], $                                                        ;errorbar overplot
      ;      ERRORBAR_COLOR = 'Orange', $
      ;      ERRORBAR_CAPSIZE = 0.0, $
      ;      errorbar_thick = 1, $
      ;      /OVERPLOT)
          ; Calculate the number of species and then overplot them individually
          n_species = N_ELEMENTS(PLS_Currents)/(n_channels*4)
          FOR jjj = 0, n_species-1 DO BEGIN
            pplot= PLOT(xvals[iii,*], PLS_Currents(iii,jjj,*), title = titles[iii], xtitle = xtitles[iii], color=colors[jjj+1],$              ; LOGAN 5/5 - Need to make this adhere to the colors rule
              ytitle = ytitles[iii,*], yrange = yrange[iii,*], /ylog, xlog = xlogs[iii], /overplot)
          ENDFOR
          dummy_tot = PLOT(xvals[iii,*], Modified_Array[iii,*], dimensions=[1000,720],  $
            xrange = xranges[iii,*], BUFFER=PLOT_BUFFER,POSITION = positionss, yrange = YRANGE[iii,*], Color = 'Deep Pink', /OVERPLOT, thick=1.5)
        ENDFOR
      ENDIF

      ; Print the relevant information on the graphic
      
      ; Determine Response, Radial Distance, and Time of plot - Not programmed well for multiple electron modes
      CASE Structure.Response OF
        0 : Response_Name = ' Sittler Response'
        1 : Response_Name = ' Kappa Fit'
      ENDCASE
      
      ;------------------ GX CHANGE FROM HERE
      IF Structure.Spacecraft EQ 2 THEN BEGIN
        COMMON VGR2Data_NEPTUNE
        SSEDR_Data = VGR2SSEDR_NEPTUNE
      ENDIF

      ; Time in SSEDR file
      SSEDR_Time = SSEDR_Data[1,*] + (SSEDR_Data[2,*] + SSEDR_Data[3,*]/60d + SSEDR_Data[4,*]/3600d)/24d
      
      ; Time of Input
      Time_Requested = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
      ; Time of Data
      ;      IF Structure.L_or_M_Mode EQ 0 THEN Data_Time = Dum_Struct.Ltimes_Array ELSE $
      ;      IF Structure.L_or_M_Mode EQ 1 THEN Data_Time = Dum_Struct.Mtimes_Array
      
      Data_Time = Dum_Struct.times_Array

      ; Create DOY, Hour, Minute, and Second from Time_Data
      ; LOGAN 1/7/15 - Check to see if this works
      Data_DOY  = FIX(Data_Time)
      Data_Hour = FIX((Data_Time - Data_DOY)*24d)
      Data_Min  = FIX(((Data_Time - Data_DOY)*24d - Data_Hour)*60d)
      Data_Sec  = FIX((((Data_Time - Data_DOY)*24d - Data_Hour)*60d - Data_Min)*60d)
      
      ; Find SSEDR index of data to find S/C location
      SSEDR_ind = CLOSEST(TRANSPOSE(SSEDR_Time[0,*]), Data_Time[0])
      ; Get X, Y, and Z locations of spacecraft and compute radial distance (r)
      xSC_V2_NEPTUNE_PLS = SSEDR_Data[6,SSEDR_ind] ; V2_NEPTUNE_PLS cartesian coordinates
      ySC_V2_NEPTUNE_PLS = SSEDR_Data[7,SSEDR_ind]
      zSC_V2_NEPTUNE_PLS = SSEDR_Data[8,SSEDR_ind]
      r_RU = SQRT(xSC_V2_NEPTUNE_PLS*xSC_V2_NEPTUNE_PLS+ySC_V2_NEPTUNE_PLS*ySC_V2_NEPTUNE_PLS+zSC_V2_NEPTUNE_PLS*zSC_V2_NEPTUNE_PLS)/r_Neptune

;      ; Print information to graphic
;      t = TEXT(0.5,0.98, 'Mode:' + Response_Name + '           Year:' + STRING(ULONG(Structure.Year)) + $
;        '           Date (DOY, Hour, Min, Second):' + $
;        STRING(ULONG(Data_DOY)) + STRING(ULONG(Data_Hour)) + STRING(ULONG(Data_Min)) + STRING(ULONG(Data_Sec)) + $
;        '           r:    ' + STRING(r_RU, format='(F0.3)'),  font_size = 10, color = 'black', $
;        ALIGNMENT = .5, vertical_ALIGNMENT = .5, /normal)
      ; Print integration time:
      IF Structure.SHORT_OR_LONG EQ 0 THEN integr_time = 'SHORT' ELSE $
      IF Structure.SHORT_OR_LONG EQ 1 THEN integr_time = 'LONG'
      t = TEXT(0.03,0.96, 'Year-DOY Hour:Min:Second: ' + STRTRIM(ULONG(Structure.Year),2) + '-'  + STRING(ROUND(Data_DOY), format='(I03)') + '  ' $
                 + STRING(ROUND(Data_Hour), format='(I02)') + ':' + STRING(ROUND(Data_Min), format='(I02)') + ':' + STRING(ROUND(Data_Sec), format='(I02)') + $
        '          r: ' + STRING(r_RU, format='(F0.3)')+' $R_U$' + '          Method: ' + Response_Name + '          Integr. time: ' + integr_time,  $
        font_size = 11, color = 'black', /normal)
      

      ; Print remaining argument tags
      u = TEXT(.02, .12, 'Electron Mass, Z (q):', /normal)
      u = TEXT(.02, .09, '$n (cm^{-3}):$',  /normal)
      u = TEXT(.02, .06, 'T (eV):',  /normal)
      
      
      IF Structure.Response EQ 1 THEN u = TEXT(.02, .03, 'Kappa:',  /normal)
      
      IF Structure.Response EQ 1 THEN kappatext = Structure.Species_1_Kappa ; A filler line for now since we are only fitting one Kappa species and probably will only ever fit one

      ; Use Generate PLS Currents to get the plasma properties
      Plasma_Properties = Generate_PLS_Current_Neptune(Structure, /RETURN_PLASMA_PROPERTIES)
      Masses       = Plasma_Properties[0,*]
      Charges      = Plasma_Properties[1,*]
      Densities    = Plasma_Properties[2,*]
      Temperatures = Plasma_Properties[3,*]
      IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
        Densities_Unc = TRANSPOSE(PLS_UNC[INDGEN(N_ELEMENTS(Masses)) * 2])
        Temperatures_Unc = TRANSPOSE(PLS_UNC[INDGEN(N_ELEMENTS(Masses)) * 2 + 1])
      ENDIF
      
      ; Determine seperation based on the number of species
      seps = [0.29,0.26,0.23,0.20,0.17,0.15,0.12,0.10,0.08]
      seperation = seps[N_ELEMENTS(Masses)-1]
      
      ; Create number format arrays. Format different for hot and cold component. Needs to change if using more than 2 components
      density_format= STRARR(N_ELEMENTS(masses)) & temperature_format = density_format
      density_format[*,0] = '(F0.3)' & temperature_format[*,0] = '(F5.1)'
      IF N_ELEMENTS(masses) GE 2 THEN BEGIN ; this is specifically for the hot component
         density_format[1,0] = '(F0.4)' & temperature_format[1,0] = '(F6.0)'
      ENDIF
      
      ; PLOT
      FOR k = 0, N_ELEMENTS(Masses)-1 DO BEGIN
        v = TEXT( (1+k)*seperation+ .05,0.12, STRCOMPRESS( STRING(UINT(Masses(k))) + ', -' + STRTRIM(UINT(Charges(k)),2)), color=colors[k+1], /norm) ; the +1 at color is to get the right colour. k = 0 is dark pink, which is for the simulated total current
        IF KEYWORD_SET(UNCERTAINTIES) THEN BEGIN
          v = TEXT((1+k)*seperation+ .05,0.09, STRCOMPRESS(STRING(Densities[k],format = density_format[k]) + ' $\pm$ ' + STRING(Densities_Unc[k], format = density_format[k])), color=colors[k+1], /norm)
          v = TEXT((1+k)*seperation+ .05,0.06, STRCOMPRESS(STRING(Temperatures[k], format = temperature_format[k]) + ' $\pm$ ' + STRING(Temperatures_Unc[k], format = temperature_format[k])), color=colors[k+1], /norm)
          IF Structure.Response EQ 1 THEN v = TEXT((1+k)*seperation+ .05,0.03, STRING(Kappatext[k], format='(F0.2)'), /norm)
        ENDIF ELSE BEGIN
          v = TEXT((1+k)*seperation+ .05,0.09, STRING(Densities[k],format = density_format[k]), color=colors[k+1], /norm)
          v = TEXT((1+k)*seperation+ .05,0.06, STRING(Temperatures[k], format = temperature_format[k]), color=colors[k+1], /norm)
          IF Structure.Response EQ 1 THEN v = TEXT((1+k)*seperation+ .05,0.03, STRING(Kappatext[k], format='(F0.2)'), /norm)
        ENDELSE
        ;STOP
      ENDFOR
            

      ; Save a plot if requested
      IF Structure.Save_Plot GE 1 THEN BEGIN


        ; KEEP THIS PART FOR NOW: IT SAVES FILES WITH THE CURRENT DATE
        ; Rob: this is convenient when you're still in the process of playign with the data, as it doesn't overwrites
        ; the existent files. Would a flag that gives the option to chose the format be useful?
        ;
        ;        ; Create name for plot based on current UTC time
        ;        YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
        ;        DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
        ;        IF ULONG(DayOfRun) LT 10 THEN BEGIN
        ;          DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
        ;        ENDIF
        ;        MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
        ;        HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
        ;        MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
        ;        SecondOfRun = STRMID(SYSTIME(/UTC),17,2)
        ;
        ;
        ;        ; Name the file
        ;        IF ((Structure.Save_Plot EQ 1) OR (Structure.Save_Plot EQ 101)) THEN $
        ;        OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
        ;                                HourOfRun, MinuteOfRun, '_', SecondOfRun, '.jpg', $
        ;                                FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
        ;        ELSE IF ((Structure.Save_Plot EQ 2) OR (Structure.Save_Plot EQ 102)) THEN $
        ;        OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
        ;                                HourOfRun, MinuteOfRun, '_', SecondOfRun, '.png', $
        ;                                FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
        ;        ELSE IF ((Structure.Save_Plot EQ 3) OR (Structure.Save_Plot EQ 103)) THEN $
        ;        OutputFilename = STRING(MonthOfRun, DayOfRun, '_', YearOfRun, '_', $
        ;                                HourOfRun, MinuteOfRun, '_', SecondOfRun, '.eps', $
        ;                                FORMAT = '(A,A,A,A,A,A,A,A,A,A)') $
        ;        ELSE PRINT, 'Unable to save selected image file type.'


        ; Create as filename the time the meaurement was taken
        ; Mode for the filename
        
        
        ; Step 1 - Create name for file that uses "Structure" parameters to indicate measurement
        ; Create as filename the time the meaurement was taken
        ; Mode for the filename
        CASE Structure.Short_or_Long OF
          0: Mode_for_filename = 'E_SHORT'
          1: Mode_for_filename = 'E_LONG'
        ENDCASE
        MeasurementFilename = STRING(ULONG(Structure.Year), FORMAT='(I0)') + '-'  + STRING(ULONG(Structure.Doy), FORMAT='(I03)') + $
          '-'  + STRING(ULONG(Structure.Hour), FORMAT='(I02)') + '-'  + STRING(ULONG(Structure.Minute), FORMAT='(I02)') + $
          '-' + STRING(ULONG(Structure.Second), FORMAT='(I02)') + '_' + Mode_for_filename


        ; Step 2 - Create Timestring of current run to append to the file, so that we create unique filenames
        ; Create name for plot based on current UTC time
        YearOfRun   = STRMID(SYSTIME(/UTC),20,4)
        DayOfRun    = STRMID(SYSTIME(/UTC),8,2)
        IF ULONG(DayOfRun) LT 10 THEN BEGIN
          DayOfRun  = String('0',STRMID(SYSTIME(/UTC),9,1),format='(A,A)')
        ENDIF
        MonthOfRun  = STRMID(SYSTIME(/UTC),4,3)
        HourOfRun   = STRMID(SYSTIME(/UTC),11,2)
        MinuteOfRun = STRMID(SYSTIME(/UTC),14,2)
        SecondOfRun = STRMID(SYSTIME(/UTC),17,2)

        UniqueTimeFilename = STRING('_T-', MonthOfRun, DayOfRun, '-', YearOfRun, '-', $
          HourOfRun, MinuteOfRun, '-', SecondOfRun, FORMAT = '(A,A,A,A,A,A,A,A,A,A)')

        ; Name the file
        IF ((Structure.Save_Plot EQ 1) OR (Structure.Save_Plot EQ 101)) THEN $
          OutputFilename = STRING(MeasurementFilename, UniqueTimeFilename, '.jpg', FORMAT = '(A,A,A)') $
        ELSE IF ((Structure.Save_Plot EQ 2) OR (Structure.Save_Plot EQ 102)) THEN $
          OutputFilename = STRING(MeasurementFilename, UniqueTimeFilename, '.png', FORMAT = '(A,A,A)') $
        ELSE IF ((Structure.Save_Plot EQ 3) OR (Structure.Save_Plot EQ 103)) THEN $
          OutputFilename = STRING(MeasurementFilename, UniqueTimeFilename, '.eps', FORMAT = '(A,A,A)') $
        ELSE PRINT, 'Unable to save selected image file type.'
      
        ; Save the file if of the 3 types - Can add different image types if desired at a later date
        IF ISA(OutputFilename) NE 0 THEN BEGIN
          ; Change directories for output of plots
          CD, CURRENT = RESTORE_DIR
          CD, POD
          dummy.SAVE, OutputFilename, Border = 10, RESOLUTION = 150
          ; Change back directories
          CD, RESTORE_DIR
        ENDIF
      ENDIF      
      
    ENDIF ;that's the ending for the "L/M EQ -1" if
    
    ; Plot the data, the errorbars, and the model
    ; If Structure.Save_Plot is on then just go ahead and save the plot automatically, with a predetermined name
    
  ENDIF ; that's the "if keyword for MAKE_PLOT was set"
  
  
  ; Produce the electron data product (EDP)
  IF EDP EQ 1 THEN BEGIN
    
    ; Calculate Latitudes and Longitudes
    xSC_RN = xSC_V2_NEPTUNE_PLS/r_Neptune
    ySC_RN = ySC_V2_NEPTUNE_PLS/r_Neptune
    zSC_RN = zSC_V2_NEPTUNE_PLS/r_Neptune
    SC_lat = asin(zSC_RN/r_RN)
    SC_lon = asin(ySC_RN/(r_RN*cos(SC_lat)))
    
    ; Calculate L-Shell                         
    ; These 2 csv files are from Fran's website
    V1Mag = READ_CSV('Spacecraft_trajectory_vip4_can.csv')
    V2Mag = READ_CSV('Spacecraft_trajectory_vip4_can_V2.csv')
    MagTimes = [V1Mag.Field01,V2Mag.Field01]
    LShellVals = [V1Mag.Field05,V2Mag.Field05]
    Data_Dec_Time = DOUBLE(Structure.DOY) + DOUBLE(Structure.Hour)/24. + DOUBLE(Structure.Minute)/(24.*60.) + DOUBLE(Structure.Second)/(24.*3600.)
    indMag = closest(MagTimes[*],Data_Dec_Time)
    LShell = LShellVals[indmag]
    
    ; Calculate Electron Fluxes
    PDF = PLS_DATA_FLUX
    
    EDP_Structure = {Year: 1979, DOY: Structure.DOY, Hour: Structure.Hour, Minute: Structure.Minute, Second: Structure.Second, $
                     x_SC:xSC_RN, y_SC:ySC_RN, z_SC:zSC_RN, r: r_RN, Latitude: SC_lat*180d/(!Pi), Longitude: SC_lon*180d/(!Pi), $
                     LShell: LShell, Cold_DN: Densities[0], Cold_T: Temperatures[0], Hot_DN: Densities[1], Hot_T: Temperatures[1], $
                     Channel01:PDF[0], Channel02:PDF[1], Channel03:PDF[2], Channel04:PDF[3], Channel05:PDF[4], Channel06:PDF[5], $
                     Channel07:PDF[6], Channel08:PDF[7], Channel09:PDF[8], Channel10:PDF[9], Channel11:PDF[10], Channel12:PDF[11], $
                     Channel13:PDF[12], Channel14:PDF[13], Channel15:PDF[14], Channel16:PDF[15], Channel17:PDF[16], Channel18:PDF[17], $
                     Channel19:PDF[18], Channel20:PDF[19], Channel21:PDF[20], Channel22:PDF[21], Channel23:PDF[22], Channel24:PDF[23], $
                     Channel25:PDF[24], Channel26:PDF[25], Channel27:PDF[26], Channel28:PDF[27]}
;------------------ TO HERE
  
    RETURN, EDP_Structure
  ENDIF
  

  ; If uncertainties are requested, return the uncertainty structure. Otherwise return the best parameters.
  ; This works only if the keyword WRITE_UNCERTAINTIES_IN_FILE is not set. If it is the uncertainties will 
  ; be in the same initial file.
  IF KEYWORD_SET(UNCERTAINTIES) EQ 1 THEN RETURN, Structure_Unc ELSE BEGIN
    ; Overwrite with data timestamps in case of changes between ions and electrons
    IF KEYWORD_SET(MAKE_PLOT) THEN BEGIN
      Structure.DOY = Data_DOY
      Structure.Hour = Data_Hour
      Structure.Minute = Data_Min
      Structure.Second = Data_Sec
    ENDIF
    RETURN, Structure
  ENDELSE
  
END

FUNCTION FIT_PARAMETERS_NEPTUNE, X, P, err=err
  ; This function is essential for MPFIT operations. MPFIT calls on it and the parameters are then
  ; passed into the structure accordingly. Once they are in the structure, the function calls on the
  ; model function to return the currents. This function sums up the currents from each species
  ; so that it can return just a total current to MPFIT and MPFIT can iteratively change the parameters
  ; 
  ; Inputs
  ; X - Array of X values to be used by MPFIT
  ; P - Array of parameter values to use
  ; err - Error values on the Data
  ; 
  ; Additional Inputs:
  ; These are passed in via a common block so that the function can have the proper
  ; structure format
  ; 
  ; Structure    - A structure with the values to use - has to be augmented with P that is passed in by MPFIT
  ; P_Locations  - A structure with the locations of the parameters so that they can easily be put into the structure
  ; 
  ; Outputs:
  ; Total_Current - A 512 (M-Mode) or 64 (L-Mode) array that has the total current
  
  ; Read in Structure and P_Locations from a common block that was saved before this function was called
  COMMON Structure_Block
  Structure   = varA
  P_Locations = varB
  ; Put parameters correctly into the Structure variable by using P_Locations
  dum_arr = MAKE_ARRAY(N_TAGS(Structure))
  FOR k = 0, N_TAGS(Structure) - 1 DO dum_arr(k) = P_locations.(k)
  ind_fit = WHERE(dum_arr EQ 1)
  FOR j = 0, N_ELEMENTS(ind_fit) - 1 DO Structure.(ind_fit[j]) = P[j]
  ; Use the PLS model to generate the current
  
  FORWARD_FUNCTION Generate_PLS_Current_Neptune
  Total_Current = Generate_PLS_Current_Neptune(Structure, /RETURN_TOTAL_CURRENT, /SET_FLOAT)
  
  ; Return Total_Current to MPFIT so it can find best fit parameters
  RETURN, Total_Current
END

FUNCTION Generate_PLS_Current_Neptune, Structure, RETURN_TOTAL_CURRENT=RETURN_TOTAL_CURRENT, SET_FLOAT=SET_FLOAT, $
  RETURN_PLASMA_PROPERTIES=RETURN_PLASMA_PROPERTIES, calculate_chi = calculate_chi
  ; This function calculates the current given the structure input. It will read elements from the structure
  ; and use the given composition to produce the simulated current that the program will read. It depends on
  ; both the planet and the type of response (Maxwellian, Kappa, etc.) that is specified
  ;
  ;
  ; Currently works with the following responses:
  ; CUPINT/DCPINT - Convected Maxwellian
  ; LABCUR/LDCUR  - Convected Maxwellian (accounts better for a warm plasma distribution)
  ;
  ; Inputs:
  ; Structure - A structure that holds all of the pertinent information for the model. Information about
  ;             contents of the structure can be found close to the top of VIPER.pro

  ; Modifications
  ; ---------------
  ; LOGAN 1/5/15
  ;   Currently, both CUPINT and LABCUR use essentially the same inputs so they will both be
  ;   set up the same. There may need to be an IF loop here to setup a Kappa distribution differently
  ;   though likely not



  ; Determine Spacecraft (for Neptune is only Voyager 2) and pull appropriate SSEDR information for trajectory and location
  ; All lookup tables are included in CONRD in common blocks
  ;------------------ GX CHANGE FROM HERE
  COMMON VGR2DATA_NEPTUNE
  ;------------------ TO HERE

  ; Setup array for the number of species that current is being calculated for
  ; If the density of a species is 0, then it will not be accounted for because that just
  ; adds zero to the current
  ;
  ; Store arrays of the important arguments (A, Z, n, T, and how the ion is being processed)

  ; Define the number of different ions - all the data on the csv file you don't need should be 0
  n_species = TOTAL(STREGEX(TAG_NAMES(Structure),'^Species_[123456789][0123456789]?$', /BOOLEAN, /FOLD_CASE ))

  ; Establish arrays to put elements in
  Arguments    = MAKE_ARRAY(n_species)
  Masses       = Arguments
  Charges      = Arguments
  Densities    = Arguments
  Temperatures = Arguments
  Kappa_Values = Arguments
    
    
  ; Fill all arrays
  FOR k = 0, n_species - 1 DO BEGIN
    str0 = STRCOMPRESS(STRING('Species_', k+1), /REMOVE_ALL)
    IF k+1 LT 10 THEN str0cmp = STRCMP(str0, TAG_NAMES(Structure), 10, /FOLD_CASE) ELSE $
                        str0cmp = STRCMP(str0, TAG_NAMES(Structure), 11, /FOLD_CASE)
    ind0 = WHERE(str0cmp EQ 1)
    Arguments(k) = Structure.(ind0)

    str1 = STRCOMPRESS(STRING('Species_', k+1, '_A'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str1cmp = STRCMP(str1, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str1cmp = STRCMP(str1, TAG_NAMES(Structure), 13, /FOLD_CASE)
    ind1 = WHERE(str1cmp EQ 1)
    Masses(k) = Structure.(ind1)

    str2 = STRCOMPRESS(STRING('Species_', k+1, '_Z'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str2cmp = STRCMP(str2, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str2cmp = STRCMP(str2, TAG_NAMES(Structure), 13, /FOLD_CASE)
    ind2 = WHERE(str2cmp EQ 1)
    Charges(k) = Structure.(ind2)

    str3 = STRCOMPRESS(STRING('Species_', k+1, '_N'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str3cmp = STRCMP(str3, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str3cmp = STRCMP(str3, TAG_NAMES(Structure), 13, /FOLD_CASE)
    ind3 = WHERE(str3cmp EQ 1)
    Densities(k) = Structure.(ind3)

    str4 = STRCOMPRESS(STRING('Species_', k+1, '_T'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str4cmp = STRCMP(str4, TAG_NAMES(Structure), 12, /FOLD_CASE) ELSE $
                        str4cmp = STRCMP(str4, TAG_NAMES(Structure), 13, /FOLD_CASE)
    ind4 = WHERE(str4cmp EQ 1)
    Temperatures(k) = Structure.(ind4)
      
    str5 = STRCOMPRESS(STRING('Species_', k+1, '_KAPPA'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str5cmp = STRCMP(str5, TAG_NAMES(Structure), 16, /FOLD_CASE) ELSE $
                        str5cmp = STRCMP(str5, TAG_NAMES(Structure), 17, /FOLD_CASE)
    ind5 = WHERE(str5cmp EQ 1)
    Kappa_Values(k) = Structure.(ind5)
  ENDFOR
  
  ; Check densities - If a species has 0 density, then it will not have a contribution to the current
  ind_good = WHERE(Densities GT 0, nind_good)

  ; Augment arrays to only contain relevant ion species
  Arguments    = Arguments(ind_good)
  Masses       = Masses(ind_good)
  Charges      = Charges(ind_good)
  Densities    = Densities(ind_good)
  Temperatures = Temperatures(ind_good)
  Kappa_Values = Kappa_Values(ind_good)

  ; Fill in temperature and density arrays
  n_tied = 0
  FOR k = 0, nind_good - 1 DO BEGIN
    ; Change temperatures to a common temperature if the parameters are indicated to have those values
    IF Structure.Common_Temperature GT 0d THEN BEGIN
      IF ((Arguments(k) EQ 0) OR (Arguments(k) EQ 2) OR (Arguments(k) EQ 4)) THEN $
        Temperatures(k) = Structure.Common_Temperature
    ENDIF
      
    IF Structure.Common_Kappa GT 0D THEN BEGIN
      IF ((Arguments(k) NE 6)) THEN Kappa_Values(k) = Structure.Common_Kappa
    ENDIF
      
    ; Check if density values should be tied together and set them accordingly
    IF ((Arguments(k) EQ 4) OR (Arguments(k) EQ 5)) THEN BEGIN
      n_tied = n_tied + 1
      IF n_tied EQ 1 THEN BEGIN
        Master_Density = Densities(k)
      ENDIF ELSE BEGIN
        Densities(k) = Master_Density*Densities(k)
      ENDELSE
    ENDIF
  ENDFOR

  ; Do a quick check that no temperatures are set to zero. If they are, that ion species will be ignored as
  ; it will lead to non finite results
  ind_good2 = WHERE(Temperatures GT 0d, nind_good2)
  IF nind_good2 NE nind_good THEN BEGIN
    Arguments    = Arguments(ind_good2)
    Masses       = Masses(ind_good2)
    Charges      = Charges(ind_good2)
    Densities    = Densities(ind_good2)
    Temperatures = Temperatures(ind_good2)
    Kappa_Values = Kappa_Values(ind_good2)
  ENDIF

  ; Return plasma properties if requested - Useful for plotting procedure
  IF KEYWORD_SET(RETURN_PLASMA_PROPERTIES) THEN BEGIN
    new_arr = MAKE_ARRAY(4, nind_good2, /DOUBLE)
    new_arr[0,*] = Masses[*]
    new_arr[1,*] = Charges[*]
    new_arr[2,*] = Densities[*]
    new_arr[3,*] = Temperatures[*]
    RETURN, new_arr
  ENDIF


  IF Structure.L_or_M_Mode EQ 1  THEN n_channels = 128 ; M-Mode
  IF Structure.L_or_M_Mode EQ 0  THEN n_channels = 16  ; L-Mode
  IF Structure.L_or_M_Mode EQ -1 THEN n_channels = 28  ; E-Mode
  ; Set up array to store current information
  ; First element is the Cup (A, B, C, D)
  ; Second element is the current/ion species
  ; Third element is the channel number of the PLS instrument
  PLS_Currents = MAKE_ARRAY(4, nind_good2, n_channels, /DOUBLE)

  
  ; Get some of the constants for the arguments of the CUPINT/DCPINT and LABCUR/LDCUR programs
  ; Contains the following variables:
  ;   SS    - Step size for numerical integration
  ;   LDIST - Logical flag to return current values ('FALSE') or reduced current distribution ('TRUE')
  COMMON Generate 

  ; Now that all species should be correctly accounted for, begin calculating the currents
  FOR j = 0, 3 DO BEGIN ; Begin the loop over the Faraday cup
    Cups = ['A','B','C','D']

    ; Put flow speeds into a vector
    VelocityCylindrical = [Structure.V1, Structure.V2, Structure.V3]

    ; Calculate decimal date of the time being used
    Time = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d

    ; Calculate flowspeed into Faraday cup
    U = VGR_CupVelocity(VelocityCylindrical, Cups[j], Structure) ; FILLER - NEED THIS FUNCTION


    ; Convert the temperatures into a thermal speed (km/s)
    ;FORWARD_FUNCTION eV_ThermalSpeed ;
    ; LOGAN 1/6/15 - eV_ThermalSpeed isn't being recognized for some reason
    Thermals = eV_ThermalSpeed(Temperatures, Masses, Kappa_Values) ; FILLER - NEED THIS FUNCTION - MAY NEED TO ADJUST OR CALL ON Electron Thermals or Something?

    ; Begin loop over relevant species
    FOR k = 0, nind_good2 - 1 DO BEGIN

      ; Calculate response of the instrument
      ; CUPINT/DCPINT response
      IF Structure.Response EQ 0 AND Structure.L_or_M_Mode GE 0 THEN BEGIN
          
        ; Calculate currents for A, B, and C cups
        IF j NE 3 THEN BEGIN
          PLS_Currents(j,k,*) = CUPINT(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode,[1,n_channels],SS,LDIST)
          ; Calculate response for D cup
        ENDIF ELSE BEGIN
          PLS_Currents(j,k,*) = DCPINT(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode,[1,n_channels],SS,LDIST)
        ENDELSE

      ; LABCUR/LDCUR response
      ENDIF ELSE IF Structure.Response EQ 1 AND Structure.L_or_M_Mode GE 0 THEN BEGIN
         
        ; Calculate Currents for A, B, and C cups
        IF j NE 3 THEN BEGIN
          PLS_Currents(j,k,*) = LABCUR(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode, $
            [1,n_channels],LDIST,Cups[j])
          ; Calculate response for D cup
        ENDIF ELSE BEGIN
          PLS_Currents(j,k,*) = LDCUR(Densities[k],Thermals[k],U,Masses[k],Charges[k],Structure.L_or_M_Mode, $
            [1,n_channels],LDIST,Cups[j])
        ENDELSE
        
      ; CUPINT/DCPINT with a Kappa distribution addition
      ENDIF ELSE IF Structure.Response EQ 2 THEN BEGIN
          
        ; Not programmed to work
        PRINT, 'No Kappa distribution has been programmed yet. Returning currents of 0'
        PLS_Currents(j,k,*) = 0
          
      ENDIF ELSE IF Structure.Response EQ 0 AND Structure.L_or_M_Mode EQ -1 THEN BEGIN
        IF j NE 3 THEN BEGIN
          PLS_Currents(j,k,*) = 0d ; Don't fill in A-C cups because those didn't measure electrons
        ENDIF ELSE BEGIN
          ;STOP
          PLS_Currents(j,k,*) = SITTLER_RESPONSE(Densities[k],Temperatures[k],U,1d,1d,-1,[1,n_channels])
        ENDELSE
      ENDIF ELSE IF Structure.Response EQ 1 AND Structure.L_or_M_Mode EQ -1 THEN BEGIN
        IF j NE 3 THEN BEGIN
          PLS_Currents(j,k,*) = 0d ; Don't fill in A-C cups because those didn't measure electrons
        ENDIF ELSE BEGIN
          PLS_Currents(j,k,*) = ELECTRON_KAPPA_RESPONSE(Densities[k],Temperatures[k],U,1d,1d,-1,[1,n_channels],Kappa_Values[k])
        ENDELSE
      ENDIF

    ENDFOR
  ENDFOR

  ; STOP
  ; Have the option to return the current already summed
  ; If ions it is current, if electrons it is the RDF
  IF KEYWORD_SET(RETURN_TOTAL_CURRENT) THEN BEGIN
    ; It always goes in the order of Cup A, Cup B, Cup C, and then Cup D
    ; PLS_Currents is 4 x Number of Species x Number of Channels
    CASE Structure.L_or_M_Mode OF ; -1:E, 0:L, 1:M
      1  : n_channels = 128 ; Number of channels per cup
      0  : n_channels = 16
      -1 : n_channels = 28
    ENDCASE
    Total_Current = MAKE_ARRAY(n_channels*4, /DOUBLE)

    Modified_Array = TOTAL(PLS_Currents, 2) ; Reduce it to a 4 x Number of Channel Array, summing together the current of each species
    ; Modify Total_Current to have 512 (M-Mode), 64 (L-Mode), or 112 (E-Mode) successive channels
    Total_Current(0*n_channels:1*n_channels-1) = Modified_Array(0, 0:n_channels-1)
    Total_Current(1*n_channels:2*n_channels-1) = Modified_Array(1, 0:n_channels-1)
    Total_Current(2*n_channels:3*n_channels-1) = Modified_Array(2, 0:n_channels-1)
    Total_Current(3*n_channels:4*n_channels-1) = Modified_Array(3, 0:n_channels-1)
    IF KEYWORD_SET(SET_FLOAT) THEN RETURN, FLOAT(TOTAL_CURRENT) ELSE RETURN, Total_Current
  ENDIF ELSE BEGIN
    ; Return array of currents per individual species
    IF KEYWORD_SET(SET_FLOAT) THEN RETURN, FLOAT(PLS_Currents) ELSE RETURN, PLS_Currents
  ENDELSE

END

FUNCTION GENERATE_ELECTRON_FLUXES_NEPTUNE, DOY, OUTPUTNAME
  ; This procedure extracts data for an individual day and creates a .sav file of a structure with flux information

;------------------ GX CHANGE FROM HERE
  ; Generate stucture of the appropriate tags and modes - Nearly all tags are unnecessary but keeping with the format to prevent crashes
  Structure = {Year: 1979, DOY: 64, Hour: 10, Minute: 16, Second: 0, Spacecraft: 2, Planet_Number: 5, Response: 0, $
               L_or_M_Mode: -1, Save_Plot: 0,  Fit: 0, Iterations: 25, Channels_CupA: 001128 , Channels_CupB: 001128, $
               Channels_CupC: 001128, Channels_CupD: 001128, Vary_V1: 0, Vary_V2: 0, Vary_V3: 0, Vary_Density: 0, $
               Vary_Temperature: 0, Analytical_Vphi: 0d, V1: 0d, V2: 0d, V3: 0d, Common_Temperature: 0d, Delamere_Composition: 0d, $
               Species_1: 0, Species_1_A: 1d, Species_1_Z: 1d, Species_1_n: 50d, Species_1_T: 3d, $
               Species_2: 0, Species_2_A: 1d, Species_2_Z: 1d, Species_2_n: 50d, Species_2_T: 3d, $
               Species_3: 0, Species_3_A: 1d, Species_3_Z: 1d, Species_3_n: 50d, Species_3_T: 3d, $
               Species_4: 0, Species_4_A: 1d, Species_4_Z: 1d, Species_4_n: 50d, Species_4_T: 3d, $
               Species_5: 0, Species_5_A: 1d, Species_5_Z: 1d, Species_5_n: 50d, Species_5_T: 3d, $
               Species_6: 0, Species_6_A: 1d, Species_6_Z: 1d, Species_6_n: 50d, Species_6_T: 3d, $
               Species_7: 0, Species_7_A: 1d, Species_7_Z: 1d, Species_7_n: 50d, Species_7_T: 3d}

  Structure.DOY = DOY
  PLS_E_Data = PLS_DATA_READ(Structure, /ELECTRON_FLUX)

  ; Analytical equations for the energy widths and energy values

  Delta_E = MAKE_ARRAY(28)
  Eenergy = MAKE_ARRAY(28)
  FOR k = 0, 27 DO BEGIN
    if k LE 15 then begin
      Delta_E(k) = ((60.*(10.^((double(k)+1.)/32.))-50) - (60.*(10.^((double(k))/32.))-50))
      Eenergy(k) = (((60.*(10.^((double(k)+1.)/32.))-50) + (60.*(10.^((double(k))/32.))-50)))/2.
    endif else if k EQ 16 then begin
      Delta_E(k) = (60.*(10.^(5./8.))-50.)-((60.*(10.^(16./32.))-50.))
      Eenergy(k) = ((60.*(10.^(5./8.))-50.)+((60.*(10.^(16./32.))-50.)))/2.
    endif else if k GE 17 then begin
      Delta_E(k) = ((60.*(10.^((double(k)+1.-12.)/8.))-50) - (60.*(10.^((double(k)-12.)/8.))-50))
      Eenergy(k) = (((60.*(10.^((double(k)+1.-12.)/8.))-50) + (60.*(10.^((double(k)-12.)/8.))-50)))/2.
    endif
  endfor
  
  ; Determine number of entries
  nel = SIZE(PLS_E_Data.timearray)
  nel = nel(1)
  
  ; Constants
  Tn = 0.56
  dOmega = 2.8
  Area = 100. ; cm^2
  q = 1.6d-19
  
  Fluxes = MAKE_ARRAY(nel, 28)
  FOR j = 0, nel-1 DO BEGIN
    Fluxes(j,*) = PLS_E_Data.currentarray(j,*)*1d-15 /(Area*delta_E(*)*Tn*q*dOmega)
  ENDFOR
  
  ; Calculate radial distances
  ; Read in SSEDR data/basically Spice kernel
 IF Structure.Spacecraft EQ 2 THEN BEGIN
    COMMON VGR2Data_NEPTUNE
    SSEDR_Data = VGR2SSEDR_NEPTUNE
  ENDIF

  ; Time in SSEDR file
  SSEDR_Time = SSEDR_Data[1,*] + (SSEDR_Data[2,*] + SSEDR_Data[3,*]/60d + SSEDR_Data[4,*]/3600d)/24d

  R = PLS_E_Data.timearray ; Generate array of appropriate size
  ; Find SSEDR index of data to find S/C location for each value by looping through
  FOR k = 0, nel-1 DO BEGIN
    SSEDR_ind = CLOSEST(TRANSPOSE(SSEDR_Time[0,*]), PLS_E_Data.timearray[k])
    ; Get X, Y, and Z locations of spacecraft and compute radial distance (R)
    xSC_V2_NEPTUNE_PLS = SSEDR_Data[6,SSEDR_ind] ; V2_NEPTUNE_PLS cartesian coordinates
    ySC_V2_NEPTUNE_PLS = SSEDR_Data[7,SSEDR_ind]
    zSC_V2_NEPTUNE_PLS = SSEDR_Data[8,SSEDR_ind]
    
    ; Calculate radial distance of S/C from planet core [r_Neptune]
    R[k] = SQRT(xSC_V2_NEPTUNE_PLS*xSC_V2_NEPTUNE_PLS+ySC_V2_NEPTUNE_PLS*ySC_V2_NEPTUNE_PLS+zSC_V2_NEPTUNE_PLS*zSC_V2_NEPTUNE_PLS)/r_Neptune
  ENDFOR
;------------------ TO HERE

  EFlux = {DeltaE:Delta_E, EnergyValue:Eenergy, Times:PLS_E_Data.timearray, $
           R:R, DeltaI:PLS_E_Data.currentarray, Flux:Fluxes}
  SAVE, EFlux, FILENAME = OUTPUTNAME
;stop
END




FUNCTION read_file_Neptune, Data_Filename, CHANGE_DIRECTORY=CHANGE_DIRECTORY
  IF KEYWORD_SET(CHANGE_DIRECTORY) THEN CD, CHANGE_DIRECTORY, CURRENT = RESTORE_DIR
  openr, lun, Data_Filename, /get_lun       ; open file
  ; make lists for each type of data that will be read in from the file
  EcurrentArrays = list()         ; list for all arrays of E currents
  cupA_LcurrentArrays = list()    ; list for all arrays of L currents in cup A
  cupA_McurrentArrays = list()    ; list for all arrays of M currents in cup A
  cupB_LcurrentArrays = list()    ; list for all arrays of L currents in cup B
  cupB_McurrentArrays = list()    ; list for all arrays of M currents in cup B
  cupC_LcurrentArrays = list()    ; list for all arrays of L currents in cup C
  cupC_McurrentArrays = list()    ; list for all arrays of M currents in cup C
  cupD_LcurrentArrays = list()    ; list for all arrays of L currents in cup D
  cupD_McurrentArrays = list()    ; list for all arrays of M currents in cup D
  Etimes = list()        ; list for all times at which E data were taken
  Ltimes = list()         ; list for all times at which L data were taken
  Mtimes = list()         ; list for all times at which M data were taken

  ; Now actually start reading the file.
  WHILE ~EOF(lun) DO BEGIN
    ; read in time & garbage lines
    firstline = ""               ; predefine to be a string
    readf, lun, firstline        ; read in the first line of data
    ; This first line contains the time and status (which tells us the mode).
    ; Now there are 2 garbage lines:
    garbageLine0 = ""            ; predefine lines of garbage to be strings
    garbageLine1 = ""
    readf, lun, garbageLine0     ; read in garbage lines (will be ignored)
    readf, lun, garbageLine1
    
    ; Alter the firstline string to split into time and status parts
    newFirstLine = strmid(firstline, 16)       ; get rid of beginning unnecessary characters
    splitFirstLine = strsplit(newFirstLine, 'status', /extract)  ; separate time & status
    time = splitFirstLine[0]               ; time is the first part of the split string
    status = strmid(splitFirstLine[1],1)  ; status is second part of split string, minus extra spaces
    ; read status to determine next steps
    mode = strmid(status, 0, 2)  ; This tells us the mode of this record
    ;                               ES    Electrons Short
    ;                               EL    Electrons Long
    ;                               LS    L-mode Short
    ;                               LL    L-mode Long
    ;                               ML    M-mode Long (there is no MS for Neptune)
    ;                               
    ; Now we read the data into different lists based on what the mode is (using an IF test).
    ; The IF test below basically does the following for each section of data:
    ; puts time in whatever times list it belongs in
    ; makes an array to read currents into
    ; reads currents into said array
    ; tacks that array onto the end of whatever list it belongs in
    ; for L or M modes, repeats for every cup
    IF ((mode EQ 'LS') OR (mode EQ 'LL')) THEN BEGIN     ; L modes
        Ltimes.add, time                     ; add the time string to the L-modes times list
      ; L modes have 4 cups that measure 16 currents, with each array of currents separated by garbage in the file.
      ; predefine the 4 current arrays (of doubles) and 3 garbage lines:
      garbage0 = ""
      dummyCurrentA = dblarr(16)
      garbage1 = ""
      dummyCurrentB = dblarr(16)
      garbage2 = ""
      dummyCurrentC = dblarr(16)
      garbage3 = ""
      dummyCurrentD = dblarr(16)
      ; now read in the data and assign it to the appropriate lists
      readf, lun, garbage0                      ; will be ignored
      readf, lun, dummyCurrentA                 ; read in cup A currents
      readf, lun, garbage1                      ; will be ignored
      readf, lun, dummyCurrentB                 ; read in cup B currents
      readf, lun, garbage2                      ; will be ignored
      readf, lun, dummyCurrentC                 ; read in cup C currents
      readf, lun, garbage3                      ; will be ignored
      readf, lun, dummyCurrentD                 ; read in cup D currents
      ; weed out negative (garbage) values to prevent issues with log plots later
      FOR i = 0, 15 DO BEGIN
        IF (dummyCurrentA[i] LE 1.) THEN dummyCurrentA[i] = 1.
        IF (dummyCurrentB[i] LE 1.) THEN dummyCurrentB[i] = 1.
        IF (dummyCurrentC[i] LE 1.) THEN dummyCurrentC[i] = 1.
        IF (dummyCurrentD[i] LE 1.) THEN dummyCurrentD[i] = 1.
      ENDFOR   ; every current in each of the four cups
      cupA_LcurrentArrays.add, dummyCurrentA    ; add this array to the end of the list of cup A's L-mode currents
      cupB_LcurrentArrays.add, dummyCurrentB    ; add this array to the end of the list of cup B's L-mode currents
      cupC_LcurrentArrays.add, dummyCurrentC    ; add this array to the end of the list of cup C's L-mode currents
      cupD_LcurrentArrays.add, dummyCurrentD    ; add this array to the end of the list of cup D's L-mode currents
  
  
    ENDIF ELSE IF (mode EQ 'ML') THEN BEGIN   ; M modes
      Mtimes.add, time                     ; add the time string to the M-modes times list
      ; M modes have 4 cups that measure 128 currents, with each array of currents separated by garbage in the file.
      ; predefine the 4 current arrays (of doubles) and 3 garbage lines:
      garbage0 = ""
      dummyCurrentA = dblarr(128)
      garbage1 = ""
      dummyCurrentB = dblarr(128)
      garbage2 = ""
      dummyCurrentC = dblarr(128)
      garbage3 = ""
      dummyCurrentD = dblarr(128)
      ; now read in the data and assign it to the appropriate lists
      readf, lun, garbage0                      ; will be ignored
      readf, lun, dummyCurrentA                 ; read in cup A currents
      readf, lun, garbage1                      ; will be ignored
      readf, lun, dummyCurrentB                 ; read in cup B currents
      readf, lun, garbage2                      ; will be ignored
      readf, lun, dummyCurrentC                 ; read in cup C currents
      readf, lun, garbage3                      ; will be ignored
      readf, lun, dummyCurrentD                 ; read in cup D currents
      ; weed out negative (garbage) values to prevent issues with log plots later
      FOR i = 0, 127 DO BEGIN
        IF (dummyCurrentA[i] LE 1.) THEN dummyCurrentA[i] = 1.
        IF (dummyCurrentB[i] LE 1.) THEN dummyCurrentB[i] = 1.
        IF (dummyCurrentC[i] LE 1.) THEN dummyCurrentC[i] = 1.
        IF (dummyCurrentD[i] LE 1.) THEN dummyCurrentD[i] = 1.
      ENDFOR   ; every current in each of the four cups
      cupA_McurrentArrays.add, dummyCurrentA    ; add this array to the end of the list of cup A's M-mode currents
      cupB_McurrentArrays.add, dummyCurrentB    ; add this array to the end of the list of cup B's M-mode currents
      cupC_McurrentArrays.add, dummyCurrentC    ; add this array to the end of the list of cup C's M-mode currents
      cupD_McurrentArrays.add, dummyCurrentD    ; add this array to the end of the list of cup D's M-mode currents
  
  
    ENDIF ELSE IF ((mode EQ 'ES') OR (mode EQ 'EL')) THEN BEGIN   ; E modes
      Etimes.add, time                    ; add the time string to the E-modes times list
      ; E1 modes have 1 cup that measures 16 currents
      dummyCurrent = dblarr(28)            ; predefine array of 28 currents (doubles)
      readf, lun, dummyCurrent             ; read in currents
      ; weed out negative (garbage) values to prevent issues with log plots later
      FOR i = 0, 27 DO BEGIN
        IF (dummyCurrent[i] LE 1.) THEN dummyCurrent[i] = 1.
      ENDFOR   ; every current in the array
      EcurrentArrays.add, dummyCurrent    ; add this array to the end of the list of the E1 currents
    ENDIF        ; every mode has been accounted for
  ENDWHILE       ; we haven't encountered the end of the file
  

  ; convert the lists of times into string arrays for use with convert_times function later
  LtimesArray = Ltimes.toarray(type=STRING)
  MtimesArray = Mtimes.toarray(type=STRING)
  EtimesArray = Etimes.toarray(type=STRING)

  ; Now, make anonymous structures to hold the data.
  cupA = {L: cupA_LcurrentArrays, M: cupA_McurrentArrays}    ; structure to hold L & M mode data for cup A
  cupB = {L: cupB_LcurrentArrays, M: cupB_McurrentArrays}    ; structure to hold L & M mode data for cup B
  cupC = {L: cupC_LcurrentArrays, M: cupC_McurrentArrays}    ; structure to hold L & M mode data for cup C
  cupD = {L: cupD_LcurrentArrays, M: cupD_McurrentArrays}    ; structure to hold L & M mode data for cup D
  Emodes = {E: EcurrentArrays}                              ; structure to hold E mode data
  times = {L: LtimesArray, M: MtimesArray, E: EtimesArray}  ; structure to hold all time arrays
  ; now make one big anonymous structure to hold the other structures, so all info can be returned at once
  data = {cupA:cupA, cupB:cupB, cupC:cupC, cupD:cupD, Emodes:Emodes, times:times}
  free_lun, lun
  IF KEYWORD_SET(CHANGE_DIRECTORY) THEN CD, RESTORE_DIR
  return, data       ; return the structure full of data
END         ; read_file



FUNCTION Extrapolate_To_Equator, Structure, SC_Position, r_planet, ISOTROPIC=ISOTROPIC, LATDIST=LATDIST  
  ; GX: NOT WORKING FOR NEPTUNE

  ; The purpose of this function is to take the measurements from the Voyager spacecraft locations
  ; and then determine the plasma properties at the centrifugal equator
  ;
  ; It currently uses a magnetic field line tracing code from Drake Rainquist, using the Khurana model
  ; to determine where the location of the centrifugal equator is based off of the initial location of
  ; the Voyager spacecraft
  ;
  ; The code then uses the latdistVoy code written by Fran Bagenal (Oct 2016) to determine the plasma properties
  ; at the centrifugal equator, given the location of the spacecraft and centrifugal equator
  ;
  ; From these calculations, the plasma properties at the centrifugal equator are saved as a csv file and returned
  ;
  ; Inputs
  ; Structure     - The same structure that is used in the rest of the VIPER code
  ; SC_Position   - Location of the spacecraft in the coordination system we need (V2_NEPTUNE_PLS for Neptune)
  ; r_planet      - The radius of the planet
  ;
  ; Keyword
  ; CSV_SAVE_NAME - Name of csv file to save. If set, the next keyword is automatically set.
  ; CSV_SAVE      - Automatically set by previous keyword. If this one is set and the previous keyword is not, then
  ;                 the csv will automatically be named as filename_CentEq.csv, where filename is what the previous
  ;                 csv was called
  ; ISOTROPIC     - Keyword argument to use the isotropic calculation that doesn't take into account the magnetic field
  ;                 for calculating nmax
  ;
  ; Outputs
  ; Structure_CentEq - Structure of different format than input structure with the time of the measurement, and the plasma
  ;                    properties at the equator

  ; Determine time in J2000 seconds
  ;
  ; t_J2000 = UTC_J2000(Structure.Year, Structure.DOY, Structure.Hour, Structure.Minute, Structure.Second, /SECONDS)
  ;

;------------------ GX: INCOMPLETE!!!!! DO NOT USE / USE WITH YOUR OWN RESPONSIBILITY
;                       
;------------------   
;------------------ GX: CHANGE STUFF FROM HERE   
  ; Unpack Position Vector
  R_plan = r_planet
  X_SC = SC_Position[0]/R_plan
  Y_SC = SC_Position[1]/R_plan
  Z_SC = SC_Position[2]/R_plan
  RAD2DEG = 180d/!PI

  ; Convert S/C location to Spherical coordinates
  ; Need Radial reference distance and centrifugal latitude for Fran's code and r, theta, and phi for Drake's code
  r_SC = SQRT(X_SC*X_SC + Y_SC*Y_SC + Z_SC*Z_SC)
  phi_SC = ATAN(Y_SC, X_SC) * RAD2DEG
  theta_SC = ACOS(Z_SC/r_SC) * RAD2DEG

  rho_SC = SQRT(X_SC*X_SC + Y_SC*Y_SC)
  theta_cent_SC = ASIN(Z_SC/r_SC) * RAD2DEG

  ; Time arrays used on a different computer to convert from UTC to et for Khurana model
  ; Calculate UTC time string of measurement and match index to the same index of matching ephemeris
  ; time file
  y = ROUND(Structure.Year)
  D = Structure.DOY
  H = Structure.Hour
  M = Structure.Minute
  S = Structure.Second

  str2 = STRCOMPRESS(STRING(y,'-',D,'T',H,':',M,':',S,'.000', $
    FORMAT = '(A,A,I03,A,I02,A,I02,A,I02,A)'), /REMOVE_ALL)

  IF Structure.Spacecraft EQ 1 THEN BEGIN
    RESTORE, 'V1_JUP_M_UTC_Times.sav'
    RESTORE, 'V1_JUP_M_et_Times.sav'
  ENDIF ELSE BEGIN
    RESTORE, 'V2_JUP_M_UTC_Times.sav'
    RESTORE, 'V2_JUP_M_et_Times.sav'
  ENDELSE

  ind_str = WHERE(STRMATCH(str,str2))
  IF ind_str EQ -1 THEN BEGIN
    PRINT, 'No matching time for measurement, defaulting to first element in array so code does not crash.'
    ind_str = 0
  ENDIF ELSE BEGIN
    t_J2000 = et(ind_str)
  ENDELSE

  ; Determine location of Centrifugal Equator
  CELoc = footpoint_logan(r_SC, theta_SC, phi_SC, 'VIP4', [0.003, 0.1], t_J2000, /RVARY, /LOG, /CENTEQ, /REVDIR) ; Khurana model
  ;CELoc = footpoint_logan(r_SC, theta_SC, phi_SC, 'Khurana', [0.001, 0.1], /RVARY, /LOG, /CENTEQ, /REVDIR) ; VIP4 with Connerney current sheet
  
  ; To Calculate the full trace, you would want to use these lines. Then, assuming isotropic plasma, could include the calculation to every point from ref point
  ;trc     = footpoint_logan(r_SC, theta_SC, phi_SC, 'VIP4', [0.001, 0.1], /RVARY, /LOG, /FULL) ; Trace from measurement towards planet, returning every step
  ;trc_rev = footpoint_logan(r_SC, theta_SC, phi_SC, 'VIP4', [0.001, 0.1], /RVARY, /LOG, /FULL, /REVDIR) ; Do reverse direction to the planet, returning every step
  ;full_trc = [trc,trc_rev] ; The matrices don't line up perfectly - allows for line to be drawn to s/c location
  
  ; Convert centrifugal equator location to appropriate coordinates for Fran's code (if necessary)
  rho_CE = CELoc[0]
  phi_CE = CELoc[1]
  z_CE   = CELoc[2]

  ; Calculation in latdist: latitude is +90 at North pole and -90 at South. Not 0 to 180.
  ; Z = R * SIN(THETA)
  theta_cent_CE = ASIN(z_CE/rho_CE)

  ; Get properties of plasma from SC measurement
  Plasma_Properties = Generate_PLS_Current_Neptune(Structure, /RETURN_PLASMA_PROPERTIES)
  Aion = Plasma_Properties[0,*]
  Zion = Plasma_Properties[1,*]
  nion = Plasma_Properties[2,*]
  Tion = Plasma_Properties[3,*]

  ; Calculate fraction of corotation
  fcorot = Structure.V2/(12.6d*rho_SC)

  ; Fill arrays for latdist code
  scpos = [rho_SC, theta_cent_SC]
  CEpos = [rho_CE, theta_cent_CE]

  ; Use Fran's code to get the density at the centrifugal equator
  IF KEYWORD_SET(LATDIST) THEN BEGIN
    ; This program should be run from the Program "Field_Lines_VIPER.pro" using "Distribution_FL.pro"
    ; It will call this program and return a large format, including measurements of latitude along the field line
    nion_dist = latdistVoy(scpos, CEpos, nion, Tion, Aion, Zion, fcorot, /LATDIST)
    RETURN, nion_dist
  ENDIF ELSE BEGIN
    nion_CE = latdistVoy(scpos, CEpos, nion, Tion, Aion, Zion, fcorot)
  ENDELSE

  ; Extract values for the appropriate elements of the structure below
  Sn = MAKE_ARRAY(10, /DOUBLE)
  ST = MAKE_ARRAY(10, /DOUBLE)

  ; Fill total electron values (always the last from latdistVoy)
  Sn[9] = nion_CE[-1]
  ST[9] = Tion[-1] ; ???

  ; Counter for primary Oxygen species or hot oxygen species
  nO = 0

  ; Bulky way to code it - Does this part sort of manually for structure filling later
  FOR i = 0, N_ELEMENTS(Aion)-1 DO BEGIN
    A0 = ROUND(Aion(i))
    Z0 = ROUND(Zion(i))
    CASE A0 OF
      1: BEGIN
        Sn[5] = nion_CE[i]
        ST[5] = Tion[i]
      END
      16: CASE Z0 OF
      1: BEGIN
        nO = nO + 1
        IF nO EQ 1 THEN BEGIN
          Sn[0] = nion_CE[i]
          ST[0] = Tion[i]
        ENDIF
        IF nO EQ 2 THEN BEGIN
          Sn[7] = nion_CE[i]
          ST[7] = Tion[i]
        ENDIF
      END
      2: BEGIN
        Sn[1] = nion_CE[i]
        ST[1] = Tion[i]
      END
    ENDCASE
    23: BEGIN
      Sn[6] = nion_CE[i]
      ST[6] = Tion[i]
    END
    32: CASE Z0 OF
    3: BEGIN
      Sn[2] = nion_CE[i]
      ST[2] = Tion[i]
    END
    2: BEGIN
      Sn[3] = nion_CE[i]
      ST[3] = Tion[i]
    END
    1: BEGIN
      Sn[4] = nion_CE[i]
      ST[4] = Tion[i]
    END
  ENDCASE
  64: BEGIN
    Sn[8] = nion_CE[i]
    ST[8] = Tion[i]
  END
ENDCASE
ENDFOR

; Rewrite to appropriate output structure
; Set up similarly to the previous
Structure_CentEq = {Year: Structure.Year, DOY: Structure.DOY, Hour: Structure.Hour, Minute: Structure.Minute, Second: Structure.Second, $
  Spacecraft: Structure.Spacecraft, Planet_Number: Structure.Planet_Number, $
  SC_rho: rho_SC, SC_phi: phi_SC, SC_z: z_SC, CE_rho: rho_CE, CE_phi: phi_CE, CE_z: z_CE, $
  Species_1_A: 16d, Species_1_Z: 1d, Species_1_n: Sn[0], Species_1_T: ST[0], $ ; Primary O+ species
  Species_2_A: 16d, Species_2_Z: 2d, Species_2_n: Sn[1], Species_2_T: ST[1], $ ; O++
  Species_3_A: 32d, Species_3_Z: 3d, Species_3_n: Sn[2], Species_3_T: ST[2], $ ; S+++
  Species_4_A: 32d, Species_4_Z: 2d, Species_4_n: Sn[3], Species_4_T: ST[3], $ ; S++
  Species_5_A: 32d, Species_5_Z: 1d, Species_5_n: Sn[4], Species_5_T: ST[4], $ ; S+
  Species_6_A: 1d,  Species_6_Z: 1d, Species_6_n: Sn[5], Species_6_T: ST[5], $ ; H+
  Species_7_A: 23d, Species_7_Z: 1d, Species_7_n: Sn[6], Species_7_T: ST[6], $ ; Na+
  Species_8_A: 16d, Species_8_Z: 1d, Species_8_n: Sn[7], Species_8_T: ST[7], $ ; O+ hot
  Species_9_A: 64d, Species_9_Z: 1d, Species_9_n: Sn[8], Species_9_T: ST[8], $ ; SO2+
  Species_10_A: 1d, Species_10_Z: -1d, Species_10_n: Sn[9], Species_10_T: ST[9]} ; e- from field line code


; Return
RETURN, Structure_CentEq
END
;
