FUNCTION FITTING_RESIDUALS, Structure, PLS_Data_Currents, PLS_Simulated_Currents
  ; This function calculated the fitting residuals, i.e. the difference between the measured current 
  ; and the simulated current
  ; 
  ; INPUTS:
  ;   Structure               - A structure containing information on how to handle all parameters in the fit
  ;   PLS_Data_Currents       - The raw data current, without the threshold values. 
  ;   PLS_Simulated_Currents  - The fitted current. If multiple species are used, this is the sum of all species for each channel
  ; 
  ; OUTPUTS:
  ;     residual_file         - A file with the fitted current residual, at the Fitting_Residuals folder
  ;     
  ; ---------------------
  ; Notes:
  ; - The value of bad measurements in the PLS_Data_Currents (i.e. measurements with only the background current) 
  ;   is zero. The fitting residual for those measurements is set to NaN
  ; ---------------------
  ; Author: George Xystouris, Oct 2024
  ;
  
  ; Check input array dimensions and reform to a single row (some might be in a single row, so this will simply be abundant for them)
  ; ----
  dim_direction = WHERE(SIZE(PLS_Data_Currents,/DIMENSIONS) EQ 4)
  IF dim_direction EQ 0 THEN BEGIN
    PLS_Data_Currents = REFORM(TRANSPOSE(PLS_Data_Currents), N_ELEMENTS(PLS_Data_Currents))
  ENDIF ELSE PLS_Data_Currents = REFORM(PLS_Data_Currents, N_ELEMENTS(PLS_Data_Currents))
  
  dim_direction = WHERE(SIZE(PLS_Simulated_Currents,/DIMENSIONS) EQ 4)
  IF dim_direction EQ 0 THEN BEGIN
    PLS_Simulated_Currents = REFORM(TRANSPOSE(PLS_Simulated_Currents), N_ELEMENTS(PLS_Simulated_Currents))
  ENDIF ELSE PLS_Simulated_Currents = REFORM(PLS_Simulated_Currents, N_ELEMENTS(PLS_Simulated_Currents))
  
  
  ; Calculate residuals. Residuals of bad measurements are NaN
  ; ----
  ; Make residuals array
  curr_diff = MAKE_ARRAY(N_ELEMENTS(PLS_Data_Currents), /DOUBLE, VALUE = !VALUES.D_NAN)
  ; Find bad data
  non_zero_ind = WHERE(PLS_Data_Currents NE 0, good_ind)
  IF good_ind GT 0 THEN BEGIN
    curr_diff[non_zero_ind] = ABS(PLS_Data_Currents[non_zero_ind] - PLS_Simulated_Currents[non_zero_ind])
  ENDIF ELSE BEGIN
    curr_diff = ABS(PLS_Data_Currents - PLS_Simulated_Currents)
  ENDELSE
  
  ; Hardcoded: median of each energy bin
  CASE Structure.L_or_M_Mode OF
    1 :  e_bins = [ 11.1, 13.35, 15.65, 18.05, 20.55, 23.15, 25.85, 28.6, 31.45, 34.45, 37.55, 40.75, 44.1, 47.55, 51.1, 54.8, 58.65, 62.65, 66.8, 71.05, 75.45, 80.05, 84.8, 89.75, 94.9, 100.2, 105.7, 111.4, 117.3, 123.4, 129.75, 136.35, 143.2, 150.3, 157.6, 165.2, 173.1, 181.25, 189.7, 198.5, 207.65, 217.1, 226.85, 236.95, 247.45, 258.35, 269.65, 281.4, 293.55, 306.1, 319.15, 332.7, 346.7, 361.2, 376.3, 391.95, 408.1, 424.85, 442.25, 460.3, 479, 498.35, 518.45, 539.3, 560.9, 583.25, 606.45, 630.5, 655.4, 681.25, 708.05, 735.85, 764.65, 794.45, 825.35, 857.45, 890.7, 925.15, 960.85, 997.85, 1036.25, 1076.05, 1117.3, 1160.05, 1204.35, 1250.3, 1297.95, 1347.35, 1398.55, 1451.6, 1506.6, 1563.6, 1622.7, 1684, 1747.5, 1813.35, 1881.6, 1952.35, 2025.7, 2101.7, 2180.55, 2262.25, 2346.95, 2434.75, 2525.75, 2620.1, 2717.9, 2819.3, 2924.4, 3033.35, 3146.3, 3263.4, 3384.75, 3510.55, 3641, 3776.2, 3916.35, 4061.65, 4212.25, 4368.35, 4530.2, 4698, 4871.9, 5052.2, 5239.1, 5432.8, 5633.65, 5841.85]
    0 :  e_bins = [ 20,  43.35, 74.5, 116, 171.35, 245.2, 343.65, 474.95, 650.05, 883.55, 1194.9, 1610.05, 2163.7, 2902.05, 3886.65, 5199.6]
    -1:  e_bins = [ 12.25, 16.9, 21.9, 27.25, 33, 39.2, 45.85, 53, 60.7, 68.95, 77.8, 87.35, 97.6, 108.6, 120.45, 133.15, 171.35, 245.2, 343.65, 474.95, 650.05, 883.55, 1194.9, 1610.1, 2163.8, 2902.15, 3886.75, 5199.7]
  ENDCASE
  e_bins = [e_bins, e_bins, e_bins, e_bins]
  
  
  ; Create the structure
  Residuals_Structure = {E_bins_middle: e_bins, Current_raw_data: PLS_Data_Currents, Simulated_current: PLS_Simulated_Currents, Residuals: curr_diff}
  Residuals_Header = TAG_NAMES(Residuals_Structure)
  
  
  ; Name the file
  ; --------
  ; Mode
  CASE Structure.L_or_M_Mode OF
    -1 : mode = 'E'
    0  : mode = 'L'
    1  : mode = 'M'
  ENDCASE
  
  ; Integration time
  IF Structure.SHORT_OR_LONG EQ 0 THEN integr_time = ' SHORT' ELSE integr_time = ' LONG'
  
  ;; Response
  ;IF Structure.Response EQ 0 THEN Response_Name = 'CUPINT/DCPINT' ELSE $
  ;IF Structure.Response EQ 1 THEN Response_Name = 'LABCUR/LDCUR'
  
  ; Create filename based on the time of measurement
  Residuals_Filename = STRING(ULONG(Structure.Year), FORMAT='(I0)') + '_'  + STRING(ULONG(Structure.Doy), FORMAT='(I0)') + '_'  + STRING(ULONG(Structure.Hour), FORMAT='(I0)') + $
                       '_'  + STRING(ULONG(Structure.Minute), FORMAT='(I0)') + '_' + STRING(ULONG(Structure.Second), FORMAT='(I0)') + '_' + STRTRIM(mode,2) + '_' + STRTRIM(integr_time,2) + '.csv'

  
  ; Save a CSV file to  Fitting_Residuals folder
  ; ----
  CD, CURRENT = RESTORE_DIR 
  STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/Fitting_Residuals'))
  ; Check if the folder exist, if not create it
  IF ~FILE_TEST(STR_DIR, /DIRECTORY) THEN FILE_MKDIR, STR_DIR
  ; Change directory to where CSVs are located
  CD, STR_DIR
  ; Check if file exists. If so overwrite with warning
  IF FILE_TEST(Residuals_Filename) THEN PRINT, 'Warning: Residuals file ' + Residuals_Filename + ' exists. This process overwrites the original file.'
  WRITE_CSV, Residuals_Filename, Residuals_Structure, HEADER = Residuals_Header
  CD, RESTORE_DIR
  
  
  RETURN, curr_diff

END




