FUNCTION PLS_UNCERTAINTIES, Structure, P_Locations, PLS_Data_Fit, PLS_Error, Model_Name, index, Delta, Current_Threshold, RELATIVE=RELATIVE
  ; This function takes the data and calculates the Chi-Squared value of the best fit parameter, best on Wilson (2015),
  ; doi:10.1002/2014EA000090
  ; The function steps to the side of the best-fit, either absolute or relative step size, and
  ; then, assuming a parabolic shape, calculates the Hessian and Covariance matrix to find the formal
  ; 1 Sigma deviation from the best fit value. 
  ;
  ; Inputs
  ;
  ; Structure    - A structure that contains the information of the parameters and constants that are being fit
  ; P_Locations  - A structure that contains what parts of the above structure are the parameters being varied
  ; PLS_Data_Fit - The data that the model is being compared to. Units: fA
  ; PLS_Error    - The error in the data. Units: fA
  ; Model_Name   - Name of IDL function that will process the structure properly
  ; index        - Index of calculations to avoid in the Chi-Square calculation (Set to -9 if all data points are good)
  ; Delta        - Step size to take away from minimum Chi-Square value
  ; Relative     - Keyword argument if the step size should be a percentage of the value
  ; 
  ;
  ; Outputs
  ;
  ; PLS_Unc  - An array of the number of parameters that contains the uncertainties in the parameters
  ;
  ; ----------------------------------------------
  ; MODIFICATIONS:
  ;
  ; ----------------------------------------------

  ; RAW DATA
  ; ------------ 
  ; Determine number of parameters and their locations within the Structure and extract best-fit parameters
  dum_arr = MAKE_ARRAY(N_TAGS(Structure))
  FOR l = 0, N_TAGS(Structure) - 1 DO dum_arr(l) = P_locations.(l)
  ind1 = WHERE(dum_arr EQ 1)
  Best = MAKE_ARRAY(N_ELEMENTS(ind1))            ; [GX] Array with the species variables (T & n) of the "triggered" fit (=1) in the csv file.
  FOR l = 0, N_ELEMENTS(ind1) - 1 DO Best[l] = Structure.(ind1[l])    ; Best is a matrix with the input variables
  ;STOP
    
  ; Redefine index to not include data points below the Current_Threshold value
  ind_remove = WHERE(PLS_DATA_FIT LE Current_Threshold, nind_remove)
  IF nind_remove GE 1 THEN BEGIN
    IF index[0] NE -9 THEN index = [index, ind_remove] ELSE index = ind_remove
  ENDIF
  
  ;print, index ; testing line
  ;print, Structure.L_or_M_Mode ; testing line
  
;  ind_remove = WHERE(PLS_DATA_FIT LE Current_Threshold, nind_remove)
;  IF nind_remove GE 1 THEN BEGIN
;    IF index[0] NE -9 THEN index = [index, ind_remove] ELSE index = ind_remove
;  ENDIF

  ; Initialize arrays that are needed
  nu    = N_ELEMENTS(Best)              ; number of free parameters
  alpha = MAKE_ARRAY(nu, nu, /DOUBLE)   ; = Hess/2
  
  
  ; GET MODEL VALUES
  ; ------------ 
  Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT)
;  Model = Model
  
  ; STOP
  ; USE RAW DATA AND MODEL VALUES TO CALCULATE UNCERTAINTIES
  ; ------------
  IF index[0] NE -9 THEN Model(index) = 0d ; Set bad data points to be 0
  Model = DOUBLE(Model)
  Error = DOUBLE(PLS_Error)
  IF index[0] NE -9 THEN PLS_DATA_FIT(index) = 0d ; Set bad data points to be 0
  Dat   = PLS_DATA_FIT

  
  ; Exit if there is no data being fit at all
  CASE Structure.L_or_M_Mode OF ; -1:E, 0:L, 1:M
    -1 : N_CHANNELS = 112
    0  : N_CHANNELS = 64
    1  : N_CHANNELS = 512
  ENDCASE
  ;stop
  ; IF Structure.L_or_M_Mode EQ 0 THEN N_CHANNELS = 64 ELSE N_CHANNELS = 512
  ; IF N_ELEMENTS(index) EQ N_CHANNELS THEN BEGIN       ; original command
  IF N_ELEMENTS(index) EQ N_ELEMENTS(Dat) THEN BEGIN    ; replaced 1 Feb 24
    IF TOTAL(ABS(index - LINDGEN(N_ELEMENTS(index)))) EQ 0 THEN BEGIN ; checking the index that has no duplicates
      Best[*] = 0
      PRINT, 'No channel values were provided for calculation. Check the "CHANNELS_CUPA/B/C/D" in input csv file. Returning to main procedure.'
      RETURN, Best
    ENDIF ELSE MESSAGE, 'Index is not 0,1,2,...,n-1. It needs to be in ascending order.'
  ENDIF
  ;STOP
  ; Set a relative step with delta - Multiply best parameter by a certain percent (0.01 for 1 percent)
  ;help, relative ; test line
  IF KEYWORD_SET(RELATIVE) THEN dx = DOUBLE(Delta)*best
  ; Calculate the Chi-Square value at the best fit parameters (should be the minimum value)
  Z0    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error))


  ; Calculate A for the diagonal of the alpha matrix
  FOR ii = 0L, nu-1 DO BEGIN
    fp  = best
    fp(ii) = fp(ii) - dx(ii)
    ; Augment structure with new parameters
    FOR l = 0, N_ELEMENTS(ind1) - 1 DO Structure.(ind1[l]) = fp[l]

    Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT)
    IF index[0] NE -9 THEN Model(index) = 0d
    Model = DOUBLE(Model)
    Z1    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error)) ; this it the chi square, NOT the reduced chi square
    ;IF Z1 LT Z0 THEN PRINT, STRING('Z1 Warning:', i, Z1, Z0)

    fp = best
    fp(ii) = fp(ii) + dx(ii)
    ; Augment structure with new parameters
    FOR l = 0, N_ELEMENTS(ind1) - 1 DO Structure.(ind1[l]) = fp[l]

    Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT)
    IF index[0] NE -9 THEN Model(index) = 0d
    Model = DOUBLE(Model)
    Z2    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error)) ; this it the chi square, NOT the reduced chi square
    ;IF Z2 LT Z0 THEN PRINT, STRING('Z2 Warning:', i, Z2, Z0)

    alpha(ii,ii) = (Z1 + Z2 - 2d*Z0)/(2d*dx(ii)*dx(ii)) ; this it the chi square, NOT the reduced chi square
  ENDFOR

  ; Calculate L/2 for the off-diagonal terms of the alpha matrix
  FOR ii = 0, nu-2 DO BEGIN
    FOR jj = ii+1, nu-1 DO BEGIN
      fp = best
      fp(ii) = fp(ii) + dx(ii)
      fp(jj) = fp(jj) + dx(jj)
      ; Augment structure with new parameters
      FOR l = 0, N_ELEMENTS(ind1) - 1 DO Structure.(ind1[l]) = fp[l]

      Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT)
      IF index[0] NE -9 THEN Model(index) = 0d
      Model = DOUBLE(Model)
      Z1    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error))
      ;IF Z1 LT Z0 THEN PRINT, STRING('Z1 Warning:', i, j, Z1, Z0)

      fp = best
      fp(ii) = fp(ii) - dx(ii)
      fp(jj) = fp(jj) - dx(jj)
      ; Augment structure with new parameters
      FOR l = 0, N_ELEMENTS(ind1) - 1 DO Structure.(ind1[l]) = fp[l]

      Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT)
      IF index[0] NE -9 THEN Model(index) = 0d
      Model = DOUBLE(Model)
      Z2    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error))
      ;IF Z2 LT Z0 THEN PRINT, STRING('Z2 Warning:', i, j, Z2, Z0)

      fp = best
      fp(ii) = fp(ii) - dx(ii)
      fp(jj) = fp(jj) + dx(jj)
      ; Augment structure with new parameters
      FOR l = 0, N_ELEMENTS(ind1) - 1 DO Structure.(ind1[l]) = fp[l]

      Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT)
      IF index[0] NE -9 THEN Model(index) = 0d
      Model = DOUBLE(Model)
      Z3    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error))
      ;IF Z3 LT Z0 THEN PRINT, STRING('Z3 Warning:', i, j, Z3, Z0)

      fp = best
      fp(ii) = fp(ii) + dx(ii)
      fp(jj) = fp(jj) - dx(jj)
      ; Augment structure with new parameters
      FOR l = 0, N_ELEMENTS(ind1) - 1 DO Structure.(ind1[l]) = fp[l]

      Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT)
      IF index[0] NE -9 THEN Model(index) = 0d
      Model = DOUBLE(Model)
      Z4    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error))
      ;IF Z4 LT Z0 THEN PRINT, STRING('Z4 Warning:', i, j, Z4, Z0)

      alpha(ii,jj) = (Z1 + Z2 - Z3 - Z4)/(8d*dx(ii)*dx(jj))
      alpha(jj,ii) = alpha(ii,jj) ; Fill in opposite off diagonals
    ENDFOR
  ENDFOR

  ; Take matrix inverse of Alpha matrix to get covariance matrix
  S2 = INVERT(alpha)
  S  = MAKE_ARRAY(nu,/DOUBLE)
  
  ;STOP
  
  FOR ii = 0L, nu-1 DO BEGIN
    ; LOGAN - Made it absolute value so it doesn't return NaNs in IDL. Not sure if I should do that
    IF S2(ii,ii) LT 0d THEN BEGIN
      PRINT, S2
      PRINT, STRCOMPRESS(STRING('Element: ', ii, ' is negative. Errors are not valid and are being set to zero.'))
      S(*) = 0d
      RETURN, S
    ENDIF
    S(ii) = SQRT(S2(ii,ii)) ; Calculate standard deviation of Element i
  ENDFOR
  
  
  ; Assuming the Z0 is the chi square value, we get to the reduced chi square by dividing by DOF.
  ; This converts the diagonal from a best fit of chi square to a best fit of reduced chi square. 
  ; Steps to get here: (1) Reduced chi square = chi square / DOF. (2) Inverting the matrix puts DOF on the numerator.
  ; (3) We take the square root of diagonal. => hence we get *sqrt(DOF).
  ; 
  ; Calculate Degrees of Freedom                
  ; IF index[0] NE -9 THEN DOF = N_CHANNELS - N_ELEMENTS(index) ELSE DOF = N_CHANNELS   ; original command - this is for the chi square, not the reduced chi squared
  IF index[0] NE -9 THEN DOF = N_ELEMENTS(Dat) - N_ELEMENTS(index) - nu ELSE DOF = N_ELEMENTS(Dat) - nu ; this is the reduced chi square version
  S = DOUBLE(S)*SQRT(DOF) ; Weight properly by the Degrees of Freedom
  
  PRINT, S
  RETURN, S
  
END

FUNCTION CHI_ARRAYS, Structure, Structure_Unc, P_Locations, PLS_Data_Fit, PLS_Error, Model_Name, index, N_Chi_Steps, N_Sigma, str_date
  ; This function computes the chi-squared space for the given parameters in the model based off of the central
  ; location given in the input structure. It will compute the model using the structure and return the values
  ; for comparison to the best fit
  ;
  ; Inputs
  ;
  ; Structure     - A structure that contains the information of the parameters and constants that are being fit
  ; Structure_Unc - A structure that contains the formal 1 Sigma uncertainties in the best parameters
  ; P_Locations   - A structure that contains what parts of the above structure are the parameters being varied
  ; PLS_Data_Fit  - The data that the model is being compared to
  ; PLS_Error     - The error in the data
  ; Model_Name    - Name of IDL function that will process the structure properly
  ; index         - Index of calculations to avoid in the Chi-Square calculation. Set to -9 if all data points are good
  ; N_Chi_Steps   - Number of steps to take in Chi-Square space. Should be an odd value so that it is symmetric about best fit
  ; N_Sigma       - Number of sigma to deviate from the starting location
  ; str_date      - String Name for text file saving
  ;
  ; Outputs
  ;
  ; Chi_Arrs      - An array that contains the Chi-Square values (Matlab is used for plotting quality contour plots given this output)
  
  dum_arr = MAKE_ARRAY(N_TAGS(Structure))
  FOR l = 0, N_TAGS(Structure) - 1 DO dum_arr(l) = P_locations.(l)
  ind1 = WHERE(dum_arr EQ 1)
  Best = MAKE_ARRAY(N_ELEMENTS(ind1))
  One_Sig = Best
  FOR l = 0, N_ELEMENTS(ind1) - 1 DO Best[l]    = Structure.(ind1[l])
  FOR l = 0, N_ELEMENTS(ind1) - 1 DO One_Sig[l] = Structure_Unc.(ind1[l])

  ; Initialize arrays that are needed
  
  nu        = N_ELEMENTS(Best)
  nn        = ((nu*nu)-nu)/2.0
  Chi_Array = MAKE_ARRAY(N_Chi_Steps, N_Chi_Steps, nn, /DOUBLE)
  
  ; Redefine index to not include data points below the 100 fA value
  ; Probably we will need to see if it's a short or long integration time, as they have different threshold - GX 16/Feb/24
  ;
  ; ind_remove = WHERE(PLS_DATA_FIT LE 100, nind_remove)
  ind_remove = WHERE(PLS_DATA_FIT LE 102, nind_remove, /NULL)
  IF nind_remove GE 1 THEN BEGIN
    IF index[0] NE -9 THEN index = [index, ind_remove] ELSE index = ind_remove
  ENDIF
  
  ; Run through the model calculation
  Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT, /CALCULATE_CHI)
  IF index[0] NE -9 THEN Model(index) = 0d ; Set bad data points to be 0
  Model = DOUBLE(Model)
  Error = DOUBLE(PLS_Error)
  IF index[0] NE -9 THEN PLS_DATA_FIT(index) = 0d ; Set bad data points to be 0
  Dat   = PLS_DATA_FIT
  
  IF Structure.L_or_M_Mode EQ 1 THEN N_Channels = 512 ELSE $
                                     N_Channels = 64
                                     
  ; Exit if there is no data being fit at all
  IF N_ELEMENTS(index) EQ N_Channels THEN BEGIN
    Best[*] = 0
    PRINT, 'All data points in Chi-Square calculation were infinite. Returning to main procedure.'
    RETURN, Best
  ENDIF

  ; Adjust esig by the appropriate number of degrees of freedom
  IF index[0] NE -9 THEN DOF = N_Channels - N_ELEMENTS(index) - N_ELEMENTS(Best) ELSE $
                         DOF = N_Channels - N_ELEMENTS(Best)
  
  Esig = One_Sig
  nnn = -1
  FOR ii = 1,nu DO BEGIN
    FOR jj = 1,nu DO BEGIN
      Best_p = Best ; Redefines structure so that values are not augmented by for loops
      ; For loops would create a matrix of elements that is nn by nn. The diagonals are
      ; unimportant since it would compare one free parameter to itself. The off-diagonal
      ; elements are symmetric along the diagonal, so this takes only the upper half
      ; corner of elements into consideration.
      ;
      ; EXAMPLE
      ;
      ; [ 1,1  1,2  1,3]
      ; | 2,1  2,2  2,3|
      ; [ 3,1  3,2  3,3]
      ;
      ; For this 3x3 matrix, it would use (1,2), (1,3), and (2,3) only

      IF jj LE ii THEN CONTINUE
      print, ii, jj
      nnn = nnn + 1

      ; Use the best fit parameters from the CSV file where the fits are stored
      ; Currently allow for the program to vary between 0.1 and 1.9 times the best value since we
      ; are unsure of the sigma values.
      val_1     = Best[ii-1]
      val_2     = Best[jj-1]
      sigma1    = N_Sigma*One_Sig(ii-1)
      sigma2    = N_Sigma*One_Sig(jj-1)

      FOR i = 0, N_Chi_Steps-1 DO BEGIN
        FOR j = 0, N_Chi_Steps-1 DO BEGIN
          StepSize = DOUBLE(FIX(N_Chi_Steps/2.0))
          Best_P[ii-1] = val_1 - sigma1 + i*sigma1/StepSize
          Best_P[jj-1] = val_2 - sigma2 + j*sigma2/StepSize          
          ;print, i ; Checks that the function is working properly since it is slow
          ;print, j

          FOR l = 0, N_ELEMENTS(ind1) - 1 DO Structure.(ind1[l]) = Best_P[l]
          Model = CALL_FUNCTION(Model_Name, Structure, /RETURN_TOTAL_CURRENT, /CALCULATE_CHI)
          IF index[0] NE -9 THEN Model(index) = 0d
          Model = DOUBLE(Model)
          Z0    = TOTAL((Dat-Model)*(Dat-Model)/(Error*Error))/(DOF-1)

          ; Save values into an array
          Chi_Array(j,i,nnn) = DOUBLE(Z0)
        ENDFOR
      ENDFOR
      
      ; Check that the minimum value is in the center
      min_value = MIN(chi_array(*,*,nnn))
      ;PRINT, min_value, format = '(f0.3)'
      min_index =  WHERE(chi_array(*,*,nnn) EQ min_value)
      ;PRINT, min_index
      IF min_index NE ((N_Chi_Steps*N_Chi_Steps - 1)/2) THEN PRINT, "ERROR: MPFIT DID NOT FIND MINIMUM VALUE"
    ENDFOR
  ENDFOR
  
  ; Change directories for output of files
  CD, CURRENT = RESTORE_DIR
  STR_DIR = RESTORE_DIR + PATH_SEP() + 'Output' + PATH_SEP() + 'Chi_Space'
  CD, STR_DIR
  
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;
  ;str_date = '1979_062_10_23_58_CHI_DUMMY_ARRAYS.txt' ; TEST TESTT TEST TEST TEST
  ;
  ;;
  ;;
  ;;
  ;;
  ;;
  ;;
  ;;
  
  ; Save Chi-Space results to a text file
  OPENW, 1, str_date
  PRINTF, 1, Chi_Array, format = '(f0.3)'
  CLOSE, 1
  
  ; Save parameters to a new text file
  par_name = STRCOMPRESS(STRING(ULONG(Structure.Year), '_', ULONG(Structure.DOY), '_', ULONG(Structure.Hour), '_', $
                                ULONG(Structure.Minute), '_', ULONG(Structure.Second), '_CHI_Parameters_MATLAB.txt', $
                                FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02,A)'), /REMOVE_ALL)
  OPENW, 1, par_name
  PRINTF, 1, Best
  CLOSE, 1
  
  ; Save uncertainties to a new text file
  unc_name = STRCOMPRESS(STRING(ULONG(Structure.Year), '_', ULONG(Structure.DOY), '_', ULONG(Structure.Hour), '_', $
                                ULONG(Structure.Minute), '_', ULONG(Structure.Second), '_CHI_Uncertainties_MATLAB.txt', $
                                FORMAT='(I04,A,I03,A,I02,A,I02,A,I02,A,I02,A)'), /REMOVE_ALL)
  OPENW, 1, unc_name
  PRINTF, 1, One_Sig
  CLOSE, 1                             
  
  ; Return to proper directory
  CD, RESTORE_DIR
  
  RETURN, Chi_Array
END
