FUNCTION ADDITIONAL_PLS_FUNCTIONS
  ; FILLER - Necessary for correct compilation by RESOLVE_ROUTINE
END

FUNCTION eV_ThermalSpeed, Temperatures, Masses, Kappas
  ; This function reads in ion temperatures and their atomic masses and converts them to
  ; a thermal speed in km/s
  ;
  ; Inputs: Can be either arrays or singular values
  ; 
  ; Temperatures - The ion temperatures in eV
  ; Masses - The amu mass of the ion
  ; 
  ; Outputs
  ; 
  ; W - The thermal speeds of the ions given their temperatures and masses
  ; Constants for electric charge and proton mass
  ;
  ;Written by Kaleb Bodisch
  ;2016
  
  e = 1.6D-19
  mp = 1.67D-27

  n_T = N_ELEMENTS(Temperatures)
  n_M = N_ELEMENTS(Masses)
  IF n_T NE n_M THEN MESSAGE, 'Error: Ion temperature and mass arrays are not of the same length.'
  
  ; Create array
  W = MAKE_ARRAY(n_T)
  ; Fill array
  FOR in_T = 0, n_T-1 DO BEGIN  
    IF kappas(in_T) EQ 0 THEN W[in_T] = SQRT(2d*e*Temperatures[in_T]                    /(Masses[in_T]*mp))             /1000d ELSE $
                              W[in_T] = SQRT(2d*e*Temperatures[in_T]*(kappas(in_T)-1.5d)/(Masses[in_T]*mp*kappas(in_T)))/1000d
  ENDFOR       
  RETURN, W
END

FUNCTION Good_Channels, sixes, cup, mode
  ; Input
  ; sixes - 4-6 large integer number to indicate the cups being fit. To fit cups 16-111 as example,
  ;         the entry would be 016111
  ; cup   - The Faraday Cup for the input as 0-3 so that the error message displays where the
  ;         error in the entry is
  ;
  ; Output
  ; Both  - A two element array with the channel numbers to fit from the six digit array
  ;
  ; --------------------
  ; MODIFICATIONS:
  ; GX 26 Aug. 2024: Added electrons case
  ;

  ; Establish which cup it is in a string
  IF cup EQ 0 THEN CupString = 'CUP A: '
  IF cup EQ 1 THEN CupString = 'CUP B: '
  IF cup EQ 2 THEN CupString = 'CUP C: '
  IF cup EQ 3 THEN CupString = 'CUP D: '

  i = STRTRIM(sixes, 1)                ; Trim blank spaces off string
  j = STRLEN(i)                        ; Get size of string
  a = STRTRIM(LONG(0),1)               ; Create a string "0"

  ; Do the calculation for L-Mode (Channel range of 1-16)
  IF Mode EQ 0 THEN BEGIN
    IF j LT 4 THEN BEGIN                 ; If less than four digits cant be two numbers
      IF ((i EQ 0) OR (i EQ 1)) THEN BEGIN
        IF i EQ 0 THEN BEGIN             ; If no fitting (0) then use no range
          Both = [0,0]
          GOTO, SkipErrorStatement1
        ENDIF
        IF i EQ 1 THEN BEGIN             ; If fitting turned on (1) then use full range
          Both = [1,16]
          GOTO, SkipErrorStatement1
        ENDIF
      ENDIF
      PRINT, CupString, "Input is not long enough. Defaulting to full range."
      Both =[1, 16]                     ; Default to full range
      SkipErrorStatement1:
    ENDIF

    IF j GT 6 THEN BEGIN                 ; If greater than 6 digits cant be two numbers
      PRINT, CupString, "Input is too long. Defaulting to full range."
      Both = [1, 16]                     ; Default to full range
    ENDIF

    IF j EQ 4 THEN sixes= STRCOMPRESS(a+a+i)          ; Add zeros onto front for parsing purposed
    IF j EQ 5 THEN sixes= STRCOMPRESS(a+i)            ; Add zero onto front for parsing purposed

    IF (j GT 3) AND (j lt 7) THEN BEGIN

      last_three  = LONG(STRMID(sixes, 2, 3, /REVERSE_OFFSET))         ; Get last three elements(end channel)
      first_three = LONG(STRMID(sixes, 5, 3, /REVERSE_OFFSET))         ; Get first three elements(start channel)

      Both = [first_three, last_three]                 ; Create array to hold start and end channel

      ; Make sure that the channel numbers are both not above 16. Set to full range in this case.
      IF (Both[0] GT 16) AND (Both[1] gt 16) THEN BEGIN
        PRINT, CupString, "Both channel indexes are larger than 16. Defaulting to full range."
        Both = [1, 16]
      ENDIF

      ; Make sure that the channel numbers are both not below 1. Set to full range in this case.
      IF (Both[0] LT 1) AND (Both[1] LT 1) THEN BEGIN
        PRINT, CupString, "Both channel indexes are smaller than 1. Defaulting to full range."
        Both = [1, 16]
      ENDIF

      ; Make sure that the first channel is not below 1. Set to 1 if it is.
      ; The upper value is kept.
      IF ((Both[0] GT 16) OR (Both[0] LT 1)) THEN BEGIN
        PRINT, CupString, "The first channel index is not between 1 and 16. It will be set to to 1."
        Both = [1, Both[1]]
      ENDIF

      ; Make sure that the second channel is not above 16. Set to 16 if it is.
      ; The lower value is kept.
      IF ((Both[1] GT 16) OR (Both[1] LT 1)) THEN BEGIN
        PRINT, CupString, "The second channel index is not between 1 and 16. It will be set to to 16."
        Both = [Both[0], 16]
      ENDIF

      ; Check to make sure that the final channel value is not smaller than the first
      ; channel. Set to full range in this case.
      IF (Both[1] - Both[0]) LT 0 THEN BEGIN
        Both = [1, 16]
        PRINT, CupString, "End channel number is smaller than initial channel number. Defaulting to full range."
      ENDIF

      ; Displays cautionary message that only 1 data point is being fit from specified cup
      IF (Both[1] - Both[0]) EQ 0 THEN BEGIN
        Both = [Both[0], Both[1]]
        PRINT, CupString, "Caution - Program will only fit one value from this cup."
      ENDIF

    ENDIF
  ENDIF


  ; Do the calculation for M-Mode (Channel range of 1-128)
  IF Mode EQ 1 THEN BEGIN
    IF j LT 4 THEN BEGIN                 ; If less than four digits cant be two numbers
      IF ((i EQ 0) OR (i EQ 1)) THEN BEGIN
        IF i EQ 0 THEN BEGIN             ; If no fitting (0) then use no range
          Both = [0,0]
          GOTO, SkipErrorStatement2
        ENDIF
        IF i EQ 1 THEN BEGIN             ; If fitting turned on (1) then use full range
          Both = [1,128]
          GOTO, SkipErrorStatement2
        ENDIF
      ENDIF
      PRINT, CupString, "Input is not long enough. Defaulting to full range."
      Both =[1, 128]                     ; Default to full range
      SkipErrorStatement2:
    ENDIF

    IF j GT 6 THEN BEGIN                 ; If greater than 6 digits cant be two numbers
      PRINT, CupString, "Input is too long. Defaulting to full range."
      Both = [1, 128]                     ; Default to full range
    ENDIF

    IF j EQ 4 THEN sixes= STRCOMPRESS(a+a+i)          ; Add zeros onto front for parsing purposed
    IF j EQ 5 THEN sixes= STRCOMPRESS(a+i)            ; Add zero onto front for parsing purposed

    IF (j GT 3) AND (j lt 7) THEN BEGIN

      last_three  = LONG(STRMID(sixes, 2, 3, /REVERSE_OFFSET))         ; Get last three elements(end channel)
      first_three = LONG(STRMID(sixes, 5, 3, /REVERSE_OFFSET))         ; Get first three elements(start channel)

      Both = [first_three, last_three]                 ; Create array to hold start and end channel

      ; Make sure that the channel numbers are both not above 128. Set to full range in this case.
      IF (Both[0] GT 128) AND (Both[1] gt 128) THEN BEGIN
        PRINT, CupString, "Both channel indexes are larger than 128. Defaulting to full range."
        Both = [1, 128]
      ENDIF

      ; Make sure that the channel numbers are both not below 1. Set to full range in this case.
      IF (Both[0] LT 1) AND (Both[1] LT 1) THEN BEGIN
        PRINT, CupString, "Both channel indexes are smaller than 1. Defaulting to full range."
        Both = [1, 128]
      ENDIF

      ; Make sure that the first channel is not below 1. Set to 1 if it is.
      ; The upper value is kept.
      IF ((Both[0] GT 128) OR (Both[0] LT 1)) THEN BEGIN
        PRINT, CupString, "The first channel index is not between 1 and 128. It will be set to to 1."
        Both = [1, Both[1]]
      ENDIF

      ; Make sure that the second channel is not above 128. Set to 128 if it is.
      ; The lower value is kept.
      IF ((Both[1] GT 128) OR (Both[1] LT 1)) THEN BEGIN
        PRINT, CupString, "The second channel index is not between 1 and 128. It will be set to to 128."
        Both = [Both[0], 128]
      ENDIF

      ; Check to make sure that the final channel value is not smaller than the first
      ; channel. Set to full range in this case.
      IF (Both[1] - Both[0]) LT 0 THEN BEGIN
        Both = [1, 128]
        PRINT, CupString, "End channel number is smaller than initial channel number. Defaulting to full range."
      ENDIF

      ; Displays cautionary message that only 1 data point is being fit from specified cup
      IF (Both[1] - Both[0]) EQ 0 THEN BEGIN
        Both = [Both[0], Both[1]]
        PRINT, CupString, "Caution - Program will only fit one value from this cup."
      ENDIF

    ENDIF
  ENDIF
  
  
  ;FILLER CODE FOR ELECTRONS WHILE DEBUGGING
  ; Do the calculation for E-Mode (Channel range of 1-28) ;Both = [1,28]
  IF mode EQ -1 THEN BEGIN
    IF j LT 4 THEN BEGIN                 ; If less than four digits cant be two numbers
      IF ((i EQ 0) OR (i EQ 1)) THEN BEGIN
        IF i EQ 0 THEN BEGIN             ; If no fitting (0) then use no range
          Both = [0,0]
          GOTO, SkipErrorStatement3
        ENDIF
        IF i EQ 1 THEN BEGIN             ; If fitting turned on (1) then use full range
          Both = [1,28]
          GOTO, SkipErrorStatement3
        ENDIF
      ENDIF
      PRINT, CupString, "Input is not long enough. Defaulting to full range."
      Both =[1, 28]                     ; Default to full range
      SkipErrorStatement3:
    ENDIF

    IF j GT 6 THEN BEGIN                 ; If greater than 6 digits cant be two numbers
      PRINT, CupString, "Input is too long. Defaulting to full range."
      Both = [1, 28]                     ; Default to full range
    ENDIF

    IF j EQ 4 THEN sixes= STRCOMPRESS(a+a+i)          ; Add zeros onto front for parsing purposed
    IF j EQ 5 THEN sixes= STRCOMPRESS(a+i)            ; Add zero onto front for parsing purposed

    IF (j GT 3) AND (j lt 7) THEN BEGIN

      last_three  = LONG(STRMID(sixes, 2, 3, /REVERSE_OFFSET))         ; Get last three elements(end channel)
      first_three = LONG(STRMID(sixes, 5, 3, /REVERSE_OFFSET))         ; Get first three elements(start channel)

      Both = [first_three, last_three]                 ; Create array to hold start and end channel

      ; Make sure that the channel numbers are both not above 28. Set to full range in this case.
      IF (Both[0] GT 28) AND (Both[1] gt 16) THEN BEGIN
        PRINT, CupString, "Both channel indexes are larger than 28. Defaulting to full range."
        Both = [1, 28]
      ENDIF

      ; Make sure that the channel numbers are both not below 1. Set to full range in this case.
      IF (Both[0] LT 1) AND (Both[1] LT 1) THEN BEGIN
        PRINT, CupString, "Both channel indexes are smaller than 1. Defaulting to full range."
        Both = [1, 28]
      ENDIF

      ; Make sure that the first channel is not below 1. Set to 1 if it is.
      ; The upper value is kept.
      IF ((Both[0] GT 28) OR (Both[0] LT 1)) THEN BEGIN
        PRINT, CupString, "The first channel index is not between 1 and 28. It will be set to to 1."
        Both = [1, Both[1]]
      ENDIF

      ; Make sure that the second channel is not above 16. Set to 16 if it is.
      ; The lower value is kept.
      IF ((Both[1] GT 28) OR (Both[1] LT 1)) THEN BEGIN
        PRINT, CupString, "The second channel index is not between 1 and 28. It will be set to to 16."
        Both = [Both[0], 28]
      ENDIF

      ; Check to make sure that the final channel value is not smaller than the first
      ; channel. Set to full range in this case.
      IF (Both[1] - Both[0]) LT 0 THEN BEGIN
        Both = [1, 28]
        PRINT, CupString, "End channel number is smaller than initial channel number. Defaulting to full range."
      ENDIF

      ; Displays cautionary message that only 1 data point is being fit from specified cup
      IF (Both[1] - Both[0]) EQ 0 THEN BEGIN
        Both = [Both[0], Both[1]]
        PRINT, CupString, "Caution - Program will only fit one value from this cup."
      ENDIF

    ENDIF
  ENDIF
  
  
  RETURN, Both
END

FUNCTION closest, array, value
  ;+
  ; NAME:
  ; CLOSEST
  ;
  ; PURPOSE:
  ; Find the element of ARRAY that is the closest in value to VALUE
  ;
  ; CATEGORY:
  ; utilities
  ;
  ; CALLING SEQUENCE:
  ; index = CLOSEST(array,value)
  ;
  ; INPUTS:
  ; ARRAY = the array to search
  ; VALUE = the value we want to find the closest approach to
  ;
  ; OUTPUTS:
  ; INDEX = the index into ARRAY which is the element closest to VALUE
  ;
  ;   OPTIONAL PARAMETERS:
  ; none
  ;
  ; COMMON BLOCKS:
  ; none.
  ; SIDE EFFECTS:
  ; none.
  ; MODIFICATION HISTORY:
  ; Written by: Trevor Harris, Physics Dept., University of Adelaide,
  ;   July, 1990.
  ;
  ;-

  IF (n_elements(array) LE 0) OR (n_elements(value) LE 0) THEN index=-1 $
  ELSE IF (n_elements(array) EQ 1) THEN index=0 $
  ELSE BEGIN
    abdiff = abs(array-value)           ; Form absolute difference
    mindiff = min(abdiff,index)         ; Find smallest difference
  ENDELSE

  RETURN, index
END

FUNCTION AREAOV, X, Y, P1, P2

  ; Pg 112 of Alan Barnett's: "The Response Function of the Voyager Plasma Science Experiment"
  ; Constants dervied from geometrical analysis of the instrument

  ; This function computes the sensitive area for the main sensor of
  ; the PLS sensor using the trapezoidal approximation.
  ;
  ; Inputs
  ; X    Flow Velocity (X component)/Particle Velocity (Z Component)
  ; Y    Flow Velocity (Y component)/Particle Velocity (Z Component)
  ; P1   (Threshold Speed/Particle Velocity(Z Component))^2
  ; P2   (Suppressor Speed/Particle Velocity(Z Component))^2
  ;
  ; Outputs
  ; R    Sensitive Area

;  a = 0.338 ; Constants that help to define the sensitive area of the sensor
;  b = 0.197
;  c = 0.093
;  d = 0.372
  sqrt1 = SQRT(1d - P1)

  ; Compute Shift scalar and shift vector

;  Shift = a + 2.0*b*(Sqrt(1.0+P2)-1.0)/P2 + c/sqrt1 + $
;    2.0*d*(1.0-sqrt1)/P1

  Shift = 0.338 + 2.0*0.197*(Sqrt(1.0+P2)-1.0)/P2 + 0.093/sqrt1 + $
          2.0*0.372*(1.0-sqrt1)/P1

  X1 = Shift * X
  Y1 = Shift * Y

  ; Compute Area Overlap

  IF X1 GT 4.94 THEN RETURN, 0d
  R = (4.94-X1)/3.84
  IF X1 LT 1.1 THEN R = 1d

  ; Compute locations of upper corners of trapezoid for area overlap

  Y3 = 0.762 * Cos(1.018*X1+0.247)/(1.0+0.25*X1)
  Y4 = 2.5 - 0.125*(X1-1.0)*(X1-1.0)

  IF (Y1 LE Y3) THEN RETURN, R

  IF (Y1 GT Y4) THEN RETURN, 0d
  R = R*(Y4-Y1)/(Y4-Y3)
  
  IF Y1 LT -2.02 THEN BEGIN
    IF Y1 LT -3.63 THEN RETURN, 0d
    R = R*(3.63+Y1)/1.61
  ENDIF

  RETURN, R
END

FUNCTION CLAROV, R2, RL, RS, D

  ; This function computes the common area of two overlapping circles
  ;
  ; Inputs
  ; R2   RL^(2.0) - RS^(2.0)
  ; RL   Radius of Larger Circle
  ; RS   Radius of Smaller Circle
  ; D    Distance between the centers of the circles
  ;
  ; Outputs
  ; A    Common area of the circles

  IF (D-RL+RS) LE 0 THEN A = !PI*RS*RS
  IF (D-RL+RS) GT 0 THEN BEGIN
    IF (D-RL-RS) GE 0 THEN A = 0.0
    IF (D-RL-RS) LT 0 THEN BEGIN
      QS = (R2 - D*D)/(2.0*D*RS)
      QL = (R2 + D*D)/(2.0*D*RL)
      ;IF qs GT 1 THEN qs = 1d
      ;IF ql GT 1 THEN ql = 1d
      A = RS*RS*(!PI/2.0 + ASIN(QS) + QS*SQRT(1.0-QS*QS))+$
        RL*RL*(!PI/2.0 - ASIN(QL) - QL*SQRT(1.0-QL*QL))
    ENDIF
  ENDIF

  IF finite(A) EQ 0 THEN A = RS*RS*!Pi
  RETURN, A
END

FUNCTION TRANSP, X, Y, P1, P2

  ; This function computes the transparency for the main sensor of the PLS instrument
  ;
  ; Inputs
  ; X    Flow Velocity (X component)/Particle velocity (Z Component)
  ; Y    Flow Velocity (Y component)/Particle Velocity (Z Component)
  ; P1   (Threshold Speed/Particle Velocity(Z Component))^2
  ; P2   (Suppressor Speed/Particle Velocity(Z Component))^2
  ;
  ; Outputs
  ; T    Grid Transparency

  ; Because there are multiple grids between the Plasma environment and the sensor,
  ; the grids will absorb some particles as they can not make it through some of the holes.
  ; Naturally, particles at a higher incident angle will end up being absorbed with a higher
  ; frequency.


  ;C = 0.02380952381 ; Constant useful for grid transparency c = 1/42


  ; Original version of the code

  ;X2 = X*X
  ;Y2 = Y*Y
  ;T = ((1.0 - C*Sqrt(1.0+X2))*(1.0-C*Sqrt(1.0+Y2)))^(5.0)
  ;T = T*((1.0 - C*Sqrt(1.0 + X2/(1.0-P1)))*$
  ;    (1.0-C*Sqrt(1.0 + Y2/(1.0-P1))))^(3.0)
  ;T = T*((1.0 - C*Sqrt(1.0 + X2/(1.0+P2)))*$
  ;    (1.0-C*Sqrt(1.0 + Y2/(1.0+P2))))

  ; Below is the same as the above commented out section but uses less computer time
  ; because of how IDL processes calculations. May be a more optimal method for computing
  ; languages other than IDL.

  X2 = X*X
  Y2 = Y*Y
  T = ((1.0 - 0.02380952381*SQRT(1.0+X2))*(1.0 - 0.02380952381*SQRT(1.0+Y2)))
  T = T*T*T*T*T
  M = ((1.0 - 0.02380952381*SQRT(1.0 + X2/(1.0-P1)))*$
    (1.0 - 0.02380952381*SQRT(1.0 + Y2/(1.0-P1))))
  T = T*M*M*M
  T = T*((1.0 - 0.02380952381*SQRT(1.0 + X2/(1.0+P2)))*$
    (1.0 - 0.02380952381*SQRT(1.0 + Y2/(1.0+P2))))

  RETURN, T
END

FUNCTION TRNSPD, PSI, PHI, P1, P2
  ; This function computes the transparency for the side sensor
  ;
  ; Inputs
  ; PSI  Flow Velocity (Tangent)/Particle Velocity (Z Component)
  ; PHI  Flow Velocity Polar Angle
  ; P1   (Threshold Speed/Particle Velocity(Z Component))^2
  ; P2   (Suppressor Speed/Particle Velocity(Z Component))^2
  ;
  ; Outputs
  ; T    Grid Transparency

  ; C = 0.02381 ; Constants for Grid Transparency c = 1/42

  ; Original Code

  ;A1 = 0
  ;A2 = 1.0821
  ;A3 = 0.8727
  ;A4 = 0.9599
  ;A5 = 0.698
  ;A6 = 0.0698
  ;AM = 1.1868
  ;AS = 0

  ;T = (1.0 - C*Sqrt(1.0+PSI*COS(PHI-A1)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-A1)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*COS(PHI-A2)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-A2)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*COS(PHI-A3)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-A3)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*COS(PHI-A4)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-A4)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*COS(PHI-A5)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-A5)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*COS(PHI-A6)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-A6)^(2.0)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*COS(PHI-AM)^(2.0)/(1.0-P1)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-AM)^(2.0)/(1.0-P1)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*COS(PHI-AS)^(2.0)/(1.0+P2)))*$
  ;    (1.0 - C*Sqrt(1.0+PSI*SIN(PHI-AS)^(2.0)/(1.0+P2)))

  ; A slightly quicker and condensed version of the above code

  ;  T = [ Cos(Phi), Sin(Phi), Cos(Phi-1.0821), Sin(Phi-1.0821), Cos(Phi-0.8727), Sin(Phi-0.8727), Cos(Phi-0.9599), Sin(Phi-0.9599), $
  ;    Cos(Phi-0.698), Sin(Phi-0.698), Cos(Phi-0.0698), Sin(Phi-0.0698), Cos(Phi-1.1868)/Sqrt(1.0-P1), Sin(Phi-1.1868)/Sqrt(1.0-P1), Cos(Phi)/Sqrt(1.0+P2), Sin(Phi)/Sqrt(1.0+P2)]
  ;  T = T*T
  ;  T = 1.0 + Psi*T
  ;  T = Sqrt(T)
  ;  T = 1.0 - 0.02381*T
  ;  T = Product(T)

  ; Below is the same as the above commented out section but uses less computer time
  ; because of how IDL processes calculations. May be a more optimal method for computing
  ; languages other than IDL.
  ;
  ; IDL has a tendency to run the Sin function a little bit faster than the Cos function, thus
  ; cosines are expressed in terms of sines.

  S = sin(Phi + [ 0, -1.0821, -0.8727, -0.9599, -0.698, -0.0698, -1.1868, 0])
  S = S*S
  C = 1-S ; using 1 = cos^2 + sin^2

  T = [C[0], S[0], C[1], S[1], C[2], S[2], C[3], S[3], C[4], S[4], C[5], S[5], C[6], S[6], C[7], S[7]]

  T = T
  T[[12,13]] = T[[12,13]]/(1.0-P1)
  T[[14,15]] = T[[14,15]]/(1.0+P2)


  T = Product(1.0- 0.02381* Sqrt(1.0 + Psi*T)) ; Transparency of the grids

  RETURN, T
END

PRO SHIFT, Psi, P1, P2, DCA = DCA, DCG = DCG, DAG = DAG
  ; This function computes the shift functions for the side sensor of the PLS instrument.
  ;
  ; Inputs
  ; Psi          (Flow Velocity(tang)/Particle Velocity(Z-comp))**2
  ; P1           (Threshold Speed/Particle Velocity(Z-comp))**2
  ; P2           (Suppressor speed/Particle Velocity(Z-comp))**2
  ;
  ; Outputs
  ; DCA          Shift between aperture and collector
  ; DCG          Shift between guard ring and collector
  ; DAG          Shift between aperture and guard ring

  arg = SQRT(1.0+P2)-1.0
  S = 0.495 + 0.720*(1.0 - Sqrt(1.0 - P1))/P1 + $
    0.290*(arg)/P2
  B = S/(0.166 + 0.139*(arg)/P2)
  DCA = Sqrt(Psi)*6.0*S
  DCG = DCA/B
  DAG = DCA - DCG

  RETURN
END

FUNCTION CUPINT, DN, W, U, A, Z, MODE, NCHN, SS, LDIST

  ; This function simulates the response of the Voyager PLS experiment main sensor to a cold
  ; (M>5) convected maxwellian plasma. The common blocks are filled by the program CONRD which must
  ; first be run in IDL.
  ;
  ; This response is more complicated than a simple Maxwellian but less complicated than
  ; the LABCUR program that does a similar calculation with better accuracy.
  ;
  ; Inputs
  ; DN          Number density in number per CC
  ; W           Thermal Speed in (km/sec)
  ; U(3)        Plasma velocity in cup coordinates
  ; A           Atomic Mass of Ion species in AMU
  ; Z           Charge state of ions in proton charges
  ; MODE        Flag; 1 = L-Mode, 2 = M-Mode
  ; NCHN(2)     Lower and Upper limits on channel number
  ; SS          Step Size for numerical integration in threshold speeds. Suggested Value: 0.003-0.018
  ; LDIST       Flag; True = Reduced Distribution Function
  ;                   False = Currents
  ; Outputs
  ; DIST(128)   Reduced Distribution Function or Current

  ; Some of the inputs are located in the common block GENERATE. The remaining inputs are specified depending
  ; on the mode and the parameters that are being fit.

  ; CONSTANTS

  VZSUP  = 134.9d                  ; Proton velocity equivalent to suppressor voltage
  R0     = 100d                    ; Aperture Area
  ECHRGE = 1.6d-19                 ; Proton charge in coulombs
  PI12   = DOUBLE(1.772453851)     ; Square root of pi
  CRR    = MAKE_ARRAY(129,/DOUBLE) ; Array for holding current
  DISTR  = MAKE_ARRAY(128,/DOUBLE) ; Array for current/distribution function at the end
  Z      = Double(Z)               ; Convert to double precision for calculations
  A      = Double(A)
  DN     = DOUBLE(DN)
  W      = DOUBLE(W)
  U      = DOUBLE(U)
  SS     = DOUBLE(SS)
  SqZA   = SQRT(Z/A)               ; Square root of Z/A to calculate it once and save IDL time

  ; Overwrite for L-mode and make a smaller output
  IF MODE EQ 0 THEN DISTR = MAKE_ARRAY(16)

  COMMON VZ                ; Read in the Common Blocks
  COMMON LCOMP

  CONST = Z*R0*ECHRGE*DN/(W*PI12)*1.0d20*(2d/3d) ; A constant for converting to the correct units (femptoamps for current)
  N1 = NCHN[0] ; Lowest Channel Number
  N2 = NCHN[1] ; Highest Channel Number
  N3 = N2 + 1  ; Highest Channel Number + 1 for integration purposes
  UU = U

  ; Begin loop over desired channels
  FOR NCN = n1,n3 DO BEGIN ; Goes to Label 200 in Fortran source code
    ; Select Proper value of NCHAN for desired MODE
    NCHAN = NCN  ; Index the values in NCHAN
    IF (MODE EQ 0) THEN NCHAN = 1 + 8*(NCN-1) ; Computes L-Mode of response
    
    NCHAN1 = NCHAN-1
    ; Initialize the current and choose proper value for velocity threshold

    Cur = 0d
    VZT = DOUBLE(AVZT[NCHAN1])*SqZA ; Read in the lowest velocity measurement from the Common Blocks
    DVZ = DOUBLE(AVZT[NCHAN1]*SS)*SqZA
    VZ = VZT ; Velocity in Z Direction

    ; Evaluate arguments of Maxwellian and exit if too large

    X1 = (VZ - UU[2])/W ; Argument of Exponential in Maxwellian Distribution
    
    N = 0 ; Initialize the loop variable
    Bound = 4d
    WHILE X1 LT Bound DO BEGIN
      VZ = VZ + DVZ ; Increment velocity according to stepsize
      X1 = (VZ - UU[2])/W
      IF (X1 LT -Bound) THEN CONTINUE ; Skip contribution to integration if below bound
      N = N + 1 ; Increment the variable for Simpson's integral weight

      ; Change variables and evaluate response function
      ; Use the function Areaov and transp above

      ; Definitions for P1, P2, X, and Y can be found in the TRANSP and AREAOV functions above
      
      VZ2 = VZ*VZ
      P1 = (VZT*VZT)/VZ2
      P2 = (VZSUP*VZSUP)/VZ2*Z/A
      X = ABS(UU[0]/VZ)
      Y = UU[1]/VZ
      T = TRANSP(X,Y,P1,P2) ; Transparency of the grids
      R = AREAOV(X,Y,P1,P2) ; Area Overlap


      ; Add contribution to sum with proper weight for simpson's rule
      CUR = CUR + R*T*EXP(-X1*X1)*VZ*DOUBLE(N MOD 2 + 1)
    ENDWHILE

    ; Place sum in temporary array
    CRR[NCHAN1] = CUR*DVZ  ; Current for the channel
  ENDFOR ; Loop over the next channel
  CRR = CRR*CONST  ; Multiply by constant
  
  IF (MODE EQ 1) THEN BEGIN

    ; For M-mode, take differences for currents. If Ldist = 'TRUE' divide by channel width for
    ; the reduced distribution function
    FOR NCHAN = n1, n2 DO BEGIN ; Store the currents in the proper output array
      DISTR[NCHAN-1] = CRR[NCHAN-1] - CRR[NCHAN] ; DIST changed to DISTR from Fortran source code since
      ; DIST is an IDL function already
      IF (LDIST EQ 'TRUE') THEN DISTR[NCHAN-1] = DISTR[NCHAN-1]/DVOLT[NCHAN-1] ; Computes the Reduced Distribution Function
      ; by dividing by Voltage width of channel
    ENDFOR
  
    ; Check that all results are finite. If not, set value to zero
    inddddd = WHERE(FINITE(DISTR) NE 1, nfo)
    IF nfo GT 0 THEN DISTR(inddddd) = 0
    RETURN, DISTR
  ENDIF ELSE IF (MODE EQ 0) THEN BEGIN

    FOR NCN = N1, N2 DO BEGIN
      NCHAN = 1 + 8*(NCN-1)
      DISTR[NCN-1] = CRR[NCHAN-1] - CRR[NCHAN+7] ; Computes Current
      IF (LDIST EQ 'TRUE') THEN DISTR[NCN-1] = DISTR[NCN-1]/DVOLTL[NCN-1] ; Computes the Reduced Distribution Function
      ; by dividing by Voltage width of channel
    ENDFOR

    ; Check that all results are finite. If not, set value to zero
    inddddd = WHERE(FINITE(DISTR) NE 1, nfo)
    IF nfo GT 0 THEN DISTR(inddddd) = 0
    RETURN, DISTR
  ENDIF
END

FUNCTION DCPINT, DN, W, U, A, Z, MODE, NCHN, StepS, LDIST

  ; This function simulates the response of the Voyager PLS experiment side sensor to a cold
  ; (M>5) convected maxwellian plasma. The common blocks are filled by the program CONRD which must
  ; first be run.
  ;
  ; This response is more complicated than a simple Maxwellian but less complicated than
  ; the LDCUR program that does a similar calculation with better accuracy.
  ;
  ; Inputs
  ; DN          Number density in number per CC
  ; W           Thermal Speed in (km/sec)
  ; U(3)        Plasma velocity in cup coordinates
  ; A           Atomic Mass of Ion species in AMU
  ; Z           Charge state of ions in proton charges
  ; MODE        Flag; 1 = L-Mode, 2 = M-Mode
  ; NCHN(2)     Lower and Upper limits on channel number
  ; SS          Step Size for numerical integration in threshold speeds. Suggested Value: 0.003-0.018
  ; LDIST       Flag; True = Reduced Distribution Function
  ;                   False = Currents
  ; Outputs
  ; DIST(128)   Reduced Distribution Function or Current
  ;
  ; Some of the inputs are located in the common block GENERATE while others are called on in the argument
  ; of the function based on the parameters that are being fit.
  ; CONSTANTS
  ;

  VZSUP  = 134.9d           ; Proton velocity equivalent to suppressor voltage
  R0     = 82.7d            ; Aperture Area
  ECHRGE = 1.6d-19          ; Proton charge in coulombs
  PI12   = 1.772453851d     ; Square root of pi
  RC     = 6.35d            ; Radius of Collector [cm]
  RG     = 5.1308d          ; Radius of Guard Ring [cm]
  RA     = 5.6438d          ; Radius of Aperture [cm]
  RCG2   = 13.99739136d     ; RC**2 - RG**2
  RCA2   = 8.47002156d      ; RC**2 - RA**2
  RAG2   = 5.5273698d       ; RA**2 - RG**2
  CRR    = Make_array(129,/DOUBLE) ; Array for holding current
  DISTR  = Make_array(128,/DOUBLE) ; Array for current/distribution function at the end
  Z = Double(Z)
  A = Double(A)
  SqZA = Sqrt(Z/A)         ; Square root of charge over mass


  ; Overwrite for L-mode and make a smaller output
  IF MODE EQ 0 THEN DISTR = Make_array(16,/DOUBLE)

  Common VZ                ; Read in the Common Blocks
  Common LCOMP

  CONST = Z*ECHRGE*DN/(W*PI12)*1.0D20*(2d/3d) ; Constant for converting current to femptoamps
  N1 = NCHN[0] ; Lowest Channel Number
  N2 = NCHN[1] ; Highest Channel Number
  N3 = N2 + 1  ; Highest Channel Number + 1 for integration purposes

  ; Begin loop over desired channels

  FOR NCN = N1, N3 DO BEGIN ; Goes to Fortran 600 statement
    ; Select Proper Value of NCHAN for desired mode

    NCHAN = NCN ; Reindex the values
    IF MODE EQ 0 THEN NCHAN = 1 + 8*(NCN-1) ; For L-Mode
    
    NCHAN1 = NCHAN-1
    
    ; Initialize current, then change particle velocity to polar coordinates

    CUR = 0d ; Current is intially 0 and adds current for each iteration
    IF ((U[0] EQ 0d) AND (U[1] EQ 0d)) THEN Phi = 0d ELSE Phi = ATAN(U[1], U[0]) ; Then there is no change in the Phi angle 
    ; into the sensor because the flow is purely into the cup, Z direction
    VT2 = U[0]*U[0] + U[1]*U[1]   ; Magnitude of velocity relative to normal of cup

    ; Choose proper values for threshold velocity and step size
    ; Evaluate argument of Maxwellian and exit if too large
    VZT = DOUBLE(AVZT[NCHAN1])*SqZA ; Read in the initial values of the velocity
    VZ = VZT
    X1 = (VZ - U[2])/W               ; Argument in exponential of Maxwellian distribution
    DVZ = DOUBLE(AVZT[NCHAN1]*StepS)*SqZA ; Calculate delta Vz for integration purposes
    ; Loop for numerical integration
    
    N = 0
    Bound = 4d
    WHILE X1 LT Bound DO BEGIN
      VZ = VZ + DVZ ; Increment velocity according to stepsize
      X1 = (VZ - U[2])/W
      IF (X1 LT -Bound) THEN CONTINUE ; Skip contribution to integration if below bound
      N = N + 1 ; Increment the variable for Simpson's integral weight
      
      VZ2 = VZ*VZ
      P1 = (VZT*VZT)/VZ2
      P2 = (VZSUP*VZSUP)/VZ2*Z/A
      PSI = VT2/VZ2
      
      ; Change variables and evalute response function;
      ; TRNSPD evaluates the grid transparency
      ; SHIFT evaluates the appropriate shift functions
      ; CLAROV calculates the common area of two circles
      ; P1, P2, and PSI are all defined in TRNSPD
      
      T = TRNSPD(PSI,PHI,P1,P2)
      SHIFT, PSI, P1, P2, DCA = DCA, DCG = DCG, DAG = DAG ; Compute the Shift function
      IF DCA EQ 0.0 THEN R = R0                           ; Some overlap differently so there are different labels
      XCG = (RCG2 + DCG*DCG)/(2.0*DCG)                    ; to quickly define the values
      XCA = (RCA2 + DCA*DCA)/(2.0*DCA)
      IF (XCA-XCG) LT 0 THEN BEGIN
        A1 = CLAROV(RCG2,RC,RG,DCG) ; Determine the area overlap of the circular parts of the D-Cup Sensor
        A2 = CLAROV(RAG2,RA,RG,DAG)
        R = A1 + A2 - R0
      ENDIF ELSE IF (XCA-XCG) GE 0 THEN BEGIN
        R = CLAROV(RCA2,RC,RA,DCA)
      ENDIF ELSE IF DCA EQ 0 THEN R = R0
      
      CUR = CUR + R*T*EXP(-X1*X1)*VZ*DOUBLE(N MOD 2 + 1)
    ENDWHILE

    CRR[NCHAN1] = CUR*DVZ
  ENDFOR
  CRR = CRR*CONST
  
  IF MODE EQ 1 THEN BEGIN
    
    ; For M-MODE, take differeence for currents.
    ; Divide by channel width for reduced distribution function
    FOR NCHAN = N1,N2 DO BEGIN
      DISTR[NCHAN-1] = CRR[NCHAN-1] - CRR[NCHAN] ; Distribution is difference in Currents between the two channels
      IF (LDIST EQ 'TRUE') THEN DISTR[NCHAN-1] = DISTR[NCHAN-1]/DVOLT[NCHAN-1] ; Reduced Distribution Function is divided
      ; by the channel voltage width
    ENDFOR
    
    ; Check that all elements are finite. Set non-finite ones to zero
    inddddd = WHERE(FINITE(DISTR) NE 1, nfo)
    IF nfo GT 0 THEN DISTR(inddddd) = 0
    RETURN, DISTR
  ENDIF ELSE IF (MODE EQ 0) THEN BEGIN
    FOR NCN = N1, N2 DO BEGIN
      NCHAN = 1 + 8*(NCN-1)
      DISTR[NCN-1] = CRR[NCHAN-1] - CRR[NCHAN+7] ; Distribution is difference in currents between the two channels
      IF (LDIST EQ 'TRUE') THEN DISTR[NCN-1] = DISTR[NCN-1]/DVOLTL[NCN-1] ; Reduced Distribution Function is divided
      ; by the channel voltage width
      
    ENDFOR

    ; Check that all elements are finite. Set non-finite ones to zero
    inddddd = WHERE(FINITE(DISTR) NE 1, nfo)
    IF nfo GT 0 THEN DISTR(inddddd) = 0
    RETURN, DISTR
  ENDIF
END

FUNCTION SITTLER_RESPONSE, DN, Te, U, A, Z, MODE, NCHN
  ; Uses the response from Sittler's 1983 technical manual for a bi-maxwellian fit (Equation 19)
  ; Bi-maxwellian because of a cold and a hot component (although they argue that hot component might not be maxwellian)
  ; 
  ; Inputs
  ; DN          Number density in particles per cubic centimeter
  ; Te          Electron temperatures
  ; U(3)        Plasma velocity in cup coordinates in km/s
  ; A           1 - Unless running a superposition of electrons
  ; Z           -1 - Charge state of electrons
  ; MODE        Flag; 1 = L-Mode, 2 = M-Mode, -1 = E-Mode (Only programmed for E-Mode at the moment)
  ; NCHN(2)     Lower and Upper limits on channel number
  ;
  ; Outputs
  ; DIST(28)   Reduced Distribution Function or Current
  ;
  ; ------------------------------------
  ; Versions / modifications:
  ; 18-Aug-24: CURRENT_OUTPUT added as output option for when the function is used during PLS_UNCERTAINTIES.pro
  ;            During that procedure the output is being compared with the current (A or fA), instead of the RDF [GX]
  ; 

  ;
  ;STOP
  ; Convert input variables using lenght to meters
  DN = DN*1d6       ; convert cm^3 to m^3
  V = U*1d3         ; convert km/s to m/s
  
  
  ; Constants
  echrge = 1.6d-19                 ; Proton charge in Coulombs
  PI12   = DOUBLE(1.772453851)     ; Square root of pi
  me = 9.11d-31                    ; electron mass in kg
  Aeff = 100d*1d-4                 ; Area of sensor (cm^2), converted to m^2
  we = SQRT(2d*echrge*Te/(A*me))   ; Electron thermal speed in m/s
  A1 = 0.929d                      ; A1, a12, A2, a22: fitting parameters for the R response (eq. 15)
  a12 = 1.11d
  A2 = 0.0813d
  a22 = 19.09d
  TN = 0.56                        ; normal transmision of sensor (see section 2)
  Vperp2 = (V[0]*V[0] + V[1]*V[1]) ; perpendicular electron speed squared
  Wperp2 = Vperp2/(we*we)          ; perpendicular speed to total thermal speed squared
  
  
  ; Create lookup table of potentials (grid edges) in V
  phij = MAKE_ARRAY(29)
  FOR i = 0, 28 DO BEGIN
    IF i LE 16 THEN BEGIN
      phij(i) = 60.*10.^(double(i)/32.) - 50.
    ENDIF ELSE IF i GE 17 THEN BEGIN
      phij(i) = 60.*10.^(double(i-12)/8.) - 50.
    ENDIF
  ENDFOR
  
  ; Equivalent electron speed for each grid edge in m/s
  vj = SQRT(2. * ECHRGE * double(Z) * phij / (double(A)*me)) ; Tried scaling it real quick to see if that works better - 2/7 CHECK LATER
;  Fvj = double(dn)/(double(w) * PI12)  * EXP(- (vj - double(u[2]))^2./ (double(w)^2.))
  ;STOP
  
  ; Create appropriately sized arrays
  dvj = MAKE_ARRAY(28)
  vjavg = dvj
  curr = dvj
  reduced_dist = dvj
  alphaj = dvj ; Split response into an alpha, gamma, and beta component for easy of writing down analytical equations
  gammaj = dvj
  betaj  = dvj
  
  ;STOP
  FOR i = 0, 27 DO BEGIN
    dvj(i) = vj(i+1) - vj(i)
    vjavg(i) = (vj(i+1) + vj(i))/2
    alphaj(i) = ECHRGE * Aeff * TN * !Pi * double(Z) * vjavg(i) * vjavg(i) *vjavg(i) * dvj(i)
    gammaj(i) = A1 * exp(-(a12*Wperp2) / (a12 + (vjavg(i)*vjavg(i)/(we*we)))) / (a12 + (vjavg(i)*vjavg(i)/(we*we))) + $ ; in memo: assumed that we^3 in series denominator is typo and it should be we^2
                A2 * exp(-(a22*Wperp2) / (a22 + (vjavg(i)*vjavg(i)/(we*we)))) / (a22 + (vjavg(i)*vjavg(i)/(we*we)))
    betaj(i) = ( DN/(PI12*PI12*PI12*we*we*we) * EXP(- ((vjavg(i) - V[2])^2.0)/(we*we)) ); V[2] = u_z (u_parallel)
    curr(i) = alphaj(i) * gammaj(i) * betaj(i)
    reduced_dist(i) = gammaj(i) * betaj(i)
  ENDFOR

;  This commented-out part is the part that was used for Jupiter. If Jupiter data are used they need to be used carefully! Check the units at the current 
;  version for the inputs
;  FOR i = 0, 27 DO BEGIN
;    dvj(i) = vj(i+1) - vj(i)
;    vjavg(i) = 0.5d*(vj(i+1) + vj(i))
;    alphaj(i) = ECHRGE * Aeff * 1d-4 * TN * !Pi * double(Z) * vjavg(i) * vjavg(i) *vjavg(i) * dvj(i) ; why 10^-4? Sittler says it's just the area
;    gammaj(i) = A1 * exp(-(a12*Wperp2) / (a12 + (vjavg(i)*vjavg(i)/(we*we)))) / (a12 + vjavg(i)*vjavg(i)/(we*we)) + $ ; Memo: assumed that we^3 in series denominator should be we^2
;                A2 * exp(-(a22*Wperp2) / (a22 + (vjavg(i)*vjavg(i)/(we*we)))) / (a22 + vjavg(i)*vjavg(i)/(we*we))
;    betaj(i) = DN/(PI12*PI12*PI12*we*we*we) * EXP(- ((vjavg(i) - V[2])^2.0)/(we*we)) ; V[2] = u_z ?  DN*1d6: number density in particles per cubic meter
;    curr(i) = alphaj(i) * gammaj(i) * betaj(i)
;    reduced_dist(i) = gammaj(i) * betaj(i)
;  ENDFOR
;  ;p = plot(sqrt(phij), reduced_dist, /ylog)
;
;   RETURN, reduced_dist*1d-12 ; At some point I determined there was a 1d12
;  STOP
 
  ; Outputs
  RETURN, reduced_dist


END

FUNCTION ELECTRON_KAPPA_RESPONSE, DN, Te, U, A, Z, MODE, NCHN, Kappa_val   ; GX: THIS NEEDS TO BE FIXED!
  Kappa = Kappa_val
   
  ; Uses a Kappa response to give the distribution function
  ;
  ; Inputs
  ; DN          Number density in number per CC
  ; Te          Electron temperatures
  ; U(3)        Plasma velocity in cup coordinates in km/s
  ; A           1 - Unless running a superposition of electrons
  ; Z           -1 - Charge state of electrons
  ; MODE        Flag; 1 = L-Mode, 2 = M-Mode, -1 = E-Mode (Only programmed for E-Mode at the moment)
  ; NCHN(2)     Lower and Upper limits on channel number
  ;
  ; Outputs
  ; DIST(28)   Reduced Distribution Function or Current
  ;
  ; Constants
  ECHRGE = 1.6d-19                 ; Proton charge in coulombs
  PI12   = DOUBLE(1.772453851)     ; Square root of pi
  me = 9.11d-31                    ; electron mass
  Aeff = 100.                      ; Area of sensor (cm^2)
  we = SQRT(2d*echrge*Te/(A*me))   ; Thermal speeds
  A1 = 0.929d
  a12 = 1.11d
  A2 = 0.0813d
  a22 = 19.09d
  TN = 0.56
  V = U*1000d
  Vperp2 = (U[0]*U[0] + U[1]*U[1])*1000d*1000d
  Wperp2 = Vperp2/(we*we)

  ; Create lookup table of potentials - doing manually since it's pretty easy
  phij = MAKE_ARRAY(29)
  FOR i = 0, 28 DO BEGIN
    IF i LE 16 THEN BEGIN
      phij(i) = 60.*10.^(double(i)/32.) - 50.
    ENDIF ELSE IF i GE 17 THEN BEGIN
      phij(i) = 60.*10.^(double(i-12)/8.) - 50.
    ENDIF
  ENDFOR

  vj = SQRT(2. * ECHRGE * double(Z) * phij / (double(A)*me)) ; 

  ; Create appropriately sized arrays
  dvj = MAKE_ARRAY(28)
  vjavg = dvj
  curr = dvj
  reduced_dist = dvj
  alphaj = dvj ; Split response into an alpha, gamma, and beta component for ease of writing down analytical equations
  gammaj = dvj
  betaj  = dvj

  FOR i = 0, 27 DO BEGIN
    dvj(i) = vj(i+1) - vj(i)
    vjavg(i) = 0.5d*(vj(i+1) + vj(i))
    alphaj(i) = ECHRGE * Aeff * 1d-4 * TN * !Pi * double(Z) * vjavg(i) * vjavg(i) *vjavg(i) * dvj(i)
    gammaj(i) = A1 * exp(-a12*Wperp2/(a12 + (vjavg(i)*vjavg(i)/(we*we))))/(a12 + vjavg(i)*vjavg(i)/(we*we)) + $ ; has a we cubed term in the memo but next line has a we squared term
      A2 * exp(-a22*Wperp2/(a22 + (vjavg(i)*vjavg(i)/(we*we))))/(a22 + vjavg(i)*vjavg(i)/(we*we))
      
    betaj(i) = DN*1d6/(PI12*PI12*PI12*we*we*we) * GAMMA(DOUBLE(Kappa) + 1d) / GAMMA(DOUBLE(Kappa - 0.5d))/(DOUBLE(Kappa)^(1.5d)) * (1 + ((vjavg(i) - V[2])^2.0)/(DOUBLE(Kappa)*we*we))^(-(DOUBLE(Kappa) + 1))

    curr(i) = alphaj(i) * gammaj(i) * betaj(i)
    reduced_dist(i) = gammaj(i) * betaj(i)
  ENDFOR
  
  ; Logan 2025_06_11 - 1d-12 Error
  ; Can't re-test Jupiter easily as file format has changed. Originally we needed that factor of 1d-12, but I think we removed the /1d-12 on the other end, so we need to remove it here as well
  RETURN, reduced_dist;*1d-12 ; 1d12 factor for cm^6 to m^6 ; Logan 2025_06_11 Test - Removed the 1d12 factor for right now to see if it works
END


FUNCTION COLD_RESPONSE, DN, Te, U, A, Z, MODE, NCHN
  ; This function uses the simplified response from the "FBthesisapp.pdf" on the LASP MOP Voyager website.
  ; Currently only created for the electron modes (MODE = -1). This used to be programmed in an older version of
  ; the code base but was abandoned for CUPINT and DCPINT for the ion analysis (Logan 1/12/18)
  ; 
  ; Inputs
  ; DN          Number density in number per CC
  ; Te          Electron temperatures
  ; U(3)        Plasma velocity in cup coordinates
  ; A           1 - Unless running a superposition of electrons
  ; Z           -1 - Charge state of electrons
  ; MODE        Flag; 1 = L-Mode, 2 = M-Mode, -1 = E-Mode (Only programmed for E-Mode at the moment)
  ; NCHN(2)     Lower and Upper limits on channel number
  ;
  ; Outputs
  ; DIST(128)   Reduced Distribution Function or Current
  ; General Equations
  ; Ij = e Aeff Z vj dvj F(vj)
  ; vj = SQRT((2 e Z Phij)/(A me))
  ; dvj = vj+1 - vj
  ; F(vj) = n/(w sqrt(pi)) exp( (-vj - Vin)^2/ w^2)
  ; Vin = U[2] - the z component at the moment
  ; 
  ; Constants
  ECHRGE = 1.6d-19                 ; Proton charge in coulombs
  PI12   = DOUBLE(1.772453851)     ; Square root of pi
  me = 9.11d-31                    ; electron mass
  Aeff = 100.                      ; Area of sensor
  W = SQRT(2d*echrge*Te/(A*me))/1000d   ; Thermal speeds
  
  ; Create lookup table of potentials - doing manually since it's pretty easy
  phij = MAKE_ARRAY(29)
  FOR i = 0, 28 DO BEGIN
    IF i LE 16 THEN BEGIN
      phij(i) = 60.*10.^(double(i)/32.) - 50.
    ENDIF ELSE IF i GE 17 THEN BEGIN
      phij(i) = 60.*10.^(double(i-12)/8.) - 50.
    ENDIF
  ENDFOR
  
  ; Very simple lookup table of response parameter - Very rough, just sort of checking if it works in some sense
  ThetaTable = [0d,5d,10d,15d,20d,25d,30d,35d,40d,45d,50d,55d,60d]
  Rtable     = [1d,1d,0.95d,0.9d,0.82d,0.73d,0.64d,0.57d,0.42d,0.32d,0.18d,0.08d,0d]
  
  
  vj = SQRT(2. * ECHRGE * double(Z) * phij / (double(A)*me*1840.)) ; Tried scaling it real quick to see if that works better
  Fvj = double(dn)/(double(w) * PI12)  * EXP(- (vj - double(u[2]))^2./ (double(w)^2.))
  dvj = MAKE_ARRAY(28)
  curr = dvj
  
  FOR i = 0, 27 DO BEGIN
    dvj(i) = vj(i+1) - vj(i)
    Curr(i) = ECHRGE * Aeff * double(Z) * vj(i) * dvj(i) * Fvj(i)
  ENDFOR
  STOP 
END

FUNCTION TDERF, ba
  ; Calculates error function of argument ba

  ;  Original Code
  ;  p  = 0.32759
  ;  aa = [0.254829592,-0.284496736, 1.421413741,-1.453152027,1.061405429]
  ;  a  = abs(ba)
  ;  IF (a GT 9.0) THEN GOTO, Label100 ; Changed to 9 to avoid IDL underflow errors
  ;  t     = 1.0/(1.0+p*a)
  ;  ac    = 0.00
  ;  acc   = t*(aa(0)+t*(aa(1)+t*(aa(2)+t*(aa(3)+t*aa(4)))))
  ;  tderf =  1.0 - acc*exp(-a*a)
  ;  IF (ba LT 0.0d0) THEN tderf = -tderf
  ;  ttderf = tderf
  ;  GOTO, LabelEnd
  ;  Label100:
  ;  IF (ba LT 0.d0) THEN GOTO, Label200
  ;  tderf =  1.0
  ;  GOTO, LabelEnd
  ;  Label200:
  ;  tderf = -1.0d0
  ;  LabelEnd:
  ;  RETURN, tderf
  ; Sped up Version
  ; 
  ;  Modified by Logan Dougherty
  
  a  = abs(ba)
  IF (a GT 9.0) THEN GOTO, Label100 ; Changed to 9 to avoid IDL underflow errors
  t     = 1.0/(1.0+0.32759*a)
  ac    = 0.00
  acc   = t*(0.254829592+t*((-0.284496736)+t*(1.421413741+t*((-1.453152027)+t*1.061405429))))
  tderf =  1.0 - acc*exp(-a*a)
  IF (ba LT 0.0d0) THEN tderf = -tderf
  ttderf = tderf
  GOTO, LabelEnd
  Label100:
  IF (ba LT 0.d0) THEN GOTO, Label200
  tderf =  1.0
  GOTO, LabelEnd
  Label200:
  tderf = -1.0d0
  LabelEnd:
  RETURN, tderf
END

FUNCTION TDERFC, ba
  ; Calculates error function of ba and adds 1 to it.
  p  = 0.32759
  aa = [0.254829592,-0.284496736, 1.421413741,-1.453152027,1.061405429]
  a  = abs(ba)
  IF (a GT 9.0d0)  THEN GOTO, Label100
  t      = 1.0/(1.0+p*a)
  ac     = 0.00
  acc    = t*(aa(0)+t*((aa(1)+t*(aa(2)+t*(aa(3)+t*aa(4))))))
  acc   = t*(0.254829592+t*((-0.284496736)+t*(1.421413741+t*((-1.453152027)+t*1.061405429))))
  tderfc = acc*exp(-a*a)
  IF (ba LT 0.0d0) THEN tderfc = 2.0d0 - tderfc
  GOTO, LabelEnd
  Label100:
  IF (ba LT 0.0) THEN GOTO, Label200
  tderfc = 0.0
  GOTO, LabelEnd
  Label200:
  tderfc = 2.0
  LabelEnd:
  RETURN, tderfc
END

FUNCTION DPHI, Z
  ;
  ;    THIS SUBROUTINE EVALUATES A FUNCTION REQUIRED FOR THE EVALUATION
  ;        OF THE INTEGRALS OVER VX AND VY
  ;
  ;  SQRTPI = 1.772453851D0
  ;  IF (ABS(Z) GT 9.D0) THEN GOTO, Label1 ; Changed to 9 to avoid IDL underflow errors
  ;  DPHI = Z*tderf(Z)+EXP(-Z*Z)/SQRTPI
  ;  GOTO, LabelEnd
  ;  Label1:
  ;  DPHI = Z*tderf(Z)
  ;  LabelEnd:
  ;  RETURN, DPHI

  ; Sped up version
  SQRTPI = 1.772453851D0
  IF (ABS(Z) LT 9.D0) THEN $
    DPHI = Z*tderf(Z)+EXP(-Z*Z)/SQRTPI $
  ELSE $
    DPHI = Z*tderf(Z)
  RETURN, DPHI
END

FUNCTION DPSI, Z, SIGMA, ALPHA
  ;
  ;    THIS SUBROUTINE EVALUATES A FUNCTION REQUIRED FOR THE EVALUATION
  ;        OF THE INTEGRAL OVER VY
  ;
  SQRTPI = 1.772453851D0
  SQRTAL = SQRT(ALPHA)
  DPSI   = SQRTPI*tderf(Z)*(-.5D0+SIGMA*SQRTAL*Z) ; eq'n 2.44i (first half)
  IF (ABS(Z) GT 9.D0) THEN GOTO, LabelEnd
  DPSI = DPSI+SIGMA*SQRTAL*EXP(-Z*Z) ; eq'n 2.44i (second half)
  LabelEnd:
  RETURN, DPSI
END

FUNCTION DPHIc, z
  SQRTPI = 1.772453851D0
  IF (ABS(Z) GT 9.0) THEN GOTO, Label1
  DPHIc = -Z*tderfc(Z) + EXP(-Z*Z)/SQRTPI
  GOTO, LabelEnd
  Label1: DPHIc = Z*tderf(Z)
  LabelEnd:
  RETURN, DPHIc
END

FUNCTION PHI, Z
  z = double(z)
  IF Abs(Z) GT 8.7 THEN BEGIN ; Avoids Underflow errors in IDL
    PHI = Abs(Z)
    GOTO, ReturnLabel
  ENDIF
  Phi = Z*erf(Z) + 1.0/Sqrt(!Pi)*Exp(-z*z) ; eq'n 2.41b
  ReturnLabel:
  RETURN, Phi
END

FUNCTION XBARP, NCUP, SIGMA, ALPHA

  ;FORTRAN SOURCE IN COMMENTED SECTION BELOW
  ;      REAL*8 $$$DT(6) /'$$DATE$$','06/20/86','16:14:27',
  ;       | '   xbarP',' FORTRAN','MJS     '/
  ;REAL*8 x1(4),X(4),DZ(4),DE(4),DP(4)
  ;REAL*8 ZD(4)
  ;REAL*8 DEC(4),DPC(4)
  ;real*4 five
  ;INTEGER*4 NCUP
  ;DATA DEL/3.84D0/
  ;data   X1/-4.94D0,-1.10D0,1.10D0,4.94D0/
  ;data   ZD/-1.61D0,-.106D0,.106D0,1.61D0/
  ;data  five / 5.0/

  X    = Make_array(4) ; Initialize arrays
  DZ   = X
  DE   = X
  DP   = X
  DEC  = X
  DPC  = X
  X1   = [-4.94,-1.10,1.10,4.94]
  ZD   = [-1.61,-0.106,0.106,1.61]
  five = 5.0
  DEL  = 3.84

  ;
  ; BEGIN EXECUTION
  ;
  
  ; A, B, C (Cups 1, 2, 3) have different response than D (Cup 4)
  IF (ncup NE 4) THEN X = X1 ELSE X = ZD

  DSIGMA = SIGMA
  ABSDSIGMA = ABS(DSIGMA)
  DALPHA = ALPHA
  SQRTDAL = SQRT(DALPHA)
  DCHI   = DSIGMA*SQRTDAL

  ;
  ; BRANCH FOR |SIGMA| > 4.94 (FOR NCUP < 4)
  ; BRANCH FOR |SIGMA| > 1.61 (FOR NCUP = 4)
  ;

  IF(NCUP NE 4) THEN GOTO, Label22
  IF(ABSDSIGMA GT 1.61D0) THEN GOTO, Label300
  GOTO, Label23
  Label22:   IF (ABSDSIGMA GT 4.94D0) THEN GOTO, Label300
  Label23:

  ;
  ; |SIGMA| <,= 4.94 (FOR NCUP < 4)
  ; |SIGMA| <,= 1.61 (FOR NCUP = 4)
  ;

  FOR I = 0,3 DO BEGIN
    DZ(i) = SQRTDAL*(X(i)-DSIGMA)
    DE(i) = tderf(DZ(i))
    DP(i) = DPHI(DZ(i))
  ENDFOR

  DENOM  = DP(3)-DP(2)+DP(0)-DP(1)
  DPNUM  = DP(3)-DP(2)-DP(0)+DP(1)
  DENUM  = DE(3)-DE(2)-DE(0)+DE(1)
  DINT   = 2.00*DEL*SQRTDAL*DPHI(DCHI)
  DXNUMR = (DCHI*DPNUM-0.50*DENUM+DINT)
  DNUMER = DXNUMR/SQRTDAL
  dxbar  = DNUMER/DENOM
  ADENOM = DENOM
  IF (ABS(DENOM) LE 1.0D-60) THEN ADENOM=0.0
  xbar = dxbar

  IF (abs(xbar) LT five) THEN GOTO, LabelEnd

  xbar = five
  GOTO, LabelEnd

  ;
  ; |SIGMA| > 4.94 (FOR NCUP < 4)
  ; |SIGMA| > 1.61 (FOR NCUP = 4)
  ;

  Label300: DSNSIG = 1.0
  IF (SIGMA LT 0.0) THEN DSNSIG = -1.0
  DCHIC = ABSDSIGMA*SQRTDAL

  DZ(*)  = SQRTDAL*(X(*)-DSIGMA)
  DEC(*) = 0.D0
  FOR iii = 0,3 DO BEGIN
    IF (ABS(DZ(iii)) LE 13.3060) THEN DEC(iii) = tderfC(DZ(iii))
    DPC(iii) = DPHIC(ABS(DZ(iii)))
  ENDFOR

  DENOMC = DPC(3)-DPC(2)+DPC(0)-DPC(1)
  DPNUMC = DPC(3)-DPC(2)-DPC(0)+DPC(1)
  DENUMC = DEC(3)-DEC(2)-DEC(0)+DEC(1)
  DINTC  = 2.00*DEL*SQRTDAL*DPHIC(DCHIC)
  DXNUMC = 1.8*(DCHI*DPNUMC+DINTC)-1.8*0.50*DSNSIG*DENUMC
  IF (ABS(DXNUMC) LE 1.D-69) THEN DXNUMC = 0.0
  IF (ABS(DXNUMC) GT 1.D-69) THEN DXNUMC = 1.D-8*DXNUMC
  DNUMEC = DXNUMC/SQRTDAL
  dxbar  = 4.94
  IF (NCUP EQ 4) THEN dxbar = 1.61
  IF (DENOMC NE 0.D0 AND DNUMEC NE 0.0) THEN dxbar = DNUMEC/DENOMC
  ADENOM=DENOMC
  IF (ABS(DENOMC) LE 1.0D-60) THEN ADENOM = 0.0
  xbar = dxbar
  IF (ABS(xbar) LT five) THEN GOTO, LabelEnd
  xbar = five
  GOTO, LabelEnd

  LabelEnd:
  xbar_denom = [xbar,ADENOM]

  RETURN, xbar_denom
END

FUNCTION LABCUR, DN, W, U1, A, Z, MODE, NCHN, LDIST, CUP

  ; This function is the full response for the main sensor of the PLS instrument. It
  ; is valid for both cold and warm plasma flows. Taken from Alan Barnett's thesis, the
  ; code was originally in Fortran.

  ;     ** ORIGINAL FORTRAN CODE HEADER FOR REFERENCE **
  ;
  ;        THIS SUBROUTINE COMPUTES THE MAIN SENSOR RESPONSE OF THE
  ;     VOYAGER PLS EXPERIMENT OPERATING IN EITHER THE M-MODE OR THE
  ;     L-MODE TO A PLASMA WHOSE DISTRIBUTION FUNCTION CAN BE MODELED AS A
  ;     CONVECTED MAXWELLIAN.  IT USES AN ALGORITHM IN WHICH THE
  ;     INTEGRATIONS OVER THE COMPONENTS OF VELOCITY PERPENDICULAR TO
  ;     THE CUP NORMAL ARE PERFORMED ANALYTICALLY, WHILE THE INTEGRATION
  ;     OVER THE COMPONENT OF VELOCITY ALONG THE CUP NORMAL IS DONE
  ;     NUMERICALLY.  THE ALGORITHM IS DESCRIBED IN DETAIL IN
  ;     ALAN BARNETT'S DOCTORAL THESIS.
  ;        FOR EACH MODULATOR STEP, THE INTEGRATION OVER THE CHANNEL
  ;     WIDTH IS DONE USING SIMPSON'S RULE, WHILE A RECTANGULAR RULE
  ;     IS USED FOR THE INTEGRATION OVER THE "TAIL".   THE CURRENT IN A
  ;     GIVEN CHANNEL IS COMPUTED FROM THE DIFFERENCE BETWEEN THE
  ;     INTEGRATED CURRENT FOR ADJACENT MODULATOR GRID STEPS.  THE VZ
  ;     INTEGRATION IS CUT OFF AT A PREDETERMINED VALUE FOR EACH CHANNEL.
  ;
  ;     THE INPUT VARIABLES ARE:
  ;        DN        DENSITY OF SPECIES  (NUMBER/CC)
  ;        W         THERMAL SPEED (KM/SEC)
  ;        U1        BULK FLOW VELOCITY IN CUP COORDINATES (KM/SEC)
  ;        A         SPECIES MASS (AMU)
  ;        Z         SPECIES CHARGE (PROTON CHARGES)
  ;        MODE      1=L,2=M
  ;        NCHN(2)   FIRST AND LAST CHANNEL TO BE SIMULATED
  ;        LDIST     LOGICAL FLAG; T-RETURNS REDUCED DISTRIBUTION FCN
  ;                                F-RETURNS CURRENT
  ;        CUP       Input string of the cup in order to convert it for xbarp procedure
  ;
  ;    THE OUTPUT VARIABLE IS
  ;        DIST(128)  REDUCED DISTRIBUTION FUNCTION OR CURRENT
  ;
  ;   THE INTEGER CONSTANTS HAVE THE FOLLOWING MEANING
  ;        NSTEP(I)       NUMBER OF STEPS IN THE NUMERICAL INTEGRATION
  ;                  BETWEEN THE THESHOLD VOLTAGES FOR CHANNEL I AND I+1
  ;        NFIRST(I)      NUMBER OF THE FIRST ELEMENT OF TVZ INCLUDED
  ;                  IN THE INTEGRATION OF THE NEAR TAIL FOR CHANNEL I
  ;
  ;     ** END OF FORTRAN HEADER **

  Z = Double(Z); Force charge and mass to be double to avoid calculational error
  A = Double(A)
  SQZA = SQRT(Z/A)

  ; Define arrays to be used in the function

  U     = Make_array(3)
  ZX    = Make_array(4,2)
  ZY    = Make_array(4,2,2)
  AA    = Make_array(2)
  SS    = AA
  CC    = AA
  A1    = AA
  XBAR  = AA
  A2    = AA
  B     = Make_array(2,2)
  D     = B
  AC    = AA
  CRR   = Make_array(129)
  X     = Make_array(4)
  CRNT  = CRR
  Y     = ZX
  DISTR = Make_array(128)

  ; Overwrite for L-mode and make a smaller output
  IF MODE EQ 0 THEN DISTR = Make_array(16)

  ;    THE REAL CONSTANTS HAVE THE FOLLOWING MEANINGS
  ;        X,Y       PARAMETERS FOR TRAPEZOIDAL APPROXIMATION TO SENSITIVE
  ;                       AREA
  ;        R         NOMINAL APERTURE AREA
  ;        T         NOMINAL GRID TRANSPARENCY AT NORMAL INCIDENCE
  ;        ECHRGE    PROTON CHARGE  IN COULOMBS

  COMMON TRNPAR
  COMMON VZ
  COMMON TAIL
  COMMON LCOMP


  ;NSTEP = [22,16,10,6,41*4,84*2] ; Fortran version - Defined equivalently below in IDL
  NSTEP         = Make_array(129)
  NSTEP(0)      = 22
  NSTEP(1)      = 16
  NSTEP(2)      = 10
  NSTEP(3)      = 6
  NSTEP(4:44)   = 4
  NSTEP(45:128) = 2
  X             = [-4.94,-1.10,1.10,4.94]
  Y             = [[-3.63,-2.02,0.0,0.0],[-3.63,-2.02,0.0,0.0]]
  R             = 102.0
  T             = 0.64807
  ECHRGE        = 1.6022D-19
  SQRTPI        = 1.77245385


  ;    COMMON BLOCKS TRNPAR, VZ, TAIL, AND LCOMP ARE FILLED BY SUBROUTINE
  ;         CONRD, WHICH MUST BE CALLED BEFORE THE FIRST CALL TO LABCUR.
  ;         THE VARIABLES HAVE THE FOLLOWING MEANINGS
  ;        AMA1                      CONSTANTS
  ;        AMA2                 WHICH DESCRIBE THE CUP
  ;        AMC1                     TRANSPARENCY
  ;        AMSFT               THE SHIFT FUNCTION
  ;     THE FIRST INDEX OF THE PRECEDING FOUR VARIABLES REFERS TO THE
  ;          STEP NUMBER IN THE VZ INTEGRATION;  THE SECOND INDEX
  ;          REFERS TO THE MODULATOR STEP NUMBER.
  ;        AAAVZ               VALUES OF VZ FOR INTEGRATION ACROSS
  ;                                      CHANNELS 1-5
  ;        AAVZ                VALUES OF VZ FOR INTEGRATION ACROSS
  ;                                      CHANNELS 6-45
  ;        AVZ                 VALUES OF VZ FOR INTEGRATION ACROSS
  ;                                      CHANNELS 46-129
  ;        AAADVZ              WEIGHTS FOR SIMPSON INTEGRATION ACROSS
  ;                                      CHANNELS 1-5
  ;        AADVZ              WEIGHTS FOR SIMPSON INTEGRATION ACROSS
  ;                                      CHANNELS 6-45
  ;        ADVZ               WEIGHTS FOR SIMPSON INTEGRATION ACROSS
  ;                                      CHANNELS 46-129
  ;        DVOLT               VOLTAGE WIDTHS OF THE M-MODE CHANNELS
  ;        DVOLTL              VOLTAGE WIDTHS OF THE L-MODE CHANNELS
  ;        NM(I)               NUMBER OF INTEGRATION STEPS INCLUDED
  ;                      WHEN STEP I IS UPPER THRESHOLD  (M-MODE)
  ;        NL(I)               NUMBER OF INTEGRATION STEPS INCLUDED
  ;                      WHEN STEP I IS UPPER THRESHOLD  (L-MODE)
  ;        TVZ                 TABLE OF VALUES OF VZ FOR INTEGRATION OF
  ;                                      THE TAIL
  ;        TDVZ                 TABLE OF WEIGHTS FOR NUMERICAL INTEGRATION
  ;                                      OF THE TAIL

  ;     DEFINE CONST AND COMPUTE NORMALIZED BULK VELOCITY

  If (CUP EQ 'A') THEN ncup = 1 ; Defines the cup as a numerical value for use in XBARP
  If (CUP EQ 'B') THEN ncup = 2
  If (CUP EQ 'C') THEN ncup = 3
  If (CUP EQ 'D') THEN ncup = 4

  CONST=R*T*Z*ECHRGE*1.0E20/SQRTPI*DN
  NC1=NCHN(0)
  NC2=NCHN(1)
  NC3=NC2+1
  U(*) = U1(*)/W

  ;
  ;    ENTER OUTER LOOP OVER MODULATOR STEPS
  ;

  FOR NCN = NC1,NC3 DO BEGIN

    NCHAN=NCN
    IF(MODE EQ 0) THEN NCHAN=1+8*(NCN-1)

    ;
    ;    INITIALIZE CURENT
    ;

    CURENT=0.0

    FOR JJ = 0,99 DO BEGIN


      ;
      ;   LAST INTEGRATION POINT NOT NEEDED; SKIP LOOP IF APPROPRIATE
      ;

      IF ( ncn ne nc3) THEN GOTO, Label8
      IF (Mode EQ 0) THEN GOTO, Label4
      IF (Mode EQ 1) THEN GOTO, Label5
      GOTO, Label8
      Label4:  IF (jj+1 gt nl(ncn-1)) THEN GOTO, Label2599
      GOTO, Label8
      Label5:  IF (jj+1 gt nm(ncn-1)) THEN GOTO, Label2599
      Label8:

      ; Condensed form of Logical statement included for reference in FORTRAN language.
      ;      IF (NCN .EQ. NC3 .AND. JJ .GT. NL(NCN) .AND. MODE .EQ. 1)
      ;     +   GO TO 2600
      ;      IF (NCN .EQ. NC3 .AND. JJ .GT. NM(NCN) .AND. MODE .EQ. 2)
      ;     +   GO TO 2600
      ;
      ;    INITIALIZE G AND CHOOSE PROPER VALUE OF VZ, DVZ
      ;

      G=0.0
      IF(NCHAN GE 6 AND NCHAN LT 46) THEN GOTO, Label100
      IF(NCHAN GE 46) THEN GOTO, Label200
      IF(JJ+1 GE NSTEP(NCHAN-1)) THEN GOTO, Label10
      VZ=AAAVZ(NCHAN-1,JJ)*SQZA/W
      DVZ=AAADVZ(NCHAN-1,JJ)*SQZA
      GOTO, Label1000
      Label10: N=JJ-NSTEP(NCHAN-1)+NFIRST(NCHAN-1)-1
      VZ=TVZ(N)*SQZA/W
      DVZ=TDVZ(N)*SQZA
      GOTO, Label1000
      Label100: IF(JJ GT 3) THEN GOTO, Label110
      VZ=AAVZ(NCHAN-6,JJ)*SQZA/W
      DVZ=AADVZ(NCHAN-6,JJ)*SQZA
      GOTO, Label1000
      Label110: N=JJ-5+NFIRST(NCHAN-1);-1
      VZ=TVZ(N)*SQZA/W
      DVZ=TDVZ(N)*SQZA
      GOTO, Label1000
      Label200: IF(JJ GT 1) THEN GOTO, Label210
      VZ=AVZ(NCHAN-46,JJ)*SQZA/W
      DVZ=ADVZ(NCHAN-46,JJ)*SQZA
      GOTO, Label1000
      Label210: N=JJ-3+NFIRST(NCHAN-1)
      VZ=TVZ(N)*SQZA/W
      DVZ=TDVZ(N)*SQZA
      Label1000: X1=(VZ-U(2))^2.0

      ;
      ;    TEST FOR CUTOFF DUE TO MAXWELLIAN DEPENDENCE ON VZ
      ;

      IF (X1 GT 16.0) THEN GOTO, Label2590 ; If argument of Maxwellian is too large skip numerical integration

      ;
      ;    ASSIGN VALUES OF RESPONSE PARAMETERS
      ;

      S=AMSFT(NCHAN-1,JJ)
      CC(0)=AMC1(NCHAN-1,JJ)
      CC(1)=1.0-CC(0)
      AC(0)=AMA1(NCHAN-1,JJ)
      AC(1)=AMA2(NCHAN-1,JJ)
      AA(0)=(AC(0)+VZ^2.0)/S^2.0
      AA(1)=(AC(1)+VZ^2.0)/S^2.0
      SS(0)=S*VZ*U(0)/(VZ*VZ+AC(0))
      SS(1)=S*VZ*U(0)/(VZ*VZ+AC(1))

      ;
      ;    CALCULATE INTEGRAND FOR NUMERICAL INTEGRATION OVER VZ;
      ;        FIRST, CHANGE VARIABLES FOR ANALYTIC EVALUATION OF VX INTEGRAL
      ;

      For L = 1,4 DO BEGIN
        ZX(L-1,0)=SQRT(AA(0))*(X(L-1)-SS(0))
        ZX(L-1,1)=SQRT(AA(1))*(X(L-1)-SS(1))
      ENDFOR

      ;
      ;    COMPUTE XBAR FOR USE IN SADDLE POINT METHOD INTEGRATION OVZER VY
      ;

      ; IDL functions only return one variable. Therefore XBAR and DENOM are combined into one
      ; variable in the function and then uncompressed after XBARP is called.

      xbar1_denom1 = XBARP(NCUP,SS(0),AA(0))
      XBAR(0) = xbar1_denom1(0)
      DENOM1  = xbar1_denom1(1)
      xbar2_denom2 = XBARP(NCUP,SS(1),AA(1))
      XBAR(1) = xbar2_denom2(0)
      DENOM2  = xbar2_denom2(1)
      IF ((ABS(XBAR(1))-1.1) LE 0.0) THEN GOTO, Label2021
      IF ((ABS(XBAR(1))-1.1) GT 0.0) THEN GOTO, Label2022

      ;
      ;    COMPUTE A2
      ;     TO IMPROVE ACCURACY OF TRAPEZOIDAL APPROXIMATION
      ;

      Label2021: A2(0)=1.0E0
      GOTO, Label2023
      Label2022: A2(0)=1.257-0.0630*ABS(XBAR(0))-0.126*SQRT(XBAR(0)^2-5.10*ABS(XBAR(0))+6.612)
      Label2023:
      IF ((ABS(XBAR(1))-1.1) LE 0.0) THEN GOTO, Label2026
      IF ((ABS(XBAR(1))-1.1) GT 0.0) THEN GOTO, Label2027
      Label2026: A2(1)=1.0
      GOTO, Label2028
      Label2027: A2(1)=1.257-0.063*ABS(XBAR(1))-0.126*SQRT(XBAR(1)^2-5.10*ABS(XBAR(1))+6.612)
      Label2028:

      ;
      ;    CHANGE VARIABLES FOR ANALYTIC EVALUATION OF VY INTEGRAL
      ;

      Y(2,0)=.762*COS(1.018*XBAR(1)+.247)/(1.0+.25*XBAR(0))
      Y(2,1)=.762*COS(1.018*XBAR(1)+.247)/(1.0+.25*XBAR(1))
      Y(3,0)=2.5-0.125*(XBAR(0)-1.0)^2.0
      Y(3,1)=2.5-0.125*(XBAR(1)-1.0)^2.0
      SS(0)=S*VZ*U(1)/(VZ*VZ+AC(0))
      SS(1)=S*VZ*U(1)/(VZ*VZ+AC(1))
      FOR L = 1,4 DO BEGIN
        ZY(L-1,0,0)=SQRT(AA(0))*(Y(L-1,0)-SS(0))
        ZY(L-1,1,0)=SQRT(AA(0))*(Y(L-1,1)-SS(0))
        ZY(L-1,0,1)=SQRT(AA(1))*(Y(L-1,0)-SS(1))
        ZY(L-1,1,1)=SQRT(AA(1))*(Y(L-1,1)-SS(1))
      ENDFOR

      ;
      ;    EVALUATE INTEGRATION OVER VX
      ;

      A1(0)=DENOM1/(ZX(3,0)-ZX(2,0))
      A1(1)=DENOM2/(ZX(3,1)-ZX(2,1))

      ;
      ;    SUM OVER TERMS IN GRID TRANSPARENCY APPROXIMATION
      ;

      FOR I = 1,2 DO BEGIN
        FOR J = 1,2 DO BEGIN

          ;
          ;    CALCULATE RESULT OF ANALYTICAL INTEGRATION OVER VY
          ;

          B(I-1,J-1)=(PHI(ZY(3,I-1,J-1))-PHI(ZY(2,I-1,J-1)))/(ZY(3,I-1,J-1)-ZY(2,I-1,J-1))-$
            (PHI(ZY(0,I-1,J-1))-PHI(ZY(1,I-1,J-1)))/(ZY(0,I-1,J-1)-ZY(1,I-1,J-1))
          D(I-1,J-1)=-U(0)^2.0*AC(I-1)/(VZ^2.0+AC(I-1))-U(1)^2.0*AC(J-1)/(VZ^2.0+AC(J-1))
          IF (D(I-1,J-1) LT -30.0) THEN GOTO, Label2050
          D(I-1,J-1)=CC(I-1)*CC(J-1)*EXP(D(I-1,J-1))/SQRT(AA(I-1)*AA(J-1))
          GOTO, Label2051
          Label2050: D(I-1,J-1)=0.0
          Label2051:

          ;
          ;    ADD TO SUM FOR INTEGRAND (EXCEPT MAXWELLIAN FACTOR)
          ;
          IF(ABS(A1(I-1)*B(I-1,J-1)) LE 1.0E-30) THEN GOTO, Label169
          GOTO, Label369
          Label169:  G1=0.0
          GOTO,  Label269
          Label369:  G1=(VZ/S)^2.0*D(I-1,J-1)*A1(I-1)*A2(I-1)*B(I-1,J-1)/4.0
          Label269:  G=G1+G
        ENDFOR
      ENDFOR

      ;
      ;    ADD TO SUM OF CONTRIBUTIONS FOR NUMERICAL INTEGRATION OVER VZ
      ;              (INCLUDING MAXWELLIAN FACTOR)
      ;

      CURENT=CURENT+VZ*EXP(-X1)*G*DVZ

      ;
      ;    FILL ARRAY FOR UPPER MODULATOR STEP INTEGRATION, IF APPROPRIATE
      ;

      Label2590: IF ((mode - 1) LT 0) THEN GOTO, Label2591
      IF ((mode - 1) EQ 0) THEN GOTO, Label2592
      IF ((mode - 1) GT 0) THEN GOTO, Label2595
      Label2591: IF  (jj+1 EQ nl(ncn-1)) THEN crr(ncn-1) = curent*const
      GOTO, Label2595
      Label2592: IF ( jj+1 EQ nm(ncn-1)) THEN crr(ncn-1) = curent*const
      Label2595:

      ; Condensed logical statements from above included for reference in FORTRAN language.
      ; 2590 IF (JJ .EQ. NL(NCN) .AND. MODE .EQ. 1) CRR(NCN)=CURENT*CONST ; FORTRAN version
      ;      IF (JJ .EQ. NM(NCN) .AND. MODE .EQ. 2) CRR(NCN)=CURENT*CONST

      Label2599:
      Label2600:
    ENDFOR

    ;
    ;    FILL ARRAY FOR LOWER MODULATOR STEP INTEGRATION
    ;

    CRNT(NCN-1)=CURENT*CONST

    IF (abs(curent) lt 10.0E10) THEN GOTO, Label5000

    Label5000:
  ENDFOR

  ;
  ;    COMPUTE CURRENTS BY TAKING DIFFERENCE BETWEEN ADJACENT MODULATOR
  ;        STEPS; IF LDIST= .TRUE., DIVIDE BY CHANNEL WIDTH
  ;

  FOR JJ = NC1, NC2 DO BEGIN
    IF(MODE EQ 1) THEN GOTO, Label5100

    ;
    ;    L-MODE, LDIST =  .TRUE.
    ;

    IF(LDIST EQ 'TRUE') THEN DISTR(JJ-1)=(CRNT(JJ-1)-CRR(JJ))/(DVOLTL(JJ-1))
    GOTO, Label5110

    ;
    ;    M-MODE, LDIST =  .TRUE.
    ;

    Label5100: IF(LDIST EQ 'TRUE') THEN DISTR(JJ-1)=(CRNT(JJ-1)-CRR(JJ))/(DVOLT(JJ-1))

    ;
    ;    L- OR M-MODE, LDIST =  .FALSE.
    ;

    Label5110: IF(LDIST NE 'TRUE') THEN DISTR(JJ-1)=(CRNT(JJ-1)-CRR(JJ))
  ENDFOR

  inddddd = where(finite(DISTR) NE 1, nfo)
  IF nfo GT 0 THEN DISTR(inddddd) = 0
  RETURN, DISTR
END

FUNCTION LDCUR, DN, W, U1, A, Z, MODE, NCHN, LDIST, CUP

  ; This function is the full response for the main sensor of the PLS instrument. It
  ; is valid for both cold and warm plasma flows. Taken from Alan Barnett's thesis, the
  ; code was originally in Fortran.
  ;
  ;     ** ORIGINAL FORTRAN CODE HEADER FOR REFERENCE **
  ;
  ;     ** MODIFIED TO USE SINGLE PRECISION 5/22/84 JDR ***
  ;
  ;        THIS SUBROUTINE COMPUTES THE RESPONSE OF THE D-CUP.  THE INPUT
  ;     VARIABLES ARE
  ;             DN             NUMBER DENSITY IN NUMBER PER CC
  ;             W              THERMAL SPEED IN KM/SEC
  ;             U              PLASMA BULK VELOCITY IN CUP COORDINATES
  ;             A              ATOMIC MASS OF IONS IN AMU
  ;             Z              CHARGE STATE OF IONS PROTON CHARGES
  ;             NCHN           LIMITS ON THE CHANNEL NUMBER
  ;        THE OUTPUT VARIABLE IS
  ;             DISTR          THE REDUCED DISTRIBUTION FUNCTION
  ;
  ;      REAL*8 $$$DT(6) /'$$DATE$$','01/18/85','15
  ;     | '   LDCUR',' FORTRAN','MJ
  ;      LOGICAL*1 LDIST
  ;      LOGICAL   LDIST
  ;     INTEGER*4 NSTEP(129),MODE,NCN,N
  ;     INTEGER*4 NC1,NC2,NC3,NCHN(2),NCHAN,JJ,J,NM,NFIRST,NL,I
  ;
  ;   THE INTEGER CONSTANTS HAVE THE FOLLOWING MEANING
  ;        NSTEP(I)       NUMBER OF STEPS IN THE NUMERICAL INTEGRATION
  ;                  BETWEEN THE THESHOLD VOLTAGES FOR CHANNEL I AND I+1
  ;        NFIRST(I)      NUMBER OF THE FIRST ELEMENT OF TVZ INCLUDED
  ;                  IN THE INTEGRATION OF THE NEAR TAIL FOR CHANNEL I
  ;     REAL*4 F(3)/.6488,.1475,.07844/,D(3)/.146963,.904687,.993027/
  ;
  ;     ** NOTE 5/26/2015 LPD **
  ;
  ; FORTRAN CODE defines the above variables for F and D. However below they are defined with a conflicting
  ; definition. It appears that the one below is actually called in the FORTRAN code so that is what is used.
  ;
  ; Here is the definition of the data with the different F and D included
  ;      data      NSTEP/22,16,10,6,41*4,84*2/
  ;      data   F/.6876,.2044,.09153/,D/-.22131,.851562,.990971/
  ;      data   A0/82.7/,T0/.68007/
  ;      data   Z0/.12/,Z1/1.82/,Z10/1.70/
  ;      data   ECHRGE/1.6022E-19/,PI12/1.7725/
  ;
  ;     ** END OF FORTRAN HEADER **


  F = [0.6488,0.1475,0.7844]
  D = [0.146963,0.904687,0.993027]
  ;NSTEP = [22,16,10,6,41*4,84*2] ; Fortran version - Defined equivalently below in IDL
  NSTEP         = Make_Array(129)
  NSTEP(0)      = 22
  NSTEP(1)      = 16
  NSTEP(2)      = 10
  NSTEP(3)      = 6
  NSTEP(4:44)   = 4
  NSTEP(45:128) = 2
  F = [0.6876,0.2044,0.09153]
  D = [-0.22131,0.851562,0.990971]
  A0 = [82.7]
  T0 = [0.68007]
  Z0 = [0.12]
  Z1 = [1.82]
  Z10 = [1.70]
  ECHRGE = [1.6022E-19]
  PI12 = [1.7725]
  U = U1
  CRNT = Make_array(129)
  CRR = CRNT
  DISTR = Make_array(128)

  Z = Double(Z)
  A = Double(A)
  SQZA = SQRT(Z/A)


  ;     REAL*4 F(3), d(3)
  ;     REAL*4 DN,W,U(3),A,Z,DIST(128),X0,X,CONST,A0,T0
  ;     REAL*4 DVZ,G,ZZ,SS,Z0,Z1,Z10
  ;     REAL*4 ECHRGE,PI12,VT2,X2
  ;     REAL*4 GAMMA,KAPPA,CURENT,VT,G12,G32,KAPPA2,VZ,CRR(129),CRNT(129)
  ;
  ;    THE VELOCITY SPACE INTEGRATION IS DONE IN CYLINDRICAL COORDINATES,
  ;        WITH THE CUP NORMAL POINTING IN THE Z- DIRECTION.  THE
  ;        INTEGRATION OVER ANGLE IN THE VX-VY PLANE RESULTS IN THE
  ;        BESSEL FUNCTION I0; IN ORDER TO PERFORM THE NEXT INTEGRATION
  ;        ANALYTICALLY, IT IS NECESSARY TO APPROXIMATE THE BESSEL
  ;        FUNCTION BY A SUM OF EXPONENTIALS.  D AND F ARE THE PARAMETERS
  ;        IN THIS APPROXIMATION. THE CARD WHICH IS COMMENTED OUT CONTAINS
  ;        ALTERNATIVE VALUES FOR THESE PARAMETERS.
  ;
  ;        SOME OTHER REAL VARIABLES ARE:
  ;             A0        NOMINAL COLLECTING AREA OF DETECTOR
  ;             T0        NOMINAL GRID TRANSPARENCY AT NORMAL INCIDENCE
  ;             ECHRGE    PROTON CHARGE IN COULOMBS
  ;             PI12           SQUARE ROOT OF PI
  ;             PRESET         VALUE FOR MAXWELLIAN CUTOFF FLAG
  ;             C1,C2,C3       PARAMETERS FOR TRAPEZOIDAL APPROXIMATION TO
  ;                                SENSITIVE AREA

  ; Overwrite for L-mode and make a smaller output
  IF MODE EQ 0 THEN DISTR = Make_array(16)

  COMMON TRNPRD
  aa1 = Make_array(129,100)
  aa2 = aa1
  aa1(*,*) = AA(*,*,0)
  aa2(*,*) = AA(*,*,1)
  COMMON VZ
  COMMON TAIL
  COMMON LCOMP



  ;
  ;   DEFINE CONST AND LIMITS FOR LOOP OVER MODULATOR STEPS
  ;

  CONST=DN*2.0*Z*ECHRGE*A0*T0/W^3.0/PI12*1.0E20
  NC1=NCHN(0)
  NC2=NCHN(1)
  NC3=NCHN(1)+1
  VT2=U(0)^2+U(1)^2.0
  VT=SQRT(VT2)

  ;
  ;   ENTER OUTER LOOP OVER MODULATOR STEPS
  ;
  
  FOR NCN = NC1, NC3 DO BEGIN
    NCHAN=NCN
    IF(MODE EQ 0) THEN NCHAN=1+8*(NCN-1)

    ;
    ;   INITIALIZE CURRENT
    ;

    CURENT=0.0
    FOR JJ = 0,99 DO BEGIN

      ;
      ;   LAST INTEGRATION LOOP NOT NECCESSAY; SKIP LOOP IF APPROPRIATE
      ;

      IF ( ncn ne nc3) THEN GOTO, Label8
      IF (Mode EQ 0) THEN GOTO, Label4
      IF (Mode EQ 1) THEN GOTO, Label5

      Label4:
      IF (jj+1 GT nl(ncn-1)) THEN GOTO, Label2600
      GOTO, Label8
      Label5:
      IF (jj+1 GT nm(ncn-1)) THEN GOTO, Label2600
      Label8:

      ; Logical statements from above condensed and included for reference in FORTRAN language
      ;      IF (NCN .EQ. NC3 .AND. JJ .GT. NL(NCN) .AND. MODE .EQ. 1)
      ;     +   GO TO 2600
      ;      IF (NCN .EQ. NC3 .AND. JJ .GT. NM(NCN) .AND. MODE .EQ. 2)
      ;     +   GO TO 2600
      ;
      ;   INTIALIZE G AND CHOOSE PROPER VALUE OF VZ,DVZ
      ;


      G=0.0
      IF(NCHAN GE 6 AND NCHAN LT 46) THEN GOTO, Label100
      IF(NCHAN GE 46) THEN GOTO, Label200
      IF(JJ+1 GT NSTEP(NCHAN-1)) THEN GOTO, Label10
      VZ=AAAVZ(NCHAN-1,JJ)*SQZA
      DVZ=AAADVZ(NCHAN-1,JJ)*SQZA
      GOTO, Label1000
      Label10:
      N=JJ-NSTEP(NCHAN-1)+NFIRST(NCHAN-1)-1
      VZ=TVZ(N)*SQZA
      DVZ=TDVZ(N)*SQZA
      GOTO, Label1000
      Label100:
      IF(JJ GT 3) THEN GOTO, Label110
      VZ=AAVZ(NCHAN-6,JJ)*SQZA
      DVZ=AADVZ(NCHAN-6,JJ)*SQZA
      GOTO, Label1000
      Label110:
      N=JJ-5+NFIRST(NCHAN-1)
      VZ=TVZ(N)*SQZA
      DVZ=TDVZ(N)*SQZA
      GOTO, Label1000
      Label200:
      IF(JJ GT 1) THEN GOTO, Label210
      VZ=AVZ(NCHAN-46,JJ)*SQZA
      DVZ=ADVZ(NCHAN-46,JJ)*SQZA
      GOTO, Label1000
      Label210:
      N=JJ-3+NFIRST(NCHAN-1)
      VZ=TVZ(N)*SQZA
      DVZ=TDVZ(N)*SQZA
      Label1000:

      X0=((VZ-U(2))/W)^2.0
      IF(X0 GT 16.0) THEN GOTO, Label2650 ; Skips integration if argument of Maxwellian is large
      X0=X0+VT2/W^2.0


      ;
      ;    EVALUATE INTEGRATION OVER VT ANALYTICALLY;
      ;       SUM OVER TERMS IN GRID TRANSPARENCY APPROXIMATION
      ;

      FOR I = 0,1 DO BEGIN
        GAMMA=((VZ/W)^2.0+(AA(NCHAN-1,JJ,I)))/S(NCHAN-1,JJ)^2.0
        G12=SQRT(GAMMA)

        ;
        ;    SUM OVER TERMS IN APPROXIMATION TO BESSEL FUNCTION
        ;

        FOR J = 0,2 DO BEGIN
          KAPPA=(D(J)/S(NCHAN-1,JJ))*VZ*VT/W^2.0/G12
          KAPPA2=KAPPA^2.0
          X=X0-KAPPA2
          ;IF(X GT 10) THEN GOTO, Label2500 ; Original from Source Code
          IF(X GT 100) THEN GOTO, Label2500 ; Changed threshold to 100 to allow it to run properly for Voyager 1 data
          G32=G12^3.0
          ZZ=2.0*G32*Z10
          SS=(PI12*(KAPPA2-G12*KAPPA*Z0+0.5)*tderf(G12*Z0-KAPPA)-$
            PI12*(KAPPA2-G12*KAPPA*Z1+0.5)*tderf(G12*Z1-KAPPA))/ZZ+$
            PI12*KAPPA*tderf(KAPPA)/2.0/GAMMA
          X2=(G12*Z0-KAPPA)^2.0
          IF(X2 LT 10)THEN SS=SS-KAPPA*EXP(-X2)/ZZ
          X2=(G12*Z1-KAPPA)^2.0
          IF(X2 LT 10) THEN SS=SS+KAPPA*EXP(-X2)/ZZ
          IF(KAPPA2 LT 10) THEN SS=SS+EXP(-KAPPA2)/2.0/GAMMA

          ;
          ;   ADD TO SUM FOR INTEGRAND (MAXWELLIAN FACTOR NOT INCLUDED)
          ;

          G=G+SS*(F(J)*CC(NCHAN-1,JJ,I))*EXP(-X)
        ENDFOR
      ENDFOR
      Label2500:

      ;
      ;   ADD TO SUM FOR NUMERICAL INTEGRATION OVER VZ
      ;

      CURENT=CURENT+G*(VZ/S(NCHAN-1,JJ))^2.0*DVZ*VZ

      ;
      ;   FILL ARRAY FOR UPPER MODULATOR STEP INTEGRATION, IF APPROPRIATE
      ;
      Label2650:
      IF ((mode - 1) LT 0) THEN GOTO, Label2591
      IF ((mode - 1) EQ 0) THEN GOTO, Label2592
      IF ((mode - 1) GT 0) THEN GOTO, Label2595
      Label2591:   IF  (jj+1 EQ nl(ncn-1)) THEN crr(ncn-1) = curent*const
      GOTO, Label2595
      Label2592:   IF ( jj+1 EQ nm(ncn-1)) THEN crr(ncn-1) = curent*const
      Label2595:

      ; Condensed Logic statements from above in FORTRAN language
      ;      IF (JJ .EQ. NL(NCN) .AND. MODE .EQ. 1) CRR(NCN)=CURENT*CONST
      ;      IF (JJ .EQ. NM(NCN) .AND. MODE .EQ. 2) CRR(NCN)=CURENT*CONST

    ENDFOR
    Label2600:

    ;
    ;    FILL ARRAY FOR LOWER MODULATOR STEP INTEGRATION
    ;

    CRNT(NCN-1)=CURENT*CONST
  ENDFOR
  Label5000:

  ;
  ;    COMPUTE CURRENTS BY TAKING DIFFERENCE BETWEEN ADJACENT MODULATOR
  ;        STEPS; IF LDIST= .TRUE., DIVIDE BY CHANNEL WIDTH
  ;
  FOR JJ = NC1,NC2 DO BEGIN
    IF(MODE EQ 1) THEN GOTO, Label5100

    ;
    ;    L-MODE;  LDIST = .TRUE.
    ;

    IF(LDIST EQ 'TRUE') THEN DISTR(JJ-1)=(CRNT(JJ-1)-CRR(JJ))/(DVOLTL(JJ-1))
    GOTO, Label5110

    ;
    ;    L-MODE;  LDIST = .TRUE.
    ;

    Label5100:
    IF(LDIST EQ 'TRUE') THEN DISTR(JJ-1)=(CRNT(JJ-1)-CRR(JJ))/(DVOLT(JJ-1))

    ;
    ;   L- OF M-MODE;  LDIST = .FALSE.
    ;

    Label5110: IF(LDIST NE 'TRUE') THEN DISTR(JJ-1)=(CRNT(JJ-1)-CRR(JJ))
  ENDFOR

  inddddd = where(finite(DISTR) NE 1, nfo)
  IF nfo GT 0 THEN DISTR(inddddd) = 0
  RETURN, DISTR
END

FUNCTION VGR_CupVelocity, VelocityIn, Cup, Structure, RETURN_r_pl=RETURN_r_pl
  ; (GX: VELOCITY IN WHICH THE PLASMA IS HITTING EACH CUP)
  ;
  ; Input
  ;
  ; VelocityIn : The cylindrical flowspeed of the Plasma in a Jupiter centered frame.
  ;              Corotation must be input. [0,0,0] for the velocity means the plasma
  ;              is completely stationary.
  ;
  ; Cup : The PLS sensor cup that is being evaluated.
  ;
  ; Structure : The VIPER 2.0 structure format.
  ;
  ; RETURN_r_pl : Keyword argument to have the function just return the s/c radial distance in planetary radii
  ;
  ;
  ; Output
  ;
  ; U : A 3 element array that holds the velocity in X, Y, and Z coordinates of the specified cup
  ;     These velocities are further used for the main fitting procedure in CUPINT/DCPINT or in
  ;     LABCUR/LDCUR.
  ;
  ;
  ; Necessary Files
  ;
  ; v1jup.ssedr.62-65.79_rob.txt    [48 columns, 2925 rows]
  ;   This file cont  ame
  ;
  ; SpiceRotations-vgr1.txt              [36 columns, 2925 rows]
  ;   This file contains the 6x6 transformation matrices for the state vectors to transform
  ;   between the System III coordinate frame and the ECL50 frame. This file is not explicitly
  ;   needed if one has access to the Spice kernels. Using the time units in the first few columns
  ;   of the previous text file it is possible to establish the Ephemeris time and use the Cspice
  ;   command of xform to get the proper transformation matrices.
  ;
  ;
  ; Explanation of Code
  ; 
  ; Takes in velocity of the Plasma in a non-rotating Jupiter centered frame. It then converts
  ; that velocity to the System III velocity of Jupiter, which accounts for the corotation of the
  ; planet. From here, it takes the cylindrical velocity and converts it to a Cartesian frame.
  ;
  ; It then subtracts the spacecraft (S/C) velocity at that given time to give the flux of plasma
  ; hitting the S/C in the System III frame. Then, using the SpiceRotations file, it transforms
  ; the System III state vector [x,y,z,vx,vy,vz] into a state vector in ECL50.
  ;
  ; From the ECL50 vector, a velocity vector is obtained for the velocity of the plasma in S/C
  ; coordinates. Then transformation matrices to the A, B, C, and D cups of the PLS instrument
  ; on the Voyager 1 spacecraft are utilized to give the Vx, Vy, and Vz velocities of the plasma
  ; in the frame of the cup so that the fitting procedure can be run correctly.
  ; 
  ; URANUS:
  ; GX: For Uranus things are more complicated, as we don't have a system equivalent to 
  ; System III for Uranus, and the described system in the SEDR files does not match to a known
  ; one. To match the SEDR file a new system was created named URANUS_FIT. It was created in SPICE, 
  ; based on past IAU reports: it started with the values of the 'Old Uranus PM' in 
  ; https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/a_old_versions/pck00003.tpc and then
  ; the POLE_RA and PM were adjusted (both the first and second component: first component shows
  ; the azimuth of the PM and the second shows the rotation rate, i.e. degrees per 24h). So for 
  ; Uranus System III is replaced with URANUS_FIT.
  ;


  ; ----------
  ; JUPITER
  ; ----------
  
  IF Structure.Planet_Number EQ 5 THEN BEGIN
    VelocityCylindrical = VelocityIn ; Rename variable so IDL doesn't overwrite
    
    ; Compute time from given structure
    Time = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
    
    IF Structure.Spacecraft EQ 1 THEN BEGIN                                                                                     ; LENGTH OF THE DATASET DEFINED 
      n = 2925 ; Number of Data points for analysis
      COMMON VGR1DATA
      Data = VGR1SSEDR ; Read in SSEDR                                                                                          
      rrr  = VGR1SPICE ; Read in rotation matrix from System III to ECL50 coordinates from the SPICE kernels
    ENDIF ELSE IF Structure.Spacecraft EQ 2 THEN BEGIN
      n = 15653 ; Number of Data points for analysis
      COMMON VGR2DATA
      Data = VGR2SSEDR ; Read in SSEDR
      rrr  = VGR2SPICE ; Read in rotation matrix from System III to ECL50 coordinates from the SPICE kernels
    ENDIF
  
    ; Compute the Decimal Date
    Decimal_Date = fltarr(n)
    FOR i = 0, n-1 DO BEGIN
      Decimal_Date[i] = (Data[2,i]*3600. + Data[3,i]*60. + Data[4,i] + Data[5,i] / 1000. ) / 86400. + Data[1,i]
    ENDFOR
    indss = closest(Decimal_Date[*], time) ; Determines index value closest to specified time
  
  
    xform = Make_array(6,6) ; Set up an array for the transformation matrix
    FOR i = 0,5 DO xform(i,0) = rrr(0+i,indss) ; Fill the elements of the array given the index value
    FOR i = 0,5 DO xform(i,1) = rrr(6+i,indss)
    FOR i = 0,5 DO xform(i,2) = rrr(12+i,indss)
    FOR i = 0,5 DO xform(i,3) = rrr(18+i,indss)
    FOR i = 0,5 DO xform(i,4) = rrr(24+i,indss)
    FOR i = 0,5 DO xform(i,5) = rrr(30+i,indss)
  
    xSC_S3 = Data[6,indss] ; System III cartesian coordinates
    ySC_S3 = Data[7,indss]
    zSC_S3 = Data[8,indss]
  
    ; The VelocityCylindrical is the velocity of the plasma in r, phi, and z (Cylindrical coordinates)
    ; This is not in the corotating frame of System III, so the plasma flow should be subtracted
    ; from the number to get the appropriate value
    ;
    ; The corotational flow goes as 12.572 km/s * Radius (Jupiter radii)
  
    RJupiter = 71492.0 ; Radius of Jupiter in km
    Radius = SQRT(xSC_S3*xSC_S3 + ySC_S3*ySC_S3) ; Radius is from axis of rotation -- this is rho essentially [GX]
    IF KEYWORD_SET(RETURN_r_pl) THEN RETURN, Radius/Rjupiter
    ; Neither the r or z components of the velocity will change with corotation
    VelocityCylindrical[1] = VelocityCylindrical[1] - 12.572*Radius/RJupiter
  
    ; Velocities of the S/C in System III
    VxSC_S3 = Data[9,indss]
    VySC_S3 = Data[10,indss]
    VzSC_S3 = Data[11,indss]
  
    ; Cartesian velocity of S/C in System III
    Velocity_SC_Cartesian_S3 = [VxSC_S3,VySC_S3,VzSC_S3]
  
    ; Split velocities into individual coordinates to change to a Cartesian coordinate system
    Vr      = VelocityCylindrical[0]
    Vphi    = VelocityCylindrical[1]
    Vz      = VelocityCylindrical[2]
  
    ; Transform R, Phi, and Z to X, Y, and Z
    ; GX: the transformation eq. from polar to Cartesian for the velocity are:
    ;   u_x = u_r cos(theta) - u_theta sin(theta)
    ;   u_y = u_r sin(theta) + u_theta cos(theta)
    ;   u_z = u_z
    ; and we have that: cos(theta) = x/r, sin(theta) = y/r
    VelocityXS3 = xSC_S3*Vr/Radius - ySC_S3*Vphi/(Radius) ; Divide Vphi by Radius since Vphi is given in km/s
    VelocityYS3 = ySC_S3*Vr/Radius + xSC_S3*Vphi/(Radius)
    VelocityZS3 = Vz ; By definition the Z axes are the same
  
    ; Cartesian Velocity of Plasma in System III
    Velocity_Plasma_Cartesian_S3 = [VelocityXS3,VelocityYS3,VelocityZS3]
  
  
    ; Create State Vectors for the Plasma and the S/C in System III and convert them to ECL50
    ; coordinates from the SPICE rotation matrix for this time.
    StateVecS3_Plasma = [xSC_S3,ySC_S3,zSC_S3,Velocity_Plasma_Cartesian_S3[0],Velocity_Plasma_Cartesian_S3[1],Velocity_Plasma_Cartesian_S3[2]]
    StateVecS3_SC = [xSC_S3,ySC_S3,zSC_S3,Velocity_SC_Cartesian_S3[0],Velocity_SC_Cartesian_S3[1],Velocity_SC_Cartesian_S3[2]]
    StateVecECL50_Plasma = xform ## StateVecS3_Plasma
    StateVecECL50_SC = xform ## StateVecS3_SC
  
    ; Extract the velocities in the ECL50 frame, since the positions are not needed for the flux.
    ; Then compute flux of plasma hitting the S/C in ECL50 frame.
    Velocity_ECL50_Plasma = [StatevecECL50_plasma[3],StatevecECL50_plasma[4],StatevecECL50_plasma[5]]
    Velocity_ECL50_SC = [StatevecECL50_SC[3],StatevecECL50_SC[4],StatevecECL50_SC[5]]
    VelocityFlux_ECL50 = Velocity_ECL50_Plasma - Velocity_ECL50_SC
  
    ; Read in rotation matrix from ECL50 to S/C pointing coordinate frame. Transpose the rotation
    ; matrix so that it is in the correct order for IDL matrix multiplication.
    ; The rotation matrix in SEDR is row-major!
    Rotation_ECL50_SC = [[Data[39,indss],Data[40,indss],Data[41,indss]], $
                         [Data[42,indss],Data[43,indss],Data[44,indss]], $
                         [Data[45,indss],Data[46,indss],Data[47,indss]]]
    Rotation_ECL50_SC = Transpose(Rotation_ECL50_SC)
  
    ; Compute flux in the S/C pointing frame from transformation in SSEDR
    VelocitySC = Rotation_ECL50_SC ## VelocityFlux_ECL50
  
    ;  Rotations from SpaceCraft to A, B, C, and D cups respectively. These are numbers from Bagenal 2 May, 2012 in VoyagerPLScoords-long.pdf on the LASP website
    ;  ARotation         = [[0.5,-0.865,0.0],[0.813,0.47,-0.342],[0.296,0.171,0.94]] ; OLDER VERSION
    ;  BRotation         = [[0.5,0.865,0.0],[-0.813,0.47,-0.342],[-0.296,0.171,0.94]]
    ;  CRotation         = [[1.0,0.0,0.0],[0.0,-0.94,-0.342],[0.0,-0.342,0.94]]
    ;  DRotation         = [[-0.682,0.731,0.0],[-0.026,-0.024,-0.999],[-0.731,-0.682,0.035]]
  
    ; Rotations from SpaceCraft to A, B, C, and D cups respectively.
    
    IF Structure.Spacecraft EQ 1 THEN BEGIN
      ARotation = [[0.5105,-0.8599,0.0000], $
                   [0.8080,0.4798,-0.3419], $
                   [0.2939,0.1745,0.9398]]
      BRotation = [[0.4939,0.8695,0.0000],   $
                   [-0.8171,0.4642,-0.3420], $
                   [-0.2974,0.1689,0.9397]]
      CRotation = [[-1.0000,0.0052,0.0000],   $
                   [-0.0049,-0.9428,-0.3333], $
                   [-0.0017,-0.3333,0.9428]]
      DRotation = [[-0.6800,0.7333,0.0000],   $
                   [-0.0216,-0.0201,-0.9996], $
                   [-0.7329,-0.6797,0.0295]]
    ENDIF ELSE IF Structure.Spacecraft EQ 2 THEN BEGIN
      ARotation = [[0.5105,-0.8599,0.0000], $
                   [0.8064,0.4788,-0.3469], $
                   [0.2983,0.1771,0.9379]]
      BRotation = [[0.5180,0.8554,0.0000],   $
                   [-0.8049,0.4875,-0.3384], $
                   [-0.2895,0.1753,0.9410]]
      CRotation = [[-1.0000,0.0098,0.0000],   $
                   [-0.0092,-0.9420,-0.3355], $
                   [-0.0033,-0.3354,0.9421]]
      DRotation = [[-0.6800,0.7333,0.0000],   $
                   [-0.0230,-0.0214,-0.9995], $
                   [-0.7329,-0.6796,0.0314]]
    ENDIF ELSE BEGIN
      PRINT, PRINT, STRCOMPRESS(STRING('There was not a Voyager ', LONG(Structure.Spacecraft), ' spacecraft. Returning a flowspeed of 0 km/s.'))
      RETURN, [0d,0d,0d]
    ENDELSE
  
  
  
    ; Change velocity flux into velocity in the correct cup's coordinates (X, Y, and Z)
    IF CUP EQ 'A' THEN U = ARotation ## VelocitySC
    IF CUP EQ 'B' THEN U = BRotation ## VelocitySC
    IF CUP EQ 'C' THEN U = CRotation ## VelocitySC
    IF CUP EQ 'D' THEN U = DRotation ## VelocitySC
  
    RETURN, U
    
  ; ----------
  ; URANUS
  ; ----------

  ENDIF ELSE IF Structure.Planet_Number EQ 7 THEN BEGIN
  
    VelocityCylindrical = VelocityIn ; Rename variable so IDL doesn't overwrite
  
    ; Compute time from given structure
    Time = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
  
      n = 10274 ; Number of SEDR file Data points for analysis    HARDCODED               GX: MAYBE READ THAT FROM THE FILE IF POSSIBLE?
      COMMON VGR2DATA_URANUS
      Data = VGR2SSEDR_URANUS ; Read in SSEDR
      rrr  = VGR2SPICE_URANUS ; Read in rotation matrix from URANUS_FIT to ECL50 coordinates from the SPICE kernels
  
    ; Compute the Decimal Date
    Decimal_Date = fltarr(n)
    FOR i = 0, n-1 DO BEGIN
      Decimal_Date[i] = (Data[2,i]*3600. + Data[3,i]*60. + Data[4,i] + Data[5,i] / 1000. ) / 86400. + Data[1,i]
    ENDFOR
    indss = closest(Decimal_Date[*], time) ; Determines index value closest to specified time
  
  
    xform = Make_array(6,6) ; Set up an array for the transformation matrix
    FOR i = 0,5 DO xform(i,0) = rrr(0+i,indss) ; Fill the elements of the array given the index value
    FOR i = 0,5 DO xform(i,1) = rrr(6+i,indss)
    FOR i = 0,5 DO xform(i,2) = rrr(12+i,indss)
    FOR i = 0,5 DO xform(i,3) = rrr(18+i,indss)
    FOR i = 0,5 DO xform(i,4) = rrr(24+i,indss)
    FOR i = 0,5 DO xform(i,5) = rrr(30+i,indss)
  
    xSC_URANUS_FIT = Data[6,indss] ; URANUS_FIT     ( GX: apparently it is from URANUS_FIT to ECL50 and not vice-versa )
    ySC_URANUS_FIT = Data[7,indss]
    zSC_URANUS_FIT = Data[8,indss]
  
    ; The VelocityCylindrical is the velocity of the plasma in r, phi, and z (Cylindrical coordinates)
    ; This is not in the corotating frame of URANUS_FIT, so the plasma flow should be subtracted
    ; from the number to get the appropriate value
    ;
    ; The corotational flow goes as -2.586 km/s * Radius (Uranus radii)
    u_rotation = -2.586 ; km/s
    RUranus = 25559.0 ; Equatorial Radius of Uranus [km]
    Radius = SQRT(xSC_URANUS_FIT*xSC_URANUS_FIT + ySC_URANUS_FIT*ySC_URANUS_FIT) ; Radius is from axis of rotation
    IF KEYWORD_SET(RETURN_r_pl) THEN RETURN, Radius/RUranus
    ; Neither the r or z components of the velocity will change with corotation
    VelocityCylindrical[1] = VelocityCylindrical[1] - u_rotation*Radius/RUranus
  
    ; Split velocities into individual coordinates to change to a Cartesian coordinate system
    Vr      = VelocityCylindrical[0]
    Vphi    = VelocityCylindrical[1]
    Vz      = VelocityCylindrical[2]
  
    ; Transform R, Phi, and Z to X, Y, and Z
    ; GX: the transformation eq. from polar to Cartesian for the velocity are:
    ;   u_x = u_r cos(theta) - u_theta sin(theta)
    ;   u_y = u_r sin(theta) + u_theta cos(theta)
    ;   u_z = u_z
    ; and we have that: cos(theta) = x/r, sin(theta) = y/r
    VelocityX_URANUS_FIT = xSC_URANUS_FIT*Vr/Radius - ySC_URANUS_FIT*Vphi/(Radius) ; Divide Vphi by Radius since Vphi is given in km/s
    VelocityY_URANUS_FIT = ySC_URANUS_FIT*Vr/Radius + xSC_URANUS_FIT*Vphi/(Radius)
    VelocityZ_URANUS_FIT = Vz ; By definition the Z axes are the same
  
    ; Cartesian Velocity of Plasma in URANUS_FIT
    Velocity_Plasma_Cartesian_URANUS_FIT = [VelocityX_URANUS_FIT,VelocityY_URANUS_FIT,VelocityZ_URANUS_FIT]
    
   ; Get the velocities of the S/C in URANUS_FIT
    VxSC_URANUS_FIT = Data[9,indss]
    VySC_URANUS_FIT = Data[10,indss]
    VzSC_URANUS_FIT = Data[11,indss]

    ; Cartesian velocity of S/C in URANUS_FIT
    Velocity_SC_Cartesian_URANUS_FIT = [VxSC_URANUS_FIT,VySC_URANUS_FIT,VzSC_URANUS_FIT] 
  
    ; Create State Vectors for the Plasma and the S/C in URANUS_FIT and convert them to ECL50
    ; coordinates from the SPICE rotation matrix for this time.
    StateVec_URANUS_FIT_Plasma = [xSC_URANUS_FIT,ySC_URANUS_FIT,zSC_URANUS_FIT,Velocity_Plasma_Cartesian_URANUS_FIT[0],Velocity_Plasma_Cartesian_URANUS_FIT[1],Velocity_Plasma_Cartesian_URANUS_FIT[2]]
    StateVec_URANUS_FIT_SC = [xSC_URANUS_FIT,ySC_URANUS_FIT,zSC_URANUS_FIT,Velocity_SC_Cartesian_URANUS_FIT[0],Velocity_SC_Cartesian_URANUS_FIT[1],Velocity_SC_Cartesian_URANUS_FIT[2]]
    StateVec_ECL50_Plasma = xform ## StateVec_URANUS_FIT_Plasma
    StateVec_ECL50_SC = xform ## StateVec_URANUS_FIT_SC
  
    ; Extract the velocities in the ECL50 frame, since the positions are not needed for the flux.
    ; Then compute flux of plasma hitting the S/C in ECL50 frame.
    Velocity_ECL50_Plasma = [Statevec_ECL50_plasma[3],Statevec_ECL50_plasma[4],Statevec_ECL50_plasma[5]]
    Velocity_ECL50_SC = [Statevec_ECL50_SC[3],Statevec_ECL50_SC[4],Statevec_ECL50_SC[5]]
    VelocityFlux_ECL50 = Velocity_ECL50_Plasma - Velocity_ECL50_SC
  
    ; Read in rotation matrix from ECL50 to S/C pointing coordinate frame. Transpose the rotation
    ; matrix so that it is in the correct order for IDL matrix multiplication.
    ; The rotation matrix in SEDR is row-major!
    Rotation_ECL50_SC = [[Data[39,indss],Data[40,indss],Data[41,indss]], $
      [Data[42,indss],Data[43,indss],Data[44,indss]], $
      [Data[45,indss],Data[46,indss],Data[47,indss]]]
    Rotation_ECL50_SC = Transpose(Rotation_ECL50_SC) ; by transposing we go from A->B to B->A
  ;STOP
    ; Compute flux in the S/C pointing frame from transformation in SSEDR
    VelocitySC = Rotation_ECL50_SC ## VelocityFlux_ECL50
  
    ;  Rotations from SpaceCraft to A, B, C, and D cups respectively. These are numbers from Bagenal 2 May, 2012 in VoyagerPLScoords-long.pdf on the LASP website
    ;  ARotation         = [[0.5,-0.865,0.0],[0.813,0.47,-0.342],[0.296,0.171,0.94]] ; OLDER VERSION
    ;  BRotation         = [[0.5,0.865,0.0],[-0.813,0.47,-0.342],[-0.296,0.171,0.94]]
    ;  CRotation         = [[1.0,0.0,0.0],[0.0,-0.94,-0.342],[0.0,-0.342,0.94]]
    ;  DRotation         = [[-0.682,0.731,0.0],[-0.026,-0.024,-0.999],[-0.731,-0.682,0.035]]
  
    ; Rotations from SpaceCraft to A, B, C, and D cups respectively.
  
    IF Structure.Spacecraft EQ 2 THEN BEGIN
      ARotation = [[.5105,-0.8599,0.0000], $
        [0.8064,0.4788,-0.3469], $
        [0.2983,0.1771,0.9379]]
      BRotation = [[0.5180,0.8554,0.0000],   $
        [-0.8049,0.4875,-0.3384], $
        [-0.2895,0.1753,0.9410]]
      CRotation = [[-1.0000,0.0098,0.0000],   $
        [-0.0092,-0.9420,-0.3355], $
        [-0.0033,-0.3354,0.9421]]
      DRotation = [[-0.6800,0.7333,0.0000],   $
        [-0.0230,-0.0214,-0.9995], $
        [-0.7329,-0.6796,0.0314]]
    ENDIF ELSE BEGIN
      PRINT, PRINT, STRCOMPRESS(STRING('There was not a Voyager ', LONG(Structure.Spacecraft), ' spacecraft. Returning a flowspeed of 0 km/s.'))
      RETURN, [0d,0d,0d]
    ENDELSE
  
    ; Change velocity flux into velocity in the correct cup's coordinates (X, Y, and Z)
    IF CUP EQ 'A' THEN U = ARotation ## VelocitySC
    IF CUP EQ 'B' THEN U = BRotation ## VelocitySC
    IF CUP EQ 'C' THEN U = CRotation ## VelocitySC
    IF CUP EQ 'D' THEN U = DRotation ## VelocitySC
  
    RETURN, U
  
  
  ENDIF ELSE IF Structure.Planet_Number EQ 8 THEN BEGIN
  ; ----------
  ; NEPTUNE
  ; ----------

  VelocityCylindrical = VelocityIn ; Rename variable so IDL doesn't overwrite

  ; Compute time from given structure
  Time = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d

  n = 28325 ; Number of SEDR file Data points for analysis    HARDCODED               GX: MAYBE READ THAT FROM THE FILE IF POSSIBLE?
  COMMON VGR2DATA_NEPTUNE
  Data = VGR2SSEDR_NEPTUNE ; Read in SSEDR
  rrr  = VGR2SPICE_NEPTUNE ; Read in rotation matrix from V2_NEPTUNE_PLS to ECL50 coordinates from the SPICE kernels
  
  ; Compute the Decimal Date
  Decimal_Date = fltarr(n)
  FOR i = 0, n-1 DO BEGIN
    Decimal_Date[i] = (Data[2,i]*3600. + Data[3,i]*60. + Data[4,i] + Data[5,i] / 1000. ) / 86400. + Data[1,i]
  ENDFOR
  indss = closest(Decimal_Date[*], time) ; Determines index value closest to specified time


  xform = Make_array(6,6) ; Set up an array for the transformation matrix
  FOR i = 0,5 DO xform(i,0) = rrr(0+i,indss) ; Fill the elements of the array given the index value
  FOR i = 0,5 DO xform(i,1) = rrr(6+i,indss)
  FOR i = 0,5 DO xform(i,2) = rrr(12+i,indss)
  FOR i = 0,5 DO xform(i,3) = rrr(18+i,indss)
  FOR i = 0,5 DO xform(i,4) = rrr(24+i,indss)
  FOR i = 0,5 DO xform(i,5) = rrr(30+i,indss)

  xSC_V2_NEPTUNE_PLS = Data[6,indss] ; V2_NEPTUNE_PLS  ( GX: it is from V2_NEPTUNE_PLS to ECL50 and not vice-versa )
  ySC_V2_NEPTUNE_PLS = Data[7,indss]
  zSC_V2_NEPTUNE_PLS = Data[8,indss]

  ; The VelocityCylindrical is the velocity of the plasma in r, phi, and z (Cylindrical coordinates)
  ; This is not in the corotating frame of V2_NEPTUNE_PLS, so the plasma flow should be subtracted
  ; from the number to get the appropriate value
  ;
  ; The corotational flow goes as 2.683 km/s * Radius (Neptune radii)
  u_rotation = 2.683d ; km/s
  RNeptune = 24764d ; Equatorial Radius of Neptune [km]
  Radius = SQRT(xSC_V2_NEPTUNE_PLS*xSC_V2_NEPTUNE_PLS + ySC_V2_NEPTUNE_PLS*ySC_V2_NEPTUNE_PLS) ; Radius is from axis of rotation
  IF KEYWORD_SET(RETURN_r_pl) THEN RETURN, Radius/RNeptune
  ; Neither the r or z components of the velocity will change with corotation
  VelocityCylindrical[1] = VelocityCylindrical[1] - u_rotation*Radius/RNeptune

  ; Split velocities into individual coordinates to change to a Cartesian coordinate system
  Vr      = VelocityCylindrical[0]
  Vphi    = VelocityCylindrical[1]
  Vz      = VelocityCylindrical[2]

  ; Transform R, Phi, and Z to X, Y, and Z
  ; GX: the transformation eq. from polar to Cartesian for the velocity are:
  ;   u_x = u_r cos(theta) - u_theta sin(theta)
  ;   u_y = u_r sin(theta) + u_theta cos(theta)
  ;   u_z = u_z
  ; and we have that: cos(theta) = x/r, sin(theta) = y/r
  VelocityX_V2_NEPTUNE_PLS = xSC_V2_NEPTUNE_PLS*Vr/Radius - ySC_V2_NEPTUNE_PLS*Vphi/(Radius) ; Divide Vphi by Radius since Vphi is given in km/s
  VelocityY_V2_NEPTUNE_PLS = ySC_V2_NEPTUNE_PLS*Vr/Radius + xSC_V2_NEPTUNE_PLS*Vphi/(Radius)
  VelocityZ_V2_NEPTUNE_PLS = Vz ; By definition the Z axes are the same

  ; Cartesian Velocity of Plasma in V2_NEPTUNE_PLS
  Velocity_Plasma_Cartesian_V2_NEPTUNE_PLS = [VelocityX_V2_NEPTUNE_PLS, VelocityY_V2_NEPTUNE_PLS, VelocityZ_V2_NEPTUNE_PLS]

  ; Get the velocities of the S/C in V2_NEPTUNE_PLS
  VxSC_V2_NEPTUNE_PLS = Data[9,indss]
  VySC_V2_NEPTUNE_PLS = Data[10,indss]
  VzSC_V2_NEPTUNE_PLS = Data[11,indss]

  ; Cartesian velocity of S/C in V2_NEPTUNE_PLS
  Velocity_SC_Cartesian_V2_NEPTUNE_PLS = [VxSC_V2_NEPTUNE_PLS, VySC_V2_NEPTUNE_PLS, VzSC_V2_NEPTUNE_PLS]

  ; Create State Vectors for the Plasma and the S/C in V2_NEPTUNE_PLS and convert them to ECL50
  ; coordinates from the SPICE rotation matrix for this time.
  StateVec_V2_NEPTUNE_PLS_Plasma = [xSC_V2_NEPTUNE_PLS,ySC_V2_NEPTUNE_PLS,zSC_V2_NEPTUNE_PLS,Velocity_Plasma_Cartesian_V2_NEPTUNE_PLS[0],Velocity_Plasma_Cartesian_V2_NEPTUNE_PLS[1],Velocity_Plasma_Cartesian_V2_NEPTUNE_PLS[2]]
  StateVec_V2_NEPTUNE_PLS_SC = [xSC_V2_NEPTUNE_PLS,ySC_V2_NEPTUNE_PLS,zSC_V2_NEPTUNE_PLS,Velocity_SC_Cartesian_V2_NEPTUNE_PLS[0],Velocity_SC_Cartesian_V2_NEPTUNE_PLS[1],Velocity_SC_Cartesian_V2_NEPTUNE_PLS[2]]
  StateVec_ECL50_Plasma = xform ## StateVec_V2_NEPTUNE_PLS_Plasma
  StateVec_ECL50_SC = xform ## StateVec_V2_NEPTUNE_PLS_SC

  ; Extract the velocities in the ECL50 frame, since the positions are not needed for the flux.
  ; Then compute flux of plasma hitting the S/C in ECL50 frame.
  Velocity_ECL50_Plasma = [Statevec_ECL50_plasma[3],Statevec_ECL50_plasma[4],Statevec_ECL50_plasma[5]]
  Velocity_ECL50_SC = [Statevec_ECL50_SC[3],Statevec_ECL50_SC[4],Statevec_ECL50_SC[5]]
  VelocityFlux_ECL50 = Velocity_ECL50_Plasma - Velocity_ECL50_SC

  ; Read in rotation matrix from ECL50 to S/C pointing coordinate frame. Transpose the rotation
  ; matrix so that it is in the correct order for IDL matrix multiplication.
  ; The indices need to agree to the SEDR - each SEDR had different number of columns.
  ; The rotation matrix in SEDR is row-major! !!!!!!!!!! NEED TO CHANGE INDICES WHEN CHANGING PLANETS!!!
  Rotation_ECL50_SC = [[Data[30,indss],Data[31,indss],Data[32,indss]], $
    [Data[33,indss],Data[34,indss],Data[35,indss]], $
    [Data[36,indss],Data[37,indss],Data[38,indss]]]
  Rotation_ECL50_SC = Transpose(Rotation_ECL50_SC) ; by transposing we go from A->B to B->A

  ; Compute flux in the S/C pointing frame from transformation in SSEDR
  VelocitySC = Rotation_ECL50_SC ## VelocityFlux_ECL50

  ;  Rotations from SpaceCraft to A, B, C, and D cups respectively. These are numbers from Bagenal 2 May, 2012 in VoyagerPLScoords-long.pdf on the LASP website
  ;  ARotation         = [[0.5,-0.865,0.0],[0.813,0.47,-0.342],[0.296,0.171,0.94]] ; OLDER VERSION
  ;  BRotation         = [[0.5,0.865,0.0],[-0.813,0.47,-0.342],[-0.296,0.171,0.94]]
  ;  CRotation         = [[1.0,0.0,0.0],[0.0,-0.94,-0.342],[0.0,-0.342,0.94]]
  ;  DRotation         = [[-0.682,0.731,0.0],[-0.026,-0.024,-0.999],[-0.731,-0.682,0.035]]

  ; Rotations from SpaceCraft to A, B, C, and D cups respectively.

  IF Structure.Spacecraft EQ 2 THEN BEGIN
    ARotation = [[.5105,-0.8599,0.0000], $
      [0.8064,0.4788,-0.3469], $
      [0.2983,0.1771,0.9379]]
    BRotation = [[0.5180,0.8554,0.0000],   $
      [-0.8049,0.4875,-0.3384], $
      [-0.2895,0.1753,0.9410]]
    CRotation = [[-1.0000,0.0098,0.0000],   $
      [-0.0092,-0.9420,-0.3355], $
      [-0.0033,-0.3354,0.9421]]
    DRotation = [[-0.6800,0.7333,0.0000],   $
      [-0.0230,-0.0214,-0.9995], $
      [-0.7329,-0.6796,0.0314]]
  ENDIF ELSE BEGIN
    PRINT, PRINT, STRCOMPRESS(STRING('There was not a Voyager ', LONG(Structure.Spacecraft), ' spacecraft. Returning a flowspeed of 0 km/s.'))
    RETURN, [0d,0d,0d]
  ENDELSE

  ; Change velocity flux into velocity in the correct cup's coordinates (X, Y, and Z)
  IF CUP EQ 'A' THEN U = ARotation ## VelocitySC
  IF CUP EQ 'B' THEN U = BRotation ## VelocitySC
  IF CUP EQ 'C' THEN U = CRotation ## VelocitySC
  IF CUP EQ 'D' THEN U = DRotation ## VelocitySC

  RETURN, U

ENDIF
        
END


FUNCTION Euler, Alpha, Beta, Gamma
  ; Calculates a rotation matrix given the Euler angles Alpha, Beta, and Gamma using a
  ; ZYZ calculation. Angles input in degrees and converted to radians for IDL calculation

  ; Inputs
  ; Alpha   - First rotation, around z, to give x'y'z'
  ; Beta    - Second rotation, around y', to give x''y''z''
  ; Gamma   - Third rotation, around z'', to give XYZ
  ;
  ; Outputs
  ; R       - Standard 3x3 matrix for rotations

  a = Double(Alpha)*!Pi/180.0
  b = Double(Beta)*!Pi/180.0
  g = Double(Gamma)*!Pi/180.0

  R = Make_array(3,3)
  e1 = Cos(g)*Cos(b)*Cos(a) - Sin(g)*Sin(a)
  e2 = Cos(g)*Cos(B)*Sin(a) + Sin(g)*Cos(a)
  e3 = -Cos(g)*Sin(b)
  e4 = -Sin(g)*Cos(b)*Cos(a) - Cos(g)*Sin(a)
  e5 = -Sin(g)*Cos(b)*Sin(a) + Cos(g)*Cos(a)
  e6 = Sin(g)*Sin(b)
  e7 = Sin(b)*Cos(a)
  e8 = Sin(b)*Sin(a)
  e9 = Cos(b)
  R = [[e1,e2,e3],[e4,e5,e6],[e7,e8,e9]]

  RETURN, R
END

FUNCTION CupRotations
  ; Copy and Paste to command line or run from here, either works
  ; It generates a text file with the data for the cup rotations
  ; It is kept so that it can be run if the Voyager Memo #29 has wrong angles

  ; Voyager 1 Matrices
  RAV1 = Euler(30.70,19.99,270)
  RBV1 = Euler(150.4,20,270)
  RCV1 = Euler(269.7,19.47,270)
  RDV1 = Euler(222.84,88.31,270)

  ; Voyager 2 Matrices
  RAV2 = Euler(30.70,20.3,270)
  RBV2 = Euler(148.8,19.78,270)
  RCV2 = Euler(269.44,19.6,270)
  RDV2 = Euler(222.84,88.2,270)

  OPENW, 1, "Voyager 1 and 2 Cup Rotations From New XSLX Data.txt"
  PRINTF, 1, "Voyager 1"
  PRINTF, 1, "Cup A"
  PRINTF, 1, RAV1
  PRINTF, 1, "Cup B"
  PRINTF, 1, RBV1
  PRINTF, 1, "Cup C"
  PRINTF, 1, RCV1
  PRINTF, 1, "Cup D"
  PRINTF, 1, RDV1
  PRINTF, 1, " "
  PRINTF, 1, "Voyager 2"
  PRINTF, 1, "Cup A"
  PRINTF, 1, RAV2
  PRINTF, 1, "Cup B"
  PRINTF, 1, RBV2
  PRINTF, 1, "Cup C"
  PRINTF, 1, RCV2
  PRINTF, 1, "Cup D"
  PRINTF, 1, RDV2
  CLOSE, 1
END

FUNCTION PLOT_RANGES_LOG, Array

  ; Input
  ; Array    - An array of values to find a good logarithmic plot range for. It will return
  ;            the nearest power of 10 above the maximum value and 4 orders of magnitude below
  ;            that value.
  ;
  ;            Example: If the maximum value is 4.3D6, then it will return 1D3 to 1D7
  ;
  ; Output
  ; YRANGES   - A 2 element array consisting of ymin and ymax for the plot ranges.

  ; Find ranges for plots
  MaxArray = Max(Array)
  Result = 0 ; Initialize a value
  FOR JJJ = 0, 10 DO BEGIN
    Result_Old = Result
    Result = MaxArray MOD 10.0^JJJ
    IF Result GT Result_Old THEN ymax = JJJ
  ENDFOR

  ymin = ymax - 4
  ; Force the lower limit to be 10^2 for the plotting procedure
  IF ymin LT 2 THEN ymin = 2
  ymax = 10.0^(ymax)
  ; Make sure that the data isn't at the very top of the plot and hard to see
  IF ymax/MaxArray LT 1.4 THEN ymax = ymax*10D
  ymin = 10.0^(ymin)

  YRANGES = [ymin,ymax]
  RETURN, YRANGES
END

FUNCTION Convert_Eddie_Ratios, Input_Ratios                                                                           ; MAYBE THIS IS NOT NEEDED FOR URANUS
  ; Currently just a quick function to read in whatever Eddie gives me and output it
  ;
  ; For the VIPER analysis around Jupiter, we presume O+ to be the dominant species when we are
  ; fixing other species to that composition. Therefore, we will need to convert the density ratios
  ; that Eddie gives us so we can easily fill the tables that we do have.
  ;
  ; Eddie gives ratios of S+/S++, S+++/S++, O+/S++, O++/S++
  ;
  ; But we want ratios in terms of ion/O+
  ;
  ; Assuming the format of his array has not changed, they are in the order stated as above:
  
  
  Dummy = MAKE_ARRAY(4, /DOUBLE)
  ; Read in Eddies 4 ratios that he provides
  Dummy = Input_Ratios
  
  ; Re-order it from low to high
  Ratios_Eddie = MAKE_ARRAY(4, /DOUBLE)
  Ratios_Eddie[0] = Dummy[3]
  Ratios_Eddie[1] = Dummy[1]
  Ratios_Eddie[2] = Dummy[2]
  Ratios_Eddie[3] = Dummy[0]
  
  Ratios = MAKE_ARRAY(5, /DOUBLE)
  Ratios[0] = 1d ; O+/O+ is 1 - I will now go from lowest Mass/Charge to highest so O++, S+++, S++, S+
  FOR i = 0,3 DO BEGIN
    Ratios[i+1] = 1d/(Ratios_Eddie[2]/Ratios_Eddie[i])
    IF i EQ 2 THEN Ratios[i+1] = 1d/Ratios_Eddie[2]
  ENDFOR
  

  RETURN, Ratios
END

FUNCTION CASSINI_RATIOS, Structure_In, Uncertainty_Flag
  ; Reads in the structure used by the VIPER 2.0 program and if Ratios are requested it accordingly obtains the correct ratios for the timestamp
  ; From Cassini UVIS observations, there are ratios for 6.5 to 8.5 Rj
  ; 
  ; Inputs
  ; Structure_In - Structure from the VIPER program
  ; Uncertainty_Flag - A 0 (doesn't exist) or 1 (does exist) flag to denote if uncertainties are on.
  ;                    If uncertainties are being fit, then don't print the errors - just to clean up command line
  ; 
  ; Outputs
  ; Structure_Out - Structure with the properly tied density parameters using Cassini UVIS data
  
  ; Count number of species
  StringCompare = STRCMP('Species_', TAG_NAMES(Structure_In), 3, /FOLD_CASE)
  ind = WHERE(StringCompare EQ 1, nind)
  n_species = nind/5
  
  Arguments = MAKE_ARRAY(n_species)
  Masses    = Arguments
  Charges   = Arguments
  Densities = Arguments
  
  ; Determine what elements in the structure should be tied together
  FOR k = 1, n_species DO BEGIN
    str0 = STRCOMPRESS(STRING('Species_', k), /REMOVE_ALL)
    IF k+1 LT 10 THEN str0cmp = STRCMP(str0, TAG_NAMES(Structure_In), 10, /FOLD_CASE) ELSE $
      str0cmp = STRCMP(str0, TAG_NAMES(Structure_In), 11, /FOLD_CASE)
    ind0 = WHERE(str0cmp EQ 1)
    Arguments(k-1) = Structure_In.(ind0)

    str1 = STRCOMPRESS(STRING('Species_', k, '_A'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str1cmp = STRCMP(str1, TAG_NAMES(Structure_In), 12, /FOLD_CASE) ELSE $
      str1cmp = STRCMP(str1, TAG_NAMES(Structure_In), 13, /FOLD_CASE)
    ind1 = WHERE(str1cmp EQ 1)
    Masses(k-1) = Structure_In.(ind1)

    str2 = STRCOMPRESS(STRING('Species_', k, '_Z'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str2cmp = STRCMP(str2, TAG_NAMES(Structure_In), 12, /FOLD_CASE) ELSE $
      str2cmp = STRCMP(str2, TAG_NAMES(Structure_In), 13, /FOLD_CASE)
    ind2 = WHERE(str2cmp EQ 1)
    Charges(k-1) = Structure_In.(ind2)
    
    ; Additionally store locations of densities
    str3 = STRCOMPRESS(STRING('Species_', k, '_N'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str3cmp = STRCMP(str3, TAG_NAMES(Structure_In), 12, /FOLD_CASE) ELSE $
      str3cmp = STRCMP(str3, TAG_NAMES(Structure_In), 13, /FOLD_CASE)
    ind3 = WHERE(str3cmp EQ 1)
    Densities(k-1) = ind3 ; Store array locations of densities for later indexing use
  ENDFOR
  
  ; Create a location array
  Locations = Arguments
  ind = WHERE((Arguments EQ 4) OR (Arguments EQ 5), nind)
  Locations[*] = 0L
  IF nind GT 0 THEN Locations[ind] = 1L
  
  ; Read in Cassini Ratio text file
  ; Columns are as follows: Radial value, O+/O+, O++/O+, S+++/O+, S++/O+, S+/O+
  Cass_UVIS = MAKE_ARRAY(6, 10001)
  
  ; Change to proper directory and read in the text file
  CD, CURRENT = RESTORE_DIR
  STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/Data/Jupiter/Chemistry_Ratios'))
  CD, STR_DIR
  OPENR, 1, "Cassini_UVIS_Ratios.txt"
  READF, 1, Cass_UVIS
  CLOSE, 1
  CD, RESTORE_DIR
  
  ; Determine radial distance of spacecraft and get matching value in Cassini UVIS data - Use Vgr_CupVelocity function and return Rj value
  Rj = VGR_CupVelocity([0d,0d,0d], 'A', Structure_In, /RETURN_RJ) ; Can use dummy arguments for the first 2 arguments of the function
  Cass_UVIS_RJ = Cass_UVIS[0,*]
  ind_Rj = CLOSEST(Cass_UVIS_RJ[*], Rj)
  
  ; Print some warning messages if chemistry table hasn't been run for this region yet
  IF Uncertainty_Flag EQ 0 THEN BEGIN
    IF Rj LT MIN(Cass_UVIS_RJ) THEN PRINT, STRCOMPRESS(STRING('Distance less than ', STRING(MIN(Cass_UVIS_RJ)), $
                                                              ' Rj. The chemistry model has not been run for this region. Assuming composition from ', $ 
                                                               STRING(MIN(Cass_UVIS_RJ)),' Rj.', FORMAT = '(A,F0.2,A,F0.2,A)'))
    IF Rj GT MAX(Cass_UVIS_RJ) THEN PRINT, STRCOMPRESS(STRING('Distance greater than ', STRING(MAX(Cass_UVIS_RJ)), $
                                                               ' Rj. Assuming that ions have stopped mixing and using the fixed composition from ', $
                                                               STRING(MAX(Cass_UVIS_RJ)),' Rj.', FORMAT = '(A,F0.2,A,F0.2,A)'))                                                           
  ENDIF
  
  ; Set output Structure
  Structure_Out = Structure_In
  
  ; Loop over Locations and determine where to replace with appropriate values
  FOR i = 0, N_ELEMENTS(ind)-1 DO BEGIN
    A0 = Masses[ind[i]]
    Z0 = Charges[ind[i]]
    
    DELVAR, Dens_Ratio ; Reset the variable
    ; Begin conditional statements
    CASE A0 OF
      1:  Dens_Ratio = Structure_In.(Densities[ind[i]]) ; Don't change the setting on the H+ species - must be specified currently
      16: CASE Z0 OF
            1: Dens_Ratio = Structure_In.(Densities[ind[i]]) ; Don't change the setting on the O+ species - It is the dominant species
            2: Dens_Ratio = Cass_UVIS[2, ind_Rj]
          ENDCASE
      32: CASE Z0 OF
            3: Dens_Ratio = Cass_UVIS[3, ind_Rj]
            2: Dens_Ratio = Cass_UVIS[4, ind_Rj]
            1: Dens_Ratio = Cass_UVIS[5, ind_Rj]
          ENDCASE
      ELSE: Dens_Ratio = Structure_In.(Densities[ind[i]])
    ENDCASE
    
    ; Check that variable has been defined - If not, use input by default
    IF ISA(Dens_Ratio) EQ 0 THEN BEGIN
      PRINT, STRCOMPRESS(STRING('Do not have a ratio to O+ for Charge: ', LONG(Z0), ' and Mass: ', LONG(A0), '.'))
      Dens_Ratio = Structure_In.(Densities[ind[i]])
    ENDIF
    
    Structure_Out.(Densities[ind[i]]) = Dens_Ratio  
  ENDFOR
  
  Return, Structure_Out
END

FUNCTION DELAMERE_RATIOS, Structure_In, Uncertainty_Flag
  ; Reads in the structure used by the VIPER 2.0 program and if Ratios are requested it accordingly obtains the correct ratios for the timestamp
  ; From Cassini UVIS observations, there are ratios for 6.5 to 8.5 Rj
  ;
  ; Inputs
  ; Structure_In - Structure from the VIPER program
  ; Uncertainty_Flag - A 0 (doesn't exist) or 1 (does exist) flag to denote if uncertainties are on.
  ;                    If uncertainties are being fit, then don't print the errors - just to clean up command line
  ;
  ; Outputs
  ; Structure_Out - Structure with the properly tied density parameters using Cassini UVIS data

  ; Count number of species
  StringCompare = STRCMP('Species_', TAG_NAMES(Structure_In), 3, /FOLD_CASE)
  ind = WHERE(StringCompare EQ 1, nind)
  n_species = nind/6

  Arguments = MAKE_ARRAY(n_species)
  Masses    = Arguments
  Charges   = Arguments
  Densities = Arguments
  
  ; Determine what elements in the structure should be tied together
  FOR k = 1, n_species DO BEGIN
    str0 = STRCOMPRESS(STRING('Species_', k), /REMOVE_ALL)
    IF k+1 LT 10 THEN str0cmp = STRCMP(str0, TAG_NAMES(Structure_In), 10, /FOLD_CASE) ELSE $
      str0cmp = STRCMP(str0, TAG_NAMES(Structure_In), 11, /FOLD_CASE)
    ind0 = WHERE(str0cmp EQ 1)
    Arguments(k-1) = Structure_In.(ind0)

    str1 = STRCOMPRESS(STRING('Species_', k, '_A'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str1cmp = STRCMP(str1, TAG_NAMES(Structure_In), 12, /FOLD_CASE) ELSE $
      str1cmp = STRCMP(str1, TAG_NAMES(Structure_In), 13, /FOLD_CASE)
    ind1 = WHERE(str1cmp EQ 1)
    Masses(k-1) = Structure_In.(ind1)

    str2 = STRCOMPRESS(STRING('Species_', k, '_Z'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str2cmp = STRCMP(str2, TAG_NAMES(Structure_In), 12, /FOLD_CASE) ELSE $
      str2cmp = STRCMP(str2, TAG_NAMES(Structure_In), 13, /FOLD_CASE)
    ind2 = WHERE(str2cmp EQ 1)
    Charges(k-1) = Structure_In.(ind2)

    ; Additionally store locations of densities
    str3 = STRCOMPRESS(STRING('Species_', k, '_N'), /REMOVE_ALL)
    IF k+1 LT 10 THEN str3cmp = STRCMP(str3, TAG_NAMES(Structure_In), 12, /FOLD_CASE) ELSE $
      str3cmp = STRCMP(str3, TAG_NAMES(Structure_In), 13, /FOLD_CASE)
    ind3 = WHERE(str3cmp EQ 1)
    Densities(k-1) = ind3 ; Store array locations of densities for later indexing use
  ENDFOR

  ; Create a location array
  Locations = Arguments
  ind = WHERE((Arguments EQ 4) OR (Arguments EQ 5), nind)
  Locations[*] = 0L
  IF nind GT 0 THEN Locations[ind] = 1L

  ; Read in Delamere Ratio text file
  ; Columns are as follows: Radial value, O+/O+, O++/O+, S+++/O+, S++/O+, S+/O+
  DEL_UVIS = MAKE_ARRAY(6, 2969)
  
  ; Change to directory with chemistry ratios
  CD, CURRENT = RESTORE_DIR
  STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/Data/Jupiter/Chemistry_Ratios'))
  CD, STR_DIR
  OPENR, 1, "Delamere_Ratios.txt"
  READF, 1, DEL_UVIS
  CLOSE, 1
  CD, RESTORE_DIR

  ; Determine radial distance of spacecraft and get matching value in Cassini UVIS data - Use Vgr_CupVelocity function and return Rj value
  Rj = VGR_CupVelocity([0d,0d,0d], 'A', Structure_In, /RETURN_RJ) ; Can use dummy arguments for the first 2 arguments of the function
  DEL_UVIS_RJ = DEL_UVIS[0,*]
  ind_Rj = CLOSEST(DEL_UVIS_RJ[*], Rj)

  ; Print some warning messages if chemistry table hasn't been run for this region yet
  IF Uncertainty_Flag EQ 0 THEN BEGIN
    IF Rj LT MIN(DEL_UVIS_RJ) THEN PRINT, STRCOMPRESS(STRING('Distance less than ', STRING(MIN(DEL_UVIS_RJ)), $
                                          ' Rj. The chemistry model has not been run for this region. Assuming composition from ', $
                                          STRING(MIN(DEL_UVIS_RJ)),' Rj.', FORMAT = '(A,F0.2,A,F0.2,A)'))
    IF Rj GT MAX(DEL_UVIS_RJ) THEN PRINT, STRCOMPRESS(STRING('Distance greater than ', STRING(MAX(DEL_UVIS_RJ)), $
                                          ' Rj. Assuming that ions have stopped mixing and using the fixed composition from ', $
                                          STRING(MAX(DEL_UVIS_RJ)),' Rj.', FORMAT = '(A,F0.2,A,F0.2,A)'))
  ENDIF

  ; Set output Structure
  Structure_Out = Structure_In

  ; Loop over Locations and determine where to replace with appropriate values
  FOR i = 0, N_ELEMENTS(ind)-1 DO BEGIN
    A0 = Masses[ind[i]]
    Z0 = Charges[ind[i]]

    DELVAR, Dens_Ratio ; Reset the variable
    ; Begin conditional statements
    CASE A0 OF
      1:  Dens_Ratio = Structure_In.(Densities[ind[i]]) ; Don't change the setting on the H+ species - must be specified currently
      16: CASE Z0 OF
      1: Dens_Ratio = Structure_In.(Densities[ind[i]]) ; Don't change the setting on the O+ species - It is the dominant species
      2: Dens_Ratio = DEL_UVIS[2, ind_Rj]
    ENDCASE
    32: CASE Z0 OF
    3: Dens_Ratio = DEL_UVIS[3, ind_Rj]
    2: Dens_Ratio = DEL_UVIS[4, ind_Rj]
    1: Dens_Ratio = DEL_UVIS[5, ind_Rj]
  ENDCASE
  ELSE: Dens_Ratio = Structure_In.(Densities[ind[i]])
ENDCASE

; Check that variable has been defined - If not, use input by default
IF ISA(Dens_Ratio) EQ 0 THEN BEGIN
  PRINT, STRCOMPRESS(STRING('Do not have a ratio to O+ for Charge: ', LONG(Z0), ' and Mass: ', LONG(A0), '.'))
  Dens_Ratio = Structure_In.(Densities[ind[i]])
ENDIF

Structure_Out.(Densities[ind[i]]) = Dens_Ratio
ENDFOR

Return, Structure_Out
END

FUNCTION VIPER_ST_JOIN, st1_AAAAAAAAAAAAAAAAAAAAA, st2_BBBBBBBBBBBBBBBBBBBBB, NOTIMECHECK = NOTIMECHECK
  ; Code originally from Rob Wilson (JADE_ST_JOIN.pro). Received May 18th, 2016 and modified for the purposes of the VIPER combination file - Some parts of this code are likely ignored

  ON_ERROR,2

  st1 = st1_AAAAAAAAAAAAAAAAAAAAA ; hack not to overwrite input
  st2 = st2_BBBBBBBBBBBBBBBBBBBBB ; hack not to overwrite input
  IF KEYWORD_SET(NOTIMECHECK) THEN NOTIMECHECK = NOTIMECHECK ELSE NOTIMECHECK = 0 ; default of off - i.e. do a time check

  ; If a jade_read_sci can not find filename, it returns J = -1 (type INT)
  ; Nothing to add!  Return other one.
  IF (ISA(st1,'Struct') + ISA(st2,'Struct') NE 2) THEN BEGIN
    ; I expect -1 for fill, but what if [] or !Null
    IF ISA(st1,/Null) THEN st1 = -1
    IF ISA(st2,/Null) THEN st2 = -1
    ; now continue
    ; if st1 is -1 but st2 is a struct
    IF ISA(st1,/NUMBER) + ISA(st2,'Struct') EQ 2 THEN IF (ISA(st1,/SCALAR) EQ 1) AND (st1[0] EQ -1) THEN RETURN, st2
    ; if st2 is -1 but st1 is a struct
    IF ISA(st2,/NUMBER) + ISA(st1,'Struct') EQ 2 THEN IF (ISA(st2,/SCALAR) EQ 1) AND (st2[0] EQ -1) THEN RETURN, st1
    ; what if both are -1
    IF ISA(st1,/NUMBER) + ISA(st2,/NUMBER) EQ 2 THEN IF (ISA(st1,/SCALAR) EQ 1) AND (st1[0] EQ -1) AND (ISA(st2,/SCALAR) EQ 1) AND (st2[0] EQ -1) THEN RETURN, st1
    ; Should not get here otherwise
    MESSAGE,'ERROR, Must feed jade_st_join two structures (or -1 and a structure)'
  ENDIF

  ; original way
  IF ((N_ELEMENTS(st1) EQ 0 ) OR (N_ELEMENTS(st2) EQ 0)) THEN MESSAGE,'STRUCTURES MUST BE DEFINED'
  IF SIZE(st1,/TYPE) NE 8 THEN MESSAGE,'FIRST ARGUMENT MUST BE A *STRUCTURE* YOU WISH TO APPEND'
  IF SIZE(st2,/TYPE) NE 8 THEN MESSAGE,'SECOND ARGUMENT MUST BE A *STRUCTURE* YOU WISH TO PREPEND'

  TNAMES1 = TAG_NAMES(st1)
  TNAMES2 = TAG_NAMES(st2)

  IF (N_ELEMENTS(TNAMES1) NE N_ELEMENTS(TNAMES2)) THEN MESSAGE,'BOTH STRUCTURES MUST HAVE SAME SIZE!'

  FOR J=0L,N_ELEMENTS(TNAMES1)-1L DO BEGIN
    IF (STRCMP(TNAMES1[J],TNAMES2[J],/FOLD_CASE) NE 1) THEN MESSAGE, 'FIELD NAMES ARE NOT THE SAME, NOR IN THE SAME ORDER'
    ;IF (SIZE(st1.(J),/TYPE) NE SIZE(st2.(J),/TYPE)) THEN PRINT, j;MESSAGE,' FIELD NAMES ARE NOT THE SAME TYPE'
  ENDFOR

  CASE 1 OF ; switched from case 0 in April 2016
    0: BEGIN
      ; RJW: THIS COULD BE FASTER!!!
      DATA = CREATE_STRUCT(TNAMES1[0], [st1.(0), st2.(0)])
      FOR J=1L,N_ELEMENTS(TNAMES1)-1 DO BEGIN ; if only 1 value, this for loop does nothing
        DATA = CREATE_STRUCT(TEMPORARY(DATA), TNAMES1[J], [st1.(J), st2.(J)]) ; THIS CAN BE MADE FASTER WITH ONLY ONE CREATE_STRUCT COMMAND
        ; IF STRING and the same (i.e. EXPLAIN_DATA*) then don't double up?, just use one...
      ENDFOR
    END

    1:BEGIN
      st_string = 'DATA=CREATE_STRUCT('
      nTNAMESm1 = N_ELEMENTS(TNAMES1)-1L
      FOR J=0L,nTNAMESm1 DO BEGIN
        IF J NE 0 THEN st_string = st_string + ',' ; add in commas, except on first one.
        jnumstr = STRTRIM(STRING(J),2)
        IF (STREGEX(TNAMES1[J],'^EXPLAIN_',/BOOLEAN)) OR (STREGEX(TNAMES1[J],'_EXPLAIN$',/BOOLEAN)) THEN BEGIN
          IF (N_ELEMENTS(st1.(J)) EQ 1) AND (N_ELEMENTS(st2.(J)) EQ 1) AND (ISA(st1.(J),/STRING)) AND (ISA(st2.(J),/STRING)) THEN BEGIN
            IF STRCMP(st1.(J), st2.(J)) EQ 1 THEN BEGIN
              st_string = st_string + ''''+TNAMES1[J]+''',st1.('+jnumstr+')'
              CONTINUE ; skip to next iteration
            ENDIF
            ; If not the same, fall through and add.
          ENDIF ELSE BEGIN
            MESSAGE,'ERROR: If using jade_st_join and there is an EXPLAIN_ field it has to be a string of size 1'
          ENDELSE
        ENDIF

        st_string = st_string + ''''+TNAMES1[J]+''',[st1.('+jnumstr+'), st2.('+jnumstr+')]'
      ENDFOR
      st_string = st_string + ')'
      ;      print,st_string
      st_result = execute(st_string)
      IF st_result NE 1 THEN MESSAGE,'ERROR: new structure failed to be made... ask Rob W. if the above error is not useful.'
    END
  ENDCASE

  RETURN, DATA
END

FUNCTION CSV_COMBINE, Input_Files, Output_File
  ; Combines csvs from the VIPER 2.0 Program
  ; Read in a string array of CSV_Filenames to combine
  ; Outputs a csv file of appropriate name

  nel = N_ELEMENTS(Input_Files)
  nloop = 0
  FOR i = 0, nel-1 DO BEGIN
    CSV_Name = Input_Files(i)
    Input_CSV = READ_CSV(CSV_Name, HEADER = Header)
    nloop = nloop + 1
    IF nloop EQ 1 THEN BEGIN
      Out_File = Input_CSV
    ENDIF ELSE BEGIN
      Out_File = VIPER_ST_JOIN(out_file, input_csv)
    ENDELSE
  ENDFOR

  WRITE_CSV, Output_File, Out_File, HEADER = Header
  RETURN, Output_File
END

FUNCTION Convert_Kaleb_Delamere_Ratios
  ; Convert what Kaleb gave me
  
  ; Lowest Value is 6.002 Rj in Odouble file
  ; Highest Value is 8.970 Rj in Splus file
  
  ; Need to go from 6.002 to 8.970 which is 2989 points
  
  R_steps = MAKE_ARRAY(2969, /DOUBLE)
  
  FOR i = 0, 2968 DO BEGIN
    R_Steps[i] = 6.002d + 0.001d*DOUBLE(i)
  ENDFOR
  
  ; Read in Kaleb's files
  ; S+++
  ; 31x2
  OPENR, 1, 'stripple.txt'
  stplus_file = DBLARR(2,31)
  READF, 1, stplus_file
  CLOSE, 1

  ; S++
  ; 24x2
  OPENR, 1, 'sdouble.txt'
  sdplus_file = DBLARR(2,24)
  READF, 1, sdplus_file
  CLOSE, 1

  ; S+
  ; 31x2
  OPENR, 1, 'splus.txt'
  splus_file = DBLARR(2,31)
  READF, 1, splus_file
  CLOSE, 1

  ; O++
  ; 38x2
  OPENR, 1, 'Odouble.txt'
  odplus_file = DBLARR(2,38)
  READF, 1, odplus_file
  CLOSE, 1
  
  ; O+
  ; 26x2
  OPENR, 1, 'oplus.txt'
  oplus_file = DBLARR(2,26)
  READF, 1, oplus_file
  CLOSE, 1

  ; Now interpolate all of these properly so that they have the same ranges of Rs and stuff

  x1 = INTERPOL(stplus_file[1,*], stplus_file[0,*], R_steps)
  x2 = INTERPOL(sdplus_file[1,*], sdplus_file[0,*], R_steps)
  x3 = INTERPOL(splus_file[1,*], splus_file[0,*], R_steps)
  x4 = INTERPOL(odplus_file[1,*], odplus_file[0,*], R_steps)
  x5 = INTERPOL(oplus_file[1,*], oplus_file[0,*], R_steps)
  
  ; The above are ratios of that density over n_e so now I just need to divide all of them by x5 to get them in the proper format
  x1 = x1/x5
  x2 = x2/x5
  x3 = x3/x5
  x4 = x4/x5
  x5 = x5/x5
  
  ; Reorder them so they are similar to the other file
  ; Rj is 0, O+/O+ is 1 - I will now go from lowest Mass/Charge to highest so O++, S+++, S++, S+
  Final_Out = MAKE_ARRAY(6,2969,/DOUBLE)
  Final_Out[0,*] = R_Steps[*]
  Final_Out[1,*] = x5[*]
  Final_Out[2,*] = x4[*]
  Final_Out[3,*] = x1[*]
  Final_Out[4,*] = x2[*]
  Final_Out[5,*] = x3[*]
  
  ; Change directories for file saving
  CD, CURRENT = RESTORE_DIR
  STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/Data/Jupiter/Chemistry_Ratios'))
  CD, STR_DIR
  OPENW, 1, "Delamere_Ratios.txt"
  PRINTF, 1, Final_Out
  CLOSE, 1
  CD, RESTORE_DIR
  
  STOP
  
  RETURN, 0
END

FUNCTION VPHI_PROFILE, R, Vphi_Eqn_no, Structure
  ; Read in an R distance and return the vphi for the location
  ; Read in the equation number as well - Currently only have one function
  ;     1: The hyperbolic tangent function - a*tanh(bx)
  ;     2: Corotation
  ;     3: Binned function from fits to V1 data
  
  IF Structure.Planet_Number EQ 5 THEN BEGIN
  R = DOUBLE(R)
    IF Vphi_Eqn_no EQ 1 THEN Vphi = 262.5d * TANH(0.05d * R)
    IF Vphi_Eqn_no EQ 2 THEN Vphi = 12.572 * R
    IF Vphi_Eqn_no EQ 3 THEN BEGIN
      ; Currently valid for the region of 19 to 22 R, using a linear fit to that data from the Voyager 1 fits
      ; I read in the data for Voyager 1 fits, extracted data from 19-22 R and then did a linfit on the data to get these numbers
      Vphi = -127.915d + 16.3563d * R
    ENDIF
;;--------- GX: Probably unnecessary?
;  ENDIF ELSE IF (Structure.Planet_Number EQ 7) AND (Vphi_Eqn_no EQ 2) THEN BEGIN
;    Vphi = 12.572 * R
;;--------- TO HERE
  ENDIF ELSE BEGIN
    Vphi = 0
    MESSAGE,'ERROR: Vphi is available fully available only for Jupiter. For other planets only the corrotation is (currently) available.'
  ENDELSE

  RETURN, Vphi
  
END


FUNCTION VGR_CUP_FLOW_DIRECTIONS, Spacecraft, Structure
  ; Inputs
  ; Structure.Planet_Number
  ; Spacecraft - Read the spacecraft number, either a 1 or 2
  ;
  ; Hardcoded: number of SEDR datapoints
  ; 
  ; 
  ; Outputs
  ; Text array that contains the following
  ;   Decimal Date, Rj, Vx_Acup, Vy_Acup ... Vz_Dcup
  ;   
  ; ----------
  ; JUPITER
  ; ----------
  IF Structure.Planet_Number EQ 5 THEN BEGIN
    ; For the desired spacecraft, determine all times that are necessary
    IF Spacecraft EQ 1 THEN BEGIN
      n = 2925 ; Number of Data points for analysis
      COMMON VGR1DATA
      Data = VGR1SSEDR ; Read in SSEDR
    ENDIF ELSE IF Spacecraft EQ 2 THEN BEGIN
      n = 15653 ; Number of Data points for analysis
      COMMON VGR2DATA
      Data = VGR2SSEDR ; Read in SSEDR
    ENDIF
    
    Decimal_Date = fltarr(n)
    FOR i = 0, n-1 DO BEGIN
      Decimal_Date[i] = (Data[2,i]*3600. + Data[3,i]*60. + Data[4,i] + Data[5,i] / 1000. ) / 86400. + Data[1,i]
    ENDFOR
    
    xSC_S3 = Data[6,*] ; System III cartesian coordinates
    ySC_S3 = Data[7,*]
    zSC_S3 = Data[8,*]
  
    ; The VelocityCylindrical is the velocity of the plasma in r, phi, and z (Cylindrical coordinates)
    ; This is not in the corotating frame of System III, so the plasma flow should be subtracted
    ; from the number to get the appropriate value
    ;
    ; The corotational flow goes as 12.572 km/s * Radius (Jupiter radii)
  
    RJupiter = 71492.0 ; Radius of Jupiter in km
    Radius   = SQRT(xSC_S3*xSC_S3 + ySC_S3*ySC_S3) ; Radius is from axis of rotation
    Vphi     = Radius/Rjupiter*12.572
    cups = ['A','B','C','D']
    Output_Array = MAKE_ARRAY(14,n)
    
    ; Fill in array with time and radial values
    Output_Array(0,*) = Decimal_Date(*)
    Output_Array(1,*) = Radius(*)/Rjupiter
    
    FOR i = 0, n-1 DO BEGIN
      imod = i MOD 100
      IF imod EQ 1 THEN PRINT, i
      
      ; Create a structure format that the function understands
      Structure = {Year: 1979, DOY: Data[1,i], Hour: Data[2,i], Minute: Data[3,i], Second: Data[4,i], Spacecraft: Spacecraft, Planet_Number: 5, Response: 0, $
                   L_or_M_Mode: 1, Save_Plot: 0,  Fit: 0, Iterations: 25, Channels_CupA: 001128 , Channels_CupB: 001128, $
                   Channels_CupC: 001128, Channels_CupD: 001128, Vary_V1: 0, Vary_V2: 0, Vary_V3: 0, Vary_Density: 0, $
                   Vary_Temperature: 0, Analytical_Vphi: 0d, V1: 0d, V2: Vphi[i], V3: 0d, Common_Temperature: 0d, Delamere_Composition: 0d, $
                   Species_1: 0, Species_1_A: 1d, Species_1_Z: 1d, Species_1_n: 50d, Species_1_T: 3d, $
                   Species_2: 0, Species_2_A: 1d, Species_2_Z: 1d, Species_2_n: 50d, Species_2_T: 3d, $
                   Species_3: 0, Species_3_A: 1d, Species_3_Z: 1d, Species_3_n: 50d, Species_3_T: 3d, $
                   Species_4: 0, Species_4_A: 1d, Species_4_Z: 1d, Species_4_n: 50d, Species_4_T: 3d, $
                   Species_5: 0, Species_5_A: 1d, Species_5_Z: 1d, Species_5_n: 50d, Species_5_T: 3d, $
                   Species_6: 0, Species_6_A: 1d, Species_6_Z: 1d, Species_6_n: 50d, Species_6_T: 3d, $
                   Species_7: 0, Species_7_A: 1d, Species_7_Z: 1d, Species_7_n: 50d, Species_7_T: 3d}
      VelocityIn = [Structure.V1, Structure.V2, Structure.V3]
      
      ; Call on appropriate function
      FOR j = 0, 3 DO BEGIN
        U = VGR_CupVelocity(VelocityIn, cups[j], Structure)
        ; Fill in array values
        Output_Array(2 + j*3    ,i) = U[0]
        Output_Array(2 + j*3 + 1,i) = U[1]
        Output_Array(2 + j*3 + 2,i) = U[2]
        
      ENDFOR
      
    ENDFOR
    
  ; ----------
  ; URANUS
  ; ----------
  ENDIF ELSE IF Structure.Planet_Number EQ 7 THEN BEGIN
    ; MESSAGE,'Function VGR_CUP_FLOW_DIRECTIONS is used in VIPER processing for Uranus' ; Logan 2025_04_09 - Commented out

      n = 10274 ; Number of Data points for analysis
      COMMON VGR2DATA_URANUS
      Data = VGR2SSEDR_URANUS ; Read in SSEDR
    
    Decimal_Date = fltarr(n)
    FOR i = 0, n-1 DO BEGIN
      Decimal_Date[i] = (Data[2,i]*3600. + Data[3,i]*60. + Data[4,i] + Data[5,i] / 1000. ) / 86400. + Data[1,i]
    ENDFOR
    
    xSC_URANUS_FIT = Data[6,*] ; URANUS_FIT
    ySC_URANUS_FIT = Data[7,*]
    zSC_URANUS_FIT = Data[8,*]
  
    ; The VelocityCylindrical is the velocity of the plasma in r, phi, and z (Cylindrical coordinates)
    ; This is not in the corotating frame of URANUS_FIT, so the plasma flow should be subtracted
    ; from the number to get the appropriate value
    ;
    ; The corotational flow goes as -2.586 km/s * Radius (Uranus radii)
    
    RUranus = 25559d ; Radius of Uranus [km]
    u_Uranus = -2.586  ; Rotation speed of Uranus [km/s]
    Radius   = SQRT(xSC_URANUS_FIT*xSC_URANUS_FIT + ySC_URANUS_FIT*ySC_URANUS_FIT) ; Radius is from axis of rotation
    Vphi     = Radius/RUranus*u_Uranus
    cups = ['A','B','C','D']
    Output_Array = MAKE_ARRAY(14,n)
    
    ; Fill in array with time and radial values
    Output_Array(0,*) = Decimal_Date(*)
    Output_Array(1,*) = Radius(*)/RUranus
    
    FOR i = 0, n-1 DO BEGIN
      imod = i MOD 100
      IF imod EQ 1 THEN PRINT, i
      
      ; Create a structure format that the function understands
      Structure = {Year: 1979, DOY: Data[1,i], Hour: Data[2,i], Minute: Data[3,i], Second: Data[4,i], Spacecraft: Spacecraft, Planet_Number: 7, Response: 0, $
                   L_or_M_Mode: 1, Save_Plot: 0,  Fit: 0, Iterations: 25, Channels_CupA: 001128 , Channels_CupB: 001128, $
                   Channels_CupC: 001128, Channels_CupD: 001128, Vary_V1: 0, Vary_V2: 0, Vary_V3: 0, Vary_Density: 0, $
                   Vary_Temperature: 0, Analytical_Vphi: 0d, V1: 0d, V2: Vphi[i], V3: 0d, Common_Temperature: 0d, Delamere_Composition: 0d, $
                   Species_1: 0, Species_1_A: 1d, Species_1_Z: 1d, Species_1_n: 50d, Species_1_T: 3d, $
                   Species_2: 0, Species_2_A: 1d, Species_2_Z: 1d, Species_2_n: 50d, Species_2_T: 3d, $
                   Species_3: 0, Species_3_A: 1d, Species_3_Z: 1d, Species_3_n: 50d, Species_3_T: 3d, $
                   Species_4: 0, Species_4_A: 1d, Species_4_Z: 1d, Species_4_n: 50d, Species_4_T: 3d, $
                   Species_5: 0, Species_5_A: 1d, Species_5_Z: 1d, Species_5_n: 50d, Species_5_T: 3d, $
                   Species_6: 0, Species_6_A: 1d, Species_6_Z: 1d, Species_6_n: 50d, Species_6_T: 3d, $
                   Species_7: 0, Species_7_A: 1d, Species_7_Z: 1d, Species_7_n: 50d, Species_7_T: 3d}
      VelocityIn = [Structure.V1, Structure.V2, Structure.V3]
      
      ; Call on appropriate function
      FOR j = 0, 3 DO BEGIN
        U = VGR_CupVelocity(VelocityIn, cups[j], Structure)
        ; Fill in array values
        Output_Array(2 + j*3    ,i) = U[0]
        Output_Array(2 + j*3 + 1,i) = U[1]
        Output_Array(2 + j*3 + 2,i) = U[2]
        
      ENDFOR
    ENDFOR
    
  ; ----------
  ; NEPTUNE
  ; ----------
  ENDIF ELSE IF Structure.Planet_Number EQ 8 THEN BEGIN

      n = 28325 ; Number of SEDR Data points for analysis
      COMMON VGR2DATA_NEPTUNE
      Data = VGR2SSEDR_NEPTUNE ; Read in SSEDR
    
    Decimal_Date = fltarr(n)
    FOR i = 0, n-1 DO BEGIN
      Decimal_Date[i] = (Data[2,i]*3600. + Data[3,i]*60. + Data[4,i] + Data[5,i] / 1000. ) / 86400. + Data[1,i]
    ENDFOR
    
    xSC_V2_NEPTUNE_PLS = Data[6,*]
    ySC_V2_NEPTUNE_PLS = Data[7,*]
    zSC_V2_NEPTUNE_PLS = Data[8,*]
  
    ; The VelocityCylindrical is the velocity of the plasma in r, phi, and z (Cylindrical coordinates)
    ; This is not in the corotating frame of V2_NEPTUNE_PLS, so the plasma flow should be subtracted
    ; from the number to get the appropriate value
    ;
    ; The corotational flow goes as 2.683 km/s * Radius (Neptune radii)
    RNeptune = 24764d ; Radius of Neptune [km]
    u_Neptune = 2.683d  ; Rotation speed of Neptune [km/s]
    Radius   = SQRT(xSC_V2_NEPTUNE_PLS*xSC_V2_NEPTUNE_PLS+ ySC_V2_NEPTUNE_PLS*ySC_V2_NEPTUNE_PLS) ; Radius is from axis of rotation
    Vphi     = Radius/RNeptune*u_Neptune
    cups = ['A','B','C','D']
    Output_Array = MAKE_ARRAY(14,n)
    
    ; Fill in array with time and radial values
    Output_Array(0,*) = Decimal_Date(*)
    Output_Array(1,*) = Radius(*)/RNeptune
    
    FOR i = 0, n-1 DO BEGIN
      imod = i MOD 100
      IF imod EQ 1 THEN PRINT, i
      
      ; Create a structure format that the function understands
      ; GX: MAYBE UNNECESSARY HERE??
      Structure = {Year: 1979, DOY: Data[1,i], Hour: Data[2,i], Minute: Data[3,i], Second: Data[4,i], Spacecraft: Spacecraft, Planet_Number: 8, Response: 0, $
                   L_or_M_Mode: 1, Save_Plot: 0,  Fit: 0, Iterations: 25, Channels_CupA: 001128 , Channels_CupB: 001128, $
                   Channels_CupC: 001128, Channels_CupD: 001128, Vary_V1: 0, Vary_V2: 0, Vary_V3: 0, Vary_Density: 0, $
                   Vary_Temperature: 0, Analytical_Vphi: 0d, V1: 0d, V2: Vphi[i], V3: 0d, Common_Temperature: 0d, Delamere_Composition: 0d, $
                   Species_1: 0, Species_1_A: 1d, Species_1_Z: 1d, Species_1_n: 50d, Species_1_T: 3d, $
                   Species_2: 0, Species_2_A: 1d, Species_2_Z: 1d, Species_2_n: 50d, Species_2_T: 3d, $
                   Species_3: 0, Species_3_A: 1d, Species_3_Z: 1d, Species_3_n: 50d, Species_3_T: 3d, $
                   Species_4: 0, Species_4_A: 1d, Species_4_Z: 1d, Species_4_n: 50d, Species_4_T: 3d, $
                   Species_5: 0, Species_5_A: 1d, Species_5_Z: 1d, Species_5_n: 50d, Species_5_T: 3d, $
                   Species_6: 0, Species_6_A: 1d, Species_6_Z: 1d, Species_6_n: 50d, Species_6_T: 3d, $
                   Species_7: 0, Species_7_A: 1d, Species_7_Z: 1d, Species_7_n: 50d, Species_7_T: 3d}
      VelocityIn = [Structure.V1, Structure.V2, Structure.V3]
      
      ; Call on appropriate function
      FOR j = 0, 3 DO BEGIN
        U = VGR_CupVelocity(VelocityIn, cups[j], Structure)
        ; Fill in array values
        Output_Array(2 + j*3    ,i) = U[0]
        Output_Array(2 + j*3 + 1,i) = U[1]
        Output_Array(2 + j*3 + 2,i) = U[2]
        
      ENDFOR
      
    ENDFOR
    
    ENDIF

 
  RETURN, Output_Array
END

FUNCTION VGR_CUP_RESPONSE_PLOT, Spacecraft
  ;
  ; NOTE!!!
  ; As of June 8, 2016, this function is no longer used. Instead a Matlab function of the name "Cup_Responses.m"
  ; is used so that the plots look better
  ;
  ; Inputs
  ; Spacecraft
  ; 
  ; Outputs
  ; Plots a graphical window to the screen with the response of the instrument to a cold rotational beam for each Faraday cup
  
  ; Read in the data for specified spacecraft
  CD, CURRENT = RESTORE_DIR
  STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/Data/Instrument_Constants'))
  CD, STR_DIR
  
  IF Spacecraft EQ 1 THEN BEGIN
    n = 2925
    a = MAKE_ARRAY(14,n)
    OPENR, 1, 'V1_Cup_Responses.txt'
    READF, 1, a
    CLOSE, 1
  ENDIF ELSE IF Spacecraft EQ 2 THEN BEGIN
    n = 15653
    a = MAKE_ARRAY(14,n)
    OPENR, 1, 'V2_Cup_Responses.txt'
    READF, 1, a
    CLOSE, 1
  ENDIF
  
  ; Get the proper response for the cold beam
  ; The response is defined as follows:
  ;
  ; R = Vz/V
  ; 
  ; Therefore it is 1 if all of the flow is directly into the cup and 0 if none is directly in
  
  ; Limit the range appropriately
  If Spacecraft EQ 2 THEN BEGIN
    ind = WHERE(a(0,*) GT 186 AND a(0,*) LT 192)
    a = a[*,ind]
  ENDIF
  
  Time   = a(0,*)
  Radius = a(1,*)
  CupA   = a(4,*)/(SQRT(a(2,*)^2.0 + a(3,*)^2.0 + a(4,*)^2.0))
  CupB   = a(7,*)/(SQRT(a(5,*)^2.0 + a(6,*)^2.0 + a(7,*)^2.0))
  CupC   = a(10,*)/(SQRT(a(8,*)^2.0 + a(9,*)^2.0 + a(10,*)^2.0))
  CupD   = a(13,*)/(SQRT(a(11,*)^2.0 + a(12,*)^2.0 + a(13,*)^2.0))
  

  
  ; Begin plotting
  IF Spacecraft EQ 1 THEN title = 'Voyager 1 Cup Response to a Cold Corotating Beam' ELSE $
  IF Spacecraft EQ 2 THEN title = 'Voyager 2 Cup Response to a Cold Corotating Beam'
  
  xtitle = 'DOY 1979'
  ytitle = '$R = V_z/|V|$'
  ; xrange = [183,192] ; Range of V2 data that seems prevalent
  
  p1 = plot(Time, CupA, linestyle = '', color = 'red', name = 'Cup A', symbol='d', yrange = [0,1], title = title, xtitle = xtitle, ytitle = ytitle)
  p2 = plot(Time, CupB, linestyle = '', color = 'blue', name = 'Cup B', symbol='d', /OVERPLOT)
  p3 = plot(Time, CupC, linestyle = '', color = 'black', name = 'Cup C', symbol='d', /OVERPLOT)
  p4 = plot(Time, CupD, linestyle = '', color = 'green', name = 'Cup D', symbol='d', /OVERPLOT)
  leg = legend(TARGET=[p1,p2,p3,p4],/NORMAL,POSITION=[0.8,0.8],/AUTO_TEXT_COLOR)
  
  ; Add radial axis
  ; Create tick values
  FOR i = 0, N_ELEMENTS(Time)-1 DO BEGIN
    IF Radius[i] MOD 5 LT 0.1 THEN BEGIN
      value = ROUND(Radius[i])
      IF ISA(indval) EQ 0 THEN BEGIN
        radii = value
        indval = i
      ENDIF ELSE BEGIN
        IF radii[N_ELEMENTS(radii)-1] EQ value THEN CONTINUE
        radii = [radii,value]
        indval = [indval,i]
      ENDELSE
    ENDIF
  ENDFOR  
  
  tickname = STRTRIM(ROUND(Radius[indval[0]]))
  tickvalue = Time[indval[0]]
  FOR k = 1, N_ELEMENTS(indval)-1 DO BEGIN
    ;ind = WHERE(Radius 
    tickname = [tickname,STRTRIM(ROUND(Radius[indval[k]]))]
    tickvalue = [tickvalue,Time[indval[k]]]
  ENDFOR
  xaxis = axis('X', LOCATION = [-0.1], title = '$Radius (R_J)$', tickname = tickname, tickvalue = tickvalue, minor = 4)
  
  ; Return to main directory
  CD, RESTORE_DIR
  
  RETURN, 0
END

FUNCTION YEAR_DOY_to_UTC_String, csv, outputname
  ; Read in an array of Year, DOY, Hour, Minute, Second and return YYYY-DOYTHH:MM:SS.SSS (For VIPER code for Khurana model)
  ; Make sure that the csv file is in the right subdirectory
  ; 
  ; Inputs
  ; csv - csv name of file with the structure that we use - This code is more of a one off code
  ; 
  ; Outputs
  ; IDL save file of strings of the time to the output filename
  
  a = READ_CSV(csv, HEADER = Structure_Names)
  
  y = a.Field01
  D = a.Field02
  H = a.Field03
  M = a.Field04
  S = a.Field05
  
  str = MAKE_ARRAY(N_ELEMENTS(y),/STRING)
  
  FOR i = 0, N_ELEMENTS(y) - 1 DO BEGIN
    str(i) = STRCOMPRESS(STRING(y(i),'-',D(i),'T',H(i),':',M(i),':',S(i),'.000', $
                                FORMAT = '(A,A,I03,A,I02,A,I02,A,I02,A)'), /REMOVE_ALL)    
  ENDFOR
  
  SAVE, str, FILENAME = outputname ; Create a save file
  STOP
END  
  
FUNCTION FIX_DELAMERE_RATIOS, csv, outputcsv
  ; In transferring over to a new csv format in January of 2016, some of the tags didn't perfectly
  ; get translated. This function just fixes the delamere column to properly reflect where delamere
  ; composition is being used.
  
  a = READ_CSV(csv, HEADER = Structure_Names)
  b = a ; Structure for output
  
  ; For Delamere the ratios should be as follows:
  ; Column AJ: 0.37443411 (36)
  ; Column AO: 0.3696681  (41)
  ; Column AT: 0.83713514 (46)
  ; Column AY: 0.12761205 (51)
  
  n1 = 0.37443411
  n2 = 0.3696681
  n3 = 0.83713514
  n4 = 0.12761205
  
  FOR i = 0, N_ELEMENTS(a.Field01)-1 DO BEGIN
    
    c1 = 0
    c2 = 0
    c3 = 0
    c4 = 0
    
    ; If already set to Delamere composition, it is correct and continue
    IF a.Field27[i] EQ 1 THEN CONTINUE
    
    ; Check if the other columns are within 0.1% of correct values - Don't check perfect equality due to rounding errors from read/write
    IF a.Field36[i] LT n1*1.001 AND a.Field36[i] GT n1*0.999 THEN c1 = 1
    IF a.Field41[i] LT n2*1.001 AND a.Field41[i] GT n2*0.999 THEN c2 = 1
    IF a.Field46[i] LT n3*1.001 AND a.Field46[i] GT n3*0.999 THEN c3 = 1
    IF a.Field51[i] LT n4*1.001 AND a.Field51[i] GT n4*0.999 THEN c4 = 1
    
    ; If all species follow Delamere, change to Delamere tag
    IF (c1 EQ 1) AND (c2 EQ 1) AND (c3 EQ 1) AND (c4 EQ 1) THEN BEGIN
      b.Field27[i] = 1
      PRINT, i+2 ; Check to see what indices were changed
    ENDIF
    
  ENDFOR
  
  WRITE_CSV, outputcsv, b, HEADER = Structure_Names
  
  RETURN, 0
END 
  
FUNCTION EMPIRICAL_TE_AND_TEHOT

; Read in Fran's Tabulated values for ease of calculation and get an empirical fit to the data of some form

Rj   = [5.7d,6d,7d,7.5d,8d,10d,13d,16d]
Te   = [2.8d,5d,7d,10d,20d,30d,40d,50d]
Teh  = [20d,60d,85d,120d,240d,360d,480d,600d]
nhnc   = [0.0001d,0.001d,0.01d,0.1d] ; Nans don't work here for some reason
Rjnhnc = [5.7d,6d,7.5d,10d]

; For the power law fits of the form as follows:
; Te = exp(b)*exp(m*alog(Rj))
; Te and Rj
; b = -3.1660668
; m = 2.6848884
;
; Teh and Rj
; b = -1.2250957
; m = 2.9082436
; 
; nhnc and Rjnhnc
; b = -27.672004
; m = 11.167562


stop

; Testing fits
p1 = plot(rjnhnc,nhnc,/xlog,/ylog,linestyle='',symbol='+')


pari = linfit(alog(rjnhnc),alog(nhnc))

dumx = DINDGEN(100)+1d
y    = exp(pari[0])*exp(pari[1]*alog(dumx))
;y    = par[0]+par[1]*dumx
op1 = plot(dumx,y,/OVERPLOT,color='red')

END 



FUNCTION PLS_DATA_READ, Structure, ELECTRON_FLUX = ELECTRON_FLUX
;------------------ GX: CHANGE STUFF HERE (add the next planet)
; Dougherty 2018 - Added keyword ELECTRON_FLUX to just return data for the electron fluxes
; Dougherty 2020 - Added Structure.Planet_Number NE 5 - Files not at Jupiter may have a different data format
; Xystouris 2024 - Added Uranus and Neptune files format
;   
  CASE Structure.Planet_Number OF
    5 : BEGIN
          CD, CURRENT = RESTORE_DIR
          STR_DIR = RESTORE_DIR + PATH_SEP() + 'Data' + PATH_SEP() + 'Jupiter' + PATH_SEP() + 'PLS_Data'
          ;STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/Data/Jupiter/PLS_Data')) ; we removed that to make the script friendly to both MAC and PC
        END
    6 : Data_Filename = 'Saturn_PLS_.csv'
    7 : BEGIN
          CD, CURRENT = RESTORE_DIR
          STR_DIR = RESTORE_DIR + PATH_SEP() + 'Data' + PATH_SEP() + 'Uranus' + PATH_SEP() + 'PLS_Data'
        END
    8 : BEGIN
          CD, CURRENT = RESTORE_DIR
          STR_DIR = RESTORE_DIR + PATH_SEP() + 'Data' + PATH_SEP() + 'Neptune' + PATH_SEP() + 'PLS_Data'
        END
  ENDCASE
;------------------ TO HERE

  ; Determine DOY and Time
  DOY = Structure.DOY
  Time = Structure.DOY + (Structure.Hour + Structure.Minute/60d + Structure.Second/3600d)/24d
    
  ; Saves information in Common Blocks to avoid excessive memory leaks from reading in the file every CSV row
  ; Jupiter files were saved as 'day___.txt' from a Fortran code and require this step
  ; GX: for Uranus and Neptune we used the same format
  COMMON DOY_COMMON_BLOCK
  IF DOY GT DOY_COMMON_BLOCK_Var THEN BEGIN
    COMMON DOY_COMMON_BLOCK, DOY_COMMON_BLOCK_Var
    DOY_COMMON_BLOCK_Var = DOY
    Data_Filename = STRCOMPRESS(STRING('day', LONG(DOY), '.txt'), /REMOVE_ALL) ; Name of text file that contains pertinent data
  ENDIF ELSE BEGIN
    Data_Filename = -1
  ENDELSE
    
  IF KEYWORD_SET(ELECTRON_FLUX) THEN E_Flux = 1 ELSE E_Flux = 0
    Old_Format = 1
  FORWARD_FUNCTION PLS_Current_Data
  RETURN, PLS_Current_Data(Structure, Data_Filename, Time, 'D', STR_DIR, E_Flux) ; Cup argument 'D' exists just to run code
END


FUNCTION PLS_Current_Data, Structure, Data_Filename, Time, Cup, STR_DIR, E_Flux, Old_Format    ; !!!!! GX: CAREFUL! I DON'T KNOW HWAT THE "OLD FORMAT" IS
  ; Inputs
  ; Structure     - General VIPER structure
  ; Data_Filename - Name of data file ('day___.txt' for Jupiter files and 'Planet_PLS.csv' for Saturn, Uranus, and Neptune)
  ; Time          - The decimal time of the measurement required
  ; Cup           - Which Cup you want to extract data for
  ; STR_DIR       - A directory where the data file is stored
  ; E_Flux        - Needed for calculating the E-Flux. This is currently only used for Jupiter files and has not been fully tested with S,U,N
  ; Old_Format    - Either 0 or 1 for the Jupiter format (1) or for the new condensed CSV (1)

  FORWARD_FUNCTION read_file_Jupiter
  FORWARD_FUNCTION read_file_Uranus
  FORWARD_FUNCTION read_file_Neptune

  ; Store the file into a common block to prevent it from being read over and over again. This prevents a memory leak

  IF ISA(Old_Format) EQ 0 THEN Old_Format = 1 ; Just a check for the initial setup
  
  IF Old_Format EQ 1 THEN BEGIN
  
    ; In the old format, the data is saved into common blocks
    IF ISA(Data_Filename, /STRING) EQ 1 THEN BEGIN
      COMMON DATA_FILE_COMMON_BLOCK, DATA_FILE_COMMON_BLOCK_Var ; Redefine the data to be of that day
      ; Each planet has a different format of data
      ;STOP
      CASE Structure.Planet_Number OF
        5 : data = read_file_Jupiter(Data_Filename, CHANGE_DIRECTORY=STR_DIR)
 ;      6 : 
        7 : data = read_file_Uranus(Data_Filename, CHANGE_DIRECTORY=STR_DIR)
        8 : data = read_file_Neptune(Data_Filename, CHANGE_DIRECTORY=STR_DIR)
      ENDCASE
      DATA_FILE_COMMON_BLOCK_Var = data ; Store the data if reading from a new data file
    ENDIF ELSE BEGIN
      data = DATA_FILE_COMMON_BLOCK_Var ; Pull from stored data if not using a new data file
    ENDELSE

  ENDIF ELSE BEGIN ; Format for S, U, and N files
    CD, STR_DIR, CURRENT = RESTORE_DIR
    ; 
    ; GX: I think this is not used, as we have the case in the code just above?
    ; 
    ; Read in all of the appropiate files for L, M, and E modes so that the format of returned structure matches with Jupiter files
    ; Difference for Uranus: the E are passed as ONE record (they are already combined - this was done due to different Uranus operating mode)
    ; Neptune follows same as Uranus (data are with same mode)
    Data_Filename = STRJOIN([Data_Filename,'L.csv'])
    Data_00 = READ_CSV(Data_Filename)
    Data_Filename = STRJOIN([Data_Filename,'M.csv'])
    Data_01 = READ_CSV(Data_Filename)
    Data_Filename = STRJOIN([Data_Filename,'E.csv'])
    Data_02 = READ_CSV(Data_Filename)


    ; May need to unpack this data into an appropriate format for the following bit of code
    CD, RESTORE_DIR
  ENDELSE


  ; This is the same for both versions (Jupiter/SUN) of the code as it is just channel numbers
  Mode = Structure.L_or_M_Mode
  ; make energy arrays for M modes, L modes, and the E modes:
  ; get array of energies using the conversionM function
  FORWARD_FUNCTION conversionM
  ; There are 128 channels for M modes - each corresponds to an energy bin
  MchannelsLow = findgen(128)             ; lower bound on bins (channel number 0-127)
  MchannelsHigh = 1. + findgen(128)       ; upper bound on bins (channel number 1-128)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  Menergies = ( conversionM(MchannelsLow) + conversionM(MchannelsHigh) ) / 2.   ; convert channel numbers to energy values for M modes

  ; get array of energies using the conversionLE2 function
  FORWARD_FUNCTION conversionLE2
  ; There are 16 channels for L modes
  LchannelsLow = findgen(16)              ; lower bound on bins (channel number 0-15)
  LchannelsHigh = 1. + findgen(16)        ; upper bound on bins (channel number 1-16)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  Lenergies = ( conversionLE2(LchannelsLow) + conversionLE2(LchannelsHigh) ) / 2.    ; convert channel numbers to energy values for L modes

  ; make energy arrays for both E modes using conversion functions
  FORWARD_FUNCTION conversionLE2
  FORWARD_FUNCTION conversionE1
  ; There are 16 channels for E1 modes
  E1channelsLow = findgen(16)             ; lower bound on bins (0-15)
  E1channelsHigh = 1. + findgen(16)       ; upper bound on bins (1-16)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  E1energies = ( conversionE1(E1channelsLow) + conversionE1(E1channelsHigh) ) / 2.  ; convert channel numbers to energy values for E1 modes
  ; There are 16 channels for E2 modes, but the first 4 are garbage, so use only the last 12
  E2channelsLow = 4. + findgen(12)        ; lower bound on bins (4-15)
  E2channelsHigh = 5. + findgen(12)       ; upper bound on bins (5-16)
  ; use the average of lower energy and upper energy of each bin as the bin's energy level
  E2energies = ( conversionLE2(E2channelsLow) + conversionLE2(E2channelsHigh) ) / 2. ; convert channel numbers to energy values for E2 modes
  ; now, combine E1 and E2 energies into one big energy array
  Eenergies = fltarr(28)                   ; make one array for all energies
  Eenergies[0:15] = E1energies             ; first half of energy array is the E1 energy levels
  Eenergies[16:27] = E2energies            ; other half of energy array is the E2 energies

  FORWARD_FUNCTION convert_times
  ; use mode chosen by user to define the important quantities for the plots (currents, energies, and times)
  ; 'Day___.txt' Jupiter files are in a vastly different format
  IF (mode EQ 1) THEN BEGIN                  ; if the user chose M modes
    IF Old_Format EQ 1 THEN BEGIN
      modetimes = data.times.M                        ; important time string array
      datatimes = convert_times(data.times.M, 1)      ; time strings converted to decimal day-of-year (to compare with user's chosen time)
      energies = Menergies                            ; important energy array
      Acurrents = data.cupA.M.toarray()               ; cup A currents are M modes
      Bcurrents = data.cupB.M.toarray()               ; cup B currents are M modes
      Ccurrents = data.cupC.M.toarray()               ; cup C currents are M modes
      Dcurrents = data.cupD.M.toarray()               ; cup D currents are M modes
    ENDIF ELSE BEGIN
      modetimes = Data_00.Field001                    ; important time string array
      datatimes = convert_times(modetimes, 0)      ; time strings converted to decimal day-of-year (to compare with user's chosen time)
      energies = Menergies                            ; important energy array

      index_val = closest(datatimes, time)
      ; I need to generalize this format to pull the correct timestamp first and match the correct index

      ; Define and fill arrays
      Acurrents = FINDGEN(128)
      Bcurrents = Acurrents
      Ccurrents = Acurrents
      Dcurrents = Acurrents
      FOR iiiii = 0  , 127 DO BEGIN
        Acurrents(iiiii) = Data_00.(iiiii + 5 + 0  )[index_val]
        Bcurrents(iiiii) = Data_00.(iiiii + 5 + 128)[index_val] ; May need to check these values and see if it should be +127, 255, and 311
        Ccurrents(iiiii) = Data_00.(iiiii + 5 + 256)[index_val]
        Dcurrents(iiiii) = Data_00.(iiiii + 5 + 312)[index_val]
      ENDFOR
    ENDELSE


    current__holder =make_array(3612)
    current_total= [Acurrents,Bcurrents,Ccurrents,Dcurrents]
    current_total = round(current_total)
    current__holder = current_total[uniq(current_total, sort(current_total))]
    current__array_size = (size(current__holder))[1]
    xvals=findgen(current__array_size)
    ;pp = plot(xvals, current__holder, /ylog, ytitle = 'unique values', xtitle = 'elements')
    ;ppl = plot(xvals, current__holder, ytitle = 'unique values', xtitle = 'elements')
    ;gradient =alog(current__holder[*])/xvals[*]
    ;print, gradient
    ;gradplot = plot(xvals, gradient, ytitle = 'gradient', xtitle = 'elements')

  ENDIF ELSE IF (mode EQ 0) THEN BEGIN       ; if the user chose L modes
    IF Old_Format EQ 1 THEN BEGIN
      modetimes = data.times.L                     ; important time string array
      datatimes = convert_times(data.times.L, 1)   ; time strings converted to decimal day-of-year (to compare with user's chosen time)
      energies = Lenergies                         ; important energy array
      Acurrents = data.cupA.L.toarray()            ; cup A currents are L modes
      Bcurrents = data.cupB.L.toarray()            ; cup B currents are L modes
      Ccurrents = data.cupC.L.toarray()            ; cup C currents are L modes
      Dcurrents = data.cupD.L.toarray()            ; cup D currents are L modes
    ENDIF ELSE BEGIN
      modetimes = Data_01.Field001                    ; important time string array
      datatimes = convert_times(modetimes, 0)      ; time strings converted to decimal day-of-year (to compare with user's chosen time)
      energies = Menergies                            ; important energy array

      index_val = closest(datatimes, time)
      ; I need to generalize this format to pull the correct timestamp first and match the correct index

      ; Define and fill arrays
      Acurrents = FINDGEN(16)
      Bcurrents = Acurrents
      Ccurrents = Acurrents
      Dcurrents = Acurrents
      FOR iiiii = 0  , 15 DO BEGIN
        Acurrents(iiiii) = Data_00.(iiiii + 5 + 0 )[index_val]
        Bcurrents(iiiii) = Data_00.(iiiii + 5 + 16)[index_val] ; May need to check these values and see if it should be +15, 31, and 47
        Ccurrents(iiiii) = Data_00.(iiiii + 5 + 32)[index_val]
        Dcurrents(iiiii) = Data_00.(iiiii + 5 + 48)[index_val]
      ENDFOR
    ENDELSE


    ; For some reason, L Mode checks for unique values? This is an older version and not sure of the reason
    Around = round(Acurrents)
    Bround = round(Bcurrents)
    Cround = round(Ccurrents)
    Dround = round(Dcurrents)

    Acurrentsunique = [uniq(Around, sort(Around))] ; sorting arrays and removing duplicate values
    Bcurrentsunique = [uniq(Bround, sort(Bround))] ; sorting arrays and removing duplicate values
    Ccurrentsunique = [uniq(Cround, sort(Cround))] ; sorting arrays and removing duplicate values
    Dcurrentsunique = [uniq(Dround, sort(Dround))] ; sorting arrays and removing duplicate values
    current__holder =make_array(3612)
    current_total= [Acurrentsunique,Bcurrentsunique,Ccurrentsunique, Dcurrentsunique]


    current_total = round(current_total)
    current__holder = current_total

    current__array_size = (size(current__holder))[1]
    xvals=findgen(current__array_size)
    

    ;  ; Dougherty 2020 - CAREFUL - I HAVE NOT PROGRAMMED THIS TO WORK FOR S,U,N YET
    ; THIS PART IS ONLY FOR JUPITER. For the rest of the planets the electrons data are prepared in advance in MATLAB
;    IF Structure.Planet_Number EQ 5 THEN BEGIN
;      IF mode EQ -1 THEN BEGIN       ; if the user chose E modes
;        modetimes = data.times.E1                    ; time string array (going with E1 modes for this since the 2 E-modes will be combined)
;        datatimes = convert_times(data.times.E1)     ; time strings converted to decimal day-of-year (to compare with user's chosen time)
;        energies = Eenergies                         ; important energy array
;        ; combine the currents
;        E1currents = data.Emodes.E1.toarray()        ; array for E1 currents
;        E2currents = data.Emodes.E2.toarray()        ; array for E2 currents
;        ; data may be missing for one mode and not for another, so ensure that the E2 mode to be plotted with each E1 mode
;        ;   is from a time close enough to the E1 mode that it must be from the same "scan"
;        t = n_elements(data.times.E1)                ; define total number of times to be number of times of E1 mode data
;        Ecurrents = dblarr(t, 28)                    ; make a 2D double array big enough to hold all E data (16 E1 channels + 12 E2 channels)
;        Ecurrents[*,0:15] = E1currents               ; first 16 elements of Ecurrents array correspond to the E1 data
;        E1times = convert_times(data.times.E1)      ; convert E1 times to decimal DOY to compare how close E2 times are
;        E2times = convert_times(data.times.E2)      ; convert E2 times to decimal DOY to compare how close they are to E1
;        ; There are 48 seconds between the time that an E1 mode is taken and the time an E2 mode is taken.
;        timespan = 4.8d1 / 8.64d4         ; this is the fraction of a day taken up by 48 seconds (the time between when E1 and E2 were taken)
;    
;        ; match up E2 currents to where their corresponding E1 currents are in the Ecurrents array
;        FOR i = 0, t-1 DO BEGIN           ; for each E1 mode time
;          FOR j = 0, n_elements(E2times)-1 DO BEGIN     ; for each E2 mode time
;            ; if an E2 mode time falls within 48 seconds of an E1 mode time, the E2 currents correspond to those E1 currents
;            IF ( (E2times[j] LE (E1times[i] + timespan)) AND (E2times[j] GT (E1times[i] - timespan)) ) THEN BEGIN  ; if E2 data is near E1 data      ; removed 1/7 no need for e mode info
;              Ecurrents[i,16:27] = E2currents[j,*]       ; add these E2 currents to the Ecurrents array in the same spot as the E1 currents
;              ; if there's an E1 mode with no matching E2 mode, the E2 part of the spectrum will be considered 0 and not plotted
;            ENDIF       ; an E2 mode is within 48 seconds of an E1 mode
;          ENDFOR      ; each element of E2times
;        ENDFOR      ; each element of E1times
;        ; now all E mode quantities are defined
;      ENDIF      ; important quantities have been defined for every mode the user could have selected
;    ENDIF
    ;TO HERE
 
  ; Logan 2025_03_03 - Edited from ".E" to ".E1" in places as code was crashing on E-Modes
  ENDIF ELSE IF (mode EQ -1) THEN BEGIN       ; if the user chose E modes
    ; Logan 2025_03_04 - Jupiter needs .E1, Uranus/Neptune need just .E (REDUNDANCY ALERT AS I COPIED A CODE BLOCK FROM ABOVE)
    IF Structure.Planet_Number EQ 5 THEN BEGIN
      IF mode EQ -1 THEN BEGIN       ; if the user chose E modes
        modetimes = data.times.E1                    ; time string array (going with E1 modes for this since the 2 E-modes will be combined)
        datatimes = convert_times(data.times.E1)     ; time strings converted to decimal day-of-year (to compare with user's chosen time)
        energies = Eenergies                         ; important energy array
        ; combine the currents
        E1currents = data.Emodes.E1.toarray()        ; array for E1 currents
        E2currents = data.Emodes.E2.toarray()        ; array for E2 currents
        ; data may be missing for one mode and not for another, so ensure that the E2 mode to be plotted with each E1 mode
        ;   is from a time close enough to the E1 mode that it must be from the same "scan"
        t = n_elements(data.times.E1)                ; define total number of times to be number of times of E1 mode data
        Ecurrents = dblarr(t, 28)                    ; make a 2D double array big enough to hold all E data (16 E1 channels + 12 E2 channels)
        Ecurrents[*,0:15] = E1currents               ; first 16 elements of Ecurrents array correspond to the E1 data
        E1times = convert_times(data.times.E1)      ; convert E1 times to decimal DOY to compare how close E2 times are
        E2times = convert_times(data.times.E2)      ; convert E2 times to decimal DOY to compare how close they are to E1
        ; There are 48 seconds between the time that an E1 mode is taken and the time an E2 mode is taken.
        timespan = 4.8d1 / 8.64d4         ; this is the fraction of a day taken up by 48 seconds (the time between when E1 and E2 were taken)
    
        ; match up E2 currents to where their corresponding E1 currents are in the Ecurrents array
        FOR i = 0, t-1 DO BEGIN           ; for each E1 mode time
          FOR j = 0, n_elements(E2times)-1 DO BEGIN     ; for each E2 mode time
            ; if an E2 mode time falls within 48 seconds of an E1 mode time, the E2 currents correspond to those E1 currents
            IF ( (E2times[j] LE (E1times[i] + timespan)) AND (E2times[j] GT (E1times[i] - timespan)) ) THEN BEGIN  ; if E2 data is near E1 data      ; removed 1/7 no need for e mode info
              Ecurrents[i,16:27] = E2currents[j,*]       ; add these E2 currents to the Ecurrents array in the same spot as the E1 currents
              ; if there's an E1 mode with no matching E2 mode, the E2 part of the spectrum will be considered 0 and not plotted
            ENDIF       ; an E2 mode is within 48 seconds of an E1 mode
          ENDFOR      ; each element of E2times
        ENDFOR      ; each element of E1times
        ; now all E mode quantities are defined
      ENDIF      ; important quantities have been defined for every mode the user could have selected
    ENDIF ELSE BEGIN
      modetimes = data.times.E                     ; time string array (going with E1 modes for this since the 2 E-modes will be combined)
      datatimes = convert_times(data.times.E)     ; time strings converted to decimal day-of-year (to compare with user's chosen time)
      energies = Eenergies                         ; important energy array
      ; combine the currents
      Ecurrents = data.Emodes.E.toarray()        ; array for E1 currents
    ENDELSE
  ENDIF      ; important quantities have been defined for every mode the user could have selected


  IF ((mode EQ 0) OR (mode EQ 1)) THEN BEGIN
    ; Get closest data time stamp
    IF Old_Format EQ 1 THEN BEGIN
      times= closest(datatimes, time)
      new_stuff={current_array_A:Acurrents[times,*], current_array_B:Bcurrents[times,*], current_array_C:Ccurrents[times,*], current_array_D:Dcurrents[times,*], times_array:datatimes[times,*]}
    ENDIF ELSE BEGIN
      ; Because of file format, the correct time has already been determined for the S,U,N files to pull the Currents
      new_stuff={current_array_A:Acurrents[*], current_array_B:Bcurrents[*], current_array_C:Ccurrents[*], current_array_D:Dcurrents[*], times_array:datatimes[index_value]}
    ENDELSE
  ENDIF ELSE IF mode EQ -1 THEN BEGIN

    times= closest(datatimes, time)
    ; For electrons there is only a current in the D cup, thus A, B, and C are multiplied by 0
    ; Dougherty 2020 - CAREFUL - I HAVE NOT PROGRAMMED THIS TO WORK FOR S,U,N YET
    IF E_Flux EQ 1 THEN BEGIN
      new_stuff = {timearray:datatimes, currentarray:Ecurrents}
    ENDIF ELSE IF E_Flux EQ 0 THEN BEGIN
      new_stuff={current_array_A:Ecurrents[times,*]*0d, current_array_B:Ecurrents[times,*]*0d, current_array_C:Ecurrents[times,*]*0d, current_array_D:Ecurrents[times,*], times_array:datatimes[times,*]}
    ENDIF
  ENDIF
  ;create a structure of current values and the corresponding time value for this measurement to be returned from the function
  return, new_stuff
END




;------------------------------------------------------------------------
; ENERGIES AND TIME CONVERSIONS 
;------------------------------------------------------------------------
;
;------------------------------------------------------------------------
; conversion functions: conversionM, conversionLE2, conversionE1
;------------------------------------------------------------------------
;
;------------------------------------------------------------------------
;
; conversionM
;   function that converts an array of channels into an array of energies
;   (to be used with M-Mode Voyager data)
;
;  INPUT: channelArray - an array of channels representing energy levels
;          (should be 128 channels ranging from 1 to 128)
;
;  OUTPUT: converts channelArray into an array of energy values, energyArray
;            and returns energyArray
;            (steps in energyArray are logarithmic between 10eV and 5950eV)
;
;------------------------------------------------------------------------
FUNCTION conversionM, channelArray
  energyArray = 60. * 10.^(channelArray / 64.) - 50.    ; conversion formula for M modes
  return, energyArray           ; return an array of the energies corresponding to the channels
END     ; conversionM

;
;------------------------------------------------------------------------
;
; conversionLE2
;   function that converts an array of channels into an array of energies
;   (to be used with L-mode or E2-mode Voyager data)
;
;  INPUT: channelArray - an array of channels representing energy levels
;          (should be 16 channels ranging from 1 to 16)
;
;  OUTPUT: converts channelArray into an array of energy values, energyArray
;            and returns energyArray
;            (steps in energyArray are logarithmic between 10eV and 5950eV)
;
;------------------------------------------------------------------------
FUNCTION conversionLE2, channelArray
  energyArray = 60. * 10.^(channelArray / 8.) - 50.    ; conversion formula for L and E2 modes
  return, energyArray           ; return an array of the energies corresponding to the channels
END     ; conversionLE2

;
;------------------------------------------------------------------------
;
; conversionE1
;   function that converts an array of channels into an array of energies
;   (to be used with E1-mode Voyager data)
;
;  INPUT: channelArray - an array of channels representing energy levels
;          (should be 16 channels ranging from 1 to 16)
;
;  OUTPUT: converts channelArray into an array of energy values, energyArray
;            and returns energyArray140eV)
;
;------------------------------------------------------------------------
FUNCTION conversionE1, channelArray
  energyArray = 60. * 10.^(channelArray / 32.) - 50.    ; conversion formula for E1 modes
  return, energyArray           ; return an array of the energies corresponding to the channels
END     ; conversionE1

;------------------------------------------------------------------------
; convert_times function
;
;  function to convert an array of strings representing times to an
;    array of doubles (decimal day of year)
;  to be used with data from Voyager
;
;  INPUT: timesStringArray - an array of strings that each give a time
;           format of string is: year, day of year, hour, minute, second
;           #### ### ## ## ##.####
;
;  OUTPUT: decimalDOYarray, the same array of times but converted
;            from strings into doubles with format of decimal DOY
;            (ie, day.fraction-of-day, ###.#######)
;
;  Note: because this is data from Voyager, the year is always 1979,
;           so the year will be ignored (not converted)
;
;------------------------------------------------------------------------
FUNCTION convert_times, timesStringArray, Old_Format

  IF ISA(Old_Format) EQ 0 THEN Old_Format = 1
  IF Old_Format EQ 1 THEN BEGIN
    decimalDOYarray = dblarr( n_elements(timesStringArray) )  ; make an array of doubles that's as big as the string array
    FOR time = 0, n_elements(timesStringArray)-1 DO BEGIN     ; for every element (every time) in the string array
      compressedTime = strcompress(timesStringArray[time])      ; get rid of extra spaces in the string
      timePieces = strsplit(compressedTime, " ", /extract)      ; split up the string into substrings wherever there are spaces
      ; Now timePieces is an array of substrings: [year, DOY, hour, minute, second]
      ; we know which element of timePieces corresponds to each kind of time unit, so convert each element into a double
      year = double(timePieces[0])               ; convert the first element of timePieces, year, into a double
      day = double(timePieces[1])                ; convert the DOY into a double
      hour = double(timePieces[2])               ; convert the hour into a double
      minute = double(timePieces[3])             ; convert the minute into a double
      second = double(timePieces[4])             ; convert the second into a double
      ; turn these into a value for decimal DOY by converting each part to fractions of a day, then adding them all together:
      ;    day / (1 day per day)
      ; +  hour / (24 hours per day)
      ; +  minute / (60 minutes per hour) / (24 hours per day)
      ; +  second / (60 seconds per minute) / (60 minutes per hour) / (24 hours per day)
      ;  = decimal DOY
      ; then fill in the corresponding element of decimalDOYarray with the result
      decimalDOYarray[time] = day + hour/24. + minute/(60.*24.) + second/(3600.*24.)
    ENDFOR         ; every time in timesStringArray, so we fill up decimalDOYarray

    ; The S,U,N files have a slightly different naming scheme for the times
  ENDIF ELSE BEGIN
    decimalDOYarray = dblarr( n_elements(timesStringArray) ) ; make an array of doubles that's as big as the string array
    FOR time = 0, n_elements(timeStringsArray)-1 DO BEGIN    ; for every element (every time) in the string array
      compressedTime = timesStringArray[time]                ; Format is YYYY-DOYTHH:MM:SS.xxx (where xxx is milliseconds)
      year   = double(compressedTime.substring(0 ,3 ))
      day    = double(compressedTime.substring(5 ,7 ))
      hour   = double(compressedTime.substring(9 ,10))
      minute = double(compressedTime.substring(12,13))
      second = double(compressedTime.substring(15,20))

      decimalDOYarray[time] = day + hour/24. + minute/(60.*24.) + second/(3600.*24.)
    ENDFOR
  ENDELSE
  return, decimalDOYarray     ; return the new array of times converted to doubles of decimal DOY
END     ; convert_times


