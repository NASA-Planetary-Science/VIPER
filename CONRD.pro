PRO CONRD
  ;; This function defines commom block variables that were taken from Alan Barnett's thesis:
  ;; Common blocks act as global variables for data that will be used throughout the VIPER analysis
  ;;   The Response Function of the Voyager Plasma Science Experiment March 1984:
  ;; File: Simconst Data N4 (As stated near the beginning of Appendix B)
  ;; A large portion of the common blocks have a read in file to get the values that they
  ;; they need

  CD, CURRENT = RESTORE_DIR
  STR_DIR = RESTORE_DIR + PATH_SEP() + 'Data' + PATH_SEP() + 'Instrument_Constants' ; works for both LINUS/Mac and Windows [GX]
  CD, STR_DIR
  COMMON TRNPAR, AMA1, AMA2, AMC1, AMSFT
  ; Read in data and get the proper arrays
  data = Make_array(121991,/DOUBLE)
  openr, 1, 'mjs.out.txt'
  readf, 1, data
  close, 1
  AMSFT_initial    = data[0:12899]
  AMA1_initial     = data[12900:25799]
  AMA2_initial     = data[25800:38699]
  AMC1_initial     = data[38700:51599]
  AMSFT      = Make_array(129,100,/DOUBLE)
  AMA1       = AMSFT
  AMA2       = AMSFT
  AMC1       = AMSFT

  For i = 0,128 Do begin
    For j = 0,99 Do begin
      AMSFT(i,j)  = DOUBLE(AMSFT_initial[j+i*100])
      AMA1(i,j)   = DOUBLE(AMA1_initial[j+i*100])
      AMA2(i,j)   = DOUBLE(AMA2_initial[j+i*100])
      AMC1(i,j)   = DOUBLE(AMC1_initial[j+i*100])
    Endfor
  Endfor



  COMMON TRNPRD, AA, CC, S
  ; Read in data and get the proper arrays
  data = Make_array(121991,/DOUBLE)
  openr, 1, 'mjs.out.txt'
  readf, 1, data
  close, 1
  S_initial        = data[51600:64499]
  aa1_initial      = data[64500:77399]
  aa2_initial      = data[77400:90299]
  cc1_initial      = data[90300:103199]
  cc2_initial      = data[103200:116099]
  S          = Make_array(129,100,/DOUBLE)
  aa1        = S
  aa2        = S
  cc1        = S
  cc2        = S
  For i = 0,128 Do begin
    For j = 0,99 Do begin
      S(i,j)      = DOUBLE(S_initial[j+i*100])
      aa1(i,j)    = DOUBLE(aa1_initial[j+i*100])
      aa2(i,j)    = DOUBLE(aa2_initial[j+i*100])
      cc1(i,j)    = DOUBLE(cc1_initial[j+i*100])
      cc2(i,j)    = DOUBLE(cc2_initial[j+i*100])
    Endfor
  Endfor
  AA = Make_array(129,100,2,/DOUBLE) ; Condense the aa1 and aa2 arrays into the proper format
  CC = AA
  AA(*,*,0) = DOUBLE(aa1(*,*))
  AA(*,*,1) = DOUBLE(aa2(*,*))
  CC(*,*,0) = DOUBLE(cc1(*,*))
  CC(*,*,1) = DOUBLE(cc2(*,*))



  COMMON ABCSET, ABCS
  ; Read in data and get the proper arrays
  data = Make_array(121991,/DOUBLE)
  openr, 1, 'mjs.out.txt'
  readf, 1, data
  close, 1
  ABCS_initial     = data[118006:119933]
  ABCS       = Make_array(241,4,2,/DOUBLE)

  For k = 0,1 Do begin
    For i = 0,3 Do begin
      For j = 0,240 Do begin
        ABCS(j,i,k) = DOUBLE(ABCS_initial[j+i*241+k*4*241])
      Endfor
    Endfor
  Endfor



  COMMON DSET, DS ; Missing Values
  ; Read in data and get the proper arrays
  data = Make_array(121991,/DOUBLE)
  openr, 1, 'mjs.out.txt'
  readf, 1, data
  close, 1
  DS_initial       = data[119934:121861]
  DS         = Make_array(241,4,2,/DOUBLE)
  For k = 0,1 Do begin
    For i = 0,3 Do begin
      For j = 0,240 Do begin
        DS(j,i,k)   = DOUBLE(DS_initial[j+i*241+k*4*241])
      Endfor
    Endfor
  Endfor


  COMMON VZ, AAAVZ, AAVZ, AVZ, AAADVZ, AADVZ, ADVZ, DVOLT, NFIRST, MFIRST
  ; Read in data and get the proper arrays
  data = Make_array(121991,/DOUBLE)
  openr, 1, 'mjs.out.txt'
  readf, 1, data
  close, 1
  AAAVZ_initial    = data[116100:116209]
  AAVZ_initial     = data[116210:116369]
  AVZ_initial      = data[116370:116537]
  AAADVZ_initial   = data[116538:116647]
  AADVZ_initial    = data[116648:116807]
  ADVZ_initial     = data[116808:116975]
  DVOLT_initial    = data[116976:117103]
  NFIRST_initial   = data[117602:117730]
  MFIRST_initial   = data[117731:117859]
  AAAVZ      = Make_array(5,22,/DOUBLE)
  AAVZ       = Make_array(40,4,/DOUBLE)
  AVZ        = Make_array(84,2,/DOUBLE)
  AAADVZ     = AAAVZ
  AADVZ      = AAVZ
  ADVZ       = AVZ
  DVOLT      = Make_array(128,/DOUBLE)
  NFIRST     = Make_array(129,/DOUBLE)
  MFIRST     = NFIRST

  For i = 0,4 Do begin
    For j = 0,21 Do begin
      AAAVZ(i,j) = DOUBLE(AAAVZ_initial[j+i*22])
      AAADVZ(i,j) = DOUBLE(AAADVZ_initial[j+i*22])
    Endfor
  Endfor

  For i = 0,39 Do begin
    For j = 0,3 Do begin
      AAVZ(i,j) = DOUBLE(AAVZ_initial[j+i*4])
      AADVZ(i,j) = DOUBLE(AADVZ_initial[j+i*4])
    Endfor
  Endfor

  For i = 0,83 Do begin
    For j = 0,1 Do begin
      AVZ(i,j) = DOUBLE(AVZ_initial[j+i*2])
      ADVZ(i,j) = DOUBLE(ADVZ_initial[j+i*2])
    Endfor
  Endfor

  DVOLT(*) = DOUBLE(DVOLT_initial(*))
  NFIRST(*) = DOUBLE(NFIRST_initial(*))
  MFIRST(*) = DOUBLE(MFIRST_initial(*))




  COMMON TAIL, TVZ, TDVZ
  ; Read in data and get the proper arrays
  data = Make_array(121991,/DOUBLE)
  openr, 1, 'mjs.out.txt'
  readf, 1, data
  close, 1
  TVZ_initial      = data[117120:117360]
  TDVZ_initial     = data[117361:117601]
  TVZ        = Make_array(241,/DOUBLE)
  TDVZ       = TVZ

  TVZ(*) = DOUBLE(TVZ_initial(*))
  TDVZ(*) = DOUBLE(TDVZ_initial(*))


  COMMON LCOMP, DVOLTL, NM, NL, AVZT, AMDVZ
  ; Read in data and get the proper arrays
  data = Make_array(121991,/DOUBLE)
  openr, 1, 'mjs.out.txt'
  readf, 1, data
  close, 1

  DVOLTL_initial   = data[117104:117119]
  NM_initial       = data[117860:117988]
  NL_initial       = data[117989:118005]
  AVZT_initial     = data[121862:121990]
  DVOLTL     = Make_Array(16,/DOUBLE)
  NM         = Make_array(129,/DOUBLE)
  NL         = Make_array(17,/DOUBLE)
  AVZT       = NM

  DVOLTL(*) = DOUBLE(DVOLTL_initial(*))
  NM(*) = DOUBLE(NM_initial(*))
  NL(*) = DOUBLE(NL_initial(*))
  AVZT(*) = DOUBLE(AVZT_initial(*))
  ; This table is taken from Barnett's 1983 thesis
  AMDVZ = DOUBLE([4.57,4.32,4.13,3.99,3.88,3.79,3.72,3.66,$
    3.62,3.58,3.56,3.54,3.53,3.52,3.52,3.53,$
    3.53,3.54,3.56,3.58,3.59,3.62,3.65,3.67,$
    3.70,3.74,3.77,3.81,3.84,3.89,3.93,3.97,$
    4.02,4.07,4.11,4.17,4.22,4.28,4.33,4.40,$
    4.45,4.51,4.58,4.64,4.70,4.78,4.85,4.92,$
    4.99,5.06,5.15,5.22,5.31,5.38,5.47,5.56,$
    5.65,5.73,5.83,5.92,6.02,6.11,6.22,6.32,$
    6.42,6.53,6.63,6.75,6.86,6.97,7.10,7.21,$
    7.33,7.46,7.59,7.71,7.85,7.98,8.12,8.25,$
    8.41,8.54,8.70,8.84,9.00,9.15,9.32,9.47,$
    9.65,9.81,9.98,10.16,10.34,10.52,10.70,10.89,$
    11.09,11.28,11.48,11.69,11.89,12.10,12.32,12.53,$
    12.76,12.99,13.22,13.45,13.69,13.94,14.18,14.44,$
    14.70,14.96,15.23,15.50,15.78,16.16,16.35,16.64,$
    16.94,17.24,17.56,17.87,18.20,18.50,18.90,19.10,19.60])
  
  ; Switch out of Instrument_Constants directory
  CD, RESTORE_DIR
  



  ; ------- JUPITER STUFF ----------
  ; 
  ; Switch to trajectories directory where Voyager 1 trajectory around Jupiter is held
  CD, CURRENT = RESTORE_DIR
  STR_DIR = RESTORE_DIR + PATH_SEP() + 'Data' + PATH_SEP() + 'Jupiter' + PATH_SEP() + 'Trajectories'
  CD, STR_DIR
  
  COMMON VGR1DATA, VGR1SSEDR, VGR1SPICE
  ; Read in the Data file that contains spacecraft locations/velocites in ECL50 and S3
  Data = Make_array(48,2925,/DOUBLE)
  openr, 1, 'v1jup.ssedr.62-65.79_rob.txt'
  readf, 1, Data
  close, 1
  VGR1SSEDR = Data
  openr, 1, 'SpiceRotations-vgr1.txt'
  rrr = Make_array(36,2925,/DOUBLE) ; Size of the rotation matrices [6x6=36]
  readf, 1, rrr
  close, 1
  VGR1SPICE = rrr

  COMMON VGR2DATA, VGR2SSEDR, VGR2SPICE
  ; Read in the Data file that contains spacecraft locations/velocites in ECL50 and S3
  Data = Make_array(48,15653,/DOUBLE)
  openr, 1, 'v2jup.ssedr.rob.txt'
  readf, 1, Data
  close, 1
  VGR2SSEDR = Data
  ; Read in rotation matrix from System III to ECL50 coordinates from the SPICE kernels
  openr, 1, 'SpiceRotations-vgr2.txt'
  rrr = Make_array(36,15653,/DOUBLE) ; Size of the rotation matrices [6x6=36]
  readf, 1, rrr
  close, 1
  VGR2SPICE = rrr
  ; Switch out of trajectories directory and return to main level
  CD, RESTORE_DIR
  
  
  ; ------- URANUS STUFF ----------
  ; Switch to trajectories directory of the planet
  CD, CURRENT = RESTORE_DIR
  STR_DIR = RESTORE_DIR + PATH_SEP() + 'Data' + PATH_SEP() + 'Uranus' + PATH_SEP() + 'Trajectories'
  CD, STR_DIR

  COMMON VGR2DATA_URANUS, VGR2SSEDR_URANUS, VGR2SPICE_URANUS
  ; Read in the Data file that contains spacecraft locations/velocites in URANUS_FIT and S3
  Data = Make_array(48,10274,/DOUBLE)
  openr, 1, 'ssedr.v2.ura.txt'
  readf, 1, Data
  close, 1
  VGR2SSEDR_URANUS = Data
  
  ;stop ; Logan 2025-04-10
  ; Read in rotation matrix from URANUS_FIT to ECL50 coordinates from the SPICE kernels
  openr, 1, 'SPICE_rotation_matrices_URANUS_FIT_to_ECLIPB1950.txt'
  rrr = Make_array(36,10274,/DOUBLE) ; Size of the rotation matrices [6x6=36]
  readf, 1, rrr
  close, 1
  VGR2SPICE_URANUS = rrr
  ; Switch out of trajectories directory and return to main level
  CD, RESTORE_DIR
  
  
  ; ------- NEPTUNE STUFF ----------
  ; Switch to trajectories directory of the planet
  CD, CURRENT = RESTORE_DIR
  STR_DIR = RESTORE_DIR + PATH_SEP() + 'Data' + PATH_SEP() + 'Neptune' + PATH_SEP() + 'Trajectories'
  CD, STR_DIR

  COMMON VGR2DATA_NEPTUNE, VGR2SSEDR_NEPTUNE, VGR2SPICE_NEPTUNE
  ; Read in the Data file that contains spacecraft locations/velocites in ECL50 and V2_NEPTUNE_PLS
  Data = Make_array(39,28325,/DOUBLE)  ; CHANGE ARRAY SIZE
  openr, 1, 'ssedr.v2.nep.modified.txt'
  readf, 1, Data
  close, 1
  VGR2SSEDR_NEPTUNE = Data
  ; Read in rotation matrix from V2_NEPTUNE_PLS to ECL50 coordinates from the SPICE kernels
  openr, 1, 'SPICE_rotation_matrices_V2_NEPTUNE_PLS_to_ECLIPB1950.txt'
  rrr = Make_array(36,28325,/DOUBLE) ; Size of the rotation matrices [6x6=36]    ; CHANGE ARRAY SIZE
  readf, 1, rrr
  close, 1
  VGR2SPICE_NEPTUNE = rrr
  ; Switch out of trajectories directory and return to main level
  CD, RESTORE_DIR
  
  
  
  
  
  COMMON GENERATE, SS, LDIST
  SS = 0.018d     ; Numerical Integration Step Size
  LDIST = 'FALSE' ; Logical flag. FALSE for current. TRUE for Reduced Distribution Function
  
  COMMON DOY_COMMON_BLOCK, DOY_COMMON_BLOCK_Var ; Define a DOY that will be less than any DOY for initial
  DOY_COMMON_BLOCK_Var = -1
  
;  COMMON DATA_FILE_COMMON_BLOCK, DATA
;  DATA = 0
  
END

; A quick common block that is used to save structures for the program so that they can be accessed
; globally. This is required because of how MPFIT handles parameters and thus cannot pass external
; arguments, like a structure, through the MPFIT code. If it were able to, or if it were recoded, then
; these common blocks could be removed
PRO Save_Structures, Struct_In_A=Struct_In_A, Struct_In_B=Struct_In_B
  COMMON Structure_Block, varA, varB
  varA = Struct_In_A
  varB = Struct_In_B
END