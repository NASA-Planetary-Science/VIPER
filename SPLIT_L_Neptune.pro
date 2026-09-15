PRO SPLIT_L_Neptune 

;; A quick function to split the Neptune input data files by Short/Long for E and L. Uranus was already split.

Input_Filename  = 'Neptune_R4_E.csv'
OutputFilenameS = 'Neptune_R4_E_Short.csv'
OutputFilenameL = 'Neptune_R4_E_Long.csv'

CD, CURRENT = RESTORE_DIR
STR_DIR = STRCOMPRESS(STRING(RESTORE_DIR, '/CSV_Files'))
CD, STR_DIR
Input_CSV = READ_CSV(Input_Filename, HEADER = Structure_Names)

    
; If Full CSV keyword is set, use every row of CSV
CSV_Indices = LINDGEN(N_ELEMENTS(Input_CSV.Field01))


FOR i = 0, N_ELEMENTS(CSV_Indices) - 1 DO BEGIN                       ;start loop
  ;nloop = nloop + 1 ; Augment loop for next row of CSV file
  CSV_Index = CSV_Indices(i)

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
  
  ; Short (Short_or_Long = 0)
  if LONG(Structure.Short_or_Long) EQ LONG(0) THEN BEGIN     
      ; If variable exists, append, otherwise create
      IF isa(OutputStructS) THEN BEGIN
        OutputStructS = [OutputStructS, Structure]
      ENDIF ELSE BEGIN
        OutputStructS = CREATE_STRUCT(Structure)
        Header = TAG_NAMES(Structure) ; Grab tag names
      ENDELSE 
  ; Same for Long (Short_or_Long = 1) 
  ENDIF ELSE BEGIN
      ; If variable exists, append, otherwise create
      IF isa(OutputStructL) THEN BEGIN
        OutputStructL = [OutputStructL, Structure]
      ENDIF ELSE BEGIN
        OutputStructL = CREATE_STRUCT(Structure)
        Header = TAG_NAMES(Structure) ; Grab tag names
      ENDELSE      
  ENDELSE
     
ENDFOR

; Write the 2 CSV files, remember to check manually - one extra row that exsisted
WRITE_CSV, OutputFilenameS, OutputStructS, HEADER = Header
WRITE_CSV, OutputFilenameL, OutputStructL, HEADER = Header

; Return to main directory
CD, RESTORE_DIR
END