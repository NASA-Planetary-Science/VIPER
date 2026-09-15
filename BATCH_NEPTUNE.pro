PRO BATCH_NEPTUNE

; This runs through the Neptune files that are in the CSV_Files directory in VIPER
; We performed final post-pruning to the data for PDS upload.
; Efforts were made to align the input files with our final output files, but due to time constraints, a final check of these files has not been performed
;   and a few "bad" fits may be present if you recreate with this file


;;;; ELECTRONS ;;;;
VIPER, 'n', Input_Filename = 'Neptune_R2_E_Long.csv', OUTPUT_FILENAME='Neptune_E_Final.csv', POD='/00_Fits/Neptune/EFinal/R2', $
  /INITIAL, /MAKE_PLOT, /PLOT_QUIETLY, /FIT, /FULL_CSV, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE
  
  
;;;; IONS ;;;;
  
; B0 - Hotter signal prevalent in C-Cup  
VIPER, 'n', Input_Filename = 'Neptune_R2_L_Long_B0.csv', OUTPUT_FILENAME='Neptune_B0_Final.csv', POD='/00_Fits/Neptune/LlFinal/B0', $
  /INITIAL, /MAKE_PLOT, /PLOT_QUIETLY, /FIT, /FULL_CSV, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE
  
; B1 - Colder signal D-Cup
VIPER, 'n', Input_Filename = 'Neptune_R2_L_Long_B1.csv', OUTPUT_FILENAME='Neptune_B1_Final.csv', POD='/00_Fits/Neptune/LlFinal/B1', $
  /INITIAL, /MAKE_PLOT, /PLOT_QUIETLY, /FIT, /FULL_CSV, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE

; Warm/Hot signal - Potentially H+, N+, and then a warm proxy N+ (some combo of H/N, but impossible to discern, like how we used a single warm O+ at Jupiter)   
VIPER, 'n', Input_Filename = 'Neptune_R2_L_Long_B2.csv', OUTPUT_FILENAME='Neptune_B2_Final.csv', POD='/00_Fits/Neptune/LlFinal/B2', $
  /INITIAL, /MAKE_PLOT, /PLOT_QUIETLY, /FIT, /FULL_CSV, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE

; B3 - Like B1 - Colder signal in D-Cup   
VIPER, 'n', Input_Filename = 'Neptune_R2_L_Long_B3.csv', OUTPUT_FILENAME='Neptune_B3_Final.csv', POD='/00_Fits/Neptune/LlFinal/B3', $
  /INITIAL, /MAKE_PLOT, /PLOT_QUIETLY, /FIT, /FULL_CSV, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE


STOP

;;;; IONS ;;;;
; Single M-Mode Spectra - Not included in PDS
; M-Mode did not have enough signal:noise, and was not included for Neptune fits. Included is a fit of the "best", and it still isn't great
VIPER, 'n', Input_Filename = 'Neptune_R2_M.csv', OUTPUT_FILENAME='Delete.csv', POD='/00_Fits/Delete', $
  /INITIAL, /MAKE_PLOT, /PLOT_QUIETLY, /FIT, CSV_INDICES = [73], /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE

END