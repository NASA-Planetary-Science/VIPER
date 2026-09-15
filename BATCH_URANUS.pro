PRO Batch_Uranus

; This runs through the Uranus files that are in the CSV_Files directory in VIPER
; We performed final post-pruning to the data for PDS upload. 
; Efforts were made to align the input files with our final output files, but due to time constraints, a final check of these files has not been performed
;   and a few "bad" fits may be present if you recreate with this file

; PROTONS - M MODE
VIPER, 'u', Input_Filename = 'U_Full_M.csv', OUTPUT_FILENAME='U_Full_M_Fits.csv', POD='/00_Fits/Uranus/2025_07_23/Ions/M/Full', $
  /INITIAL, /MAKE_PLOT, /FULL_CSV, /PLOT_QUIETLY, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE


; PROTONS - L MODE
VIPER, 'u', Input_Filename = 'U_B1_0604_Ll.csv', OUTPUT_FILENAME='U_B1_Ions_Ll_Fit.csv', POD='/00_Fits/Uranus/2025_07_23/Ions/L/B1', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY
VIPER, 'u', Input_Filename = 'U_B2_0604_Ll.csv', OUTPUT_FILENAME='U_B2_Ions_Ll_Fit.csv', POD='/00_Fits/Uranus/2025_07_23/Ions/L/B2', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY
VIPER, 'u', Input_Filename = 'U_B3_0604_Ll.csv', OUTPUT_FILENAME='U_B3_Ions_Ll_Fit.csv', POD='/00_Fits/Uranus/2025_07_23/Ions/L/B3', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY
VIPER, 'u', Input_Filename = 'U_B6_0604_Ll.csv', OUTPUT_FILENAME='U_B6_Ions_Ll_Fit.csv', POD='/00_Fits/Uranus/2025_07_23/Ions/L/B6', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY
; This is a combination of all 4 above files
VIPER, 'u', Input_Filename = 'U_Full_Ll.csv', OUTPUT_FILENAME='U_Full_Ions_Ll_Fit.csv', POD='/00_Fits/Uranus/2025_07_23/Ions/L/Full', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY


; ELECTRONS
VIPER, 'u', Input_Filename = 'U_R2_El_Final.csv', OUTPUT_FILENAME='U_R2_El_FinalUnc_Fits.csv', POD='/00_Fits/Uranus/2025_07_23/Electrons/R2', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY
VIPER, 'u', Input_Filename = 'U_R3_El_Final.csv', OUTPUT_FILENAME='U_R3_El_FinalUnc_Fits.csv', POD='/00_Fits/Uranus/2025_07_23/Electrons/R3', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY
VIPER, 'u', Input_Filename = 'U_R4_El_Final.csv', OUTPUT_FILENAME='U_R4_El_FinalUnc_Fits.csv', POD='/00_Fits/Uranus/2025_07_23/Electrons/R4', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY
; This is a combination of all 3 above files
VIPER, 'u', Input_Filename = 'U_Full_El_Final.csv', OUTPUT_FILENAME='U_Full_El_FinalUnc_Fits.csv', POD='/00_Fits/Uranus/2025_07_23/Electrons/Full', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV, /PLOT_QUIETLY





STOP
; Example Fit to fit a specific row of data
; CSV_Indices is what is displayed in Excel, with the header row = 1, first row of data = 2, etc.

VIPER, 'u', Input_Filename = 'U_B6_0604_Ll.csv', OUTPUT_FILENAME='Delete.csv', POD='/00_Fits/Uranus/Delete', $
  /INITIAL, /MAKE_PLOT, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, CSV_INDICES = [50]
 
; Solar Wind Fits - Not included in PDS
VIPER, 'u', Input_Filename = 'Uranus_R1_M.csv', OUTPUT_FILENAME='Uranus_R1_M_Test.csv', POD='/00_Fits/Uranus/SW_Inbound/R1/M', $
  /INITIAL, /MAKE_PLOT, /PLOT_QUIETLY, /UNCERTAINTIES, /WRITE_UNCERTAINTIES_IN_FILE, /FULL_CSV

END