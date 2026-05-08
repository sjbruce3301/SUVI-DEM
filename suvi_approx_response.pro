function suvi_approx_response, wvln, fw1, fw2, ccd_temp

;; +
;; SUVI_APPROX_RESPONSE:
;;
;; Returns coefficients required to convert SUVI L1b radiance to various other systems
;;    SUVI L1b has units W/m^2 sr-1.
;;
;; INPUTS:
;;    wvln: Wavelength (either as a string ('171 A') or integer (171))
;;    fwl: String with name of the focal plane filter 1:
;;         Valid options are: 'Open', 'Thin_Al', 'Thin_Zr', 'Thick_Zr', 'Thick_Al'
;;    fwl: String with name of the focal plane filter 2:
;;         Valid options are: 'Glass', 'Open', 'Thin_Al', 'Thin_Zr', Thick_Al'
;;    ccd_temp: Detector temperature (needed to compute Gain)
;;         A typical value is -50. Use FITS Keyword CCD_TMP1 if 
;;         working with a specific observation.
;;
;; OUTPUT:
;;    response: A five-element array containing the following elements:
;;        0: Instrument response in J/DN m^-2 sr^-1. Divide L1b by this number to get original DNs.
;;        1: Instrumental effective area in m^2 (Includes optics and QE)
;;        2: CCD Quantum Yield (electrons per photon for a given channel)
;;        3: Energy per photon used in radiance calculation
;;        4: CCD Gain (electrons per DN)
;;
;; HISTORY:
;;    Created by Dan Seaton (dseaton@boulder.swri.edu) March 2022
;;
;; CAVEATS:
;;    This code works for original L1b GOES-16 SUVI observations obtained up to at least March 2022.
;;    Future compatibility may depend on updates to calibration codes, changes made during
;;    data reprocessing, or changes made to account for instrument degradation. An update
;;    will probably be required at some point to support both reprocessed and future observations.
;; -



  if size(wvln, /type) eq 7 then $
     wvlns = ['94 A', '131 A', '171 A', '195 A', '284 A', '304 A'] $
   else $
     wvlns = [94, 131, 171, 195, 284, 304] 
  fw1_names = ['Open',  'Thin_Al', 'Thin_Zr', 'Thick_Zr', 'Thick_Al'] 
  fw2_names = ['Glass', 'Open',    'Thin_Al', 'Thin_Zr',  'Thick_Al']
  flat_name_fw1 = ['OPEN',  'THINAL', 'THINZR', 'THICKZR', 'THICKAL']
  flat_name_fw2 = ['GLASS',  'OPEN',  'THINAL', 'THINZR',  'THICKAL']
  flat_name_wvln = ['94', '131', '171', '195', '284', '304']
 
  ;; Define some constants.
  geom_area = 0.001936D ; m^2 or 19.36 cm^2
  solid_angle = 1.46903000D-10 ; steradian

  ;; Mirror reflectivity arrays for primary and secondary, by channel (using indexing above)
  m1_reflectivity = [0.3249D,  0.6615D,  0.4833D,  0.3615D,  0.2415D,  0.2264D]
  m2_reflectivity = [0.2452D,  0.6171D,  0.4459D,  0.3302D,  0.2193D,  0.2103D]
 
  ;; Entrance filter transmissivities by channel
  ent_trans = [0.3655D,  0.3289D,  0.5392D,  0.5323D,  0.3907D,  0.3645D] 

  ;; Transmission for each filter in FW1 as a function of wavelength. Wavelength is indexed
  ;; vertically, different filters horizontally.
  fw1_trans = $
     [[  1.00000000D00, 5.95142100D-03, 3.67894000D-01, 2.64501197D-01, 4.31837000D-04], $ ;  94
      [  1.00000000D00, 8.20062500D-03, 3.30119966D-01, 2.35383536D-01, 7.26064000D-04], $ ; 131
      [  1.00000000D00, 5.34282080D-01, 1.66932636D-01, 9.02533750D-02, 4.67286276D-01], $ ; 171
      [  1.00000000D00, 5.29553176D-01, 4.24509360D-02, 1.23050400D-02, 4.71508674D-01], $ ; 195
      [  1.00000000D00, 3.85627244D-01, 1.58321000D-10, 4.64234000D-15, 3.18311852D-01], $ ; 284
      [  1.00000000D00, 3.63653942D-01, 1.96980000D-12, 7.21305000D-18, 2.97853032D-01]]   ; 304
  ;;     Open           Thin_Al         Thin_Zr         Thick_Zr        Thick_Al

  ;; Same as above but for FW2
  fw2_trans = $
     [[  0.00000000D00, 1.00000000D00,  7.03444200D-03, 3.66690547D-01, 4.03902000D-04], $ ;  94
      [  0.00000000D00, 1.00000000D00,  9.57078800D-03, 3.29012990D-01, 6.82550000D-04], $ ; 131
      [  0.00000000D00, 1.00000000D00,  5.38863752D-01, 1.65916294D-01, 4.65692969D-01], $ ; 171
      [  0.00000000D00, 1.00000000D00,  5.33485838D-01, 4.19321000D-02, 4.70115272D-01], $ ; 195
      [  0.00000000D00, 1.00000000D00,  3.90370926D-01, 1.42733000D-10, 3.16758959D-01], $ ; 284
      [  0.00000000D00, 1.00000000D00,  3.68309497D-01, 1.73955000D-12, 2.96341297D-01]]   ; 304
  ;;     Glass          Open            Thin_Al         Thin_Zr         Thick_Al

  ;; CCD Quantum Efficiency as a function of wavelength
  ccd_qe_list = [0.5255D,  0.872D,   0.8256D,  0.8093D,  0.7676D,  0.7575D] 

  ;; This array contains the temperatures against which CCD gain is indexed
  ccd_gain_temp_array = $ ; Deg C
     [-88.93333333D, -87.91666667D, -86.9D       , -86.02142857D, -85.14285714D, $
      -84.26428571D, -83.38571429D, -82.50714286D, -81.62857143D, -80.75D      , $
      -79.87142857D, -78.99285714D, -78.11428571D, -77.23571429D, -76.35714286D, $
      -75.05555556D, -74.6D,        -73.65384615D, -72.70769231D, -71.76153846D, $
      -70.81538462D, -69.86923077D, -68.92307692D, -67.97692308D, -67.03076923D, $
      -66.08461538D, -65.13846154D, -64.19230769D, -63.2D       , -62.3D       , $
      -61.34615385D, -60.39230769D, -59.43846154D, -58.48461538D, -57.53076923D, $
      -56.57692308D, -55.62307692D, -54.66923077D, -53.71538462D, -52.76153846D, $
      -51.388D     , -50.396D     , -49.9D,        -48.86666667D, -47.83333333D, $
      -46.8D       , -45.76666667D, -44.73333333D, -43.7D,        -42.66666667D, $
      -41.63333333D, -40.6D       , -39.56666667D, -38.944D     , -37.5D,        $
      -36.53846154D, -35.57692308D, -34.61538462D, -33.65384615D, -32.69230769D, $
      -31.73076923D] 

  ;; This array contains the gain as a function of temperature, as indexed in the
  ;; array above. This one is for the left amplifier.
  ccd_gain_left = $ ; e- / DN
     [ 35.42917853D,  35.47891651D, 35.52865449D, 35.57163649D, 35.61461849D, $
       35.6576005D ,  35.7005825D , 35.7435645D , 35.7865465D , 35.8295285D, $
       35.8725105D ,  35.9154925D , 35.95847451D, 36.00145651D, 36.04443851D, $
       36.10811555D,  36.13040251D, 36.17669082D, 36.22297913D, 36.26926744D, $
       36.31555575D,  36.36184406D, 36.40813237D, 36.45442068D, 36.50070899D, $
       36.5469973D ,  36.5932856D , 36.63957391D, 36.68812019D, 36.73215053D, $
       36.77881517D,  36.82547981D, 36.87214444D, 36.91880908D, 36.96547372D, $
       37.01213835D,  37.05880299D, 37.10546763D, 37.15213226D, 37.1987969D, $
       37.26599398D,  37.3145252D , 37.33879081D, 37.38934417D, 37.43989753D, $
       37.49045088D,  37.54100424D, 37.5915576D , 37.64211095D, 37.69266431D, $
       37.74321767D,  37.79377102D, 37.84432438D, 37.87478685D, 37.94543109D, $
       37.99247206D,  38.03951302D, 38.08655399D, 38.13359495D, 38.18063592D, $
       38.22767688D]
  ;; Same as above but for the right amplifier. Presently these are the same
  ;; but this is because only one amplifier is used. This needs an update 
  ;; if there is a switch.
  ccd_gain_right = $ ; e- / DN
     [ 35.42917853D, 35.47891651D, 35.52865449D, 35.57163649D, 35.61461849D, $
       35.6576005D , 35.7005825D , 35.7435645D , 35.7865465D , 35.8295285D, $
       35.8725105D , 35.9154925D , 35.95847451D, 36.00145651D, 36.04443851D, $
       36.10811555D, 36.13040251D, 36.17669082D, 36.22297913D, 36.26926744D, $
       36.31555575D, 36.36184406D, 36.40813237D, 36.45442068D, 36.50070899D, $
       36.5469973D , 36.5932856D , 36.63957391D, 36.68812019D, 36.73215053D, $
       36.77881517D, 36.82547981D, 36.87214444D, 36.91880908D, 36.96547372D, $
       37.01213835D, 37.05880299D, 37.10546763D, 37.15213226D, 37.1987969D, $
       37.26599398D, 37.3145252D , 37.33879081D, 37.38934417D, 37.43989753D, $
       37.49045088D, 37.54100424D, 37.5915576D , 37.64211095D, 37.69266431D, $
       37.74321767D, 37.79377102D, 37.84432438D, 37.87478685D, 37.94543109D, $
       37.99247206D, 38.03951302D, 38.08655399D, 38.13359495D, 38.18063592D, $
       38.22767688D] 
  
     wvln_loc = where(wvln eq wvlns)
     fw1_loc = where(fw1 eq fw1_names)
     fw2_loc = where(fw2 eq fw2_names)


  ;; These can be computed from physical constants, but we use the LUT values
  ;; to ensure we handle this the same as the LUT.
  phot_elec_conversion_list = Double([ 36.17,  25.89,  19.85,  17.41,  11.95,  11.18]) ; e- / photon
  phot_energy_list =   [ 2.12000000D-17, 1.51000000D-17, 1.16000000D-17, 1.02000000D-17, $
                    6.99000000D-18, 6.54000000e-18] ; Joules / photon


  eff_area = (reform(geom_area * ent_trans[wvln_loc] * fw1_trans[fw1_loc, wvln_loc] * fw2_trans[fw2_loc, wvln_loc] * $
                        m1_reflectivity[wvln_loc] * m2_reflectivity[wvln_loc] * ccd_qe_list[wvln_loc]))[0]
  ccd_gain = interpol(ccd_gain_right, ccd_gain_temp_array, ccd_temp)
        

response = ccd_gain ; e-/DN
response = response / (phot_elec_conversion_list[wvln_loc])[0] ; phot/DN
response = response / eff_area ; phot/DN 1//m^2
response = response * (phot_energy_list[wvln_loc])[0] ; J/DN m^-2
response = response / solid_angle ; J/DN m^-2 sr^-1 

return, [response, eff_area, (phot_elec_conversion_list[wvln_loc])[0], (phot_energy_list[wvln_loc])[0], ccd_gain]
end