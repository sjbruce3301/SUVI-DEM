function suvi_get_tresp, spacecraft = spacecraft

;; +
;; SUVI_GET_TRESP:
;;
;; Returns the temperature response functions for the GOES-16 Solar Ultraviolet Imager.
;;
;; INPUTS:
;;    None
;;
;; OUTPUT:
;;    response: A six element data structure with the following for each passband:
;;        Name: String Identifying what the structure is for
;;        Wavelength: Applicable passband (angstroms)
;;        Units: Units of the temperature response
;;        LogTe: Array with Log Te for which the response is computed
;;        Response: Instrumental Temperature Response
;;
;; REQUIREMENTS:
;;    Must have GOES-R/SUVI and SDO/AIA packages installed.
;;
;; HISTORY:
;;    Created by Dan Seaton (dseaton@boulder.swri.edu) March 2022
;; -

	if ~keyword_set(spacecraft) then spacecraft = 16


	;; get coronal emissivities
	emiss = aia_get_response(/emiss)

	wavelength_list = [94, 131, 171, 195, 284, 304]

	for w = 0, 5 do begin 

		wav = wavelength_list[w]

		;; get the spectral response
		suvi_resp = suvi_get_response(wav, spacecraft)

		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		;; interpolate the response to match the wavelenghts in the emissivities ;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

		;; wavelength to be interpolated
		wav_int = suvi_resp.wavelength

		;; response function to be interpolated
		cal_int_new = suvi_resp.thin_open_resp * suvi_resp.platescale

		;; region of overlap in spectral response and emissivities
		known_spec = where(emiss.wave ge min(wav_int) AND emiss.wave le max(wav_int))

		;; we only interpolate where we have overlap
		new_wav = emiss.wave[known_spec]

		;; do the interpolation
		new_resp_init = spl_init(wav_int, cal_int_new)
		new_resp = spl_interp(wav_int, cal_int_new, new_resp_init, new_wav)

		;; get emissivities where we know the spectrum
		useful_emiss = emiss.emissivity[known_spec, *]
		wavestep = new_wav[1] - new_wav[0]

		; Temperature response is matrix multiply of interpolated response and useful part of emissivty matrix
		; Have to add conversion factors for wavelength bin size
		tresp_out = REFORM( (new_resp>0) # useful_emiss)
		tresp_out = tresp_out * wavestep

		if w eq 0 then $
			response = create_struct('NAME', 'GOES-R SUVI Temperature Response', $
								'WAVELENGTH', wav, $
								'UNITS', 'DN cm^5 s^-1 pix^-1', $
								'LOGTE', emiss.logte, $
								'RESPONSE', tresp_out) $
		else $
			response = [ response, create_struct('NAME', 'GOES-R SUVI Temperature Response', $
								'WAVELENGTH', wav, $
								'UNITS', 'DN cm^5 s^-1 pix^-1', $
								'LOGTE', emiss.logte, $
								'RESPONSE', tresp_out) ]
	endfor

	return, response

end
