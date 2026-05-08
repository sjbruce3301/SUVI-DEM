pro DEM_work, index, data, index_array, dem_cube, logte

; scratch work procedure for making a DEM cube
;
; fetch and load data using the following:
; 
;result = vso_search('2023-04-20 02:58:00', '2023-04-20 03:54:00', inst = 'suvi', source = 'goes16')
;status = vso_get(result)

print,SYSTIME()
files = file_search('/Users/sbruce/SUVI_400', '*.fits.gz') ;/Desktop
;files = file_search('/Users/sbruce/Desktop/Work/AIA_data', '*.fits') ;/Desktop
mreadfits, files, index
;sorted_files = files[sort(index.date_d$obs)]
sorted_files = files
read_suvi, sorted_files, index, data;, /despike
;stop

; INPUTS:
;    index - FITS header data structure
;    data - data cube
;
; OUTPUTS:
;    index_array - array of indexes for long exposures used in DEM calculation
;    dem_cube - 4 D array with DEMs (Dimensions: [x, y, LogT, time])
;    logte - list of Log T's used in calculation
;
; REQUIREMENTS:
;    This requires the following SSW packages:
;      PROBA2/SWAP - used for image coalignment before computing DEMs
;      GOES-R/SUVI - used to generate response functions
;      SDO/AIA - used to generate emissivities needed for response functions
;      SIMPLE_REG_DEM - used to compute DEMS
;          see: https://ui.adsabs.harvard.edu/abs/2020ApJ...905...17P/abstract


;;;;;;;;;;;;;;;;;;;;;
;;; Initial Setup ;;;
;;;;;;;;;;;;;;;;;;;;;

;; get the SUVI temperature response
suvi_tresp = suvi_get_tresp()
; logte = suvi_tresp[0].logte[25:65]
logte = suvi_tresp[0].logte[30:52:1] ; cut off between 4 - 5, check difference, try one at 2.8/3. 5.7 LOGTE - 6.3. 34 - 46
;stop

resp = fltarr(n_elements(logte), 5)
; for n = 0, 4 do resp[*, n] = suvi_tresp[n].response[25:65]
for n = 0, 4 do resp[*, n] = suvi_tresp[n].response[30:52:1] ; change here too


;; all images are normalized so this is just a dummy value
exptimes = replicate(1., 5)

;; need this to generate responses
filter_list = [replicate('Thin_Zr', 2), replicate('Thin_Al', 4)]

;; estimate camera noise in DN
approx_camera_noise = 1.25

;; set up images and timing data
;wvln_list = [94, 131, 171, 195, 284, 304]
wvln_list = [131, 171, 195, 284, 304]
long_ims = where(index.exptime gt 0.9)
short_ims = where(index.exptime lt 0.9)

time_utc = anytim2utc(index.date_d$obs)
time_sec = (time_utc.time - time_utc[0].time)/1000.

;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Make some HDR images ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; NOTE: skip this whole step if you are not working with flare data and just use long exposure images ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; We'll try to find the best correspondence between 131 images and the others
target_times = time_sec[long_ims[where(index[long_ims].wavelnth eq 131)]]
data_size = size(data)

;; this cube receives the DEM output
dem_cube = fltarr( data_size[1], data_size[2], n_elements(logte), n_elements(target_times) )

;; loop over all the target times
for q = 0, n_elements(target_times) - 1 do begin

	;; here we find the long exposures that best match the target time image
	subdata = fltarr(1280, 1280, 6)
	subimages = lonarr(6)
		for n = 0, 4 do begin & $ ;5 instead of 4****
		candidate_ims = where( index[long_ims].wavelnth eq wvln_list[n] ) & $
		delta_t = time_sec[long_ims[candidate_ims]] - target_times[q] & $
		target_im = where(min(abs(delta_t)) eq abs(delta_t)) & $
		subdata[*, *, n] = data[*, *, long_ims[candidate_ims[target_im]]] & $
		subimages[n] = long_ims[ candidate_ims[ target_im ]] & $
	endfor
	subindex = index[subimages]
	
	hdr_data = subdata
;
;	;; here we find the corresponding short exposures, needed to deal with saturated pixels
;	shortdata = subdata * 0.
;	shortimages = subimages * 0.
;	for n = 0, 5 do begin & $
;		candidate_ims = where( index[short_ims].wavelnth eq wvln_list[n] ) & $
;		image_time_sec = time_sec[subimages[n]] & $
;		delta_t = time_sec[short_ims[candidate_ims]] - image_time_sec & $
;		target_im = where(min(abs(delta_t)) eq abs(delta_t)) & $
;		shortdata[*, *, n] = data[*, *, short_ims[candidate_ims[target_im]]] & $
;		shortimages[n] = short_ims[ candidate_ims[ target_im ]] & $
;	endfor
;	shortindex = index[shortimages]
;
;	;; this ad hoc array keeps track of the approximate theoretical maximum
;	;; value in a SUVI long exposure for each channel. This information is
;	;; supposed to be in the FITS headers, but is currently not implemented correctly
;	sat_thr = [218.486, 31.6662, 23.9244, 45.1470, 200.399, 256.768]
;
;	;; set some values to determine how we want to roll over from long to short image
;	threshold = 0.9   ; 
;	range = 0.05
;
;	;; merge the long and short images using the parameters above to roll smoothly
;	;; from one to the other
;	hdr_data = subdata
;	for n = 0, 5 do begin
;		image = subdata[* , *, n]
;		short_image = shortdata[*, *, n]
;
;		xvals = findgen(16000.)/16000. * sat_thr[n]
;
;		lower_limit = (threshold - range) * sat_thr[n] 
;		upper_limit = (threshold + range) * sat_thr[n] 
;		
;		weights = (1./(lower_limit - upper_limit) * (xvals - upper_limit) > 0) < 1
;		weights = gauss_smooth(weights, 3, /edge_truncate)
;		mask = weights[round((image > 0)/max(image) * 16000)]
;
;		hdr_data[*, *, n] = image * mask + short_image * (1. - mask)
;	endfor


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Align channels and attempt to separate 304/284 contamination ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;; array to receive aligned images
	aldata = hdr_data

	;; this uses the PROBA2/SWAP limb finder to find the pointing for each channel
	c_param_arr = fltarr(3, 6)
	for n = 0, 4 do begin & $ ;*******5 instead of 4
		c_param = [index[n].crpix1 - 1, index[n].crpix2 - 1, index[n].diam_sun/2.] & $
		if n ne 4 then $;*****5 4
			result = p2sw_fit_limb(hdr_data[*, *, n], c_param, 'max', 10000, 40, [0.99, 1.01], roll = -index[n].crota, tol = 1e-3, /holes) $
		else $
			result = p2sw_fit_limb(hdr_data[*, *, n], c_param, 'min', 10000, 40, [0.99, 1.01], roll = -index[n].crota, tol = 1e-3, /holes) & $
		c_param_arr[*, n] = c_param & $
	endfor

	;; calculate pointing offsets for each channel
	disp = fltarr(2, 6)
	disp[0, *] = c_param_arr[0, 0] - c_param_arr[0, *]
	disp[1, *] = c_param_arr[1, 0] - c_param_arr[1, *]

	;; align the images
	aldata = shift_img(aldata, disp)


	;; Initial setup to determine approximate response levels for specific lines
	r284 = suvi_get_response(284) 
	r304 = suvi_get_response(304)

	;; need to fix the hard coding below
	;; to do this use that Fe XV line is at 284.1630 Å
	;; and He II line is at 303.785 Å

	r284at284 = r284.thin_open_resp[2742]
	r284at304 = r284.thin_open_resp[2938]
	r304at284 = r304.thin_open_resp[2742]
	r304at304 = r304.thin_open_resp[2938]

	;; it's not really right to assume that the ratio of spectra is the same in every pixel
	;; but we don't have a good way to adjust as a function of pixel (yet -- maybe a DEM-based forward model can help)
	;; so we adjust to get a clean result in the bulk of the pixels - 0.8 works well for quiet regions, but is not good for flares
	empirical_factor = 0.8
	r284at304 = r284at304 * empirical_factor

	;; contaminated images
	image284_raw = aldata[*, *, 4]
	image304_raw = aldata[*, *, 5]

	;; decontamination calculation
	image284 = ((image304_raw * r284at304 - image284_raw * r304at304) / ( r284at304 * r304at284 - r284at284 * r304at304)) * r284at284
	image304 = ((image304_raw * r284at284 - image284_raw * r304at284) / (-r284at304 * r304at284 + r284at284 * r304at304)) * r304at304

	;; put the corrected images into the final repaired data array
	;; this is working poorly for flares, so we just skip this step for now
	fixed_data = aldata
	fixed_data[*, *, 4] = image284
	fixed_data[*, *, 5] = image304
	fixed_data = float(fixed_data)


	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Revert data from radiance to counts ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

	;; this routine allows us to get the data in counts and also roughly estimate the error level in every pixel
	data_dn = fixed_data
	errors = data_dn * 0.
	for n = 0, 5 do begin & $
		data_dn[*, *, n] = suvi_restore_dn(fixed_data[*, *, n], subindex[n], approx_resp = resp_conversion, error_estimate = error_image ) & $
		errors[*, *, n] = error_image & $
	endfor

	;; negative values screw up the DEM calculation
	;; we map them to small positive numbers with a huge error so they get ignored
	dem_input_data = data_dn
	mask = where(dem_input_data lt 0)
	errors[mask] = 10 * max(errors)
	dem_input_data[mask] = 1e-5

	;;;;;;;;;;;;;;;;;;;;;;;;;;
	;;; Calculate the DEMs ;;;
	;;;;;;;;;;;;;;;;;;;;;;;;;;
	print, '---- CHECKS ----'
	print, 'errors min/max: ', min(errors), max(errors)
	print, 'exptimes: ', exptimes
	print, 'resp min/max: ', min(resp), max(resp)
	print, 'finite resp bad: ', total(finite(resp) eq 0)
	print, 'finite data bad: ', total(finite(dem_input_data) eq 0)
	print, 'logte min/max: ', min(logte), max(logte)
	print, 'min/max data: ', min(dem_input_data), max(dem_input_data)
	;stop

  print,'pre-dem'
  print,SYSTIME()
	dems = simple_reg_dem(dem_input_data, errors, exptimes, logte, resp, chi2)
	print,'dem done'
	print,SYSTIME()

	;; put all the output data into a data cube
	dem_cube[*, *, *, q] = dems
	if q eq 0 then index_array = subindex else index_array = [[index_array], [subindex]]
	;stop
	print,max(dem_cube),min(dem_cube)
	;print,dem_cube
	;stop
	;print,dem_cube
	SAVE, dem_cube, FILENAME = '/Users/sbruce/Desktop/Work/DEM_no94_new1.sav'
	print,'file saved.'

  endfor

print,'for_end'

end


