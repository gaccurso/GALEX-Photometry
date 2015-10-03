FUNCTION galex_imaging_new_routine, obsid, fits_image_to_read_fuv, fits_image_to_read_nuv 

READCOL, 'LOWMASS_catalog.dat', id, rac, decl, redshift, logM, SFR, r50, r90, D25, colour, format='L,F,F,F,F,F,F,F,F,F', comment='#'
AD = WHERE(id EQ obsid)
ra = rac[AD]
dec = decl[AD]
redshift = redshift[AD]
optical_uv_sfr = SFR[AD]
beamsize = 1.5*D25[AD]

img = mrdfits(fits_image_to_read_fuv, 0, hdr)
d=size(img,/dimension)
nx=d[0]
ny=d[1]

ra_sexagesimal = sixty(ra/15)   ;converting degree decimals into degree longtitude and latitude coordinates
dec_sexagesimal = sixty(dec) ;converting degree decimals into degree longtitude and latitude coordinates
degree_radius = ((beamsize/2.0)/3600.0) 
new_ra = ra + degree_radius 
new_dec = dec + degree_radius
ra_radius_sexagesimal = sixty(new_ra/15) 
dec_radius_sexagesimal = sixty(new_dec) 
coordinates = string([ra_radius_sexagesimal, dec_radius_sexagesimal], /PRINT)   ;co-ordinates of radius point
stringad, coordinates, x, y           ; cordinates of radius point in decimals
ADXY, hdr, ra, dec, x_center, y_center  ; coordinates of center in pixel coordinates
ADXY, hdr, x, y, x_radius, y_radius     ; coordinates of radius point in pixel coordinates
radius_in_pixels = abs(long(x_radius)-long(x_center))
print, 'radius in pixels' , radius_in_pixels
;radius_for_photometry = 3.0*radius_in_pixels		;watch with this, could introduce errors
print, x_center, y_center, 'centers of the image'

radin=1.0*radius_in_pixels
radout=5.0*radius_in_pixels
;areacircle = 3.14159*radin*radin
;areashell = (3.14159*radout*radout) - areacircle

totalflux=0.0
background=0.0
error = 0.0 

;counter_shell = 0.0

;for i=0, nx-1 do begin
;	for j=0, ny-1 do begin
;		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin
;			background = background + (img[i, j])
;			counter_shell = counter_shell + 1.0
;		endif else begin 
;			totalflux = totalflux
;			background = background
;		endelse
;	endfor
;endfor

;skylevel = ((background)/(counter_shell))
;print, totalflux, background, skylevel, 'answers'

skylevel = 0.0
counter_shell = 0.0
for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin			
			error = error + (img[i, j] - skylevel)^2.0  ;^2.0
			counter_shell = counter_shell + 1.0
		endif else begin 
			error = error 
		endelse
	endfor
endfor

skylevel = 0.0
counter_circle = 0.0
for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radin) then begin
			radius = sqrt(((abs(i-x_center))*(abs(i-x_center)))+((abs(j-y_center))*(abs(j-y_center))))
			weighting = 1.0 ;exp(-(radius/radius_in_pixels)*(radius/radius_in_pixels)*alog(2))
			totalflux = totalflux + ((img[i, j]-skylevel)*weighting)
			counter_circle = counter_circle + 1.0
		endif else begin 
			totalflux = totalflux
			background = background	
		endelse
	endfor
endfor

answer_fuv_error_flux =  (sqrt(error)/sqrt(counter_shell))*sqrt(counter_circle)    ;;(error*counter_circle)/(counter_shell)     ; ;sqrt((error)/(counter_shell))*sqrt(counter_circle)

print, totalflux, background

x =  float(strsplit(hdr[52], 'EXPTIME=', /EXTRACT))
print, 'EXPTIME=' , x[1]

real_flux = totalflux			;this has units of CPS
print, real_flux
;the above flux is in poisson counts per second

;To convert from GALEX counts per second (CPS) to magnitudes in the AB system (Oke 1990):
;http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html

FUV_mag = (-2.5*alog10(real_flux)) + 18.82 ;- FUV_attenuation     ;*exposure_time_fuv
FUV_lambda = 1516.0				;taken from http://galexgi.gsfc.nasa.gov/docs/galex/Documents/ERO_data_description_2.htm
FUV_flux = ((1.40E-15)*real_flux*FUV_lambda)		;this is a correction. It now has units of erg sec-1 cm-2 ;(10^((-48.60 - FUV_mag)/2.5))*1.978891820580474E+15     ; *5.07E+14  or 6.429E+14          ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2/hz range taken from the GALEX website;  



distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758e+18)

FUV_luminosity = (FUV_flux*4.0*3.1415*distance_in_cm*distance_in_cm)/(3.846E+33)		;this now has units of solar luminosity 

; from the jarret paper http://www.ast.uct.ac.za/~jarrett/papers/jarrett.II.preprint.pdf
 
answer_fuv = (10^(-9.69))*FUV_luminosity



answer_fuv_error_flux_mag =  (0.434*answer_fuv_error_flux)/(real_flux)            ; at first done wrongly (-2.5*alog10(answer_fuv_error_flux)) + 18.82 ;- FUV_attenuation     ;*exposure_time_fuv
answer_fuv_error_flux_flux = FUV_flux*2.303*answer_fuv_error_flux_mag             ; at first done wrongly (10^((-48.60 - answer_fuv_error_flux_mag)/2.5))*1.978891820580474E+15     ; *5.07E+14  or 6.429E+14          ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2/hz range taken from the GALEX website
distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758e+18)

answer_fuv_error_flux_luminosity = (FUV_luminosity*answer_fuv_error_flux_flux)/(FUV_flux)        ; this was originally done wrongly (answer_fuv_error_flux_flux*4*3.1415*distance_in_cm*distance_in_cm)/(3.846E+33)		;this now has units of solar luminosity 
answer_fuv_error = (answer_fuv*answer_fuv_error_flux_luminosity)/(FUV_luminosity)      ; this was originally done wrongly (10^(-9.69))*answer_fuv_error_flux_luminosity

signal_to_noise_fuv = (totalflux/(answer_fuv_error_flux))

;###########################################################################################################################
;##################################The NUV image photometry#################################################################
;###########################################################################################################################

READCOL, 'LOWMASS_catalog.dat', id, rac, decl, redshift, logM, SFR, r50, r90, D25, colour, format='L,F,F,F,F,F,F,F,F,F', comment='#'
AD = WHERE(id EQ obsid)
ra = rac[AD]
dec = decl[AD]
redshift = redshift[AD]
optical_uv_sfr = SFR[AD]
beamsize = 1.5*D25[AD]

img = mrdfits(fits_image_to_read_nuv, 0, hdr)
d=size(img,/dimension)
nx=d[0]
ny=d[1]

ra_sexagesimal = sixty(ra/15)   ;converting degree decimals into degree longtitude and latitude coordinates
dec_sexagesimal = sixty(dec) ;converting degree decimals into degree longtitude and latitude coordinates
degree_radius = ((beamsize/2.0)/3600.0) 
new_ra = ra + degree_radius 
new_dec = dec + degree_radius
ra_radius_sexagesimal = sixty(new_ra/15) 
dec_radius_sexagesimal = sixty(new_dec) 
coordinates = string([ra_radius_sexagesimal, dec_radius_sexagesimal], /PRINT)   ;co-ordinates of radius point
stringad, coordinates, x, y           ; cordinates of radius point in decimals
ADXY, hdr, ra, dec, x_center, y_center  ; coordinates of center in pixel coordinates
ADXY, hdr, x, y, x_radius, y_radius     ; coordinates of radius point in pixel coordinates
radius_in_pixels = abs(long(x_radius)-long(x_center))
print, 'radius in pixels' , radius_in_pixels
;radius_for_photometry = 3.0*radius_in_pixels		;watch with this, could introduce errors

radin=1.0*radius_in_pixels
radout=5.0*radius_in_pixels
;areacircle = 3.14159*radin*radin
;areashell = (3.14159*radout*radout) - areacircle

totalflux=0.0
background=0.0
error = 0.0  

;for i=0, nx-1 do begin
;	for j=0, ny-1 do begin
;		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin
;			background = background + (img[i, j])
;		endif else begin 
;			totalflux = totalflux
;			background = background
;		endelse
;	endfor
;endfor

skylevel = 0.0
;print, totalflux, background, skylevel, 'answers'

counter_shell = 0.0
for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radout) and (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) GT radin) then begin
			error = error + (img[i, j] - skylevel)^2.0  ;^2.0
			counter_shell = counter_shell + 1.0
		endif else begin 
			error = error 
		endelse
	endfor
endfor



skylevel = 0.0
counter_circle = 0.0

for i=0, nx-1 do begin
	for j=0, ny-1 do begin
		if (sqrt((abs(i-x_center))^2.0 + (abs(j-y_center))^2.0) LT radin) then begin
			radius = sqrt(((abs(i-x_center))*(abs(i-x_center)))+((abs(j-y_center))*(abs(j-y_center))))
			weighting = 1.0 ;exp(-(radius/radius_in_pixels)*(radius/radius_in_pixels)*alog(2))
			totalflux = totalflux + ((img[i, j]-skylevel)*weighting)
			counter_circle = counter_circle + 1.0
		endif else begin 
			totalflux = totalflux
			background = background	
		endelse
	endfor
endfor

answer_nuv_error_flux = (sqrt(error)/sqrt(counter_shell))*sqrt(counter_circle)  ;(error*counter_circle)/(counter_shell)   ; sqrt((error)/(counter_shell))*sqrt(counter_circle)

print, totalflux, background

x =  float(strsplit(hdr[52], 'EXPTIME=', /EXTRACT))
print, 'EXPTIME=' , x[1]

real_flux = totalflux

;the above flux is in poisson counts per second

;To convert from GALEX counts per second (CPS) to magnitudes in the AB system (Oke 1990):
;http://galexgi.gsfc.nasa.gov/docs/galex/FAQ/counts_background.html

NUV_mag = (-2.5*alog10(real_flux)) + 20.08 ;- FUV_attenuation     ;*exposure_time_fuv
NUV_lambda = 2267.0		

NUV_flux = ((2.06E-16)*real_flux*NUV_lambda*NUV_lambda*1E-10)/(3.00E+08)		;this is a correction. It now has units of erg sec-1 cm-2 Hz-1	;(10^((-48.60 - NUV_mag)/2.5))     ;;;;;*1.978891820580474E+15     ; *5.55E+14            ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2/hz range taken from the GALEX website
print, NUV_flux					             ; this is the flux seen at our observational point.

distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758e+18)

NUV_luminosity = (NUV_flux*4.0*3.1415*distance_in_cm*distance_in_cm)    ;/((3.00E+08)/(NUV_lambda*1E-10))		;this now has units of ergs/s/hz
 
answer_nuv = (10^(-28.165))*NUV_luminosity



answer_nuv_error_flux_mag = (0.434*answer_nuv_error_flux)/(real_flux)            ; at first done wrongly (-2.5*alog10(answer_nuv_error_flux)) + 20.08 ;- FUV_attenuation     ;*exposure_time_fuv

answer_nuv_error_flux_flux = NUV_flux*2.303*answer_nuv_error_flux_mag             ; at first done wrongly (10^((-48.60 - answer_nuv_error_flux_mag)/2.5))   ; *5.07E+14  or 6.429E+14          ; conversion see http://en.wikipedia.org/wiki/AB_magnitude units of F_UV ergs/s/cm2/hz range taken from the GALEX website
distance_in_cm = lumdist(redshift)*(1E+06)*(3.08567758e+18)
answer_nuv_error_flux_luminosity = (NUV_luminosity*answer_nuv_error_flux_flux)/(NUV_flux)        ; this was originally done wrongly (answer_nuv_error_flux_flux*4*3.1415*distance_in_cm*distance_in_cm)		;this now has units of solar luminosity 
answer_nuv_error = (answer_nuv*answer_nuv_error_flux_luminosity)/(NUV_luminosity)      ; this was originally done wrongly (10^(-28.165))*answer_nuv_error_flux_luminosity

signal_to_noise_nuv = (totalflux/(answer_nuv_error_flux))

phot = mrdfits('SED_SFR_170113.fits', 1)
proper_nuv = phot[WHERE(phot.id EQ obsid)].mag[5]
proper_fuv = phot[WHERE(phot.id EQ obsid)].mag[6]

A = fltarr(14)
A = [obsid, answer_fuv, answer_nuv, optical_uv_sfr, answer_fuv_error, answer_nuv_error, NUV_mag, FUV_mag, answer_nuv_error_flux_mag, answer_fuv_error_flux_mag, signal_to_noise_fuv, signal_to_noise_nuv, proper_nuv, proper_fuv]   
return, A


end




