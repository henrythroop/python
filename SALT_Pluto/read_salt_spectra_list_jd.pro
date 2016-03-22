pro read_salt_spectra_list_jd

; Read in JD of an observation. For each JD, compute the Pluto sub-solar longitude.
;
; NB: This is an IDL program. But it is in the ~/python directory because both the input and output are 
; python programs.
;
; HBT 14-Jun-2015

common units

lines = read_txt_file('salt_spectra_list_jd.txt')				; Created by salt_read.py

naifinit, kernelfile = expand_path('~/gv/dev/gv_kernels_new_horizons.txt')  	; PCK0008 (ie, LHR, old-style)

name_observer = 'Earth'
name_target   = 'Pluto'

method = 'Near Point: Ellipsoid'
fixref = 'IAU_PLUTO'
abcorr = 'LT'

print, 'JD, UTC, Subsollon_Pluto [deg], Subsollat_Pluto [deg]'

for i = 0, sizex(lines)-1 do begin
  line = lines[i]
  if grep(line, 'meanjd') then begin
    s = strsplit(line, ' ', /EXT)
    jd = s[5]
    CSPICE_UTC2ET, 'JD' + string(jd), et
    CSPICE_ET2UTC, et, 'C', 3, utc


    CSPICE_SUBSLR, method, name_target, et, fixref, abcorr, name_observer, $
      spoint, trgepc, srfvec

    cspice_reclat, spoint, radius, lon, lat

    print, jd + ', ' + utc + ', ' + st(lon*r2d) + ', ' + st(lat*r2d)

  end

end


stop

end
