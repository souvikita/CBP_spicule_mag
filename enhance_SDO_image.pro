dpath_SDO = '/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/edvarda/sdo2sst/'
sdo_cube = lp_read(dpath_SDO+'sdo_crispex_for_px_nw11_nt96.bcube')
;IDL> lp_header, sdo_cube
;stokes=[I], ns=1 :  datatype=1 (byte), dims=3, nx=1796, ny=1138, nt=1056, endian=l
sdo_cube_rearr = reform(sdo_cube,1796,1138,11,96)

enhanced_cube = sdo_cube_rearr*0
for wav=0,10 do begin
	temp_avg = mean(sdo_cube_rearr[*,*,wav,*],dimension=4) ; This gives temporal average for a given wavelength channel (2D image)
	for time =0,95 do begin
		unshp_msk = unsharp_mask(sdo_cube_rearr[*,*,wav,time],radius=50,amount=2.) ;amount is the strength of the high pass
		enhanced_cube[*,*,wav,time] = unshp_msk - temp_avg ; subtracting the temporal average to enhance the short-lived features
	endfor
	print,string(13b)+'Enhancing the images:  % finished: ',$
           float(wav)*100./(11),format='(a,f4.0,$)'
endfor
enhanced_cube = reform(enhanced_cube,1796,1138,1056)

lp_write, enhanced_cube, '/mn/stornext/d9/souvikb/CBP_spicules_SDO/enhanced_sdo_crispex_for_px_nw11_nt96.bcube'
end



