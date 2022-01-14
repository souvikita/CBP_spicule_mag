dpath ='/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/souvikb_CBP/SDO_only_aia/'
files = file_search(dpath+'*.fits', count=count)
for f=0,count-1 do begin
	str_channel = strmid(files[f],92,3)
	aia_temp = readfits(files[f], h1)
	dimension = size(aia_temp)
	temp_avg = mean(aia_temp,dimension=3)
	enhanced_aia_cube = fltarr(dimension[1], dimension[2], dimension[3])
	for time =0,dimension[3]-1 do begin
		unnshp_mask = unsharp_mask(aia_temp[*,*,time], radius = 30, amount=2)
		enhanced_aia_cube[*,*,time] = unnshp_mask - temp_avg
		plot_image, reform(enhanced_aia_cube[*,*,time]), min=-3500, max=4000
		;stop
		;save, enhanced_aia_cube, filename =dpath +'aia'+str_channel+'sav'
		
	endfor
	print,string(13b)+'Enhancing the images:  % finished: ',$
           float(f)*100./(count),format='(a,f4.0,$)'
;	save, enhanced_aia_cube, filename =dpath +'aia'+str_channel+'_enhanced_small_radius'+'.sav'
	;stop
endfor
end


