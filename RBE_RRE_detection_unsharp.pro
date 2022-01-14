dpath_sst = '/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/'
nx_H =1796
ny_H = 1138
nw_H = 27
nt_H = 96

openr, lun, dpath_sst+'nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube', /get_lun

dat_nb_H = assoc(lun, intarr(nx_H,ny_H,nw_H,/nozer),512)
closing_RBE = intarr(nx_H,ny_H,nt_H)
closing_RRE = closing_RBE*0
dopplergram_enhanced = closing_RBE*0
RBE_RRE_combined = closing_RBE*0

for t=0, nt_H-1 do begin
	image_at_time= dat_nb_H[t]
	blue_wing_avg = mean(image_at_time[*,*,8:10],dimension=3)
	red_wing_avg = mean(image_at_time[*,*,16:18],dimension=3)
	mod_dopp_image = blue_wing_avg - red_wing_avg
	;Dopp_image = image_at_time[*,*,9]-image_at_time[*,*,17]
	result_unsharped = UNSHARP_MASK(mod_dopp_image, RADIUS=20,amount=2)
	mask_RBE = result_unsharped*0
	mask_RRE = result_unsharped*0
	w_RBE = where(result_unsharped lt -1.4*stddev(result_unsharped)) ;lt -250
	w_RRE = where(result_unsharped gt 2.5*stddev(result_unsharped)) ; gt 450
	mask_RBE[w_RBE] = 1
	mask_RRE[w_RRE] = 1
	opening_RBE = morph_open(mask_RBE, REPLICATE(1,3,3))
	closing_RBE[*,*,t] = morph_close(opening_RBE, REPLICATE(1,3,3))
	opening_RRE = morph_open(mask_RRE, REPLICATE(1,3,3))
	closing_RRE[*,*,t] = morph_close(opening_RRE, REPLICATE(1,3,3))
	dopplergram_enhanced[*,*,t] = result_unsharped
	RBE_RRE_combined[*,*,t] = closing_RBE[*,*,t]+closing_RRE[*,*,t]
	print,string(13b)+'Detecting spicules. % finished: ',float(t)*100./(nt_H-1),format='(a,f4.0,$)'
endfor

print, 'Creating CRISPEX compatible cubes'
write_data_to ='/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/souvikb_CBP/'
lp_write, closing_RBE, write_data_to+'Unshrped_RBEs_nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube'
lp_write, closing_RRE, write_data_to+'Unsrhped_RREs_nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube'
lp_write, dopplergram_enhanced, write_data_to+'Dopplergram_nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube'
lp_write, RBE_RRE_combined, write_data_to+'Unshrped_RBEs_RREs_comb_nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube'
end

