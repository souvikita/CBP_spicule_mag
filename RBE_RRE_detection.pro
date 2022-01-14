dpath_sst = '/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/'
nx_H =1796
ny_H = 1138
nw_H = 27
nt_H = 96

openr, lun, dpath_sst+'nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube', /get_lun

dat_nb_H = assoc(lun, intarr(nx_H,ny_H,nw_H,/nozer),512)
closing_RBE = intarr(nx_H,ny_H,nt_H)
closing_RRE = closing_RBE*0

for t=0, nt_H-1 do begin
	image_at_time= dat_nb_H[t]
	Dopp_image = image_at_time[*,*,9]-image_at_time[*,*,17]
	mask_RBE = Dopp_image*0
	mask_RRE = Dopp_image*0
	w_RBE = where(Dopp_image lt -130)
	w_RRE = where(Dopp_image gt 330)
	mask_RBE[w_RBE] = 1
	mask_RRE[w_RRE] = 1
	opening_RBE = morph_open(mask_RBE, REPLICATE(1,3,3))
	closing_RBE[*,*,t] = morph_close(opening_RBE, REPLICATE(1,3,3))
	opening_RRE = morph_open(mask_RRE, REPLICATE(1,3,3))
	closing_RRE[*,*,t] = morph_close(opening_RRE, REPLICATE(1,3,3))
	print,string(13b)+'Detecting spicules. % finished: ',float(t)*100./(nt_H-1),format='(a,f4.0,$)'
endfor

print, 'Creating CRISPEX compatible cubes'
;write_data_to ='/mn/stornext/d18/lapalma/reduc/2021/2021-08-04/CHROMIS/cubes_nb/souvikb_CBP/'
;lp_write, closing_RBE, write_data_to+'RBEs_nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube'
;lp_write, closing_RRE, write_data_to+'RREs_nb_4846_2021-08-04T09:56:50_scans=0-95_corrected_cmapcorr_im.icube'

end

