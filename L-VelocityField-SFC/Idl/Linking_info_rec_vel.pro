
path = "/u/mingli/data/work_ksz_all_data/"
base = "Halo"

Group_info = {$
              M_Crit200    : 0.0 , $
              R_Crit200    : 0.0 , $
              CM           : fltarr(3) , $
              Pos          : fltarr(3) , $
              Vel          : fltarr(3) , $
              VelDisp      : 0.0   $
             }

snum = 54                       ; Z=0.25
outbase = base + "_Z=0.25"

; -----
info      = read_ascii("./zlist.txt")
info      = info.(0)
ExpFactor = info(1,*)
afact  = ExpFactor[snum]
scale  = [afact, afact, afact]
center = [1500., 1500., 1500.]

deg2rad     = !pi/180.
;----------------------------------------------------------
sp =1500
th =sp
nns=2048L
nr=2

fix1="ns"+strcompress(string(fix(nns)),/remove_all)
fix2="r"+strcompress(string(fix(sp)),/remove_all)
fix3="d"+strcompress(string(fix(th)),/remove_all)

ICASE = 3
case ICASE of
0:fix4="nr"+strcompress(string(fix(nr)),/remove_all)
1:fix4="nr_l2"
2:fix4="maxbcg"
;3:fix4 = 'sdss_m1.5e3_m2e3'
;3:fix4 = 'sdss_m2e3_m3e3'
;3:fix4 = 'sdss_m3e3_m5e3'
3:fix4 = 'sdss_m5e3_mmax'
endcase

velbase = "_"

;goto,analysis
goto,projecting

;----------------------------------------------------------
analysis:

openr,lun,path+'dat.project.3/Project.3.Halo_Pos_Vel_mass_bin_'+fix1+'_'+fix2+'_'+fix4+'.dat',/get_lun
Ntot = 0.
readu,lun,Ntot
print, 'Total Halo Num ', Ntot

RMhalo=fltarr(Ntot)
RPhalo=fltarr(3,Ntot)
RVhalo=fltarr(3,Ntot)

readu,lun,RMhalo
readu,lun,RPhalo
readu,lun,RVhalo
free_lun,lun


;----------------------------------------------------------
BaseDir =  "/u/mingli/data/" 

nfile = 54L
fexts='000'
fexts=fexts+strcompress(string(nfile),/remove_all)
fexts=strmid(fexts,strlen(fexts)-3,3)

Fstruct = {$
			g0    : fltarr(3), $
			g1_25 : fltarr(3), $
			g2_5  : fltarr(3), $
			g5    : fltarr(3), $
			g10   : fltarr(3), $
			g20   : fltarr(3)  $
		  }

FieldData = replicate(Fstruct, Ntot)

ntask = 0l
nfiles = 0l
boxsize = 0.D
SI = 0l
nslab = 0l
fslab = 0l

rep = 0
dim = 0
fv  = BaseDir+'/halo_fields_'+fexts+'/'+str(rep)+'/'+$
	  'halo_rec_vel'+velbase+str(dim)+'_2048_'+fexts+'_'+str(rep)+'.0'
print,fv
openr,lun, fv
readu,lun, ntask, nfiles, boxsize, SI, nslab, fslab
close,lun

pindx = (RPhalo[0,*] + center[0]) / BoxSize * SI
pindy = (RPhalo[1,*] + center[1]) / BoxSize * SI
pindz = (RPhalo[2,*] + center[2]) / BoxSize * SI

pcumu = 0
for nrf=0,31 do begin
  print, '-------------------------------------------------------'
  print, 'FILE ', nrf
  print, ''
  rep = 0
  dim = 0
  fnr = strcompress(string(nrf),/remove_all)
  fv  = BaseDir+'/halo_fields_'+fexts+'/'+str(rep)+'/'+$
		'halo_rec_vel'+velbase+str(dim)+'_2048_'+fexts+'_'+str(rep)+'.'+fnr
  openr, lun, fv
  readu, lun, ntask, nfiles, boxsize, SI, nslab, fslab
  close, lun

  bdx1 = fslab
  bdx2 = fslab + nslab

  indx = where(pindx ge bdx1 and pindx lt bdx2,cntx)

  if cntx ne 0 then begin 
	print, ' - There is ', cntx, ' halos in this file'
	pcumu += cntx

	for rep = 0,5 do begin 
	  print, '--------------------------------'
	  print, " -- Smoothing scale", 1.25*2^(rep-1)
	  print, ''

	  nloop = 2

	  for dim=0,2 do begin
		print, ''
		print,' --- VEL-FIELD DIM ', dim
		print, ''
		fdim = strcompress(string(dim),/remove_all)

		velr  =fltarr(SI, SI, nloop*nslab)
		for loop = 0, nloop-1 do begin
		  print,' ---- Generating loop ', loop+1

		  nrr = (nrf + loop + ntask) mod ntask
		  fnrr = strcompress(string(nrr),/remove_all)
		  fv  = BaseDir+'/rec_vel_fields_'+fexts+'/'+str(rep)+'/'+$
				'rec_vel'+velbase+fdim+'_2048_'+fexts+'_'+str(rep)+'.'+fnrr

		  tmp  =fltarr(SI, SI, nslab)
		  openr, lun, fv
		  skip_lun, lun, 4L*2 + 8L + 4L*3
		  readu, lun,tmp
		  close, lun
		  velr[*,*,nslab*loop:nslab*(loop+1)-1] = tmp[*,*,*]
		endfor
		;field = TRANSPOSE(field, [2, 1, 0])

		for ii = 0, cntx-1 do begin
		  ipx = floor(pindx[indx[ii]]) - fslab
		  ipy = floor(pindy[indx[ii]])
		  ipz = floor(pindz[indx[ii]])

		  ipx1 = (ipx + 1)
		  ipy1 = (ipy + 1) mod SI
		  ipz1 = (ipz + 1) mod SI

		  dx  = pindx[indx[ii]] - floor(pindx[indx[ii]])
		  dy  = pindy[indx[ii]] - floor(pindy[indx[ii]])
		  dz  = pindz[indx[ii]] - floor(pindz[indx[ii]])
	  
		  if ipx lt 0 or ipx ge nslab or $
		     ipy lt 0 or ipy ge    SI or $
		     ipz lt 0 or ipz ge    SI then begin
			print, 'Something is wrong with halo position ', indx[ii]
			stop
		  endif

		  ;amp = momt[ipz,ipy,ipx] / dens[ipz,ipy,ipx]
		  ;amp = field[ipx,ipy,ipz] 
		  amp = velr[ipz ,ipy ,ipx ] * (1.0 - dx) * (1.0 - dy) * (1.0 - dz) + $
				velr[ipz ,ipy1,ipx ] * (1.0 - dx) *        dy  * (1.0 - dz) + $
				velr[ipz1,ipy ,ipx ] * (1.0 - dx) * (1.0 - dy) *        dz  + $
				velr[ipz1,ipy1,ipx ] * (1.0 - dx) *        dy  *        dz  + $
				velr[ipz ,ipy ,ipx1] *        dx  * (1.0 - dy) * (1.0 - dz) + $
				velr[ipz ,ipy1,ipx1] *        dx  *        dy  * (1.0 - dz) + $
				velr[ipz1,ipy ,ipx1] *        dx  * (1.0 - dy) *        dz  + $
				velr[ipz1,ipy1,ipx1] *        dx  *        dy  *        dz

		  case rep of
		  0:FieldData[indx[ii]].g0   [dim] = amp ;* sqrt(scale[dim])
		  1:FieldData[indx[ii]].g1_25[dim] = amp ;* sqrt(scale[dim])
		  2:FieldData[indx[ii]].g2_5 [dim] = amp ;* sqrt(scale[dim])
		  3:FieldData[indx[ii]].g5   [dim] = amp ;* sqrt(scale[dim])
		  4:FieldData[indx[ii]].g10  [dim] = amp ;* sqrt(scale[dim])
		  5:FieldData[indx[ii]].g20  [dim] = amp ;* sqrt(scale[dim])
		  endcase
		endfor
	  endfor
	endfor
	print, 'done for file', nrf
  endif else print, 'These is no need to load file ', nrf
endfor

openw,lun,'./Halo_Halo_CIC_Rec_Vel'+velbase+'mass_bin_'+fix1+'_'+fix2+'_'+fix4+'.dat',/get_lun

writeu,lun,Ntot
writeu,lun,FieldData
free_lun,lun

goto, endfun
;----------------------------------------------------------
projecting:

Fstruct = {$
			g0    : fltarr(3), $
			g1_25 : fltarr(3), $
			g2_5  : fltarr(3), $
			g5    : fltarr(3), $
			g10   : fltarr(3), $
			g20   : fltarr(3)  $
		  }

;----------------------------------------------------------
Ntot = 0.
openr,lun,path+'dat.project.3/Project.3.Halo_Pos_Vel_mass_bin_'+fix1+'_'+fix2+'_'+fix4+'.dat',/get_lun
readu,lun,Ntot

RRMhalo=fltarr(Ntot)
RRPhalo=fltarr(3,Ntot)
RRVhalo=fltarr(3,Ntot)

readu,lun,RRMhalo
readu,lun,RRPhalo
readu,lun,RRVhalo
free_lun,lun

;----------------------------------------------------------
openr,lun,'./Halo_Halo_CIC_Rec_Vel'+velbase+'mass_bin_'+fix1+'_'+fix2+'_'+fix4+'.dat',/get_lun
readu,lun,Ntot
FieldData = replicate(Fstruct, Ntot)
readu,lun,FieldData
free_lun,lun

;----------------------------------------------------------

RVField = fltarr(Ntot,6)

for ii = 0, Ntot-1 do begin
  pvec = RRPhalo[*,ii]
  for rep = 0, 5 do begin
	case rep of
	0:vvec=FieldData[ii].g0   
	1:vvec=FieldData[ii].g1_25
	2:vvec=FieldData[ii].g2_5   
	3:vvec=FieldData[ii].g5  
	4:vvec=FieldData[ii].g10  
	5:vvec=FieldData[ii].g20  
	endcase

	angle, pvec, vvec, 1, rad
	RVField[ii,rep]=sqrt(total(vvec^2))*rad
  endfor
endfor

openw,lun,fix4+'_'+'halo_rec_cic_vel'+velbase+'field.dat',/get_lun
writeu,lun,Ntot
writeu,lun,RVField
free_lun,lun

goto, endfun
;-----------
endfun:
stop

end
