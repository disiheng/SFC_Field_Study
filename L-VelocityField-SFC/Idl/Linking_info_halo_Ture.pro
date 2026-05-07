
path = "/u/mingli/data/work_ksz_all_data/"
base = "Halo"

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
NSide = nns

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

RVhalo = fltarr(Ntot)

for ii = 0, Ntot-1 do begin
  pvec = RRPhalo[*,ii]
  vvec = RRVhalo[*,ii]
  angle, pvec, vvec, 1, rad
  RVhalo[ii]=sqrt(total(vvec^2))*rad
endfor

openw,lun,fix4+'_'+'halo_true_vel.dat',/get_lun
writeu,lun,Ntot
writeu,lun,RVhalo
free_lun,lun

goto, endfun
;-----------
endfun:
stop

end
