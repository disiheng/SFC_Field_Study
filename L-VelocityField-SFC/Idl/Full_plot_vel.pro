
path = "/u/mingli/data/work_ksz_all_data/"
base = "Halo"

snum = 54                       ; Z=0.25
;----------------------------------------------------------

; -----
info      = read_ascii("./zlist.txt")
info      = info.(0)
ExpFactor = info(1,*)

afact  = ExpFactor[snum]

sp =1500

;----------------------------------------------------------
; -----
;nns=1024L
nns=2048L
nnpix=nside2npix(nns, ERROR=error)
pixels=4.*!pi/nnpix

th=sp
nr=2
nt=2.

fix1="ns"+strcompress(string(fix(nns)),/remove_all)
fix2="r"+strcompress(string(fix(sp)),/remove_all)
fix3="d"+strcompress(string(fix(th)),/remove_all)


hnum = [7071, 6853, 5342, 5263]
prefix = ["sdss_m1.5e3_m2e3", "sdss_m2e3_m3e3", "sdss_m3e3_m5e3", "sdss_m5e3_mmax"]


ki=3 & kf=3
fac_tot = fltarr(total(hnum[ki:kf]))
RVhalo_tot = fltarr(total(hnum[ki:kf]))

;----------------------------------------------------------
rv_00_tot = fltarr(total(hnum[ki:kf]))
rv_0125_tot = fltarr(total(hnum[ki:kf]))
rv_025_tot = fltarr(total(hnum[ki:kf]))
rv_05_tot = fltarr(total(hnum[ki:kf]))
rv_10_tot = fltarr(total(hnum[ki:kf]))

sv_00_tot = fltarr(total(hnum[ki:kf]))
sv_0125_tot = fltarr(total(hnum[ki:kf]))
sv_025_tot = fltarr(total(hnum[ki:kf]))
sv_05_tot = fltarr(total(hnum[ki:kf]))
sv_10_tot = fltarr(total(hnum[ki:kf]))

;----------------------------------------------------------
skip = 0l
for kk = ki, kf do begin
fix4=prefix[kk]

Nbins=0
;----------------------------------------------------------
openr,lun,path+'dat.project.3/Project.3.Halo_R_vel_mass_bin_'+fix1+'_'+fix2+'_'+fix4+'.dat',/get_lun
readu,lun,Nbins

bincnt=lonarr(Nbins)
readu,lun,bincnt

Ntot=total(bincnt)
mbin=intarr(Ntot)
RVhalo=fltarr(Ntot)
readu,lun,mbin
readu,lun,RVhalo
free_lun,lun

print,Ntot
print,"read done ..."

;----------------------------------------------------------
RVhalo_tot[skip:skip+hnum[kk]-1] = RVhalo


openr,lun,'./'+fix4+'_'+'subhalo_smooth_cic_vel_field.dat',/get_lun

readu,lun,Ntot
RVField=fltarr(Ntot,6)
readu,lun,RVField
free_lun,lun

rv_00   = RVField[*,0] / afact
rv_0125 = RVField[*,1] / afact

rv_025  = RVField[*,2] / afact
rv_05   = RVField[*,3] / afact
rv_10   = RVField[*,4] / afact

openr,lun,'./'+fix4+'_'+'halo_smooth_cic_vel_field.dat',/get_lun

readu,lun,Ntot
SVField=fltarr(Ntot,6)
readu,lun,SVField
free_lun,lun

sv_00   = SVField[*,0] / afact
sv_0125 = SVField[*,1] / afact

sv_025  = SVField[*,2] / afact
sv_05   = SVField[*,3] / afact
sv_10   = SVField[*,4] / afact

; -----
rv_00_tot[skip:skip+hnum[kk]-1] = rv_00
rv_0125_tot[skip:skip+hnum[kk]-1] = rv_0125

rv_025_tot[skip:skip+hnum[kk]-1] = rv_025
rv_05_tot[skip:skip+hnum[kk]-1] =  rv_05
rv_10_tot[skip:skip+hnum[kk]-1] =  rv_10

sv_00_tot[skip:skip+hnum[kk]-1] = sv_00
sv_0125_tot[skip:skip+hnum[kk]-1] = sv_0125

sv_025_tot[skip:skip+hnum[kk]-1] = sv_025
sv_05_tot[skip:skip+hnum[kk]-1] =  sv_05
sv_10_tot[skip:skip+hnum[kk]-1] =  sv_10

skip += hnum[kk]
endfor

stop
;----------------------------------------------

plotps, /win, xsize=1000, ysize=600
;plotps, /ps, filename="relation_1.ps", xsize=800, ysize=600

xt1='\tex[c][B][1.2]{$v_{\mathrm{halo}}\ [\mathrm{km/s}]$}'
xt2='\tex[c][B][1.2]{$v_{\mathrm{CIC},\ r_{\mathrm{s}}=2.5\ h^{-1}\mathrm{Mpc}}\ [\mathrm{km/s}]$}'

yt1='\tex[c][B][1.2]{$v_{\mathrm{CIC}}\ [\mathrm{km/s}]$}'
yt2='\tex[c][B][1.2]{$v_{\mathrm{rec}}\ [\mathrm{km/s}]$}'

psc1='\tex[bl][bl][1.2]{$r_{\mathrm{S}}=2.5\ h^{-1}\mathrm{Mpc}$}'
psc2='\tex[bl][bl][1.2]{$r_{\mathrm{S}}=5\ h^{-1}\mathrm{Mpc}$}'
psc3='\tex[bl][bl][1.2]{$r_{\mathrm{S}}=10\ h^{-1}\mathrm{Mpc}$}'

;yt1='\tex[c][B][1.0]{$v_{\mathrm{CIC,smth=2.5} \mathrm{Mpc}/h}\ [\mathrm{km/s}]$}'
;yt2='\tex[c][B][1.0]{$v_{\mathrm{CIC,smth=5} \mathrm{Mpc}/h}\ [\mathrm{km/s}]$}'
;yt3='\tex[c][B][1.0]{$v_{\mathrm{CIC,smth=10} \mathrm{Mpc}/h}\ [\mathrm{km/s}]$}'
;
;yt4='\tex[c][B][1.0]{$v_{\mathrm{rec,smth=2.5} \mathrm{Mpc}/h}\ [\mathrm{km/s}]$}'
;yt5='\tex[c][B][1.0]{$v_{\mathrm{rec,smth=5} \mathrm{Mpc}/h}\ [\mathrm{km/s}]$}'
;yt6='\tex[c][B][1.0]{$v_{\mathrm{rec,smth=10} \mathrm{Mpc}/h}\ [\mathrm{km/s}]$}'

mfont=0
!p.font=mfont

cgerase &multiplot,/initialize
multiplot,[3,2],gap=0.01,/square;,/doxaxis,/doyaxis

bins=50
;----------------------
for iplot=0,5 do begin

;if iplot eq 0 or $
;   iplot eq 1 or $
;   iplot eq 2 or $
;   iplot eq 3 or $
;   iplot eq 4 or $
;   iplot eq 5 then begin vel1 = RVhalo_tot & xt=xt1 & endif
;
;if iplot eq 6 or $
;   iplot eq 7 or $
;   iplot eq 8 then begin vel1 = sv_025_tot & xt=xt2 & endif
;
;if iplot eq 0 then begin vel2 = sv_025_tot & yt=yt1 & pschar=psc1 & endif
;if iplot eq 1 then begin vel2 = sv_05_tot  & yt=yt1 & pschar=psc2 & endif
;if iplot eq 2 then begin vel2 = sv_10_tot  & yt=yt1 & pschar=psc3 & endif
;
;if iplot eq 3 or iplot eq 6 then begin vel2 = rv_025_tot & yt=yt2 & pschar=psc1 & endif
;if iplot eq 4 or iplot eq 7 then begin vel2 = rv_05_tot  & yt=yt2 & pschar=psc2 & endif
;if iplot eq 5 or iplot eq 8 then begin vel2 = rv_10_tot  & yt=yt2 & pschar=psc3 & endif

vel1 = RVhalo_tot

if iplot eq 0 then begin vel2 = sv_025_tot & xr=[-1000,1000] & yt=yt1 & pschar=psc1 & endif
if iplot eq 1 then begin vel2 = sv_05_tot  & xr=[-999.99,1000] & yt=''  & pschar=psc2 & endif
if iplot eq 2 then begin vel2 = sv_10_tot  & xr=[-999.99,1000] & yt=''  & pschar=psc3 & endif

if iplot eq 3 then begin vel2 = rv_025_tot & xr=[-1000,1000] & xt=''  & yt=yt2 & pschar=psc1 & endif
if iplot eq 4 then begin vel2 = rv_05_tot  & xr=[-999.99,1000] & xt=xt1 & yt=''  & pschar=psc2 & endif
if iplot eq 5 then begin vel2 = rv_10_tot  & xr=[-999.99,1000] & xt=''  & yt=''  & pschar=psc3 & endif

plot,[0],[0],xr=xr ,yr=[-1000,1000],xst=1,yst=1,xtitle=xt,ytitle=yt,charsize=0.9
make_2dcolor,vel1,vel2,bins,bins,/rho,/rev,ct=39,/topx
plots,[-1000,1000],[-1000,1000],thick=3.5,color=cgColor("red")

aa = mean2d(vel1, vel2, 10, xmin=-850, xmax=850)
 oplot, aa[*,0], aa[*,4], color=cgcolor("blue"), thick=3.5
 ;oplot, aa[*,0], aa[*,1], color=cgcolor("blue"), thick=3.5
 oplot, aa[*,0], aa[*,6], color=cgcolor("white"), thick=3.5, linestyle=0
 oplot, aa[*,0], aa[*,7], color=cgcolor("white"), thick=3.5, linestyle=0
 oplot, aa[*,0], aa[*,6], color=cgcolor("blue"), thick=3.5, linestyle=2
 oplot, aa[*,0], aa[*,7], color=cgcolor("blue"), thick=3.5, linestyle=2

 iis = where(aa[*,3] gt 0)
 scatter = mean((aa[*,7]-aa[*,6])/2.)
 print,scatter
 fit = linfit(aa[iis,0], aa[iis,1])
 print,fit

xsize=!p.position[2]-!p.position[0]
ysize=!p.position[3]-!p.position[1]
xp=!p.position[0]+0.1*xsize
;xp=!p.position[0]+0.025*xsize
yp=!p.position[1]+0.8*ysize

XYOUTS, xp, yp, pschar, /NORMAL, font=mfont, charsize=1.0

plc1='\tex[br][br][1.2]{$\mathrm{slope}='+strn(fit[1],LENGTH = 4)+'$}'
plc2='\tex[br][br][1.2]{$\mathrm{scatter}='+strn(scatter,LENGTH = 6)+'$}'
al_legend, [plc1,plc2],/bottom,/right,textcolor=[cgColor('black'),cgColor('black')],box=0, charsize=1.1

multiplot;,/doyaxis,/doxaxis
endfor

multiplot,/reset

;plotps,/ps,/close
;plotps,/win,/close

stop
;endfor

;------
endfun:
stop

end

