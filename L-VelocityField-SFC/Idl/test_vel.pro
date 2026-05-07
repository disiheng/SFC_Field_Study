

Ntot = 0.
fix4 = "sdss_m5e3_mmax"
openr,lun,'./'+fix4+'_'+'halo_smooth_cic_vel_field.dat',/get_lun

readu,lun,Ntot
RVField=fltarr(Ntot,6)
readu,lun,RVField
free_lun,lun

sv_00   = RVField[*,0]
sv_0125 = RVField[*,1]

sv_025  = RVField[*,2]
sv_05   = RVField[*,3]
sv_10   = RVField[*,4]
free_lun,lun


openr,lun,'./'+fix4+'_'+'subhalo_smooth_cic_vel_field.dat',/get_lun
readu,lun,Ntot
sRVField=fltarr(Ntot,6)
readu,lun,sRVField
free_lun,lun

ssv_00   = sRVField[*,0]
ssv_0125 = sRVField[*,1]

ssv_025  = sRVField[*,2]
ssv_05   = sRVField[*,3]
ssv_10   = sRVField[*,4]
free_lun,lun


openr,lun,'./'+fix4+'_'+'halo_true_vel.dat',/get_lun
readu,lun,Ntot
RVhalo=fltarr(Ntot)
readu,lun,RVhalo
free_lun,lun

openr,lun,'/u/mingli/data/Project.3.Stacking/smooth_cic_vel/'+fix4+'_'+'smooth_cic_vel_field.dat',/get_lun
readu,lun,Ntot
dRVField=fltarr(Ntot,6)
readu,lun,dRVField
free_lun,lun

dsv_00   = dRVField[*,0]
dsv_0125 = dRVField[*,1]

dsv_025  = dRVField[*,2]
dsv_05   = dRVField[*,3]
dsv_10   = dRVField[*,4]
free_lun,lun
end
