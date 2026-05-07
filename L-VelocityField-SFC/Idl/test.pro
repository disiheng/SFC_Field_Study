
ntask = 0l
nfiles = 0l
boxsize = 0.D
SI = 0l
nslab = 0l
fslab = 0l

fv  = '~/data/halo_fields_054/3/halo_density_smooth_2048_054_3.0'
openr, 1, fv
readu, 1, ntask, nfiles, boxsize, SI, nslab, fslab
print, ntask, nfiles, boxsize, SI, nslab, fslab
dens  =fltarr(SI, SI, nslab)
readu, 1, dens
free_lun, 1

fv  = '~/data/halo_fields_054/3/halo_velocity_smooth_0_2048_054_3.0'
openr, 1, fv
readu, 1, ntask, nfiles, boxsize, SI, nslab, fslab
print, ntask, nfiles, boxsize, SI, nslab, fslab
momt  =fltarr(SI, SI, nslab)
readu, 1, momt
free_lun, 1

vel = momt/dens

fv  = '~/data/smoothed_fields_054/3/field_2048_054_3.0'
openr, 1, fv
readu, 1, ntask, nfiles, boxsize, SI, nslab, fslab
print, ntask, nfiles, boxsize, SI, nslab, fslab
dens1  =fltarr(SI, SI, nslab)
readu, 1, dens1
free_lun, 1

fv  = '~/data/smoothed_fields_054/3/vel_0_2048_054_3.0'
openr, 1, fv
readu, 1, ntask, nfiles, boxsize, SI, nslab, fslab
print, ntask, nfiles, boxsize, SI, nslab, fslab
momt1  =fltarr(SI, SI, nslab)
readu, 1, momt1
free_lun, 1

vel1 = momt1 /dens1
stop

fix4 = "sdss_m5e3_mmax"
openr,lun,'./'+fix4+'_'+'halo_rec_cic_vel_field.dat',/get_lun

readu,lun,Ntot
RVField=fltarr(Ntot,6)
readu,lun,RVField
free_lun,lun

rv_00   = RVField[*,0]
rv_0125 = RVField[*,1]

rv_025  = RVField[*,2]
rv_05   = RVField[*,3]
rv_10   = RVField[*,4]
free_lun,lun

openr,lun,'/u/mingli/data/backup/SimuCode/L-VelocityField-MXXL/Idl/'+fix4+'_'+'rec_cic_vel_field.dat',/get_lun

readu,lun,Ntot
oRVField=fltarr(Ntot,6)
readu,lun,oRVField
free_lun,lun

orv_00   = oRVField[*,0]
orv_0125 = oRVField[*,1]

orv_025  = oRVField[*,2]
orv_05   = oRVField[*,3]
orv_10   = oRVField[*,4]
free_lun,lun


end
