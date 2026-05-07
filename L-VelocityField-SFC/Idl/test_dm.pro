
ntask = 0l
nfiles = 0l
boxsize = 0.D
SI = 0l
nslab = 0l
fslab = 0l

fac = 35.3279
fd  = '~/data/halo_fields_054/0/dm_density_smooth_2048_054_0.0'
openr, 1, fd
readu, 1, ntask, nfiles, boxsize, SI, nslab, fslab
print, ntask, nfiles, boxsize, SI, nslab, fslab
dens  =fltarr(SI, SI, nslab)
readu, 1, dens
free_lun, 1

;fdo  = '~/data/halo_fields_054/dm_density_054.0'
fdo  = '~/data/smoothed_fields_054/0/field_2048_054_0.0'
openr, 1, fdo
readu, 1, ntask, nfiles, boxsize, SI, nslab, fslab
print, ntask, nfiles, boxsize, SI, nslab, fslab
denso  =fltarr(SI, SI, nslab)
readu, 1, denso
free_lun, 1


end
