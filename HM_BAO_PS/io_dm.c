
void load_dm_single( int nr)
{
  int j=0, i, k, rep, thisN, num_files, ngroup, totngroup, ntrees, *galsPerTree, dummy;
  char buf[512];
  FILE *fd;
  float *halo_read;
  io_header header;
  gsl_rng *random_generator;

  gsl_rng_env_setup();
  random_generator = gsl_rng_alloc(gsl_rng_default);

  sprintf(buf, "%s/%s.%d", SimulationDir, SnapshotFileBase, nr);
  if (SnapNum != -1)
    sprintf(buf, "%s/snapdir_%03d/%s_%03d.%d", SimulationDir, SnapNum, SnapshotFileBase, SnapNum, nr);
   
  printf("%s\n",buf);
  fd = my_fopen(buf,"r");
  my_fread(&dummy, sizeof(int),1,fd);
  my_fread(&header, sizeof(io_header),1,fd);
  my_fread(&dummy, sizeof(int),1,fd);

  ngroup = header.npart[1];

  if(PP == NULL)
    PP = mymalloc_movable(&PP, "PP", (LocNgroups + ngroup) * sizeof(Particle));
  else
    PP = myrealloc_movable(PP, (LocNgroups + ngroup) * sizeof(Particle));

  my_fread(&dummy, sizeof(int),1,fd);
   for (rep=0; rep < CYCLES; rep++)
   {
      if(rep == (CYCLES - 1))
        thisN = ngroup - (CYCLES - 1) * (ngroup / CYCLES);
      else
        thisN = (ngroup / CYCLES);

     halo_read = mymalloc_movable(&halo_read,"halo_read", 3 * thisN* sizeof(float));
     my_fread(halo_read, thisN * 3 * sizeof(float),1,fd);

     for(i = 0; i < thisN; i++)
     {
       for(k=0; k<3; k++)
       {
         PP[LocNgroups].pos[k] = halo_read[3 * i + k];

         if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
         if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
  //       if (PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
  //         printf("(%f|%f|%f) Box=%f\n",PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize);

         //assert(PP[LocNgroups].pos[k] >= 0 && PP[LocNgroups].pos[k] < BoxSize);
       }
        LocNgroups++;
     }
     myfree_movable(halo_read);
   }
  my_fread(&dummy, sizeof(int),1,fd);
  
#if defined(Z_SPACE) || defined(DIVERGENCE)
   my_fread(&dummy, sizeof(int),1,fd);
   LocNgroups -= ngroup;
   for (rep=0; rep < CYCLES; rep++)
   {
      if(rep == (CYCLES - 1))
        thisN = ngroup - (CYCLES - 1) * (ngroup / CYCLES);
      else
        thisN = (ngroup / CYCLES);

     halo_read = mymalloc_movable(&halo_read,"halo_read", 3 * thisN* sizeof(float));
     my_fread(halo_read, thisN * 3 * sizeof(float),1,fd);

     for(i = 0; i < thisN; i++)
     {
        for(k=0; k<3; k++)
        {
#ifdef Z_SPACE
          if (k == 2) 
            PP[LocNgroups].pos[k] +=  halo_read[3*i + k] / Hubble / sqrt(ExpFactor); 
#endif
#ifdef DIVERGENCE
          PP[LocNgroups].vel[k] =  halo_read[3*i + k] / Hubble / sqrt(ExpFactor); 
#endif

        if (PP[LocNgroups].pos[k] >= BoxSize) PP[LocNgroups].pos[k] -= BoxSize;
        if (PP[LocNgroups].pos[k] < 0)  PP[LocNgroups].pos[k] += BoxSize;       
        //if (PP[LocNgroups].pos[k] < 0 || PP[LocNgroups].pos[k] >= BoxSize)
        //  printf("2: (%f|%f|%f) Box=%f\n",PP[LocNgroups].pos[0],PP[LocNgroups].pos[1],PP[LocNgroups].pos[2], BoxSize);
       } 
       LocNgroups++;
     }
     myfree_movable(halo_read);
   }
  my_fread(&dummy, sizeof(int),1,fd);
#endif
  fclose(fd);

  gsl_rng_free(random_generator);

}

