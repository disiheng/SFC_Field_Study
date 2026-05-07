#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>
#include <assert.h>

#include "allvars.h"
#include "proto.h"


void density(Particle * P, fftw_real *grid, long long NumPart, int mode, int value)
{
  int rep;
  long long nstart, npart;

  if( value == 0 ) MSG0(" Computing Density Field", t0, ThisTask);
  if( value == 1 ) MSG0(" Computing Vx Field", t0, ThisTask);
  if( value == 2 ) MSG0(" Computing Vy Field", t0, ThisTask);
  if( value == 3 ) MSG0(" Computing Vz Field", t0, ThisTask);

  //nstart = 0;
  //npart = NumPart;

  for(rep = 0; rep < CYCLES; rep++)
    {
      if(rep == (CYCLES - 1))
		npart = NumPart - (CYCLES - 1) * (NumPart / CYCLES);
      else
		npart = (NumPart / CYCLES);

      nstart = rep * (NumPart / CYCLES);

      if(ThisTask == 0)
		printf("CYCLE %d-%d nstart=%ld npart=%ld, %ld\n", rep, CYCLES, nstart, npart, NumPart);

      if(ThisTask == 0)
        printf("Mode 0\n");
      mapping(P, grid, nstart, npart, 0, value);

#ifdef CIC
      if(ThisTask == 0)
        printf("Mode 1\n");
      mapping(P, grid, nstart, npart, 1, value);
#endif

    }
//    if( value == 0 )
//      density_statistics(Grid);

}


void mapping(Particle * P, fftw_real *grid, long long FirstPart, long long NumPart, int mode, int value)
{   
  int i, j, k, itask, send_task, jpos;
  long long ii, kk;
  long *counter, *send_size, *rec_size;
  long **send_buffer, **rec_buffer;
  float **send_weight_buffer, **rec_weight_buffer, **send_pos_buffer, **rec_pos_buffer;
  double weight, dx,dy,dz;
  //float weight, dx,dy,dz;
  int i1,j1,k1;
  char buf[255];
  MPI_Status mpi_status;
  
  counter = mymalloc("counter", NTask * sizeof(long));
  send_size = mymalloc("send_size", NTask * sizeof(long));
  rec_size = mymalloc("rec_size", NTask * sizeof(long));
  
  for(ii = 0; ii < NTask; ii++)
    {
      rec_size[ii] = 0;
      send_size[ii] = 0;
      counter[ii] = 0;
    } 
      
  if (ThisTask == 0)
      printf("First %lld -- NumPart %lld \n", FirstPart, NumPart);
  for(ii = FirstPart; ii < FirstPart + NumPart; ii++)
    {
      i = (to_slab_fac * P[ii].pos[0]);
      if (mode == 0) 
        i = (i + Ndim) % Ndim;
      if (mode == 1)  
        i = (i + 1 + Ndim) % Ndim;
      
      if ( i < 0 || i > Ndim)
        printf("ii=%lld pos=%f grid=%d\n",ii,P[ii].pos[0], (int) Ndim);

      assert(i >= 0 && i < Ndim);
      itask = slab_to_task[i];
      
      assert(itask >= 0 && itask < NTask);
      send_size[itask]++;
    } 
      
  for(ii = 0; ii < NTask; ii++)
    if(send_size[ii] > 0)
      DBG("Task:%d sending:%d particles to task %d", t0, dbg2, ThisTask, send_size[ii], ii);
      
  DBG(" [Task=%d] We allocate a buffer to send to the processors", dbg2, t0, ThisTask);

  long test = 0; 
  for(ii = 0; ii < NTask; ii++)
      test += send_size[ii];
  assert(test = NumPart);


  rec_buffer = (long **) mymalloc("rec_buffer", NTask * sizeof(*rec_buffer));
  send_buffer = (long **) mymalloc("send_buffer", NTask * sizeof(*send_buffer));
  
  rec_weight_buffer = (float **) mymalloc("rec_weight_buffer", NTask * sizeof(*rec_weight_buffer));
  send_weight_buffer = (float **) mymalloc("send_weight_buffer", NTask * sizeof(*send_weight_buffer));
  
  rec_pos_buffer = (float **) mymalloc("rec_pos_buffer", NTask * sizeof(*rec_pos_buffer));
  send_pos_buffer = (float **) mymalloc("send_pos_buffer", NTask * sizeof(*send_pos_buffer));
  
  for(ii = 0; ii < NTask; ii++)
    send_buffer[ii] = (long *) mymalloc("send_buffer[ii]", MAX(send_size[ii], 1) * sizeof(long));
  
  for(ii = 0; ii < NTask; ii++)
    send_weight_buffer[ii] = (float *) mymalloc("send_weight_buffer[ii]", MAX(send_size[ii], 1) * sizeof(float));
  
  for(ii = 0; ii < NTask; ii++)
    send_pos_buffer[ii] = (float *) mymalloc("send_pos_buffer[ii]", 3 * MAX(send_size[ii], 1) * sizeof(float));

  DBG(" [Task=%d] Filling that buffer", dbg2, t0, ThisTask);

  for(ii = FirstPart; ii < FirstPart + NumPart; ii++)
    {
      if (value == 0 ) weight = 1.0;
      if (value == 1) weight = P[ii].vel[0];
      if (value == 2) weight = P[ii].vel[1];
      if (value == 3) weight = P[ii].vel[2];

      i = (int) (to_slab_fac * P[ii].pos[0]);
      j = (int) (to_slab_fac * P[ii].pos[1]);
      k = (int) (to_slab_fac * P[ii].pos[2]);

      i = (i + Ndim) % Ndim;
      j = (j + Ndim) % Ndim;
      k = (k + Ndim) % Ndim;

      if (mode == 0)
        i = (i + Ndim) % Ndim;
      else
        i = (i + 1 + Ndim) % Ndim;

      itask = slab_to_task[i];

      assert(itask >= 0 && itask < NTask);
      i -= first_slab_of_task[itask];
      //assert(i >= 0 && i < slabs_per_task[itask]);

      send_buffer[itask][counter[itask]] = k + Ndim * (j + Ndim * i);
      send_weight_buffer[itask][counter[itask]] = weight;
      for (jpos=0; jpos<3; jpos++)
        send_pos_buffer[itask][3*counter[itask] + jpos] = P[ii].pos[jpos];
      counter[itask]++;
    }

  for(ii = 0; ii < NTask; ii++)
    assert(counter[ii] == send_size[ii]);

  // Sending and rceiving buffers
  DBG(" [Task=%d] Now we can send the buffers", dbg2, t0, ThisTask);
  for(ii = 0; ii < NTask; ii++)
    {
      send_task = ii;
      if(send_task == ThisTask)
	  {
		for(kk = 0; kk < NTask; kk++)
		  {
			if(kk != ThisTask)
			{
			  MPI_Send(&send_size[kk], 1, MPI_LONG, kk, TAG_N, MPI_COMM_WORLD);
			  if(send_size[kk] > 0)
			  {
				MPI_Send(send_buffer[kk], send_size[kk], MPI_LONG, kk, TAG_PDATA, MPI_COMM_WORLD);
				MPI_Send(send_weight_buffer[kk],  send_size[kk], SEND_TYPE, kk, TAG_MDATA, MPI_COMM_WORLD);
				MPI_Send(send_pos_buffer[kk], 3 * send_size[kk], SEND_TYPE, kk, TAG_POSDATA, MPI_COMM_WORLD);
			  }
			}
		  }
		// to be consistent
		//rec_buffer[ThisTask] = (long*) mymalloc( "rec_buffer[ii]",1 * sizeof(long) );
	  }
      else
	  {
		MPI_Recv(&rec_size[send_task], 1, MPI_LONG, send_task, TAG_N, MPI_COMM_WORLD, &mpi_status);
		rec_buffer[ii] = (long *) mymalloc("rec_buffer[ii]", MAX(rec_size[send_task], 1) * sizeof(long));
		rec_weight_buffer[ii] = (float *) mymalloc("rec_weight_buffer[ii]", MAX(rec_size[send_task], 1) * sizeof(float));
		rec_pos_buffer[ii] = (float *) mymalloc("rec_pos_buffer[ii]", 3 * MAX(rec_size[send_task], 1) * sizeof(float));

		if(rec_size[ii] > 0)
		{
		  DBG("Task:%d received %d particles from:%d", dbg2, t0, ThisTask, rec_size[ii], ii);
		  MPI_Recv(rec_buffer[ii], rec_size[ii], MPI_LONG, send_task, TAG_PDATA, MPI_COMM_WORLD, &mpi_status);
		  MPI_Recv(rec_weight_buffer[ii], rec_size[ii], SEND_TYPE, send_task, TAG_MDATA, MPI_COMM_WORLD, &mpi_status);
		  MPI_Recv(rec_pos_buffer[ii], 3 * rec_size[ii], SEND_TYPE, send_task, TAG_POSDATA, MPI_COMM_WORLD, &mpi_status);

		  for(kk = 0; kk < rec_size[ii]; kk++)
			{
#ifdef CIC
			  i = (int) (rec_pos_buffer[ii][3*kk+0] * to_slab_fac);
			  j = (int) (rec_pos_buffer[ii][3*kk+1] * to_slab_fac);
			  k = (int) (rec_pos_buffer[ii][3*kk+2] * to_slab_fac);

			  dx = rec_pos_buffer[ii][3*kk+0] * to_slab_fac - i;
			  dy = rec_pos_buffer[ii][3*kk+1] * to_slab_fac - j;
			  dz = rec_pos_buffer[ii][3*kk+2] * to_slab_fac - k;

			  i = (i + Ndim) % Ndim;
			  j = (j + Ndim) % Ndim;
			  k = (k + Ndim) % Ndim;

			  i -= first_slab_of_task[ThisTask];

			  i1 = i + 1;
			  j1 = j + 1;
			  k1 = k + 1;

			  if(i1 >= Ndim) i1 -= Ndim;
			  if(j1 >= Ndim) j1 -= Ndim;
			  if(k1 >= Ndim) k1 -= Ndim;

			  weight = rec_weight_buffer[ii][kk];

			  if(mode == 0)
			  {
				grid[k  + Ndim * (j  + Ndim * i)] += weight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
				grid[k  + Ndim * (j1 + Ndim * i)] += weight * (1.0 - dx) * dy * (1.0 - dz);
				grid[k1 + Ndim * (j  + Ndim * i)] += weight * (1.0 - dx) * (1.0 - dy) * dz;
				grid[k1 + Ndim * (j1 + Ndim * i)] += weight * (1.0 - dx) * dy * dz;
			  }

			  if(mode == 1)
			  {
				grid[k  + Ndim * (j  + Ndim * i1)] += weight * dx * (1.0 - dy) * (1.0 - dz);
				grid[k  + Ndim * (j1 + Ndim * i1)] += weight * dx * dy * (1.0 - dz);
				grid[k1 + Ndim * (j  + Ndim * i1)] += weight * dx * (1.0 - dy) * dz;
				grid[k1 + Ndim * (j1 + Ndim * i1)] += weight * dx * dy * dz;
			  }
              //assert(rec_buffer[ii][kk] == k + Ndim * (j + Ndim * i));
#else
			  grid[rec_buffer[ii][kk]] += (float) rec_weight_buffer[ii][kk];
#endif
			}
        }

		myfree(rec_pos_buffer[ii]);
		myfree(rec_weight_buffer[ii]);
		myfree(rec_buffer[ii]);
	  }
    }

  for(kk = 0; kk < send_size[ThisTask]; kk++)
  {
	ii = ThisTask;
#ifdef CIC
    i = (int) (send_pos_buffer[ii][3*kk+0] * to_slab_fac);
    j = (int) (send_pos_buffer[ii][3*kk+1] * to_slab_fac);
    k = (int) (send_pos_buffer[ii][3*kk+2] * to_slab_fac);

    dx = send_pos_buffer[ii][3*kk+0]  * to_slab_fac - i;
    dy = send_pos_buffer[ii][3*kk+1]  * to_slab_fac - j;
    dz = send_pos_buffer[ii][3*kk+2]  * to_slab_fac - k;

    i = (i + Ndim) % Ndim;
    j = (j + Ndim) % Ndim;
    k = (k + Ndim) % Ndim;

    i -= first_slab_of_task[ThisTask];

    i1 = i + 1;
    if(i1 >= Ndim) i1 -= Ndim;
    j1 = j + 1;
    if(j1 >= Ndim) j1 -= Ndim;
    k1 = k + 1;
    if(k1 >= Ndim) k1 -= Ndim;

    weight = send_weight_buffer[ii][kk];

	if(mode == 0)
	{
	  grid[k  + Ndim * (j  + Ndim * i)] += weight * (1.0 - dx) * (1.0 - dy) * (1.0 - dz);
	  grid[k  + Ndim * (j1 + Ndim * i)] += weight * (1.0 - dx) * dy * (1.0 - dz);
	  grid[k1 + Ndim * (j  + Ndim * i)] += weight * (1.0 - dx) * (1.0 - dy) * dz;
	  grid[k1 + Ndim * (j1 + Ndim * i)] += weight * (1.0 - dx) * dy * dz;
	}

	if(mode == 1)
	{
	  grid[k  + Ndim * (j  + Ndim * i1)] += weight * dx * (1.0 - dy) * (1.0 - dz);
	  grid[k  + Ndim * (j1 + Ndim * i1)] += weight * dx * dy * (1.0 - dz);
	  grid[k1 + Ndim * (j  + Ndim * i1)] += weight * dx * (1.0 - dy) * dz;
	  grid[k1 + Ndim * (j1 + Ndim * i1)] += weight * dx * dy * dz;
	}
#else
	grid[send_buffer[ii][kk]] += (float) send_weight_buffer[ii][kk];
#endif
  }

  //report_memory_usage(&HighMark, "DENSITY");

  DBG0(" Freeing  buffers", dbg1, t0, ThisTask);

  for(ii = NTask - 1; ii >= 0; ii--)
    myfree(send_pos_buffer[ii]);

  for(ii = NTask - 1; ii >= 0; ii--)
    myfree(send_weight_buffer[ii]);

  for(ii = NTask - 1; ii >= 0; ii--)
    myfree(send_buffer[ii]);

  myfree(send_pos_buffer);
  myfree(rec_pos_buffer);

  myfree(send_weight_buffer);
  myfree(rec_weight_buffer);

  myfree(send_buffer);
  myfree(rec_buffer);

  myfree(rec_size);
  myfree(send_size);
  myfree(counter);
}


// --------------------------------------------- //
/* 初始化MRA参数 */
void init_mra_parameters()
{
    if (ThisTask == 0) {
        printf("Initializing MRA parameters...\n");
    }
    
    /* 从配置文件或全局变量获取参数 */
    int mra_base_type = 0;      /* 0: B-Spline, 1: Daubechies */
    int mra_genus = 3;          /* 尺度函数阶数 */
    SampRate = 1000;        /* 固定采样率 */
    
    /* 计算尺度函数 */
    if (mra_base_type == 0) {
        /* B-Spline */
        phiSupport = phiGenus + 1;
        phi_data = compute_b_spline(phiGenus, SampRate, &phiSupport);
    } else {
        /* Daubechies */
        phiSupport = 2 * phiGenus - 1;
        phi_data = compute_daubechies(phiGenus, SampRate, &phiSupport);
    }
    
    /* 计算缩放因子 */
    scale_factor_mra = Ndim / BoxSize;
    
    if (ThisTask == 0) {
        printf("MRA initialized: support=%d, samprate=%d\n", phiSupport, SampRate);
    }
}


/* B-Spline尺度函数计算（C语言版本） */
double compute_b_spline(int n, int sample_rate)
{
    int total_points = (n + 1) * sample_rate + 1;
    double *spline = (double *)mymalloc("b_spline", total_points * sizeof(double));
    double delta_x = 1.0 / sample_rate;
    
    /* 初始化0阶B-Spline */
    for (int i = 1; i < sample_rate; i++)
    {
        spline[i] = 1.0;
    }

    spline[0] = 0.5;
    spline[sample_rate] = 0.5;
    
    /* 递归计算高阶B-Spline */
    for (int m = 1; m <= n; m++)
    {
        for (int i = 0; i <= m * sample_rate; i++)
        {
            spline[i] *= i * delta_x / m;
        }

        for (int i = sample_rate; i <= (m + 1) * sample_rate / 2; i++)
        {
            spline[i] += spline[(m + 1) * sample_rate - i];
        }

        for (int i = 0; i <= (m + 1) * sample_rate / 2; i++)
        {
            spline[(m + 1) * sample_rate - i] = spline[i];
        }
    }
    
    return spline;
}

/* Daubechies尺度函数（需要预计算的数据） */
void get_daubechies(int genus, int sample_rate)
{
    char fname[1000], buf[500];
    FILE *fd;

    int total_points = phiSupport * SampRate + 1;
    phi = mymalloc("daubechies", total_points * sizeof(double));

    if (ThisTask == 0)
    {
        sprintf(fname, "./wavelets/binary/DaubechiesG%dPhi.bin", phiGenus);
        fd = my_fopen(fname, "r");
        my_fread(phi , sizeof(double), total_points, fd);
        fclose(fd);
    }

    MPI_Bcast(phi, sizeof(double)*total_points, MPI_BYTE, 0, MPI_COMM_WORLD);
    
    return 0;
}


/* MRA版本的密度场计算 */
void density_mra(Particle *P, fftw_real *grid, long long NumPart, int mode, int value)
{
    int rep;
    long long nstart, npart;
    
    if (value == 0) MSG0(" Computing Density Field (MRA Method)", t0, ThisTask);
    if (value == 1) MSG0(" Computing Vx Field (MRA)", t0, ThisTask);
    if (value == 2) MSG0(" Computing Vy Field (MRA)", t0, ThisTask);
    if (value == 3) MSG0(" Computing Vz Field (MRA)", t0, ThisTask);
    
    /* MRA需要更细粒度的批次处理 */
    int CYCLES_MRA = CYCLES * 4;  /* 增加批次数量 */
    
    for (rep = 0; rep < CYCLES_MRA; rep++) {
        if (rep == (CYCLES_MRA - 1))
            npart = NumPart - (CYCLES_MRA - 1) * (NumPart / CYCLES_MRA);
        else
            npart = (NumPart / CYCLES_MRA);
        
        nstart = rep * (NumPart / CYCLES_MRA);
        
        if (ThisTask == 0)
            printf("MRA CYCLE %d-%d nstart=%lld npart=%lld\n", 
                   rep, CYCLES_MRA, nstart, npart);
        
        mapping_mra(P, grid, nstart, npart, mode, value);
    }
}

/* MRA版本的映射函数 - 核心实现 */
void mapping_mra(Particle *P, fftw_real *grid, long long FirstPart, 
                 long long NumPart, int mode, int value)
{
    int i, j, k, itask, send_task, di, dj, dk;
    long long ii, kk;
    long *counter, *send_size, *rec_size;
    long   **send_buffer, **rec_buffer;
    double **send_weight_buffer, **rec_weight_buffer;
    MPI_Status mpi_status;

    int grid_x, grid_y, grid_z;
    double xx, yy, zz;
    int xxc, yyc, zzc, xxf, yyf, zzf;
    
    /* 分配计数器 */
    counter = mymalloc("counter", NTask * sizeof(long));
    send_size = mymalloc("send_size", NTask * sizeof(long));
    rec_size = mymalloc("rec_size", NTask * sizeof(long));
    
    for (ii = 0; ii < NTask; ii++) {
        rec_size[ii] = 0;
        send_size[ii] = 0;
        counter[ii] = 0;
    }
    
    if (ThisTask == 0)
       printf("First %lld -- NumPart %lld \n", FirstPart, NumPart);

    /* 计算每个粒子影响的所有slab */
    for (ii = FirstPart; ii < FirstPart + NumPart; ii++) {
        /* 坐标缩放 */
        xx  = P[ii].pos[0] * scale_factor_mra;
        xxc = (int)floor(xx);
        
        /* 确定这个粒子会影响哪些slab */
        for (di = 0; di < phiSupport; di++) {
            grid_x = (xxc - di + Ndim) % Ndim;
            itask = slab_to_task[grid_x];
            
            /* 确保slab在有效范围内 */
            assert(itask >= 0 && itask < NTask);
            send_size[itask]++;
        }
    }
    
    for(ii = 0; ii < NTask; ii++)
      if(send_size[ii] > 0)
        DBG("Task:%d sending:%d vs. %d particles to task %d", t0, dbg2, ThisTask, send_size[ii]/phiSupport, NumPart, ii);

    DBG(" [Task=%d] We allocate a buffer to send to the processors", dbg2, t0, ThisTask);

    /* 分配发送缓冲区 */
    rec_buffer         = (long   **)mymalloc("rec_buffer",         NTask * sizeof(*rec_buffer));
    send_buffer        = (long   **)mymalloc("send_buffer",        NTask * sizeof(*send_buffer));

    rec_weight_buffer  = (double **)mymalloc("rec_weight_buffer",  NTask * sizeof(*rec_weight_buffer));
    send_weight_buffer = (double **)mymalloc("send_weight_buffer", NTask * sizeof(*send_weight_buffer));

    for (ii = 0; ii < NTask; ii++) 
        send_buffer[ii]        = (long   *)mymalloc("send_buffer[ii]",        MAX(send_size[ii], 1) * sizeof(long));

    for (ii = 0; ii < NTask; ii++)
        send_weight_buffer[ii] = (double *)mymalloc("send_weight_buffer[ii]", MAX(send_size[ii], 1) * sizeof(double));
    
    DBG(" [Task=%d] Filling that buffer", dbg2, t0, ThisTask);

    /* 填充缓冲区 */
    for (ii = FirstPart; ii < FirstPart + NumPart; ii++) {
        double weight = 1.0;
        if (value == 1) weight = P[ii].vel[0];
        if (value == 2) weight = P[ii].vel[1];
        if (value == 3) weight = P[ii].vel[2];
        
        xx = P[ii].pos[0] * scale_factor_mra;
        yy = P[ii].pos[1] * scale_factor_mra;
        zz = P[ii].pos[2] * scale_factor_mra;
        
        xxc = (int)floor(xx);
        yyc = (int)floor(yy);
        zzc = (int)floor(zz);
        
        xxf = (int)(SampRate * (xx - xxc));
        yyf = (int)(SampRate * (yy - yyc));
        zzf = (int)(SampRate * (zz - zzc));
        
        /* 为每个受影响的slab计算贡献 */
        for (int di = 0; di < phiSupport; di++) {
            for (int dj = 0; dj < phiSupport; dj++) {
                for (int dk = 0; dk < phiSupport; dk++) {
                    grid_x = (xxc - di + Ndim) % Ndim;
                    grid_y = (yyc - dj + Ndim) % Ndim;
                    grid_z = (zzc - dk + Ndim) % Ndim;
                    
                    itask = slab_to_task[slab_x];
                    assert(itask >= 0 && itask < NTask);

                    slab_x = grid_x - first_slab_of_task[ThisTask]; 

                    /* 计算权重 */
                    phi_val = phi_data[xxf + di * SampRate] *
                              phi_data[yyf + dj * SampRate] *
                              phi_data[zzf + dk * SampRate];
                        
                    contrib = weight * phi_val;
                        
                    /* 存储到缓冲区 */
                    int idx = counter[itask];
                    send_buffer[itask][idx] = grid_z + Ndim * (grid_y + Ndim * slab_x);
                    send_weight_buffer[itask][idx] = contrib;
                    counter[itask]++;
                    }
                }
            }
    }
    
    /* 通信：发送和接收数据 */
    DBG(" [Task=%d] Now we can send the buffers", dbg2, t0, ThisTask);
    for (ii = 0; ii < NTask; ii++) {
        send_task = ii;
        
        if (send_task == ThisTask) {
            /* 发送数据到其他进程 */
            for (kk = 0; kk < NTask; kk++) {
                if (kk != ThisTask) 
                {
                    MPI_Send(&send_size[kk], 1, MPI_LONG, kk, TAG_N, MPI_COMM_WORLD);
                    if (send_size[kk] > 0) 
                    {
                        MPI_Send(send_buffer[kk],        send_size[kk], MPI_LONG  , kk, TAG_PDATA, MPI_COMM_WORLD);
                        MPI_Send(send_weight_buffer[kk], send_size[kk], MPI_DOUBLE, kk, TAG_MDATA, MPI_COMM_WORLD);
                    }
                }
            }
        } else {
            /* 接收其他进程的数据 */
            MPI_Recv(&rec_size[send_task], 1, MPI_LONG, send_task, TAG_N, MPI_COMM_WORLD, &mpi_status);
            rec_buffer[ii]        = (long   *)mymalloc("rec_buffer[ii]",        MAX(rec_size[send_task], 1) * sizeof(long));
            rec_weight_buffer[ii] = (double *)mymalloc("rec_weight_buffer[ii]", MAX(rec_size[send_task], 1) * sizeof(double));
            
            if (rec_size[ii] > 0) {
                DBG("Task:%d received %d particles from:%d", dbg2, t0, ThisTask, rec_size[ii], ii);
                
                MPI_Recv(rec_buffer       [ii], rec_size[ii], MPI_LONG  , send_task, TAG_PDATA, MPI_COMM_WORLD, &mpi_status);
                MPI_Recv(rec_weight_buffer[ii], rec_size[ii], MPI_DOUBLE, send_task, TAG_MDATA, MPI_COMM_WORLD, &mpi_status);
                
                /* 处理接收到的数据 */
                for (kk = 0; kk < rec_size[ii]; kk++) {
                    idx = rec_buffer[ii][kk];
                    grid[idx] += rec_weight_buffer[ii][kk];
                    }
                }
                
            myfree(rec_weight_buffer[ii]);
            myfree(rec_buffer[ii]);
            }
    }
    
    /* 处理本地数据 */
    for (kk = 0; kk < send_size[ThisTask]; kk++) {
        ii = ThisTask;
        
        grid_idx = send_buffer[ii][kk];
        grid[grid_idx] += send_weight_buffer[ii][kk];
    }
   
    //report_memory_usage(&HighMark, "DENSITY");
    DBG0(" Freeing  buffers", dbg1, t0, ThisTask);

    /* 清理内存 */
    for(ii = NTask - 1; ii >= 0; ii--)
        myfree(send_weight_buffer[ii]);

    for(ii = NTask - 1; ii >= 0; ii--)
        myfree(send_buffer[ii]);
    
    myfree(send_weight_buffer);
    myfree(rec_weight_buffer);

    myfree(send_buffer);
    myfree(rec_buffer);

    myfree(rec_size);
    myfree(send_size);
    myfree(counter);
}


