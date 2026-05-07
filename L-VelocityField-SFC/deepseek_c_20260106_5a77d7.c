/* ============================ MRA网格化函数（C语言风格） ============================ */

#include "allva
#include "proto.h"

/* 全局MRA参数 */
static double *phi_data = NULL;        /* 尺度函数数据 */
static int phi_support = 0;            /* 尺度函数支撑长度 */
static int phi_samprate = 1000;        /* 采样率 */
static double scale_factor_mra = 1.0;  /* 坐标缩放因子 */

/* 初始化MRA参数 */
void init_mra_parameters()
{
    if (ThisTask == 0) {
        printf("Initializing MRA parameters...\n");
    }
    
    /* 从配置文件或全局变量获取参数 */
    int mra_base_type = 0;      /* 0: B-Spline, 1: Daubechies */
    int mra_genus = 3;          /* 尺度函数阶数 */
    phi_samprate = 1000;        /* 固定采样率 */
    
    /* 计算尺度函数 */
    if (mra_base_type == 0) {
        /* B-Spline */
        phi_support = mra_genus + 1;
        phi_data = compute_b_spline(mra_genus, phi_samprate, &phi_support);
    } else {
        /* Daubechies */
        phi_support = 2 * mra_genus - 1;
        phi_data = compute_daubechies(mra_genus, phi_samprate, &phi_support);
    }
    
    /* 计算缩放因子 */
    scale_factor_mra = GridLen / BoxSize;
    
    if (ThisTask == 0) {
        printf("MRA initialized: support=%d, samprate=%d\n", 
               phi_support, phi_samprate);
    }
}

/* B-Spline尺度函数计算（C语言版本） */
double* compute_b_spline(int n, int sample_rate, int *support)
{
    int total_points = (n + 1) * sample_rate + 1;
    double *spline = (double *)mymalloc("b_spline", total_points * sizeof(double));
    double delta_x = 1.0 / sample_rate;
    
    /* 初始化0阶B-Spline */
    for (int i = 1; i < sample_rate; i++) {
        spline[i] = 1.0;
    }
    spline[0] = 0.5;
    spline[sample_rate] = 0.5;
    
    /* 递归计算高阶B-Spline */
    for (int m = 1; m <= n; m++) {
        for (int i = 0; i <= m * sample_rate; i++) {
            spline[i] *= i * delta_x / m;
        }
        for (int i = sample_rate; i <= (m + 1) * sample_rate / 2; i++) {
            spline[i] += spline[(m + 1) * sample_rate - i];
        }
        for (int i = 0; i <= (m + 1) * sample_rate / 2; i++) {
            spline[(m + 1) * sample_rate - i] = spline[i];
        }
    }
    
    *support = n + 1;
    return spline;
}

/* Daubechies尺度函数（需要预计算的数据） */
double* compute_daubechies(int genus, int sample_rate, int *support)
{
    /* 这里简化处理，实际应该从文件读取预计算的Daubechies函数 */
    *support = 2 * genus - 1;
    int total_points = (*support) * sample_rate + 1;
    double *phi = (double *)mymalloc("daubechies", total_points * sizeof(double));
    
    /* 简单初始化（实际应用需要正确的Daubechies系数） */
    for (int i = 0; i < total_points; i++) {
        phi[i] = 1.0; /* 占位符 */
    }
    
    return phi;
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
    int i, j, k, itask, send_task;
    long long ii, kk;
    long *counter, *send_size, *rec_size;
    double **send_weight_buffer, **rec_weight_buffer;
    float **send_pos_buffer, **rec_pos_buffer;
    MPI_Status mpi_status;
    
    /* 分配计数器 */
    counter = mymalloc("counter", NTask * sizeof(long));
    send_size = mymalloc("send_size", NTask * sizeof(long));
    rec_size = mymalloc("rec_size", NTask * sizeof(long));
    
    for (ii = 0; ii < NTask; ii++) {
        rec_size[ii] = 0;
        send_size[ii] = 0;
        counter[ii] = 0;
    }
    
    /* 计算每个粒子影响的所有slab */
    for (ii = FirstPart; ii < FirstPart + NumPart; ii++) {
        double weight = 1.0;
        if (value == 1) weight = P[ii].vel[0];
        if (value == 2) weight = P[ii].vel[1];
        if (value == 3) weight = P[ii].vel[2];
        
        /* 坐标缩放 */
        double xx = P[ii].pos[0] * scale_factor_mra;
        double yy = P[ii].pos[1] * scale_factor_mra;
        double zz = P[ii].pos[2] * scale_factor_mra;
        
        int xxc = (int)floor(xx);
        int yyc = (int)floor(yy);
        int zzc = (int)floor(zz);
        
        /* 确定这个粒子会影响哪些slab */
        for (int di = 0; di < phi_support; di++) {
            int grid_x = (xxc - di + GridLen) % GridLen;
            int slab_x = grid_x / (GridLen / NTask);
            
            /* 确保slab在有效范围内 */
            if (slab_x >= 0 && slab_x < NTask) {
                send_size[slab_x]++;
            }
        }
    }
    
    /* 分配发送缓冲区 */
    send_weight_buffer = (double **)mymalloc("send_weight_buffer", NTask * sizeof(double *));
    send_pos_buffer = (float **)mymalloc("send_pos_buffer", NTask * sizeof(float *));
    
    for (ii = 0; ii < NTask; ii++) {
        send_weight_buffer[ii] = (double *)mymalloc("send_weight_buffer[ii]", 
                                                    MAX(send_size[ii], 1) * sizeof(double));
        send_pos_buffer[ii] = (float *)mymalloc("send_pos_buffer[ii]", 
                                                3 * MAX(send_size[ii], 1) * sizeof(float));
    }
    
    /* 填充缓冲区 */
    for (ii = FirstPart; ii < FirstPart + NumPart; ii++) {
        double weight = 1.0;
        if (value == 1) weight = P[ii].vel[0];
        if (value == 2) weight = P[ii].vel[1];
        if (value == 3) weight = P[ii].vel[2];
        
        double xx = P[ii].pos[0] * scale_factor_mra;
        double yy = P[ii].pos[1] * scale_factor_mra;
        double zz = P[ii].pos[2] * scale_factor_mra;
        
        int xxc = (int)floor(xx);
        int yyc = (int)floor(yy);
        int zzc = (int)floor(zz);
        
        int xxf = (int)(phi_samprate * (xx - xxc));
        int yyf = (int)(phi_samprate * (yy - yyc));
        int zzf = (int)(phi_samprate * (zz - zzc));
        
        /* 为每个受影响的slab计算贡献 */
        for (int di = 0; di < phi_support; di++) {
            for (int dj = 0; dj < phi_support; dj++) {
                for (int dk = 0; dk < phi_support; dk++) {
                    int grid_x = (xxc - di + GridLen) % GridLen;
                    int grid_y = (yyc - dj + GridLen) % GridLen;
                    int grid_z = (zzc - dk + GridLen) % GridLen;
                    
                    int slab_x = grid_x / (GridLen / NTask);
                    
                    if (slab_x >= 0 && slab_x < NTask) {
                        /* 计算权重 */
                        double phi_val = phi_data[xxf + di * phi_samprate] *
                                         phi_data[yyf + dj * phi_samprate] *
                                         phi_data[zzf + dk * phi_samprate];
                        
                        double contrib = weight * phi_val;
                        
                        /* 存储到缓冲区 */
                        int idx = counter[slab_x];
                        send_weight_buffer[slab_x][idx] = contrib;
                        send_pos_buffer[slab_x][3 * idx] = grid_x;
                        send_pos_buffer[slab_x][3 * idx + 1] = grid_y;
                        send_pos_buffer[slab_x][3 * idx + 2] = grid_z;
                        counter[slab_x]++;
                    }
                }
            }
        }
    }
    
    /* 通信：发送和接收数据 */
    rec_weight_buffer = (double **)mymalloc("rec_weight_buffer", NTask * sizeof(double *));
    rec_pos_buffer = (float **)mymalloc("rec_pos_buffer", NTask * sizeof(float *));
    
    for (ii = 0; ii < NTask; ii++) {
        send_task = ii;
        
        if (send_task == ThisTask) {
            /* 发送数据到其他进程 */
            for (kk = 0; kk < NTask; kk++) {
                if (kk != ThisTask && send_size[kk] > 0) {
                    MPI_Send(&send_size[kk], 1, MPI_LONG, kk, TAG_N, MPI_COMM_WORLD);
                    MPI_Send(send_weight_buffer[kk], send_size[kk], MPI_DOUBLE, 
                            kk, TAG_MDATA, MPI_COMM_WORLD);
                    MPI_Send(send_pos_buffer[kk], 3 * send_size[kk], MPI_FLOAT,
                            kk, TAG_POSDATA, MPI_COMM_WORLD);
                }
            }
        } else {
            /* 接收其他进程的数据 */
            MPI_Recv(&rec_size[send_task], 1, MPI_LONG, send_task, TAG_N, 
                    MPI_COMM_WORLD, &mpi_status);
            
            if (rec_size[send_task] > 0) {
                rec_weight_buffer[send_task] = (double *)mymalloc("rec_weight_buffer[ii]",
                    rec_size[send_task] * sizeof(double));
                rec_pos_buffer[send_task] = (float *)mymalloc("rec_pos_buffer[ii]",
                    3 * rec_size[send_task] * sizeof(float));
                
                MPI_Recv(rec_weight_buffer[send_task], rec_size[send_task], MPI_DOUBLE,
                        send_task, TAG_MDATA, MPI_COMM_WORLD, &mpi_status);
                MPI_Recv(rec_pos_buffer[send_task], 3 * rec_size[send_task], MPI_FLOAT,
                        send_task, TAG_POSDATA, MPI_COMM_WORLD, &mpi_status);
                
                /* 处理接收到的数据 */
                for (kk = 0; kk < rec_size[send_task]; kk++) {
                    int gx = (int)rec_pos_buffer[send_task][3 * kk];
                    int gy = (int)rec_pos_buffer[send_task][3 * kk + 1];
                    int gz = (int)rec_pos_buffer[send_task][3 * kk + 2];
                    
                    /* 转换到本地索引 */
                    int local_x = gx - first_slab_of_task[ThisTask] * (GridLen / NTask);
                    
                    if (local_x >= 0 && local_x < slabs_per_task[ThisTask]) {
                        long grid_idx = gz + GridLen * (gy + GridLen * local_x);
                        if (grid_idx >= 0 && grid_idx < maxfftsize) {
                            grid[grid_idx] += rec_weight_buffer[send_task][kk];
                        }
                    }
                }
                
                myfree(rec_weight_buffer[send_task]);
                myfree(rec_pos_buffer[send_task]);
            }
        }
    }
    
    /* 处理本地数据 */
    for (kk = 0; kk < send_size[ThisTask]; kk++) {
        int gx = (int)send_pos_buffer[ThisTask][3 * kk];
        int gy = (int)send_pos_buffer[ThisTask][3 * kk + 1];
        int gz = (int)send_pos_buffer[ThisTask][3 * kk + 2];
        
        int local_x = gx - first_slab_of_task[ThisTask] * (GridLen / NTask);
        
        if (local_x >= 0 && local_x < slabs_per_task[ThisTask]) {
            long grid_idx = gz + GridLen * (gy + GridLen * local_x);
            if (grid_idx >= 0 && grid_idx < maxfftsize) {
                grid[grid_idx] += send_weight_buffer[ThisTask][kk];
            }
        }
    }
    
    /* 清理内存 */
    for (ii = 0; ii < NTask; ii++) {
        if (send_size[ii] > 0) {
            myfree(send_weight_buffer[ii]);
            myfree(send_pos_buffer[ii]);
        }
    }
    
    myfree(send_weight_buffer);
    myfree(send_pos_buffer);
    myfree(rec_weight_buffer);
    myfree(rec_pos_buffer);
    myfree(rec_size);
    myfree(send_size);
    myfree(counter);
}

/* 统一的密度场计算接口 */
void density_unified(Particle *P, fftw_real *grid, long long NumPart, int mode, int value)
{
#ifdef USE_MRA_METHOD
    /* 使用MRA方法 */
    if (phi_data == NULL) {
        init_mra_parameters();
    }
    density_mra(P, grid, NumPart, mode, value);
#else
    /* 使用传统的CIC方法 */
    density(P, grid, NumPart, mode, value);
#endif
}

/* 优化的批次处理 - 自适应MRA */
void density_adaptive(Particle *P, fftw_real *grid, long long NumPart, int mode, int value)
{
#ifdef USE_MRA_METHOD
    /* MRA需要更小的批次 */
    int adaptive_cycles = CYCLES;
    
    /* 根据phi_support调整批次 */
    if (phi_support > 8) {
        adaptive_cycles = CYCLES * 4;  /* 大支撑区域 */
    } else if (phi_support > 4) {
        adaptive_cycles = CYCLES * 2;  /* 中等支撑区域 */
    }
    
    if (ThisTask == 0) {
        printf("Adaptive cycles for MRA: %d (support=%d)\n", 
               adaptive_cycles, phi_support);
    }
    
    for (int rep = 0; rep < adaptive_cycles; rep++) {
        long long npart, nstart;
        
        if (rep == (adaptive_cycles - 1))
            npart = NumPart - (adaptive_cycles - 1) * (NumPart / adaptive_cycles);
        else
            npart = (NumPart / adaptive_cycles);
        
        nstart = rep * (NumPart / adaptive_cycles);
        
        mapping_mra(P, grid, nstart, npart, mode, value);
    }
#else
    /* 传统CIC方法使用原有批次 */
    density(P, grid, NumPart, mode, value);
#endif
}

/* 清理MRA资源 */
void cleanup_mra()
{
    if (phi_data != NULL) {
        myfree(phi_data);
        phi_data = NULL;
    }
    phi_support = 0;
}
