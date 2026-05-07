# SFC Field Study: L-VelocityField-SFC 重构设计

2026-05-07

## 目标

将三个子项目（L-VelocityField-SFC、mracs_cpp、HM_BAO_PS）的功能整合到一个规范化的 C/MPI 程序包中。

## 范围与优先级

- **B 策略**：功能优先（移植 → 集成 → 规整）
- **规整分两阶段**：先 A（清理死代码、统一命名、header guard），后 B（标准目录布局、统一 Makefile、运行时参数文件替代编译宏）
- mracs_cpp 全部功能移植，用**纯 C** 重写
- bispectrum 两种方法都要：传统 FFT（HM_BAO_PS）和小波（mracs_cpp 方法）
- 最小区分：本 spec 聚焦第一阶段（清理 + 新模块添加，保持现有目录结构可用）

## 模块架构

```
L-VelocityField-SFC/
├── main.c                    # 入口，CLI 解析，按 mode dispatch
├── Makefile                  # 统一构建
│
├── core/                     # [清理] 基础设施
│   ├── allvars.h             # 全局类型、宏
│   ├── proto.h               # 函数声明
│   ├── io.c / io.h           # 日志、安全 I/O
│   ├── mymalloc.c / mymalloc.h  # 内存追踪
│   └── read_parameters.c     # 参数文件解析
│
├── field/                    # [清理] 场引擎
│   ├── velocity_field.c      # FFTW 初始化、速度重建、宇宙学函数
│   ├── smooth.c              # 傅里叶空间高斯平滑
│   ├── field_assign.c        # CIC 粒子→网格赋值
│   └── mxxl_io.c             # MXXL 模拟格式 I/O
│
├── mracs/                    # [新增] MRA 引擎
│   ├── particle_io.c         # Millennium/TNG 粒子读取
│   ├── sfc.c                 # 空间填充曲线排序
│   ├── wavelet.c             # Daubechies/B-spline 尺度函数
│   ├── kernel.c              # 窗函数
│   └── mra_grid.c            # MRA 网格 + 3D 卷积
│
├── statistics/               # [清理+新增] 统计量
│   ├── power.c               # 功率谱
│   ├── correl.c              # 两点相关函数
│   ├── bispectrum_fft.c      # [新增] FFT bispectrum
│   ├── bispectrum_wavelet.c  # [新增] 小波 bispectrum
│   └── three_pcf.c           # [新增] 三点相关函数
│
└── job-scripts/              # HPC 作业脚本
```

## 核心数据结构

跨模块统一的三个基础结构体：

```c
/* 粒子 */
typedef struct {
    double x, y, z;
    double weight;
} Particle;

/* 3D 网格 */
typedef struct {
    long long Ngrid;
    long long Ntotal;
    double BoxSize;
    int local_ny;
    int local_y_start;
    fftw_real *data;
} Grid;

/* 小波配置 */
typedef struct {
    int resolution;     /* J, grid=2^J */
    int base_type;      /* 0=B-Spline, 1=Daubechies */
    int genus;          /* Daubechies 阶数 */
    int samp_rate;
    double *phi;
} WaveletConfig;

/* 窗函数类型 */
typedef enum { SHELL, SPHERE, GAUSSIAN, DUAL_RING } KernelType;

/* 统计结果 bin */
typedef struct {
    long long n;
    double x, y;
    double data, data2;
} Result;
```

## 主处理管线

通过参数文件 `mode` 字段选择，每个 mode 有 init/run/clean 生命周期：

```
main.c
  ├── mode=density          已有
  ├── mode=velocity         已有
  ├── mode=power            已有
  ├── mode=correl           已有
  ├── mode=paircount        新增, mracs_cpp
  │     read_particles → sfc_order → mra_grid → convolve → count
  ├── mode=bispec_fft       新增, HM_BAO_PS
  │     CIC → FFT → bispectrum_fft(triangles) → write
  ├── mode=bispec_wavelet   新增, mracs_cpp
  │     read_particles → sfc → mra_grid → wavelet_decomp → bispec
  └── mode=three_pcf        新增
        read_particles → sfc → mra_grid → wavelet_decomp → angular_3pcf
```

## 第一阶段：清理（现有代码）

1. 删除 `.bak` 文件：`main.c.bak`, `main.c.bak1`, `mxxl_io.c.bak1`
2. 删除 `.o` 文件
3. 统一命名风格为 snake_case（camelCase → snake_case）
4. 每个 `.h` 加 `#ifndef` guard
5. 移除无用的 Makefile SYSTYPE 路径（保留 COSM 和 ShuGuang，删除 Ubuntu/OpenSuse/Genius/HLRB2/VIP 等不活跃的）
6. 清理 `main.c` 中注释掉或 `#if 0` 的死代码
7. 提取 `density_statistics()` 和 `log_density()` 从 `main.c` 移到独立的 `statistics/stats_utils.c`

## 第一阶段：新增模块（从 mracs_cpp 和 HM_BAO_PS）

按以下顺序实现，每个子任务完成后可独立编译测试：

### Step 1: particle_io.c（移植自 mracs_cpp/src/readin.cpp）
- `read_Millennium_galaxies(char *path, Particle **p, long long *n)` — 读取二进制星系目录
- `read_TNG_particles(char *path, Particle **p, long long *n)` — 读取 IllustrisTNG 粒子
- 返回动态分配的 `Particle` 数组，调用者负责释放

### Step 2: sfc.c（移植自 mracs_cpp/src/csmain.cpp 的 SFC 部分）
- `sfc_order(Particle *p, long long n, long long grid_len)` — 空间填充曲线排序
- 依赖：粒子数据、MRA 网格边长

### Step 3: wavelet.c（移植自 mracs_cpp/src/csmain.cpp 的小波部分）
- `daubechies_phi(int genus, int samp_rate, double **phi, int *n)` — 生成 Daubechies 尺度函数
- `b_spline(int n, int samp_rate, double **phi, int *len)` — 生成 B-spline
- 依赖：无外部依赖（纯数值计算）

### Step 4: kernel.c（移植自 mracs_cpp/src/kernel.cpp）
- `window_shell(double R, double kx, double ky, double kz)` → 傅里叶空间壳函数
- `window_sphere(double R, double kx, double ky, double kz)` → 球顶帽
- `window_gaussian(double R, double kx, double ky, double kz)` → 高斯窗

### Step 5: mra_grid.c（移植自 mracs_cpp/src/csmain.cpp 的卷积部分）
- `mra_convolve(Grid *input, Grid *output, WaveletConfig *wav, KernelType kernel, double radius)` — MRA 网格上的 3D 卷积
- 使用 FFTW2（替代 MKL DFTI）

### Step 6: bispectrum_fft.c（提取自 HM_BAO_PS/power.c）
- `bispectrum_fft(Grid *grid, double k1, double k2, int nbins, Result *out)` — FFT 方法 bispectrum
- 从 `power.c` 的 `#ifdef BISPECTRUM` 区域提取，改为接受 `Grid*` 参数而非全局变量
- 依赖：`Grid` 结构体（已经做完 CIC + FFT）

### Step 7: bispectrum_wavelet.c + three_pcf.c（基于 MRA 的新实现）
- `bispectrum_wavelet(Grid *grid, WaveletConfig *wav, KernelType kernel, double radius, int nbins, Result *out)`
- `three_pcf(Particle *p, long long n, WaveletConfig *wav, KernelType kernel, double radius, Result *out)`

## 构建

- 保持单一 Makefile，所有 `.c` 编译为 `.o`，链接为可执行文件
- 编译选项通过 Makefile 顶部 `OPTIONS` 变量控制
- 依赖：MPI、FFTW 2.x、GSL、HDF5（可选）

## 后续阶段（暂不执行）

- 重组为 `src/` `include/` `lib/` 标准布局
- 所有配置走参数文件而非编译宏
- 形成 libvelocityfield.a 可复用库，公共 API：`vf_init()` / `vf_run()` / `vf_clean()`
