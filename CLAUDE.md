# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What Is SFC Field

宇宙学大尺度结构的场级别分析工具包，包含三个子项目：

- **L-VelocityField-SFC** — C/MPI 程序，生成密度场、速度场、计算功率谱。基于 FFTW 2.x + GSL
- **mracs_cpp** — C++ 程序包，pair-counting、相关函数、功率谱计算。基于 Intel MKL + FFTW3，使用 Daubechies 小波多分辨率分析
- **HM_BAO_PS** — C/MPI 程序，计算相关函数、功率谱、bipectrum，用于星系样本分析。基于 FFTW 2.x + GSL + HDF5

## Build Commands

### L-VelocityField-SFC
```bash
cd L-VelocityField-SFC
make                    # SYSTYPE=COSM, mpicc, 需要 $GSL_HOME 和 $FFTW_HOME
```
产物：`L-VelHalo-MXXL`。依赖 FFTW 2.x MPI（默认单精度）、GSL。

编译选项（Makefile 中 `OPTIONS` 行）：

| 选项 | 作用 |
|------|------|
| `GEN_HALOFIELD` | 从 halo catalogue 生成场（否则读取预计算 mesh slab） |
| `DMP` | 读取 DM 粒子 snapshot 替代 halo catalogue |
| `SMOOTH_HALOFIELD` | 在傅里叶空间高斯平滑密度场 |
| `LINEAR_VELFIELD` | 从密度场用线性理论重建速度场：**v**_k = -i(Hf) **k**/|k|² δ_k |
| `SMOOTH_HALOVELFIELD` | 平滑速度场各分量（与 `SMOOTH_HALOFIELD` 联用） |
| `CIC` | Cloud-in-Cell 质量分配（触发 k 空间反卷积） |
| `SUBHALO` | 包含 subhalo 数据 |
| `LOGDENS` | 输出 log-density 统计 |

运行：
```
mpirun -np <N> L-VelHalo-MXXL <parameterfile> <snapshot_number>
```

### HM_BAO_PS
```bash
cd HM_BAO_PS
make -f Makefile.Halo    # MACHINE=ShuGuang, mpicc, 需要 GSL, FFTW2, HDF5
```
产物：`Power.exe`。默认开启双精度 FFTW (`-DDOUBLEPRECISION_FFTW`)。

关键编译选项（Makefile.Halo 中 `OPT` 行）：

| 选项 | 作用 |
|------|------|
| `DM` | 暗物质粒子模式 |
| `CIC` | Cloud-in-Cell 分配 |
| `BISPECTRUM` | 双谱计算（目标复用的模块） |
| `SUBSAMPLE` | 子采样模式 |
| `DENSITY_DUMP/LOAD` | 密度场保存/加载 |
| `Z_SPACE` / `PHOTO_SPACE` | 红移空间 / 测光空间 |
| `CORREL` | 相关函数计算 |
| `DIVERGENCE` | 速度散度场 |
| `HOD` | Halo Occupation Distribution |
| `PRELOAD` | 预加载模式 |

### mracs_cpp
```bash
cd mracs_cpp/mracs
make                    # SYSTEM=Linux (g++), 或改为 SYSTEM=MacOS (clang++)
```
产物：`mracs`, `counting`, `2pc`, `sigma_r`, `sosta`, `orcorr`, `dipole`。
需要 Intel MKL 和 FFTW3。运行时读取 `param.txt`，无需重编译。

## Architecture

### L-VelocityField-SFC 数据流

```
main.c → read_parameter_file() → init_grid() (FFTW plan + MPI slab 分解)
  ├── [若 GEN_HALOFIELD] load_catalogue() → density() CIC赋值 → save_grids()
  └── [若已预计算 mesh] get_slabs_file() → load_grid_field()
       ├── fftwnd_mpi (forward FFT)
       ├── [若 SMOOTH_HALOFIELD] smooth() → 高斯核 exp(-k²σ²/2) → 逆FFT → save_grids()
       └── [若 LINEAR_VELFIELD] velocity_field(dim) → k空间速度重建 → 逆FFT → save_grids()
```

核心文件职责：

| 文件 | 职责 |
|------|------|
| `main.c` | 入口，CLI 解析，分支逻辑 |
| `read_parameters.c` | 参数文件解析（Smoothing, BoxSize, Time, Omega0, OutputDir 等） |
| `velocity_field.c` | `init_grid()` MPI 初始化 + `velocity_field(dim)` 线性理论速度重建 + 宇宙学函数（growth, F_Omega, hubble_function，用 GSL 积分） |
| `smooth.c` | 傅里叶空间高斯平滑 |
| `mxxl_io.c` | MXXL 模拟数据读写（slab mesh, catalogue, snapshot），`save_grids()` 输出场 |
| `field_assign.c` | CIC 粒子→网格赋值 |
| `allvars.h` | 全局变量、FFTW 类型选择、Peano-Hilbert 键宏、自定义内存宏 |
| `mymalloc.c` | 自定义内存追踪分配器（含 high-water-mark） |
| `io.c` | 日志宏（MSG/MSG0/ERR/DBG/DBG0）+ 安全 I/O 封装 |

`Grid` 缓冲区在 FFTW 就地变换中复用——前向 FFT 后复数数据存储在同一 `fftw_real` 数组中。

### HM_BAO_PS 数据流

```
main_power.c → read_parameter_file() → 密度场赋值(CIC) → FFT
  ├── power() → 功率谱计算（P(k) 分 bin）
  ├── correl() → 两点相关函数
  ├── divergence() → 速度散度场 theta = ∇·v
  └── [若 BISPECTRUM] bispectrum 模块
```

核心文件：

| 文件 | 职责 |
|------|------|
| `main_power.c` | 入口，counts-in-cells 矩、交叉相关 |
| `power.c` / `power.h` | 功率谱计算（含反卷积因子） |
| `correl.c` / `correl.h` | 两点相关函数 |
| `density.c` / `density.h` | CIC 密度场赋值 |
| `divergence.c` | 速度散度场计算 |
| `gadget.c` / `gadget.h` | Gadget snapshot 格式读写 |
| `brute_correl.c` | 暴力 pairwise 相关函数（验证用） |
| `parallel_sort.c` | 大规模并行排序 |
| `io.c` / `io.h` | I/O 和 catalogue 加载 |

### mracs_cpp 数据流

```
param.txt → read_parameter() → read_in_*() 加载粒子
  → sfc() 空间填充曲线排序到 MRA 网格
  → FFT 3D 卷积（MKL DFTI）与窗函数核
  → 结果解释（取决于分析工具：counting / 2pc / sosta 等）
```

核心文件：

| 文件 | 职责 |
|------|------|
| `src/csmain.cpp/.h` | 核心库：SFC 排序、3D 卷积、功率谱估计、小波函数生成。定义 Particle/Galaxy/Index/Offset 结构体 |
| `src/readin.cpp` | 读取 Millennium Run 星系目录 / IllustrisTNG 粒子数据 |
| `src/kernel.cpp` | 窗函数：shell, sphere, Gaussian, dual-ring |
| `src/counting.cpp` | 暴力 sphere 内 pair-counting（MRA 验证基准） |
| `src/sosta.cpp` | 小波密度估计 |
| `src/mracs.h` | 全局变量：Resolution J, BaseType, phiGenus, KernelFunc, Radius, SimBoxL, Threads |
| `param.txt` | 运行时参数（无需重编译） |
| `wavelets/` | 预计算 Daubechies 小波 phi .bin 文件 |

## 关键依赖差异

三个子项目都依赖 FFT 库但版本不同：

| 子项目 | FFT 库 | 语言 | 并行 |
|--------|--------|------|------|
| L-VelocityField-SFC | FFTW 2.x MPI | C | MPI |
| HM_BAO_PS | FFTW 2.x MPI (双精度) | C | MPI |
| mracs_cpp | FFTW 3.x + Intel MKL | C++ | OpenMP |

L-VelocityField-SFC 和 HM_BAO_PS 共享类似的代码结构（mymalloc、CIC、FFTW2 MPI slab 分解），但 HM_BAO_PS 功能更丰富（HOD、bipectrum、多种空间）。

## Reference

- mracs_cpp 背景文献：
  - <https://iopscience.iop.org/article/10.1086/511024>
  - <https://academic.oup.com/mnras/article/535/4/3500/7899967>
  - <https://academic.oup.com/mnras/article/546/1/staf2275/8405292>
- HM_BAO_PS 背景文献：
  - <https://academic.oup.com/mnras/article/442/3/2131/1045856>

## Planning & Tracking

- Plan in `docs/plan/phase_NN.md`（如 `phase_01.md`）
- Each task has ID: `PNN_NNN`（如 `P01_001` = phase 1, task 1）
- Track in `docs/TODO.md` with per-phase tables
- Record mistakes/instructions in `docs/LESSONS.md`

## Workflow

- 将 mracs_cpp 的算法移植到 L-VelocityField-SFC 中
- 结合 HM_BAO_PS 中计算 bispectrum 的模块，将 mracs_cpp 算法实现出计算 bispectrum 和三点相关函数的模块集成到 L-VelocityField-SFC 中
- 规整 L-VelocityField-SFC 中的程序，形成一个规范的程序包
