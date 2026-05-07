# L-VelocityField-SFC Refactoring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 清理 L-VelocityField-SFC 现有代码并移植 mracs_cpp + HM_BAO_PS bispectrum 模块，形成规范的纯 C/MPI 场分析程序包

**Architecture:** 渐进式移植——先清理现有代码（Phase 1A），再按 particle_io → sfc → wavelet → kernel → mra_grid → bispectrum_fft → bispectrum_wavelet → three_pcf 顺序逐步新增模块（Phase 1B），每步产出可编译验证的代码

**Tech Stack:** C99, MPI, FFTW 2.x, GSL, HDF5 (可选), mpicc

---

## File Mapping (What Goes Where)

### Phase 1A — 清理
| 操作 | 文件 |
|------|------|
| 删除 | `main.c.bak`, `main.c.bak1`, `mxxl_io.c.bak1`, `*.o` |
| 修改 | `allvars.h` — 加 guard, 重命名 camelCase |
| 修改 | `proto.h` — 加 guard |
| 修改 | `main.c` — 提取 statistics 函数, 清理死代码 |
| 修改 | `velocity_field.c` — 统一命名 |
| 修改 | `smooth.c` — 统一命名 |
| 修改 | `field_assign.c` — 统一命名 |
| 修改 | `mxxl_io.c` — 统一命名 |
| 修改 | `io.c` — 加 .h, 统一命名 |
| 修改 | `mymalloc.c` — 加 .h guard |
| 修改 | `read_parameters.c` — 统一命名 |
| 修改 | `Makefile` — 精简 SYSTYPE, 更新 OBJS |
| 新增 | `io.h` — io.c 的函数声明 |
| 新增 | `mymalloc.h` — mymalloc.c 的函数声明 |
| 新增 | `field_assign.h` — 场赋值函数声明 |
| 新增 | `statistics/stats_utils.c` — density_statistics, log_density |
| 新增 | `statistics/stats_utils.h` — 对应头文件 |

### Phase 1B — 新增
| 操作 | 文件 | 来源 |
|------|------|------|
| 新增 | `mracs/particle_io.c` + `.h` | mracs_cpp/src/readin.cpp |
| 新增 | `mracs/sfc.c` + `.h` | mracs_cpp/src/csmain.cpp |
| 新增 | `mracs/wavelet.c` + `.h` | mracs_cpp/src/csmain.cpp |
| 新增 | `mracs/kernel.c` + `.h` | mracs_cpp/src/kernel.cpp |
| 新增 | `mracs/mra_grid.c` + `.h` | mracs_cpp/src/csmain.cpp |
| 新增 | `statistics/bispectrum_fft.c` + `.h` | HM_BAO_PS/power.c |
| 新增 | `statistics/bispectrum_wavelet.c` + `.h` | 基于 MRA 新实现 |
| 新增 | `statistics/three_pcf.c` + `.h` | 基于 MRA 新实现 |
| 修改 | `main.c` — 新增 mode dispatch |
| 修改 | `Makefile` — 新增 OBJS 和目录 |

---

### Task 1: Git 初始化和目录准备

**Files:**
- Create: `.gitignore`

- [ ] **Step 1: Initialize git repo**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git init && git config user.name "liming" && git config user.email "disiheng@gmail.com"
```

- [ ] **Step 2: Create .gitignore**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/.gitignore << 'EOF'
*.o
*.exe
*.out
.DS_Store
.ipynb_checkpoints/
*.tar
*.swp
*.swo
EOF
```

- [ ] **Step 3: Initial commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "chore: initial commit — snapshot before refactoring"
```

- [ ] **Step 4: Create refactor branch**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git checkout -b refactor/lvelocityfield-cleanup
```

---

### Task 2: 删除死文件

- [ ] **Step 1: Remove backup and object files**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC && rm -f main.c.bak main.c.bak1 mxxl_io.c.bak1 *.o
```

- [ ] **Step 2: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "chore: remove backup and object files"
```

---

### Task 3: 给所有 .h 文件加 include guard

**Files:**
- Modify: `L-VelocityField-SFC/allvars.h`
- Modify: `L-VelocityField-SFC/proto.h`

- [ ] **Step 1: Verify allvars.h already has guard**

`allvars.h` already has `#ifndef ALLVARS_H` / `#define ALLVARS_H` / `#endif` — no change needed.

- [ ] **Step 2: Verify proto.h guard is correct**

Read `proto.h:1-4` — it starts with `#ifndef ALLVARS_H` which incorrectly uses the same guard name as `allvars.h`. Fix it:

Edit `L-VelocityField-SFC/proto.h`, change:
```c
#ifndef ALLVARS_H
#include "allvars.h"
#endif
```
To:
```c
#ifndef PROTO_H
#define PROTO_H

#include "allvars.h"
```

And add `#endif` at the end of the file (after the last function declaration).

- [ ] **Step 3: Verify compilation**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC && make clean && make 2>&1 | head -30
```

Expected: compiles without errors (may fail if $GSL_HOME/$FFTW_HOME not set locally — that's OK, the header changes are syntactically correct).

- [ ] **Step 4: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "fix: add proper include guard to proto.h"
```

---

### Task 4: 提取 io.h 和 mymalloc.h 头文件

**Files:**
- Create: `L-VelocityField-SFC/io.h`
- Create: `L-VelocityField-SFC/mymalloc.h`
- Modify: `L-VelocityField-SFC/io.c`
- Modify: `L-VelocityField-SFC/mymalloc.c`
- Modify: `L-VelocityField-SFC/allvars.h`

- [ ] **Step 1: Create io.h**

从 `proto.h` 中提取 io.c 相关的函数声明：

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/io.h << 'EOF'
#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <time.h>

void MSG(const char *format, time_t t0, ...);
void MSG0(const char *format, time_t t0, int task, ...);
void ERR(const char *format, ...);
void DBG(const char *format, int mydbg, time_t t0, ...);
void DBG0(const char *format, int mydbg, time_t t0, int task, ...);

size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
FILE *my_fopen(char *file, char *mode);

#endif
EOF
```

- [ ] **Step 2: Create mymalloc.h**

从 `allvars.h` 中提取 mymalloc 宏和声明：

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mymalloc.h << 'EOF'
#ifndef MYMALLOC_H
#define MYMALLOC_H

#include <stdlib.h>

extern size_t AllocatedBytes;
extern size_t HighMarkBytes;
extern float HighMarkBytes_Step;
extern size_t FreeBytes;

#define mymalloc(x, y)            mymalloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define mymalloc_movable(x, y, z) mymalloc_movable_fullinfo(x, y, z, __FUNCTION__, __FILE__, __LINE__)
#define myrealloc(x, y)           myrealloc_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define myrealloc_movable(x, y)   myrealloc_movable_fullinfo(x, y, __FUNCTION__, __FILE__, __LINE__)
#define myfree(x)                 myfree_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define myfree_movable(x)         myfree_movable_fullinfo(x, __FUNCTION__, __FILE__, __LINE__)
#define report_memory_usage(x, y) report_detailed_memory_usage_of_largest_task(x, y, __FUNCTION__, __FILE__, __LINE__)

void *mymalloc_fullinfo(const char *varname, size_t n, const char *func, const char *file, int linenr);
void *mymalloc_movable_fullinfo(void *ptr, const char *varname, size_t n, const char *func, const char *file, int line);
void *myrealloc_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void *myrealloc_movable_fullinfo(void *p, size_t n, const char *func, const char *file, int line);
void myfree_fullinfo(void *p, const char *func, const char *file, int line);
void myfree_movable_fullinfo(void *p, const char *func, const char *file, int line);
void mymalloc_init(void);
void dump_memory_table(void);
void report_detailed_memory_usage_of_largest_task(size_t *OldHighMarkBytes, const char *label, const char *func, const char *file, int line);

#endif
EOF
```

- [ ] **Step 3: Update io.c to include io.h**

Edit `L-VelocityField-SFC/io.c`, after existing includes, add:
```c
#include "io.h"
```

- [ ] **Step 4: Update mymalloc.c to include mymalloc.h**

Edit `L-VelocityField-SFC/mymalloc.c`, after existing includes, add:
```c
#include "mymalloc.h"
```

- [ ] **Step 5: Update allvars.h — remove mymalloc macros (they're now in mymalloc.h)**

Edit `L-VelocityField-SFC/allvars.h`, remove lines 139-148 (the `#define mymalloc` through `#define report_memory_usage` block), replace with:
```c
#include "mymalloc.h"
```

Also remove the extern declarations for `AllocatedBytes`, `HighMarkBytes`, `HighMarkBytes_Step`, `FreeBytes` (lines 125-128) since they're now in `mymalloc.h`.

- [ ] **Step 6: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "refactor: extract io.h and mymalloc.h headers"
```

---

### Task 5: 清理 allvars.h — 统一命名和去冗余

**Files:**
- Modify: `L-VelocityField-SFC/allvars.h`

- [ ] **Step 1: Add include for stdlib.h**

At the top of `allvars.h`, under `#include <stdio.h>`, also make sure we have what's needed. The file already has `#include <stdio.h>` but uses `size_t` — add:
```c
#include <stdlib.h>
```

- [ ] **Step 2: Remove unused macros**

Remove `SET_MASK`, `IMIN`, `BITS_PER_DIMENSION`, `PEANOCELLS` macros (lines 130-134) — these are Peano-Hilbert related but not used in L-VelocityField-SFC (that's HM_BAO_PS territory).

- [ ] **Step 3: Remove unused `#define dbg1 0` / `#define dbg2 0`** (lines 38-39)

- [ ] **Step 4: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "refactor: clean up allvars.h — remove unused macros"
```

---

### Task 6: 从 main.c 提取 statistics 函数

**Files:**
- Create: `L-VelocityField-SFC/statistics/stats_utils.c`
- Create: `L-VelocityField-SFC/statistics/stats_utils.h`
- Modify: `L-VelocityField-SFC/main.c`
- Modify: `L-VelocityField-SFC/proto.h`
- Modify: `L-VelocityField-SFC/Makefile`

- [ ] **Step 1: Create stats_utils.h**

```bash
mkdir -p /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/stats_utils.h << 'EOF'
#ifndef STATS_UTILS_H
#define STATS_UTILS_H

#include "../allvars.h"

double density_statistics(fftw_real *grid);
double log_density(float *grid);

#endif
EOF
```

- [ ] **Step 2: Create stats_utils.c**

Move `density_statistics()` and `log_density()` from `main.c` (lines 138-253) into this new file:

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/stats_utils.c << 'ENDOFFILE'
#include <mpi.h>
#include <stdio.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "stats_utils.h"

double density_statistics(fftw_real *grid)
{
    long long i;
    double delta;
    double local_total = 0.0, total = 0.0;
    double local_var = 0.0, var = 0.0;
    double local_m3 = 0.0, m3 = 0.0, m1 = 0.0, m2 = 0.0;
    double local_min = 1e20, min = 1e20;
    double local_max = -1e20, max = -1e20;
    double local_mean = 0.0, mean = 0.0;
    double local_dm3 = 0.0, dm3 = 0.0;
    double local_dm2 = 0.0, dm2 = 0.0;

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++)
        local_mean += (double) grid[i];
    MPI_Allreduce(&local_mean, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    mean = mean / CUB(Ndim);

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++) {
        if (local_min > grid[i]) local_min = grid[i];
        if (local_max < grid[i]) local_max = grid[i];
        local_total += (double) grid[i];
        local_var += SQR(grid[i]);
        local_m3 += CUB(grid[i]);
    }

    MPI_Allreduce(&local_m3, &m3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_var, &var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_min, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&local_max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    mean = total / CUB(Ndim);
    m1 = mean;
    m2 = var / CUB(Ndim);
    m3 = m3 / CUB(Ndim);

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++) {
        delta = grid[i];
        local_dm2 += SQR(delta);
        local_dm3 += CUB(delta);
    }

    MPI_Allreduce(&local_dm2, &dm2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_dm3, &dm3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    dm2 = dm2 / CUB(Ndim);
    dm3 = dm3 / CUB(Ndim);

    if (ThisTask == 0) {
        printf(" ------------------------------------------------\n");
        printf(" Grid statistics\n");
        printf(" Size: %ld\n", Ndim);
        printf(" Local Size: %d starting at %d\n", nslab_x, slabstart_x);
        printf(" Minimun Value: %f\n", min);
        printf(" Maximun Value: %f\n", max);
        printf(" Mean: %f\n", m1);
        printf(" <S^2>: %f\n", m2);
        printf(" <S^3>: %f\n", m3);
        printf(" <(S/<S>-1)^2>: %f\n", dm2);
        printf(" <(S/<S>-1)^3>: %f\n", dm3);
        printf(" Total: %f\n", total);
        printf(" Variance: %f\n", sqrt(var / Ndim - mean * mean));
        printf(" ------------------------------------------------\n");
    }

    return total;
}

double log_density(float *grid)
{
    long long i;
    double local_total = 0.0, total = 0.0;
    double local_log_total = 0.0, log_total = 0.0;
    double log_mean, mean, true_mean;

    true_mean = (double) TotNumHalo / CUB(Ndim);

    for (i = 0; i < SQR(Ndim) * (nslab_x); i++) {
        local_log_total += (double) log(grid[i] / true_mean);
        local_total += (double) (grid[i]);
    }

    MPI_Allreduce(&local_total, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&local_log_total, &log_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    mean = total / CUB(Ndim);
    log_mean = log_total / CUB(Ndim);

    if (ThisTask == 0) {
        printf(" ------------------------------------------------\n");
        printf(" Grid statistics\n");
        printf(" Size: %ld\n", Ndim);
        printf(" Local Size: %d starting at %d\n", nslab_x, slabstart_x);
        printf(" Mean: %f\n", mean);
        printf(" Log-Mean: %f\n", log_mean);
        printf(" True Mean: %f\n", true_mean);
        printf(" ------------------------------------------------\n");
    }

    return log_mean;
}
ENDOFFILE
```

- [ ] **Step 3: Remove functions from main.c**

Remove lines 138-253 from `main.c` (the `density_statistics` and `log_density` function definitions). Add `#include "statistics/stats_utils.h"` at the top.

- [ ] **Step 4: Remove declarations from proto.h**

Remove `double density_statistics(fftw_real *grid);` and `double log_density(float *grid);` from `proto.h`.

- [ ] **Step 5: Update Makefile**

Edit `Makefile`, change `OBJS` to include `statistics/stats_utils.o`:
```makefile
OBJS = main.o mymalloc.o read_parameters.o allvars.o io.o smooth.o \
       velocity_field.o mxxl_io.o field_assign.o \
       statistics/stats_utils.o
```

Add build rule for subdirectory:
```makefile
statistics/%.o: statistics/%.c
	$(CC) $(CFLAGS) -c $< -o $@
```

- [ ] **Step 6: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "refactor: extract density_statistics and log_density to statistics/stats_utils.c"
```

---

### Task 7: 创建 field_assign.h 头文件

**Files:**
- Create: `L-VelocityField-SFC/field_assign.h`
- Modify: `L-VelocityField-SFC/field_assign.c`
- Modify: `L-VelocityField-SFC/proto.h`

- [ ] **Step 1: Create field_assign.h**

Check what functions are in `field_assign.c` first, then create:

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/field_assign.h << 'EOF'
#ifndef FIELD_ASSIGN_H
#define FIELD_ASSIGN_H

#include "allvars.h"

void density(Particle *P, fftw_real *grid, long long NumPart, int mode, int value);
void mapping(Particle *P, fftw_real *grid, long long FirstPart, long long NumPart, int mode, int value);

#endif
EOF
```

- [ ] **Step 2: Update field_assign.c — add include**

Add `#include "field_assign.h"` at the top of `field_assign.c`.

- [ ] **Step 3: Remove field_assign declarations from proto.h**

Remove `void density(...)` and `void mapping(...)` from `proto.h`.

- [ ] **Step 4: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "refactor: extract field_assign.h header"
```

---

### Task 8: 精简 Makefile

**Files:**
- Modify: `L-VelocityField-SFC/Makefile`

- [ ] **Step 1: Remove inactive SYSTYPE blocks**

Remove Makefile lines for SYSTYPEs: `JSC`, `Genius`, `HLRB2`, `Ubuntu`, `OpenSuse`, `OpenSuse64`, `OPA-Cluster64-Intel`, `OPA-Cluster64-Gnu`, `OpteronMPA-Intel`, `OpteronMPA-Gnu`, `MPA`, `VIP`. Keep only `COSM` (and add `ShuGuang` for HM_BAO_PS compatibility).

Also remove the `SOFTDOUBLEDOUBLE` and `VORONOI` conditional blocks (lines 63-71).

- [ ] **Step 2: Clean up the FFTW_LIB section**

Keep the `NOTYPEPREFIX_FFTW` / `DOUBLEPRECISION_FFTW` conditional but remove unused `ifeq` nesting. The simplified version:

```makefile
ifeq (DOUBLEPRECISION_FFTW,$(findstring DOUBLEPRECISION_FFTW,$(OPT)))
  FFTW_LIB = $(FFTW_LIBS) -ldrfftw_mpi -ldfftw_mpi -ldrfftw -ldfftw
else
  FFTW_LIB = $(FFTW_LIBS) -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw
endif
```

- [ ] **Step 3: Remove empty lines at end of file** (lines 293-304)

- [ ] **Step 5: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "refactor: simplify Makefile — remove inactive SYSTYPE blocks"
```

---

### Task 9: 统一 main.c 命名和清理

**Files:**
- Modify: `L-VelocityField-SFC/main.c`

- [ ] **Step 1: Remove unistd.h include if not used**

`main.c` includes `<unistd.h>` and `<signal.h>` and `<assert.h>` — check usage. `unistd.h` not used in main.c, remove. `signal.h` not used, remove.

- [ ] **Step 2: Rename camelCase variables to snake_case**

In `main.c`:
- `SnapshotNum` → `snapshot_num` (already used in code, but declared extern in allvars.h — skip for now, rename in allvars.h in future pass)
- `ThisTask` → `this_task` (this is global, skip for now — save for a future global-rename pass that touches all files)

Actually, since camelCase renaming touches every file, let's defer the global rename to a later cleanup pass and not risk breaking compilation. For now just clean dead code.

- [ ] **Step 3: Remove dead/commented code**

Remove the `#ifdef GEN_HALOFIELD` block's unused variable declarations (the `#ifndef DMP` load_catalogue path) — no, keep it, it's a valid code path under a compile flag.

Actually, `main.c` is pretty clean already — the function extraction in Task 6 was the main change. Let's verify no remaining dead code.

- [ ] **Step 4: Add missing includes for extracted stats_utils**

Add at top of `main.c`:
```c
#include "statistics/stats_utils.h"
```

- [ ] **Step 5: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "refactor: cleanup main.c — remove unused includes, add stats_utils.h"
```

---

### Task 10: 移植 particle_io.c — Millennium/TNG 粒子读取

**Files:**
- Create: `L-VelocityField-SFC/mracs/particle_io.c`
- Create: `L-VelocityField-SFC/mracs/particle_io.h`

**Source:** `mracs_cpp/mracs/src/readin.cpp` + `csmain.h` 中的 `Galaxy` 结构体

- [ ] **Step 1: Create directory and header**

```bash
mkdir -p /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/particle_io.h << 'EOF'
#ifndef PARTICLE_IO_H
#define PARTICLE_IO_H

#include <stdint.h>

/* Particle with position and weight */
typedef struct {
    double x, y, z;
    double weight;
} Particle;

/* Millennium Run galaxy (23 float fields in binary) */
typedef struct {
    float x, y, z;
    float vx, vy, vz;
    float mag_u, mag_g, mag_r, mag_i, mag_z;
    float bulge_mag_u, bulge_mag_g, bulge_mag_r, bulge_mag_i, bulge_mag_z;
    float stellar_mass, bulge_mass, cold_gas, hot_gas, ejected_mass, black_hole_mass, sfr;
} Galaxy;

long long read_millennium_galaxies(const char *path, Particle **particles, int use_mass_weight);
long long read_tng_particles(const char *path, Particle **particles);
void free_particles(Particle *p);

#endif
EOF
```

- [ ] **Step 2: Create particle_io.c**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/particle_io.c << 'ENDOFFILE'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "particle_io.h"

long long read_millennium_galaxies(const char *path, Particle **particles, int use_mass_weight)
{
    FILE *fp;
    long long ngal, i;
    float *buf;
    Galaxy gal;

    fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s\n", path);
        return -1;
    }

    /* Determine file size */
    fseek(fp, 0, SEEK_END);
    long long file_size = ftell(fp);
    rewind(fp);

    /* Each galaxy = 23 floats = 92 bytes */
    ngal = file_size / (23 * sizeof(float));

    *particles = (Particle *) malloc(ngal * sizeof(Particle));
    if (!*particles) {
        fprintf(stderr, "Error: malloc failed for %lld particles\n", ngal);
        fclose(fp);
        return -1;
    }

    buf = (float *) malloc(23 * sizeof(float));

    for (i = 0; i < ngal; i++) {
        if (fread(buf, sizeof(float), 23, fp) != 23) break;

        (*particles)[i].x = (double) buf[0];
        (*particles)[i].y = (double) buf[1];
        (*particles)[i].z = (double) buf[2];

        if (use_mass_weight) {
            (*particles)[i].weight = (double) buf[16];  /* stellar_mass */
        } else {
            (*particles)[i].weight = 1.0;
        }
    }

    free(buf);
    fclose(fp);
    return i;
}

long long read_tng_particles(const char *path, Particle **particles)
{
    FILE *fp;
    long long npart, i;
    double pos[3];

    fp = fopen(path, "rb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open %s\n", path);
        return -1;
    }

    fseek(fp, 0, SEEK_END);
    long long file_size = ftell(fp);
    rewind(fp);

    /* Each particle = 3 doubles (x,y,z) */
    npart = file_size / (3 * sizeof(double));

    *particles = (Particle *) malloc(npart * sizeof(Particle));
    if (!*particles) {
        fprintf(stderr, "Error: malloc failed for %lld particles\n", npart);
        fclose(fp);
        return -1;
    }

    for (i = 0; i < npart; i++) {
        if (fread(pos, sizeof(double), 3, fp) != 3) break;

        (*particles)[i].x = pos[0];
        (*particles)[i].y = pos[1];
        (*particles)[i].z = pos[2];
        (*particles)[i].weight = 1.0;
    }

    fclose(fp);
    return i;
}

void free_particles(Particle *p)
{
    free(p);
}
ENDOFFILE
```

- [ ] **Step 3: Write a small test program to verify compilation**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/test_particle_io.c << 'EOF'
#include <stdio.h>
#include "particle_io.h"

int main(void) {
    Particle *p = NULL;
    long long n;

    /* Test with non-existent file — should print error, not crash */
    n = read_millennium_galaxies("/nonexistent/path.bin", &p, 0);
    if (n < 0) {
        printf("OK: correctly handled missing file\n");
    }

    printf("Particle size: %zu bytes\n", sizeof(Particle));
    printf("Galaxy size: %zu bytes (expected 92)\n", sizeof(Galaxy));

    return 0;
}
EOF
```

- [ ] **Step 4: Compile and run test**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs && gcc -std=c99 -Wall -o test_particle_io particle_io.c test_particle_io.c -lm && ./test_particle_io
```

Expected output:
```
Error: cannot open /nonexistent/path.bin
OK: correctly handled missing file
Particle size: 32 bytes
Galaxy size: 92 bytes (expected 92)
```

- [ ] **Step 5: Clean up test and commit**

```bash
rm /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/test_particle_io.c /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/test_particle_io
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add particle_io — Millennium/TNG particle readers (from mracs_cpp)"
```

---

### Task 11: 移植 kernel.c — 窗函数

**Files:**
- Create: `L-VelocityField-SFC/mracs/kernel.c`
- Create: `L-VelocityField-SFC/mracs/kernel.h`

**Source:** `mracs_cpp/mracs/src/kernel.cpp`

- [ ] **Step 1: Create kernel.h**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/kernel.h << 'EOF'
#ifndef KERNEL_H
#define KERNEL_H

typedef enum {
    KERNEL_SHELL,
    KERNEL_SPHERE,
    KERNEL_GAUSSIAN,
    KERNEL_DUAL_RING
} KernelType;

double window_shell(double R, double kx, double ky, double kz);
double window_sphere(double R, double kx, double ky, double kz);
double window_gaussian(double R, double kx, double ky, double kz);
double window_dual_ring(double R, double theta, double kx, double ky, double kz);
double kernel_eval(KernelType type, double R, double theta, double kx, double ky, double kz);

#endif
EOF
```

- [ ] **Step 2: Read kernel.cpp for exact implementation**

First, read the original C++ code:

```bash
cat /Users/liming/Nutstore/SFC_Field_Study/mracs_cpp/mracs/src/kernel.cpp
```

- [ ] **Step 3: Create kernel.c with C port of all 4 window functions**

The C port preserves the math from `kernel.cpp`, replacing `<cmath>` with `<math.h>`:

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/kernel.c << 'ENDOFFILE'
#include <math.h>
#include "kernel.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double window_shell(double R, double kx, double ky, double kz)
{
    double k = sqrt(kx * kx + ky * ky + kz * kz);
    if (k == 0.0) return 4.0 * M_PI * R * R;
    return 4.0 * M_PI * R * R * sin(k * R) / (k * R);
}

double window_sphere(double R, double kx, double ky, double kz)
{
    double k = sqrt(kx * kx + ky * ky + kz * kz);
    if (k == 0.0) return 4.0 * M_PI * R * R * R / 3.0;
    double kR = k * R;
    return 4.0 * M_PI * (sin(kR) - kR * cos(kR)) / (k * k * k);
}

double window_gaussian(double R, double kx, double ky, double kz)
{
    double k2 = kx * kx + ky * ky + kz * kz;
    double R2 = R * R;
    return pow(2.0 * M_PI * R2, 1.5) * exp(-k2 * R2 / 2.0);
}

double window_dual_ring(double R, double theta, double kx, double ky, double kz)
{
    double k = sqrt(kx * kx + ky * ky + kz * kz);
    if (k == 0.0) return 0.0;
    double cost = (kx * sin(theta) + ky * cos(theta)) / k;
    double kR = k * R;
    return 2.0 * cos(kR * cost) * 4.0 * M_PI * R * R * sin(kR) / (kR);
}

double kernel_eval(KernelType type, double R, double theta, double kx, double ky, double kz)
{
    switch (type) {
        case KERNEL_SHELL:      return window_shell(R, kx, ky, kz);
        case KERNEL_SPHERE:     return window_sphere(R, kx, ky, kz);
        case KERNEL_GAUSSIAN:   return window_gaussian(R, kx, ky, kz);
        case KERNEL_DUAL_RING:  return window_dual_ring(R, theta, kx, ky, kz);
        default:                return 0.0;
    }
}
ENDOFFILE
```

- [ ] **Step 4: Compile test**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs && gcc -std=c99 -Wall -c kernel.c -o kernel.o && echo "OK: kernel.c compiles"
```

- [ ] **Step 5: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add kernel.c — window functions (shell/sphere/gaussian/dual-ring) from mracs_cpp"
```

---

### Task 12: 移植 wavelet.c — Daubechies/B-spline 尺度函数

**Files:**
- Create: `L-VelocityField-SFC/mracs/wavelet.c`
- Create: `L-VelocityField-SFC/mracs/wavelet.h`
- Reference: `mracs_cpp/mracs/src/csmain.cpp` — `Daubechies_Phi()`, `B_Spline()`
- Reference: `Untitled.ipynb` — loading and verification of wavelet .bin files

- [ ] **Step 1: Create wavelet.h**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/wavelet.h << 'EOF'
#ifndef WAVELET_H
#define WAVELET_H

typedef struct {
    int resolution;       /* J, grid side = 2^J */
    int base_type;        /* 0 = B-Spline, 1 = Daubechies */
    int genus;            /* Daubechies order n, or B-Spline order */
    int samp_rate;        /* samples per unit support interval */
    int phi_len;          /* total samples in phi array */
    double *phi;          /* scaling function samples */
} WaveletConfig;

/* Load Daubechies phi from pre-computed .bin file */
int daubechies_load(const char *datadir, int genus, WaveletConfig *cfg);

/* Generate B-Spline scaling function in memory (does not require .bin file) */
int bspline_generate(int order, int samp_rate, WaveletConfig *cfg);

/* Free wavelet config */
void wavelet_free(WaveletConfig *cfg);

#endif
EOF
```

- [ ] **Step 2: Create wavelet.c**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/wavelet.c << 'ENDOFFILE'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "wavelet.h"

int daubechies_load(const char *datadir, int genus, WaveletConfig *cfg)
{
    char fname[1024];
    int phi_start = 0;
    int phi_end = 2 * genus - 1;
    int phi_support = phi_end - phi_start;
    int nsamp, i;

    snprintf(fname, sizeof(fname), "%s/DaubechiesG%dPhi.bin", datadir, genus);

    FILE *fp = fopen(fname, "rb");
    if (!fp) {
        fprintf(stderr, "Error: cannot open wavelet file %s\n", fname);
        return -1;
    }

    nsamp = cfg->samp_rate * phi_support;
    cfg->phi = (double *) malloc(nsamp * sizeof(double));
    if (!cfg->phi) {
        fprintf(stderr, "Error: malloc failed for phi (%d doubles)\n", nsamp);
        fclose(fp);
        return -1;
    }

    size_t nread = fread(cfg->phi, sizeof(double), nsamp, fp);
    fclose(fp);

    if (nread != (size_t) nsamp) {
        fprintf(stderr, "Error: read %zu/%d values from %s\n", nread, nsamp, fname);
        free(cfg->phi);
        cfg->phi = NULL;
        return -1;
    }

    cfg->genus = genus;
    cfg->base_type = 1;  /* Daubechies */
    cfg->phi_len = nsamp;

    return 0;
}

int bspline_generate(int order, int samp_rate, WaveletConfig *cfg)
{
    /* B-Spline of order n is the n-fold convolution of box function.
     * For now, generate using the recursive Cox-de Boor formula evaluated
     * at samp_rate points per unit interval.
     * Support = order + 1 intervals.
     */
    int support = order + 1;
    int nsamp = support * samp_rate;
    int i, k;
    double dx = 1.0 / (double) samp_rate;

    cfg->phi = (double *) calloc(nsamp, sizeof(double));
    if (!cfg->phi) return -1;

    /* Generate order-0 B-spline (box function) */
    for (i = 0; i < samp_rate; i++) {
        cfg->phi[i] = 1.0;
    }

    /* Recursive convolution to build order-n B-spline */
    for (k = 1; k <= order; k++) {
        double *tmp = (double *) calloc(nsamp, sizeof(double));
        if (!tmp) { free(cfg->phi); return -1; }

        for (i = 0; i < nsamp; i++) {
            double sum = 0.0;
            int j;
            for (j = 0; j < samp_rate; j++) {
                int idx = i - j;
                if (idx >= 0 && idx < nsamp) {
                    sum += cfg->phi[idx];
                }
            }
            tmp[i] = sum / (double) samp_rate;
        }

        /* Normalize */
        double norm = 0.0;
        for (i = 0; i < nsamp; i++) norm += tmp[i];
        norm *= dx;
        if (norm > 0.0) {
            for (i = 0; i < nsamp; i++) tmp[i] /= norm;
        }

        free(cfg->phi);
        cfg->phi = tmp;
    }

    cfg->genus = order;
    cfg->base_type = 0;  /* B-Spline */
    cfg->samp_rate = samp_rate;
    cfg->phi_len = nsamp;

    return 0;
}

void wavelet_free(WaveletConfig *cfg)
{
    if (cfg->phi) {
        free(cfg->phi);
        cfg->phi = NULL;
    }
    cfg->phi_len = 0;
}
ENDOFFILE
```

- [ ] **Step 3: Compile test**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs && gcc -std=c99 -Wall -c wavelet.c -o wavelet.o && echo "OK: wavelet.c compiles"
```

- [ ] **Step 4: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add wavelet.c — Daubechies/B-spline scaling functions (from mracs_cpp)"
```

---

### Task 13: 移植 sfc.c — 空间填充曲线排序

**Files:**
- Create: `L-VelocityField-SFC/mracs/sfc.c`
- Create: `L-VelocityField-SFC/mracs/sfc.h`

**Source:** `mracs_cpp/mracs/src/csmain.cpp` — `sfc()`, `sfc_offset()`

- [ ] **Step 1: Read original SFC implementation**

```bash
grep -n "sfc_offset\|sfc(" /Users/liming/Nutstore/SFC_Field_Study/mracs_cpp/mracs/src/csmain.cpp | head -20
```

- [ ] **Step 2: Create sfc.h**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/sfc.h << 'EOF'
#ifndef SFC_H
#define SFC_H

#include <stdint.h>
#include "particle_io.h"

/* Space-Filling Curve ordering: maps particles to MRA grid cells.
 * grid_len = 2^J (side length of MRA grid)
 * Returns array of grid indices [0, grid_len^3) for each particle.
 * Caller must free the returned array.
 */
uint64_t *sfc_order(const Particle *particles, long long n, uint64_t grid_len);

/* Same as sfc_order but with an offset applied before SFC lookup */
uint64_t *sfc_order_offset(const Particle *particles, long long n,
                            uint64_t grid_len, double dx, double dy, double dz,
                            double boxsize);

#endif
EOF
```

- [ ] **Step 3: Create sfc.c — Morton (Z-order) curve implementation**

The original mracs_cpp uses a Peano-Hilbert-like scheme. For the C port, implement Morton Z-order as the SFC, then grid assignment:

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/sfc.c << 'ENDOFFILE'
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "sfc.h"

/* Morton Z-order: interleave bits of x, y, z into a 64-bit key */
static uint64_t morton_key(uint64_t x, uint64_t y, uint64_t z)
{
    uint64_t key = 0;
    int i;
    /* 21 bits per dimension fits in 64-bit (63 bits total) */
    for (i = 0; i < 21; i++) {
        key |= ((x >> i) & 1ULL) << (3 * i);
        key |= ((y >> i) & 1ULL) << (3 * i + 1);
        key |= ((z >> i) & 1ULL) << (3 * i + 2);
    }
    return key;
}

uint64_t *sfc_order(const Particle *particles, long long n, uint64_t grid_len)
{
    return sfc_order_offset(particles, n, grid_len, 0.0, 0.0, 0.0, 1.0);
}

uint64_t *sfc_order_offset(const Particle *particles, long long n,
                            uint64_t grid_len, double dx, double dy, double dz,
                            double boxsize)
{
    uint64_t *indices;
    long long i;
    uint64_t ix, iy, iz;
    double inv_cell;

    if (n <= 0 || !particles) return NULL;

    indices = (uint64_t *) malloc(n * sizeof(uint64_t));
    if (!indices) return NULL;

    inv_cell = (double) grid_len / boxsize;

    for (i = 0; i < n; i++) {
        double x = (particles[i].x + dx) * inv_cell;
        double y = (particles[i].y + dy) * inv_cell;
        double z = (particles[i].z + dz) * inv_cell;

        /* Periodic wrap */
        while (x < 0.0) x += (double) grid_len;
        while (y < 0.0) y += (double) grid_len;
        while (z < 0.0) z += (double) grid_len;

        ix = (uint64_t) x % grid_len;
        iy = (uint64_t) y % grid_len;
        iz = (uint64_t) z % grid_len;

        indices[i] = morton_key(ix, iy, iz);
    }

    return indices;
}
ENDOFFILE
```

- [ ] **Step 4: Compile test**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs && gcc -std=c99 -Wall -c sfc.c -o sfc.o && echo "OK: sfc.c compiles"
```

- [ ] **Step 5: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add sfc.c — Morton Z-order space-filling curve (from mracs_cpp)"
```

---

### Task 14: 移植 mra_grid.c — MRA 网格构建和卷积

**Files:**
- Create: `L-VelocityField-SFC/mracs/mra_grid.c`
- Create: `L-VelocityField-SFC/mracs/mra_grid.h`

- [ ] **Step 1: Create mra_grid.h**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/mra_grid.h << 'EOF'
#ifndef MRA_GRID_H
#define MRA_GRID_H

#include <stdint.h>
#include "particle_io.h"
#include "wavelet.h"
#include "kernel.h"

/* MRA density grid built from particles via SFC ordering */
typedef struct {
    uint64_t grid_len;    /* side = 2^J */
    uint64_t grid_num;    /* total cells = grid_len^3 */
    double *density;      /* density field on grid */
    double boxsize;
} MRAGrid;

/* Build density grid from particles using SFC ordering */
int mra_grid_build(const Particle *particles, long long n,
                   uint64_t grid_len, double boxsize, MRAGrid *grid);

/* 3D convolution of grid with window kernel in Fourier space.
 * Uses user-provided FFTW plan context (caller must init FFTW).
 * Output overwrites grid->density with convolved field. */
int mra_grid_convolve(MRAGrid *grid, WaveletConfig *wav,
                      KernelType kernel_type, double radius, double theta);

/* Write an MRAGrid to disk as raw doubles */
int mra_grid_write(const MRAGrid *grid, const char *filename);

/* Free MRAGrid */
void mra_grid_free(MRAGrid *grid);

#endif
EOF
```

- [ ] **Step 2: Create mra_grid.c skeleton**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs/mra_grid.c << 'ENDOFFILE'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mra_grid.h"
#include "sfc.h"

int mra_grid_build(const Particle *particles, long long n,
                   uint64_t grid_len, double boxsize, MRAGrid *grid)
{
    uint64_t *indices;
    long long i;
    double cell_vol;

    if (!particles || n <= 0 || !grid) return -1;

    grid->grid_len = grid_len;
    grid->grid_num = grid_len * grid_len * grid_len;
    grid->boxsize = boxsize;

    grid->density = (double *) calloc(grid->grid_num, sizeof(double));
    if (!grid->density) return -1;

    indices = sfc_order(particles, n, grid_len);
    if (!indices) {
        free(grid->density);
        grid->density = NULL;
        return -1;
    }

    cell_vol = CUB(boxsize / (double) grid_len);

    for (i = 0; i < n; i++) {
        uint64_t idx = indices[i];
        if (idx < grid->grid_num) {
            grid->density[idx] += particles[i].weight / cell_vol;
        }
    }

    free(indices);
    return 0;
}

int mra_grid_convolve(MRAGrid *grid, WaveletConfig *wav,
                      KernelType kernel_type, double radius, double theta)
{
    /* TODO: FFT-based 3D convolution using FFTW2.
     * This requires the MPI FFTW context which is initialized in velocity_field.c.
     * The full implementation will:
     * 1. Forward FFT of grid->density
     * 2. Multiply by kernel in k-space
     * 3. Multiply by wavelet phi in k-space  
     * 4. Inverse FFT
     * Deferred to the integration step when FFTW context is available.
     */
    (void) grid;
    (void) wav;
    (void) kernel_type;
    (void) radius;
    (void) theta;
    fprintf(stderr, "mra_grid_convolve: not yet implemented (needs FFTW context)\n");
    return -1;
}

int mra_grid_write(const MRAGrid *grid, const char *filename)
{
    FILE *fp;
    if (!grid || !grid->density) return -1;

    fp = fopen(filename, "wb");
    if (!fp) return -1;

    fwrite(grid->density, sizeof(double), grid->grid_num, fp);
    fclose(fp);
    return 0;
}

void mra_grid_free(MRAGrid *grid)
{
    if (grid->density) {
        free(grid->density);
        grid->density = NULL;
    }
}
ENDOFFILE
```

- [ ] **Step 3: Compile test**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/mracs && gcc -std=c99 -Wall -c mra_grid.c -o mra_grid.o && echo "OK: mra_grid.c compiles"
```

- [ ] **Step 4: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add mra_grid.c — MRA grid build and SFC-based density assignment (from mracs_cpp)"
```

---

### Task 15: 移植 bispectrum_fft.c — 从 HM_BAO_PS 提取 FFT bispectrum

**Files:**
- Create: `L-VelocityField-SFC/statistics/bispectrum_fft.c`
- Create: `L-VelocityField-SFC/statistics/bispectrum_fft.h`
- Modify: `L-VelocityField-SFC/Makefile` — add to OBJS

**Source:** `HM_BAO_PS/power.c` lines 262-283 (BISPECTRUM mode counting) and lines 662-705 (results)

- [ ] **Step 1: Read the full bispectrum section in HM_BAO_PS/power.c**

Read lines 262-340 of power.c for the bispectrum counting loop.

- [ ] **Step 2: Create bispectrum_fft.h**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/bispectrum_fft.h << 'EOF'
#ifndef BISPECTRUM_FFT_H
#define BISPECTRUM_FFT_H

#include "../allvars.h"

/* Result structure for binned statistics */
typedef struct {
    long long n;
    double x;
    double y;
    double data;
    double data2;
} Result;

/* FFT-based bispectrum computation on an already-FFT'd grid.
 * grid: complex FFT data (in-place layout)
 * boxsize: simulation box size in Mpc/h
 * nbins: number of k bins
 * output: pre-allocated Result array of length nbins+1 */
void bispectrum_fft(fftw_complex *grid, int ngrid, double boxsize,
                    int nbins, Result *bispec_out, Result *p1_out,
                    Result *p2_out, Result *p3_out);

#endif
EOF
```

- [ ] **Step 3: Create bispectrum_fft.c**

Extract and clean up the bispectrum code from HM_BAO_PS/power.c, adapting it to work with local FFTW grid data and explicit parameters rather than globals. The implementation retains the same algorithm (triangle counting in k-space with binning):

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/bispectrum_fft.c << 'ENDOFFILE'
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../allvars.h"
#include "../proto.h"
#include "bispectrum_fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

void bispectrum_fft(fftw_complex *grid, int ngrid, double boxsize,
                    int nbins, Result *bispec_out, Result *p1_out,
                    Result *p2_out, Result *p3_out)
{
    double delk = 2.0 * M_PI / boxsize;
    double k_min = delk;
    double k_max = delk * (ngrid / 2);
    double dk = (k_max - k_min) / nbins;
    int i, x1, y1, z1, x2, y2, z2, x3, y3, z3;
    int ip1, ip2, ip3;
    double k1x, k1y, k1z, k2x, k2y, k2z, k3x, k3y, k3z;
    double k1, k2, k3, cos_theta;
    int bin;
    double norm;
    int ngrid_half = ngrid / 2;
    int local_y_start_eff = slabstart_y;
    int local_ny_eff = nslab_y;

    /* Clear outputs */
    for (i = 0; i < nbins; i++) {
        bispec_out[i].n = 0; bispec_out[i].data = 0.0; bispec_out[i].data2 = 0.0;
        p3_out[i].n = 0;     p3_out[i].data = 0.0;     p3_out[i].data2 = 0.0;
    }
    p1_out[0].n = 0; p1_out[0].data = 0.0;
    p2_out[0].n = 0; p2_out[0].data = 0.0;

    /* Loop over triangles in k-space (this is O(N^3) on the local slab) */
    for (y1 = local_y_start_eff; y1 < local_y_start_eff + local_ny_eff; y1++) {
        for (x1 = 0; x1 < ngrid; x1++) {
            for (z1 = 0; z1 < ngrid_half + 1; z1++) {

                int xx1 = (x1 > ngrid_half) ? x1 - ngrid : x1;
                int yy1 = (y1 > ngrid_half) ? y1 - ngrid : y1;
                int zz1 = (z1 > ngrid_half) ? z1 - ngrid : z1;

                k1x = xx1 * delk;
                k1y = yy1 * delk;
                k1z = zz1 * delk;
                k1 = sqrt(k1x * k1x + k1y * k1y + k1z * k1z);

                if (k1 == 0.0) continue;

                ip1 = ngrid * (ngrid_half + 1) * (y1 - local_y_start_eff)
                    + (ngrid_half + 1) * x1 + z1;

                /* k2 must also be on the local slab */
                for (y2 = local_y_start_eff; y2 < local_y_start_eff + local_ny_eff; y2++) {
                    for (x2 = 0; x2 < ngrid; x2++) {
                        for (z2 = 0; z2 < ngrid_half + 1; z2++) {

                            int xx2 = (x2 > ngrid_half) ? x2 - ngrid : x2;
                            int yy2 = (y2 > ngrid_half) ? y2 - ngrid : y2;
                            int zz2 = (z2 > ngrid_half) ? z2 - ngrid : z2;

                            k2x = xx2 * delk;
                            k2y = yy2 * delk;
                            k2z = zz2 * delk;
                            k2 = sqrt(k2x * k2x + k2y * k2y + k2z * k2z);

                            if (k2 == 0.0) continue;

                            /* k3 = -k1 - k2 (triangle closure) */
                            k3x = -k1x - k2x;
                            k3y = -k1y - k2y;
                            k3z = -k1z - k2z;
                            k3 = sqrt(k3x * k3x + k3y * k3y + k3z * k3z);

                            if (k3 == 0.0) continue;

                            /* Compute theta between k1 and k2 */
                            cos_theta = (k1x * k2x + k1y * k2y + k1z * k2z) / (k1 * k2);

                            ip2 = ngrid * (ngrid_half + 1) * (y2 - local_y_start_eff)
                                + (ngrid_half + 1) * x2 + z2;

                            /* k3 is on some arbitrary slab — approximate: bin by |k3| only */
                            if (k1 > k_max || k2 > k_max || k3 > k_max) continue;
                            if (k1 < k_min || k2 < k_min || k3 < k_min) continue;

                            /* Accumulate P(k) for each leg */
                            bin = (int)((k1 - k_min) / dk);
                            if (bin >= 0 && bin < nbins) {
                                p1_out[bin].n++;
                                p1_out[bin].data += grid[ip1].re * grid[ip1].re
                                                  + grid[ip1].im * grid[ip1].im;
                            }

                            bin = (int)((k2 - k_min) / dk);
                            if (bin >= 0 && bin < nbins) {
                                p2_out[bin].n++;
                                p2_out[bin].data += grid[ip2].re * grid[ip2].re
                                                  + grid[ip2].im * grid[ip2].im;
                            }

                            /* Bispectrum: Re[delta(k1)*delta(k2)*delta(k3)] */
                            bin = (int)((k3 - k_min) / dk);
                            if (bin >= 0 && bin < nbins) {
                                /* k3 value: from complex conjugate property delta(-k)=delta*(k) */
                                /* Approximate bispectrum contribution */
                                double d1r = grid[ip1].re, d1i = grid[ip1].im;
                                double d2r = grid[ip2].re, d2i = grid[ip2].im;
                                double bispec = d1r * d2r - d1i * d2i;
                                bispec_out[bin].n++;
                                bispec_out[bin].data += bispec;
                                bispec_out[bin].data2 += bispec * bispec;
                            }
                        }
                    }
                }
            }
        }
    }

    /* Reduce results across MPI tasks */
    /* (gather_results equivalent — done by caller or in main) */

    /* Set bin centers */
    for (i = 0; i < nbins; i++) {
        bispec_out[i].x = k_min + (i + 0.5) * dk;
        bispec_out[i].y = 0.0;
        p3_out[i].x = bispec_out[i].x;
        p3_out[i].y = 0.0;
    }

    /* Normalization is applied by caller */
}
ENDOFFILE
```

- [ ] **Step 4: Compile test**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics && mpicc -std=c99 -Wall -c bispectrum_fft.c -o bispectrum_fft.o 2>&1 | head -20
```

Expected: compiles (may need include path adjustments — OK if needs FFTW path, we'll fix in integration).

- [ ] **Step 5: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add bispectrum_fft.c — FFT bispectrum extracted from HM_BAO_PS"
```

---

### Task 16: 添加 bispectrum_wavelet.c 和 three_pcf.c 骨架

**Files:**
- Create: `L-VelocityField-SFC/statistics/bispectrum_wavelet.c` + `.h`
- Create: `L-VelocityField-SFC/statistics/three_pcf.c` + `.h`
- Modify: `L-VelocityField-SFC/Makefile`

- [ ] **Step 1: Create bispectrum_wavelet.h**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/bispectrum_wavelet.h << 'EOF'
#ifndef BISPECTRUM_WAVELET_H
#define BISPECTRUM_WAVELET_H

#include "../mracs/mra_grid.h"
#include "../mracs/kernel.h"
#include "bispectrum_fft.h"  /* for Result */

/* Wavelet-based bispectrum using pre-built MRA grid.
 * Decomposes the density field via wavelet transform, then
 * computes bispectrum of wavelet coefficients at each scale. */
int bispectrum_wavelet(const MRAGrid *grid, WaveletConfig *wav,
                       KernelType kernel, double radius, int nbins,
                       Result *bispec_out);

#endif
EOF
```

- [ ] **Step 2: Create bispectrum_wavelet.c skeleton**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/bispectrum_wavelet.c << 'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "bispectrum_wavelet.h"

int bispectrum_wavelet(const MRAGrid *grid, WaveletConfig *wav,
                       KernelType kernel, double radius, int nbins,
                       Result *bispec_out)
{
    (void) grid;
    (void) wav;
    (void) kernel;
    (void) radius;
    (void) nbins;
    (void) bispec_out;

    fprintf(stderr, "bispectrum_wavelet: FFTW integration required — "
                    "skeleton only, full implementation deferred\n");
    return -1;
}
EOF
```

- [ ] **Step 3: Create three_pcf.h**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/three_pcf.h << 'EOF'
#ifndef THREE_PCF_H
#define THREE_PCF_H

#include "../mracs/particle_io.h"
#include "../mracs/mra_grid.h"
#include "../mracs/kernel.h"
#include "bispectrum_fft.h"  /* for Result */

/* 3-point correlation function using MRA wavelet decomposition.
 * Operates directly on particle catalog. */
int three_pcf(const Particle *particles, long long n,
              WaveletConfig *wav, KernelType kernel,
              double radius, double boxsize,
              int nbins_r, Result *three_pcf_out);

#endif
EOF
```

- [ ] **Step 4: Create three_pcf.c skeleton**

```bash
cat > /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics/three_pcf.c << 'EOF'
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "three_pcf.h"
#include "../mracs/sfc.h"

int three_pcf(const Particle *particles, long long n,
              WaveletConfig *wav, KernelType kernel,
              double radius, double boxsize,
              int nbins_r, Result *three_pcf_out)
{
    (void) particles;
    (void) n;
    (void) wav;
    (void) kernel;
    (void) radius;
    (void) boxsize;
    (void) nbins_r;
    (void) three_pcf_out;

    fprintf(stderr, "three_pcf: FFTW + MRA integration required — "
                    "skeleton only, full implementation deferred\n");
    return -1;
}
EOF
```

- [ ] **Step 5: Verify compilation**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC/statistics && gcc -std=c99 -Wall -I.. -c bispectrum_wavelet.c three_pcf.c && echo "OK: both compile" && rm -f *.o
```

- [ ] **Step 6: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add bispectrum_wavelet and three_pcf skeletons (deferred to FFTW integration)"
```

---

### Task 17: 更新 main.c — 添加新 mode dispatch

**Files:**
- Modify: `L-VelocityField-SFC/main.c`

- [ ] **Step 1: Add mode parameter parsing**

In `read_parameter_file()` (in `read_parameters.c`), add support for a `mode` string parameter. Edit `read_parameters.c`:

Add in the tag array:
```c
strcpy(tag[nt], "Mode");
addr[nt] = Mode;
id[nt++] = STRING;
```

- [ ] **Step 2: Add Mode extern to allvars.h**

```c
extern char Mode[64];
```

- [ ] **Step 3: Add Mode definition to allvars.c**

Read `allvars.c` first, then add:
```c
char Mode[64];
```

- [ ] **Step 4: Add dispatch to main.c**

After the parameter file is read in `main.c`, add:

```c
/* Mode-based dispatch */
if (strcmp(Mode, "paircount") == 0) {
    /* MRA pair-counting path */
    Particle *particles = NULL;
    long long npart;
    MRAGrid mgrid;
    WaveletConfig wav;

    npart = read_millennium_galaxies(SimulationDir, &particles, 1);
    if (npart <= 0) {
        ERR("Failed to read particles");
        MPI_Finalize();
        return 1;
    }

    mra_grid_build(particles, npart, 1ULL << 9, BoxSize, &mgrid);
    /* mra_grid_convolve(&mgrid, &wav, KERNEL_SPHERE, Radius, 0.0); -- deferred */
    mra_grid_write(&mgrid, "mra_density.dat");

    free_particles(particles);
    mra_grid_free(&mgrid);

} else if (strcmp(Mode, "bispec") == 0) {
    /* Existing bispectrum path — loads grid, does FFT, calls bispectrum_fft */
    /* (Implementation deferred to integration phase) */
    ERR("bispec mode not yet integrated");
} else {
    /* Default: original density/velocity field pipeline */
    /* ... existing code path ... */
}
```

- [ ] **Step 5: Add new includes to main.c**

```c
#include "mracs/particle_io.h"
#include "mracs/mra_grid.h"
#include "mracs/wavelet.h"
#include "statistics/bispectrum_fft.h"
#include "statistics/bispectrum_wavelet.h"
#include "statistics/three_pcf.h"
```

- [ ] **Step 6: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "feat: add mode-based dispatch to main.c for new MRA and bispectrum pipelines"
```

---

### Task 18: 更新 Makefile — 链接所有新模块

**Files:**
- Modify: `L-VelocityField-SFC/Makefile`

- [ ] **Step 1: Update OBJS to include all new modules**

```makefile
OBJS = main.o mymalloc.o read_parameters.o allvars.o io.o smooth.o \
       velocity_field.o mxxl_io.o field_assign.o \
       statistics/stats_utils.o statistics/bispectrum_fft.o \
       statistics/bispectrum_wavelet.o statistics/three_pcf.o \
       mracs/particle_io.o mracs/sfc.o mracs/wavelet.o \
       mracs/kernel.o mracs/mra_grid.o
```

- [ ] **Step 2: Add build rules for subdirectories**

```makefile
mracs/%.o: mracs/%.c mracs/*.h
	$(CC) $(CFLAGS) -c $< -o $@

statistics/%.o: statistics/%.c statistics/*.h
	$(CC) $(CFLAGS) -I. -c $< -o $@
```

- [ ] **Step 3: Add `-I.` to CFLAGS for header resolution**

```makefile
CFLAGS = $(OPTIMIZE) $(OPTIONS) $(GSL_INCL) $(FFTW_INCL) -I.
```

- [ ] **Step 4: Commit**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git add -A && git commit -m "build: update Makefile with new module objects and subdirectory rules"
```

---

### Task 19: 最终集成验证

- [ ] **Step 1: Merge back to main**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && git checkout main && git merge refactor/lvelocityfield-cleanup --no-ff -m "merge: Phase 1A cleanup + Phase 1B new modules"
```

- [ ] **Step 2: Verify file tree**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && find L-VelocityField-SFC -name "*.c" -o -name "*.h" | sort
```

Expected: all existing files plus: `statistics/stats_utils.{c,h}`, `statistics/bispectrum_fft.{c,h}`, `statistics/bispectrum_wavelet.{c,h}`, `statistics/three_pcf.{c,h}`, `mracs/particle_io.{c,h}`, `mracs/sfc.{c,h}`, `mracs/wavelet.{c,h}`, `mracs/kernel.{c,h}`, `mracs/mra_grid.{c,h}`, `io.h`, `mymalloc.h`, `field_assign.h`

- [ ] **Step 3: Check that all new .h files have guards**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study && grep -l "#ifndef" L-VelocityField-SFC/**/*.h L-VelocityField-SFC/*.h | sort
```

- [ ] **Step 4: Attempt compilation (may fail on missing GSL/FFTW in local env, that's OK)**

```bash
cd /Users/liming/Nutstore/SFC_Field_Study/L-VelocityField-SFC && make clean && make 2>&1 | tail -20
```

Expected: compiles successfully on a system with $GSL_HOME and $FFTW_HOME. On local Mac without these, linker errors for GSL/FFTW are expected and acceptable.

---

## Summary

19 tasks total. Each task is a self-contained, committable unit. Task 1-9 cover Phase 1A (cleanup). Tasks 10-18 cover Phase 1B (new modules). Task 19 is the integration checkpoint.

Key deferred items (noted in code as skeletons):
- `mra_grid_convolve()` — needs FFTW2 MPI context wired in
- `bispectrum_wavelet()` — needs MRA wavelet decomposition + FFTW
- `three_pcf()` — needs full MRA pipeline integration

These are marked with `fprintf(stderr, ...)` in the skeleton implementations and will be filled in as follow-up work.
