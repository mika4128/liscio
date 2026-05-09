[English](README.md) | **中文**

# liscio (中文名: 莉修)

<p align="center">
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green.svg" alt="MIT"></a>
  <img src="https://img.shields.io/badge/lang-C99-blue.svg" alt="C99">
  <img src="https://img.shields.io/badge/deps-none-brightgreen.svg" alt="no deps">
</p>

> **发音**: 意大利语 /ˈliʃʃo/ ≈ **莉修** (liscio 意为 "光滑 / 平滑"), 正契合轨迹平滑预处理的主题.

**CNC 轨迹规划器的纯 C99 几何预处理库 — 微段 → arc / cubic Bezier / rational Bezier / B-spline 拟合.**

## 特点

- 纯 C99, 零第三方依赖, 可移植到裸机 / RT 内核 (LinuxCNC RTAPI).
- **9 维 (XYZ + ABC + UVW)**, 每轴独立拟合, 5 轴 CAM 友好.
- 原语: **LINE**, **ARC** (Kasa LSQ), **cubic Bezier** (9D per-axis LSQ),
  **rational Bezier** (NURBS), **B-spline** (clamped cubic).
- Composite 递归分段 + Hoschek 重参数化, 逼近 Fanuc AICC II / Siemens COMPCURV 精度.
- 角点自动切分, 切向量输出, 流式 API, 固定内存 (无堆分配).
- 吞吐 14M waypoints / 秒 (x86 单核), 典型 CAM 压缩 10–25×.

## 工业对标

| 特性 | Fanuc AICC II | Siemens 840D | liscio |
|---|---|---|---|
| 微段 → 几何压缩 | ✅ Nano Smoothing | ✅ COMPCURV | ✅ |
| 9D (含 rotary/uvw) | ✅ | ✅ | ✅ |
| Cubic Bezier LSQ | ✅ | ✅ | ✅ |
| Composite 递归 | ✅ | ✅ | ✅ |
| Hoschek 重参数化 | ✅ | ✅ | ✅ |
| 角点自动切分 | ✅ | G641/642 | ✅ |
| Rational weights | ✅ | ✅ | ✅ |
| Exact arc → NURBS | ✅ | ✅ | ✅ |
| Free-form B-spline LSQ | ✅ | ✅ | ✅ |

## 推荐配置

绝大多数用户**只需调一个 `tol_xyz`**, 其余字段保持默认.

```c
liscio_cfg_t cfg;
liscio_cfg_default(&cfg);
cfg.tol_xyz            = 0.01;     /* 10 µm — 工业主流 sweet spot */
cfg.tol_abc            = 0.1;      /* 0.1° — 旋转轴 */
cfg.arc_subseg_samples = 0;        /* G2/G3 直通, 几何无损 (默认) */
liscio_ctx_t *ctx = liscio_create(&cfg);
```

### `tol_xyz` 怎么选

跟 LinuxCNC `[TRAJ] DEFAULT_LINEAR_TOLERANCE` / G-code `G64 P` 一致最省心:

| 场景 | 推荐 tol | 实测压缩 | 备注 |
|---|---:|---:|---|
| 木工 / 塑料 / 3D 打印 | 0.05 mm | 7.49× | CAM 默认下限 |
| **金属切削主流** | **0.01 mm** | **4.43×** | **推荐默认 ✅** |
| 模具 / 光学 / 半导体 | 1e-3 mm | 2.02× | 超精密 |
| 形位公差零容忍 | ≤1e-4 mm | ≈ 1.4× | 接近 passthrough |

> liscio 在双精度浮点底线 `1e-8 mm` 仍跑得动, 无性能惩罚. 紧公差时压缩
> 比逼近 1.40× 渐近线 (≈ 几何 passthrough). 详见

### G2/G3 处理: **直通** (`arc_subseg_samples = 0`)

CAM 的 G2/G3 弧调 `liscio_add_arc()` 直接通过 → emit `LISCIO_PRIM_ARC`,
**max_dev = 0 几何无损**. 下游 TP 拿解析弧公式做 jerk 规划比沿 Bezier
数值积分快得多.

`arc_subseg_samples > 0` (切分模式) 仅在三个条件**同时**成立时才打开:
1. CAM 输出大量短弧 (例如 biarc fillet)
2. 工件公差 ≥ 10 µm
3. 下游 TP buffer 紧张, 需减 prim 数

实测: 紧公差下切分模式**反向膨胀**最高 5×.  默认 0 是工业最优.

### 进阶: 邻接弧合并 + HELIX 原语 (opt-in)

```c
liscio_cfg_set_arc_merge(&cfg, LISCIO_ARC_MERGE_HELIX, 0, 0);
```

把相邻同心 + 同径 + 同 pitch 的 G2/G3 弧合并成单一 ARC (或 HELIX 当
arc_angle > 2π).  几何无损. 用户场景:

| 适用 | 不适用 |
|---|---|
| 用户手写 G-code 多段同心圆 | CAM chord-approximation 弦近似 |
| 螺纹加工 (lathe G33) | 自由曲面切削 |
| 同心多层钻孔 | 螺旋壁加工 (R 递减) |

### 进阶: Look-ahead corner 决策 (opt-in)

```c
liscio_cfg_set_corner_detection(&cfg, LISCIO_CORNER_LOOKAHEAD,
                                 15.0 /*soft°*/, 60.0 /*hard°*/);
```

把"角度判断"和"是否 split"解耦, 两级阈值:
- **软角** (15-60°): 不立即 split, 让 Bezier G1 fit 吸收 → 切向连续性大幅提升
- **硬角** (>60°): 立即 split (保 corner 意图)

实测 (16 文件 tol=0.05): 加权 tan_mean **20.4° → 14.4° (-30%)**, 代价
压缩比 7.49× → 6.05× (-19%). 受益最大 flower (-40%) / flower-one-line (-55%).

适用: 下游做 **G64 blend** 时, 把 jerk-limited 速度规划给得更平滑.
不适用: G61 精准停 + prim count 紧张. 默认仍 IMMEDIATE.

### G-code 事件流: RAPID + STOP 直通

liscio 输出 = 完整 G-code 事件流, 不只切削原语.  非切削事件**专用 API**:

```c
liscio_add_rapid(ctx, &start, &end, line);                   /* G0 */
liscio_emit_stop(ctx, LISCIO_STOP_DWELL, 1.5, &pos, line);   /* G4 P1.5 */
liscio_emit_stop(ctx, LISCIO_STOP_TOOL_CHANGE, 0, NULL, line); /* M6 */
liscio_emit_stop(ctx, LISCIO_STOP_PROGRAM_END, 0, NULL, line); /* M2/M30 */
```

输出 `LISCIO_PRIM_RAPID` / `LISCIO_PRIM_STOP` 原语, 下游 TP 单一入口
(`on_prim`) 按 `p->type` 分发: 切削原语进 jerk-limited 队列, RAPID 走
max-velocity, STOP 按 `p->stop_reason` 处理.

完整 LinuxCNC canon (STRAIGHT_FEED/STRAIGHT_TRAVERSE/DWELL/CHANGE_TOOL/
SET_ORIGIN_OFFSETS/...) → liscio API 映射表见

### 完整推荐配置流程图

```
你的 G64 P 是多少?
  ├─ 没设 ──────────────→ tol_xyz = 0.01 (LinuxCNC 默认 sweet spot)
  ├─ 0.05 (木工/塑料) ──→ tol_xyz = 0.05
  ├─ 0.01 (金属主流) ──→ tol_xyz = 0.01
  └─ 0.001 (模具) ─────→ tol_xyz = 0.001

CAM 文件含大量 G2/G3 弧?
  ├─ 否 (95% 文件) ────→ arc_subseg_samples = 0  ✅
  └─ 是 + tol≥0.01 ──→ 可试 arc_subseg_samples = 8 (实验对比再决定)
```

## 快速上手

```bash
cmake -B build
cmake --build build -j
ctest --test-dir build --output-on-failure
```

产出 `libliscio.a` + `libliscio.so`, 测试全 PASS.

### 最小用法

```c
#include <liscio/liscio.h>

static void on_prim(const liscio_primitive_t *p, void *user) {
    /* p->type = LINE / ARC / BEZIER / RBEZIER / SPLINE
     * 用 p->start / p->end / P0..P3 / cx,cy,... 消费拟合结果 */
}

liscio_cfg_t cfg;
liscio_cfg_default(&cfg);           /* tol_xyz=0.01, arc 直通 (推荐默认) */
liscio_ctx_t *ctx = liscio_create(&cfg);
liscio_set_callback(ctx, on_prim, NULL);

for (/* 每个 G01 段 */) {
    liscio_pose_t s = {...}, e = {...};
    liscio_add_line(ctx, &s, &e, feedrate, line_no);
}
/* G2/G3 弧用 add_arc 直通: */
for (/* 每个 G02/G03 段 */) {
    liscio_primitive_t arc = {...};      /* 含 cx,cy,radius,arc_angle,normal */
    liscio_add_arc(ctx, &arc, feedrate, line_no);
}
liscio_flush(ctx);
liscio_destroy(ctx);
```

### CMake 集成

```cmake
add_subdirectory(liscio)
target_link_libraries(your_target PRIVATE liscio m)
```

或 `find_package(liscio CONFIG REQUIRED)` 安装后.

### 可视化拟合结果

```bash
cmake --build build --target dump_ngcs    # test/data/*.ngc → CSV
cmake --build build --target plot_ngcs    # CSV → PNG (需 matplotlib)
ls build/ngc_out/*.png
# 灰点=原始 G01 / 蓝=LINE / 橘=ARC / 绿=BEZIER / 粉=SPLINE
```

单文件:
```bash
./build/test_dump_csv --tol_xyz=0.05 --samples=30 foo.ngc
python3 test/plot_fit.py --3d foo.ngc
```

## 布局

```
liscio/
├── include/liscio/liscio.h   公开 API
├── src/                       拟合算法实现 (arc/bezier/rbezier/bspline)
├── test/                      单元测试 + CAM 回归 + 可视化
└── results/                   测试结果 (PNG / CSV / sweep / 跨库实测)
```

## 工具

`liscio_compress` — **G-code → G-code 压缩器**, 跨控制器兼容. 把 CAM
输出的 NGC 文件压缩, 实测一般机加工文件 **1.6-1.9× 文件尺寸压缩**:

```bash
./build/liscio_compress --tol_xyz=0.01 --output=part.ngc.gz part.ngc
```

详见 [tool/README.zh.md](tool/README.zh.md).

## 添加新测试样本

丢 `.ngc` 到 [test/data/](test/data/), 重跑 `cmake -B build` 让 CMake
`file(GLOB)` 重新扫描. 自动纳入:

- `ctest verify_geometry` — 拟合质量门控 (OVER_TOL=0, DRIFT=0)
- `make dump_ngcs` / `make plot_ngcs` — CSV + PNG 可视化
- `./test_metrics --csv=report.csv test/data/*.ngc` — 量化指标

## License

MIT (SPDX `MIT`). 详见 [LICENSE](LICENSE).
