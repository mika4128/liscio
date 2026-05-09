[English](README.md) | **中文**

# liscio 测试结果

自动生成回归归档. 刷新命令:

```bash
cmake --build build --target liscio_results
```

跑 `test/data/*.ngc` 全部, 重生 CSV/PNG 拷到这里. CMake `metrics_csv`
target 在 `ALL` 中, 所以 `data/metrics.csv` 和 `data/metrics_g1.csv`
每次普通 build 自动保持最新.

## 结构

```
results/
├── README.md            ← 英文版
├── README.zh.md         ← 你正在读
├── images/              ← 每 ngc 一张 PNG (4 视角拟合可视化)
│   ├── thomam.ngc.png   ← helix 检测示例 (紫色 = HELIX 原语)
│   └── ...
├── data/
│   ├── metrics.csv      ← tol_xyz=0.01, helix 自动检测 ON
│   ├── metrics_g1.csv   ← tol_xyz=0.01, arc_subseg=8 (chord 采样模式)
│   └── xcompare.txt     ← 5 库横评 log
                         ← 功能深入文档
```

## 默认配置

`metrics.csv` 和 PNG 用默认 tol_xyz, helix 自动检测开
(`cfg.arc_to_helix_max = 16`). 没有文件特定 flag — `test/data/` 里
每个 `.ngc` (或 `.NGC`) 走同一管线:

```
G1     → 窗口 LSQ → ARC | HELIX | BEZIER | SPLINE | LINE
G2/G3  → arc-merge / arc-to-helix LSQ / 直通
```

**helix 自动检测**: fitter 缓冲连续 16 个 G2/G3 弧, 在 buffered 几何上
跑解析 LSQ 拟合 (中心共 axis 线, pitch 一致, radius 一致). 残差
≤ tol_xyz 则 emit 单个 `LISCIO_PRIM_HELIX`; 否则 fall back, 各 ARC
单独 emit (max_dev = 0 直通).

## 当前运行亮点

(完整 18 文件表见 `data/metrics.csv`.)

### Helix 检测 (thomam.ngc — 车床螺纹)

无 helix-detection 时此文件 254 ARC 直通. 自动检测后同几何 fold 成:

```
in=276  out=80  ratio=3.45×  HELIX=14  max_dev=0.0000  tan_mean=2.6°
```

PNG (`images/thomam.ngc.png`) 显紫色多 turn HELIX 原语.

## 质量门控 (ctest verify_geometry)

`OVER_TOL` 和 `DRIFT` 每文件必须 0. CI 每次 push 跑此 gate.

```bash
ctest --test-dir build -R verify_geometry --output-on-failure
```

## 手动重现

```bash
cmake -B build && cmake --build build -j
ctest --test-dir build --output-on-failure
cmake --build build --target liscio_results   # 刷新 PNG + CSV
```

## 添加新测试样本

把 `.ngc` 丢进 `test/data/`. CMake `CONFIGURE_DEPENDS` glob 自动
扫到, 所有 target 自动纳入.

## 深度文档

| 文档 | 主题 |
|---|---|
| [helix.zh.md](helix.zh.md) | **G2/G3 → HELIX LSQ 自动检测** (默认开启, 推荐) |
| [arc_merge.zh.md](arc_merge.zh.md) | 严格等于路径 (旧, opt-in) |
| [lookahead.zh.md](lookahead.zh.md) | 软/硬两级 corner 切片 |
| [events.zh.md](events.zh.md) | G-code 事件流 + LinuxCNC canon 映射 |
| [sweep.zh.md](sweep.zh.md) | tol_xyz × G2/G3 模式 8×2 网格 (历史快照) |
