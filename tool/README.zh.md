# liscio 工具集

> 用户面向工具 (`BUILD_TOOLS=ON`).  跟 `test/` 目录的回归 / 基准不同, 这些
> 是给生产环境用的.

## liscio_compress — G-code → G-code 压缩器

把 CAM 输出的 G-code 通过 liscio 拟合压缩, 重新生成标准 G-code, **跨控制器
兼容** (LinuxCNC / Fanuc / Siemens / Haas / Mach3 等都能跑).

```
input.ngc                      output.ngc
---------                      ----------
G1 X0.001 Y0.000          →    G1 X0.001 Y0.000 F500
G1 X0.002 Y0.001                ...
G1 X0.003 Y0.003                G2 X10 Y10 I5 J0  ; 100 G1 折叠成 1 ARC
... (上千条 G1) ...              ...
G1 X10 Y10                      G3 X20 Y0 I0 J-10
G3 X20 Y0 R10                   M2
```

### Usage

```bash
./build/liscio_compress [options] <input.ngc>

Options:
  --tol_xyz=N        XYZ 拟合容差 mm (默认 0.01 = 10µm 工业 sweet spot)
  --output=PATH      输出文件 (默认 stdout)
  --precision=N      G-code 数字小数位 (默认 5)
  --arc_merge=N      G2/G3 合并: 0=off 1=ARC 2=HELIX (默认 2)
  --corner_mode=N    Corner 检测: 0=IMMEDIATE 1=LOOKAHEAD (默认 0)
  --bezier_samples=N BEZIER/SPLINE 重采样为 N 个 G1 (默认 16)
  --portable         拆 helical 多圈为多个单圈 G2/G3 (无 P<turns>);
                      默认用 LinuxCNC P 形式 (更紧凑)
  --target=MODE      portable (默认) / linuxcnc
                      linuxcnc 模式: XY-planar BEZIER 用 G5 cubic Bezier 输出
                                   (大幅减小 BEZIER-密集文件; 见下表)

Examples:
  # 标准用法 (LinuxCNC + 默认精度)
  liscio_compress --output=part.compressed.ngc part.ngc

  # 跨控制器兼容 (Fanuc / Siemens 等)
  liscio_compress --portable --output=out.ngc in.ngc

  # 高精度模式 (光学/半导体)
  liscio_compress --tol_xyz=0.001 --precision=7 --output=out.ngc in.ngc

  # 木工 / 3D 打印 (容差宽, 最大压缩)
  liscio_compress --tol_xyz=0.05 --bezier_samples=8 \
                  --output=out.ngc in.ngc
```

### 输入到输出 G-code 映射

| 输入 G-code | liscio 中间原语 | 输出 G-code |
|---|---|---|
| 多条相邻 G1 | LINE (collinear 折叠) | 单条 G1 |
| 多条 G1 (圆形分布) | ARC (Kasa LSQ 拟合) | G2/G3 X- Y- Z- I- J- |
| 多条 G1 (平滑曲线) | BEZIER (LSQ 拟合) | N 条 G1 (重采样) |
| 多条 G1 (复杂曲线) | SPLINE (B-spline LSQ) | N 条 G1 (重采样) |
| 一条 G2/G3 | ARC (直通) | G2/G3 X- Y- Z- I- J- |
| 多条同心 G2/G3 (启用 arc_merge) | HELIX (合并) | G2/G3 X- Y- Z- I- J- P\<turns\> |
| G0 rapid | RAPID | G0 X- Y- Z- |
| G4 dwell | STOP/DWELL | G4 P\<sec\> |
| M0/M1 | STOP | M0 / M1 |
| M2/M30 | STOP | M2 / M30 |
| M6 | STOP | M6 |
| G92/G10/G43/G49/G54-59 | STOP | (注释占位, 详见下) |

#### 已知限制

- **portable 模式 BEZIER → G1 重采样**: 标准 G-code 无原生 Bezier.
  对 BEZIER-密集 2D 输入, 输出反向膨胀.  解决方案:
  → `--target=linuxcnc` 启用 G5 emit (XY 平面 Bezier 直接 G5 输出,
    sine 文件实测 1.83× → **8.93×**).
- **`--target=linuxcnc` 仅适用于 XY 平面 Bezier**: G5 G-code 标准只
  支持 X/Y 两轴, 含 Z 变化的 3D Bezier (如 flower) 仍回退 G1 重采样.
  9D 旋转/UVW 也回退. 对全 9D 表达, 直接接 liscio API (绕过 G-code) 是更
  好的方案.
- **G92/G10/G43/G49 偏移事件**: 当前以注释形式占位, 不做坐标系还原
  (因为 ngc_parser 已用 rs274 把所有坐标转 absolute mm). 用户注意:
  压缩后的 NGC 是 absolute mm 形式, 不能回到原坐标系.
- **平面外弧 (tilted plane)**: 当前回退为 G1 重采样. 极少见.

### 实测 (18 CAM 文件, tol=0.05, default cfg with helix9 LSQ on)

`portable` mode 文件大小压缩 (helix9 默认开启, 复用此功能):

| 文件 | 输入 | 输出 | 压缩比 | 说明 |
|---|---:|---:|---:|---|
| bto45 | 533 B | 131 B | **4.07×** | 16 G1 → 1 LINE |
| **thomam** | 8.77 KB | 4.47 KB | **1.96×** ✅ | helix9 fold 11 HELIX, 大幅改善 |
| 3_1001 | 352 KB | 188 KB | 1.87× | |
| 4_1001 | 134 KB | 71.8 KB | 1.87× | |
| sine_test | 8.85 KB | 4.83 KB | 1.83× | 17 BEZIER |
| 5_1001 | 307 KB | 178 KB | 1.72× | |
| 6_1001 | 141 KB | 82.5 KB | 1.70× | |
| 2_1001 | 19.8 KB | 12.0 KB | 1.65× | |
| 1_1001 | 502 KB | 308 KB | 1.63× | |
| 130207LZW | 130 KB | 82.9 KB | 1.57× | |
| SER-40 | 553 KB | 402 KB | 1.38× | 雕刻文字 |
| luca-long-reverse | 213 B | 189 B | 1.13× | |
| flower-one-line | 3.76 KB | 3.68 KB | 1.02× | BEZIER 重采样 ≈ 持平 |
| rechteck-10x10 | 161 B | 282 B | 0.57× ❌ | 已是 G1 矩形, 无空间 |
| **flower** | 733 KB | 1243 KB | 0.59× ❌ | 3D BEZIER 重采样反膨胀 |

> **thomam 改善** (vs 早期文档 0.78× 反膨胀): helix9 LSQ 自动检测让
> 245 个 G3 弧 fold 成 11 个 HELIX (G2/G3 P\<turns\> 多 turn form),
> 文件大小近半. 早期 ARC 直通模式每弧重写为完整 IJK-form 故膨胀.

`--target=linuxcnc` 增量收益 (XY-plane BEZIER 用 G5 emit):
- **sine_test** 1.83× → **8.93×** (BEZIER → 单 G5)
- **1_1001** 1.63× → **4.68×** (大量 2D BEZIER)
- 3D BEZIER (flower, 4/6_1001) / ARC 主导 / 角点密集 文件无变化

**总结**:
| 数据特征 | 推荐模式 |
|---|---|
| **跨控制器需通用** (Fanuc/Siemens/...) | `portable` (默认) |
| **跑 LinuxCNC 且 2D 切削** | `--target=linuxcnc` ⭐ |
| **3D 自由曲面 (flower 类)** | 任一都不理想, G-code 标准本身限制 |
| **ARC 主导 (车铣典型)** | 任一相同, 选 portable 通用 |

### 几何精度保证

- 输出文件用 liscio 重新解析 + 拟合, max_dev 仍 ≤ tol_xyz
- ARC 原样直通, max_dev = 0
- BEZIER 重采样按 chord 长度 ≤ tol/2 控制, 实测 max_dev ≤ tol_xyz

### 自验证 (round-trip)

```bash
# 1. 压缩
./build/liscio_compress --tol_xyz=0.05 --output=/tmp/out.ngc input.ngc

# 2. 把压缩输出再过一遍 verify_geometry — 应仍 OVER_TOL=0
./build/test_verify --tol_xyz=0.05 /tmp/out.ngc
```

### 集成示例 (LinuxCNC HAL / production)

后处理流程:

```bash
#!/bin/bash
# CAM 后处理脚本: 取 CAM 原 G-code, 压缩到 LinuxCNC 用
RAW=$1
OUT=${2:-${RAW%.ngc}.compressed.ngc}

./liscio_compress \
    --tol_xyz=0.01 \
    --arc_merge=2 \
    --output=$OUT \
    $RAW

echo "Compressed: $RAW → $OUT"
ls -la $RAW $OUT | awk '{print $5, $NF}'
```

### 未来改进

- [ ] LinuxCNC G5/G5.1 (cubic Bezier / quadratic B-spline) emit, `--target=linuxcnc`
- [ ] biarc 近似 BEZIER (用 G2/G3 对替代 G1 重采样)
- [ ] R-form 弧 (短弧时比 IJK 更紧凑)
- [ ] 自适应 BEZIER 重采样密度 (按局部曲率)
- [ ] 保留原 G-code 注释 (rs274 当前丢弃)
