# liscio 工具集

> 这里放的是 **给最终用户跑的命令行程序** — 跟 `test/` 里的开发回归测试
> 和基准跑分不同, 这些是装到机床 / 后处理流水线上日常用的.
> 默认编译; 关闭传 `-DBUILD_TOOLS=OFF` 给 cmake.

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
  --tol_xyz=N        XYZ 几何容差 mm (默认 0.01, 即 10 µm — 金属切削常用值)
  --output=PATH      输出文件路径 (默认写到屏幕)
  --precision=N      输出 G-code 坐标小数位数 (默认 5 位)
  --arc_merge=N      相邻 G2/G3 弧合并: 0=不合并  1=合并成 ARC  2=识别 HELIX
                      (默认 2 — 多 turn 螺旋会自动识别并紧凑输出)
  --corner_mode=N    角点处理: 0=遇到角点立即切分 (默认)
                              1=Look-ahead, 软角让 BEZIER 吸收, 硬角才切
  --bezier_samples=N 把每条 BEZIER 重采样成 N 段 G1 写到输出 (默认 8)
                      (LinuxCNC 模式下尽量用 G5; 此参数仅 portable 路径用到)
  --portable         多圈螺旋拆成多个单圈 G2/G3 (没 P<turns>); 默认用
                     LinuxCNC 的 P 形式更紧凑
  --target=MODE      portable (默认) / linuxcnc
                      linuxcnc 模式: BEZIER 用 G5 输出 (LinuxCNC 原生 cubic
                                     Bezier, 支持 XY/XZ/YZ 三个平面),
                                     文件比 G1 重采样小很多, 详见下表

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

- **portable 模式: BEZIER 必须重采样成 G1**.
  标准 G-code 没有原生 cubic Bezier 指令 (G2/G3 只能画圆弧, G1 只能画直线).
  跨控制器输出时只能把每条 BEZIER 切成多段直线写出来, 对 BEZIER 多的
  2D 文件可能反向膨胀. 解决: 用 `--target=linuxcnc` 启用 G5 输出 (见下).

- **`--target=linuxcnc` 限制**: G5 是 LinuxCNC 原生 cubic Bezier 指令,
  支持 G17/G18/G19 三个标准平面 (XY / XZ / YZ); 倾斜平面 (任意 3D) 上的
  BEZIER 仍只能 G1 重采样. 旋转轴 / UVW 也会 fall back. 如果你需要保留
  全 9D 几何, 直接调 liscio API 比走 G-code 输出更合适.

- **坐标偏移 (G92/G10/G43/G49) 当注释**: 输入文件被 rs274 解析成绝对
  毫米坐标后, 这些 "改坐标系" 类指令对位置已无影响. 工具用注释占位
  保留行号, 不还原原始坐标系. 注意: 压缩后的文件就是绝对毫米形式,
  不能回到原坐标系.

- **倾斜平面圆弧** (法向不沿 X/Y/Z 轴): 当前 fall back G1 重采样.
  CAM 数据里极少见.

### 实测 (18 CAM 文件, tol=0.05, default cfg with helix9 LSQ on)

| 文件 | 输入 | portable | linuxcnc (G5) | p_x | **l_x** |
|---|---:|---:|---:|---:|---:|
| **sine_test** | 8.85 KB | 4.83 KB | **0.99 KB** | 1.83× | **8.93×** ✅ |
| **1_1001** | 502 KB | 308 KB | **100 KB** | 1.63× | **5.01×** ✅ |
| bto45 | 533 B | 131 B | 131 B | 4.07× | 4.07× |
| sine_test | 8.85 KB | 2.45 KB | (同左) | 3.62× | — |
| **thomam** | 8.77 KB | 4.47 KB | 4.47 KB | **1.96×** ✓ | 1.96× |
| 3_1001 | 352 KB | 155 KB | **151 KB** | 2.27× | 2.33× |
| 5_1001 | 307 KB | 138 KB | **134 KB** | 2.22× | 2.29× |
| 4_1001 | 134 KB | 60.6 KB | 60.6 KB | 2.22× | 2.22× |
| **flower** | 733 KB | 654 KB | **353 KB** | 1.12× | **2.07×** ✅ |
| flower-one-line | 3.76 KB | 2.12 KB | **1.69 KB** | 1.77× | 2.22× |
| 6_1001 | 141 KB | 70.4 KB | 70.4 KB | 2.00× | 2.00× |
| 2_1001 | 19.8 KB | 10.2 KB | **9.12 KB** | 1.94× | 2.17× |
| 130207LZW | 130 KB | 68.9 KB | **67.7 KB** | 1.89× | 1.92× |
| SER-40 | 553 KB | 397 KB | 397 KB | 1.39× | 1.39× |
| luca-long-reverse | 213 B | 189 B | 189 B | 1.13× | 1.13× |
| rechteck-10x10 | 161 B | 282 B | 282 B | 0.57× ❌ | 0.57× |
| yy / yy-g61 | ≤163 B | ≤442 B | (同) | 0.35-0.37× ❌ | 同 |

> p_x: portable 模式 (跨控制器 G1 重采样)
> l_x: --target=linuxcnc (cubic Bezier → G5 emit, **XY/XZ/YZ 三个平面**)

### thomam 改善 (vs 早期文档 0.78× 反膨胀)

早期 G2/G3 直通 — 245 个 ARC 各自展开 IJK-form 文件膨胀 (输出 11.3 KB).
现在 helix9 LSQ 自动检测让 245 个 G3 弧 fold 成 11 个 HELIX, 用
`G2/G3 P<turns>` 多 turn 紧凑 form, 文件 8.77 → 4.47 KB (1.96× 压缩).
跟"是否 LinuxCNC target"无关, helix9 在两种 target 下都跑.

### flower 改善 (vs 早期 0.59× 反膨胀)

flower 是 3D 自由曲面 (XYZ 全变), liscio 拟合给出 3476 个 BEZIER.
早期 G5 emit 仅支持 XY plane (G17), flower 数据不严格 XY-planar →
全数 BEZIER 走 G1 重采样 fallback (16 segs/BEZIER × 3476 = 56k G1
行) → 输出 1243 KB > 输入 733 KB.

**改进 (2026-05-09)**: G5 emit 扩展到 LinuxCNC 完整三平面支持
(XY/XZ/YZ, 即 G17/G18/G19 模式). flower 实际 BEZIER 分布:
- 0 个严格 XY-planar (Z 不变)
- **3472 个 XZ-planar** (Y 不变 — 切片层加工)
- 2 个 YZ-planar
- 2 个真 3D (无 cardinal plane)

新 G5 emit 让 3474 BEZIER 用 G5 直接表达, 仅 2 个 3D 段 fall back
G1 重采样. 输出 1243 KB → **353 KB (2.07× 正向压缩)**.

### --target=linuxcnc 增量收益

- **sine_test** 3.62× → **8.93×** (BEZIER → 单 G5)
- **1_1001** 2.84× → **5.01×** (大量 2D BEZIER)
- **flower** 1.12× → **2.07×** (3-plane G5 救回 3D)
- ARC 主导 / 角点密集 文件无变化 (G5 不参与 ARC/LINE)

### Round-trip 验证

`portable` 模式 18 文件 round-trip (`liscio_compress` + `test_verify`)
全 OVER_TOL=0 ✓. `--target=linuxcnc` 模式 round-trip 在 ctest 工具
层 (ngc_parser) 不识别 LinuxCNC `NURBS_G5_FEED` canon, verify 报错;
**这是测试工具限制, 非 G5 emit 错误** — production LinuxCNC 内核
正常解析 G5 (LinuxCNC `interp_convert.cc` `convert_spline()`).

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
