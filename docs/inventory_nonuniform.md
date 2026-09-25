# 均勻網格假設盤點表（Stage 0）

基準版本：git tag `pre-nonuniform`（`fdtd3d_oblique.c` 843 行，凍結副本 `ref/fdtd3d_oblique_uniform_ref.c`）。
格式：`檔案:行號 | 假設內容 | 修改方式`。
「必改」＝ nonuniform 路徑會用到，必須改成讀節點陣列或局部間距。
「保留」＝ 只屬於舊的均勻測試或 uniform 產生器；在 uniform 模式下仍正確，不在 nonuniform 路徑上使用。

## 1. C 求解器 `fdtd3d_oblique.c`

| 檔案:行號 | 假設內容 | 修改方式 | 類別 |
|---|---|---|---|
| fdtd3d_oblique.c:63 | 全域 `D`（= 1/nl），被當成三軸共同間距 | 改名語意為 Δ_ref，只用於 mesh=uniform 產生 h 陣列與 meta 相容欄位；更新式不再引用 | 必改 |
| fdtd3d_oblique.c:162–168 | `ky_disc` 用 `D` 與 `P.S`（隱含 dt = S·D） | `ky_local(h, n)`：sin²(ky h/2) = (h/2)²[n²ω̃² − K̃x² − K̃z²]，ω̃ 由 dt 求得，不再經過 S | 必改 |
| fdtd3d_oblique.c:182–184 | `amplitudes` 的 K̃x、K̃y、K̃z 全用 `D` | K̃x 用 Δx、K̃z 用 Δz、K̃y 用 aux 源區間距 Δ_a | 必改 |
| fdtd3d_oblique.c:214–219 | `inc_amp` 的 y = (j + OFF)·D | y 取 aux 節點（主或對偶）座標陣列 | 必改 |
| fdtd3d_oblique.c:226–236 | aux 長度 `half = nsteps/2 + 20`（索引因果）；`ikxD = i·K̃x·D`、`ikzD = i·K̃z·D` | 長度規則保留（Yee 模板每步最多一格，與間距無關）；改存 iK̃x、iK̃z（不乘 D） | 必改 |
| fdtd3d_oblique.c:239–248 | `aux_update_H`：`c = dt/D` 同時乘 y 差分與 ikzD/ikxD 相量項 | y 差分乘 `Δt/h_{u+½}`，相量項乘 `Δt·iK̃`；兩者拆開 | 必改 |
| fdtd3d_oblique.c:251–252 | aux 1D TF/SF 修正係數 `c` | 用 j_a 處局部 Δt/Δ_a（源區必須均勻） | 必改 |
| fdtd3d_oblique.c:255–264 | `aux_update_E`：同上共用 `c` | y 差分乘 `Δt/d_u`，相量項乘 `Δt·iK̃` | 必改 |
| fdtd3d_oblique.c:266–267 | aux 1D TF/SF（E 側）係數 `c` | 同 251–252 | 必改 |
| fdtd3d_oblique.c:272–275 | `D = 1/nl`、`dt = S·D`、Courant 檢查 S < 1/√3 | Δt = S·√3/sqrt(Σ1/Δmin²)（均勻時退化為 S·D）；檢查保留；報告 Δt 代價 | 必改 |
| fdtd3d_oblique.c:276–279 | `Nx = Lx·nl`、`Nz = Lz·nl`，要求 L 為 D 的整數倍 | mesh=file 由 grid 檔讀 Nx、Nz、Ny 與 h 陣列 | 必改 |
| fdtd3d_oblique.c:296–298 | `Ny = 2 npml + sf + tf`、`j0 = npml + sf`、`j1 = j0 + y1`（以格數描述幾何） | mesh=file 由 grid zones 給 Ny、j0、npml_lo/hi；mesh=uniform 保留舊公式 | 必改 |
| fdtd3d_oblique.c:310–316 | 材料以 `j1`、`eps2` 描述；介面主節點 ε 取算術平均 (1+ε₂)/2 | ε_t[j] 以 h 加權平均：(h_{j−½}ε₁ + h_{j+½}ε₂)/(h_{j−½}+h_{j+½})；ε_n[j] 取該對偶格材料；uniform 時等於舊值 | 必改 |
| fdtd3d_oblique.c:317–318 | `ceEx = ceEz = dt/(ε_t D)`、`ceEy = dt/(ε_n D)`（材料與幾何合併、單一間距） | 拆成 cE_y[j] = Δt/(ε_t d_j)、cE_tx[j]、cE_tz[j]、cE_nx[j]、cE_nz[j] | 必改 |
| fdtd3d_oblique.c:320 | `ch = dt/D`（所有 H 分項共用） | cH_y[j] = Δt/h_{j+½}，cH_x = Δt/Δx，cH_z = Δt/Δz | 必改 |
| fdtd3d_oblique.c:329–338 | CPML：`d = npml·D`、`y = (j + h/2)·D`、`rho` 以索引距離乘 D | 以節點座標 y_j、y_{j+½} 與 PML 內界面座標求物理距離 ρ；d = N_pml·Δ_pml | 必改 |
| fdtd3d_oblique.c:341 | `smax = sig_fac(m+1)/(D·n_loc)` | 用該端 Δ_pml,lo / Δ_pml,hi | 必改 |
| fdtd3d_oblique.c:345–346 | CPML b、c 使用 dt | 不變（dt 已由新公式決定） | 保留 |
| fdtd3d_oblique.c:358–363 | 入射相位表 `p1 = kx(i+½)D + kz k D`、`p2 = kx i D + kz(k+½)D` | 使用 x、z 主/對偶節點座標陣列 | 必改 |
| fdtd3d_oblique.c:370–374 | `analytic_val` 座標 = (idx + OFF)·D | 節點座標陣列（debug 模式；init=a 只在 uniform 有精確意義） | 必改 |
| fdtd3d_oblique.c:392–398 | `dirichlet_y` 經 analytic_val 間接使用 D | 隨 analytic_val 改 | 必改 |
| fdtd3d_oblique.c:425 | Hy 更新：`ch·((δzEx) − (δxEz))` | `cH_z·δzEx − cH_x·δxEz` | 必改 |
| fdtd3d_oblique.c:431–432 | Hx、Hz 更新：`ch` 同時乘 y 差分與 x/z 差分 | Hx：`cH_y[j]·δyEz − cH_z·δzEy`；Hz：`cH_x·δxEy − cH_y[j]·δyEx` | 必改 |
| fdtd3d_oblique.c:435–442 | CPML 區 Hx、Hz：ψ 以差分儲存，乘 `ch` | ψ 仍存差分，乘 cH_y[j]；x/z 項乘 cH_x、cH_z | 必改 |
| fdtd3d_oblique.c:464–467 | Ey 更新 `ceEy[j]·(δzHx − δxHz)`；`epsy = dt/(c·D)` 反推 ε | `cE_nz[j]·δzHx − cE_nx[j]·δxHz`；能量權重改用 epsN[j] 與 W 權重 | 必改 |
| fdtd3d_oblique.c:474–483 | Ex、Ez 更新：單一 cx、cz 同乘 y 與 x/z 差分；`epsx = dt/(cx·D)` | Ex：`cE_y[j]·δyHz − cE_tz[j]·δzHy`；Ez：`cE_tx[j]·δxHy − cE_y[j]·δyHx` | 必改 |
| fdtd3d_oblique.c:485–494 | CPML 區 Ex、Ez 同上 | ψ 乘 cE_y[j] | 必改 |
| fdtd3d_oblique.c:449–451, 468, 472, 482, 494, 501–503 | 能量 ½Σ(εE^n·E^{n+1} + \|H^{n+½}\|²)，各格權重相同 | 各點乘 W_E = Δx·d_j·Δz（E 主節點 y 向用 d_j，Ey 用 h_{j+½}）等體積權重；另加 W_mod 欄位 | 必改 |
| fdtd3d_oblique.c:512–513 | `energy_sum` 乘 `D³` | 權重已逐點包含，改乘 1 | 必改 |
| fdtd3d_oblique.c:517–535 | `incident_E/H` 取 aux 陣列 `j0 − Ajlo` 索引 | 索引規則不變（aux 與主網格節點一一對應） | 保留 |
| fdtd3d_oblique.c:538–544 | 主網格 TF/SF（H 側）係數 `ch` | cH_y[j0−1]（j0 附近局部均勻） | 必改 |
| fdtd3d_oblique.c:546–552 | 主網格 TF/SF（E 側）係數 `ceEx[j0]`、`ceEz[j0]` | cE_y[j0] | 必改 |
| fdtd3d_oblique.c:565–572 | 切面 shape `Nx×SY`、`SY×Nz` 對所有分量相同 | 形狀保留；半格分量的 j = Ny 列寫 NaN，meta 記有效範圍 | 必改 |
| fdtd3d_oblique.c:597–613 | `dft_accumulate` type 1、type 2 對 Ey、Hx、Hz 也讀 j = Ny（該處從未更新，為無效值 0） | 半格分量迴圈到 Ny−1，j = Ny 輸出 NaN | 必改（bug） |
| fdtd3d_oblique.c:625–644 | `sf_max` 以索引範圍界定 SF 區 | 索引規則不變；SF 區由 zones 定義 | 保留 |
| fdtd3d_oblique.c:648–667 | `div_max`：各向差分和再整體除以 D | 非均勻算子：x 項 /Δx、y 項 /d_j（div E）或 /h_{j+½}（div H）、z 項 /Δz；`divop=uni` 保留錯誤版供對照 | 必改 |
| fdtd3d_oblique.c:716–757 | `meta.json` 只有 `Delta`、`nl`、`S` | 追加 mesh、grid_file、grid_hash、dt_courant、dt_penalty、每軸 Delta_min/max、r_max、auxref；另寫 grid_used.json | 必改 |
| fdtd3d_oblique.c:774–780 | init=b 高斯波包 y = j·D、x = i·D、z = k·D | 節點座標陣列 | 必改 |

## 2. Python 後處理與理論

| 檔案:行號 | 假設內容 | 修改方式 | 類別 |
|---|---|---|---|
| fdtd_io.py:88–109 | `slice_coords`：座標 = (idx + offset)·Delta | 新增 `coords(outdir)` 讀 `grid_used.json`；slice_coords 在 mesh=file 時改用節點陣列（mesh=uniform 結果不變） | 必改 |
| fdtd_io.py:112–121 | `valid_j_mask` 已排除半格分量 j = Ny | 保留；新 NaN 與此遮罩一致 | 保留 |
| fdtd_io.py:20–40 | 快取 key = 參數列表 + binary mtime | 參數含 grid 路徑時，另比對 grid_hash（grid 內容變了就重跑） | 必改 |
| fdtd_theory.py:39–46 | `Setup`：Delta = λ0/nl、dt = S·Delta、Nx = Lx/Delta | 舊測試保留；新增 `ky_local(h, eps, dt, Kt)`、`theta_local` 等以局部 h、dt 為參數的函式 | 保留 + 新增 |
| fdtd_theory.py:105–112 | `ky_discrete` 以單一 Delta 計算三軸 | 新函式 `ky_local` 分軸 | 保留 + 新增 |
| fdtd_theory.py:119–120 | `ktilde(k, Delta)` 單一間距 | 呼叫端分軸傳入 Δx、Δy、Δz | 保留 |
| fdtd_theory.py:137–151 | `yee_coords`：(idx + offset)·Delta | 非均勻時改用 grid_used.json | 保留（只給舊測試） |
| fdtd_theory.py:168–181 | `interp_factors`、`group_velocity` 用 Delta | 只在均勻子區有意義；nonuniform 分析不用 | 保留 |
| fdtd_theory.py:186–223 | `discrete_fresnel`：介面兩側 ±D/2、ε_if 算術平均 | tmm.tmm_discrete 為一般化版本；T-TMM1 要求均勻特例相差 < 1e-12 | 保留 + 一般化 |
| analyze.py:24 | `setup_from_meta` 以 nl、S 重建 Setup | nonuniform run 另寫 `grid_setup(meta, grid)` | 保留 + 新增 |
| analyze.py:148–155 | `div_from_planes`：兩相鄰 y 平面差分 /D | 新增 `div_nu(fields, grid)`；舊函式保留供舊測試 | 保留 + 新增 |
| analyze.py:160–164 | `conserved_flux(F, D)`：·D² | nonuniform：Δx·Δz（x、z 均勻）；新增 `conserved_flux_nu(F, dx, dz)` | 必改（新增） |
| analyze.py:205 | `theory_flux`：cos(ky·Delta/2) 胞心內插因子 | 只給均勻舊測試；nonuniform 用守恆通量 | 保留 |
| compare.py:28, 63–133, 174, 205, 239–241, 273 | `D = 1/NL`；座標 = (idx + offset)·D；Meep 取樣位置以 D 對齊 | 舊關卡 3 保留；第 3 關 (A)(B) 另寫函式，座標讀 grid_used.json | 保留（舊）+ 新增 |
| meep_ref.py:3–10, 50–52, 66–173 | `resolution = nl`，Meep 座標 = 整數倍 Delta，dt = 0.5·D | 保留給舊關卡 3；新第 3 關用 meep_ref_uniform.py、meep_ref_transform.py | 保留 |
| oblique.py:55–216 | 物理單位換算：Delta = λ/nl、dt = S·D/c、Ly = Ny/nl、圖座標 idx·Delta | 只支援 mesh=uniform（說明寫入 README）；nonuniform 不在本工具範圍 | 保留 |
| tests/level0_exact_injection.py:43–66, 81–117 | 解析殘差用單一 D 的 curl/div | 舊測試保留；新 nu_gate0 用 analyze 的非均勻算子 | 保留 |
| tests/level0b_multistep.py:23–67 | nl 決定網格 | 保留 | 保留 |
| tests/level1_leakage.py:30–38, 63–81 | 因果時窗 y = idx·D、geom(nl) | 保留；nu_gate1 以 grid 座標重寫 | 保留 |
| tests/level1_all.py:36–54, 61–113, 220–236, 323–324, 447–451 | 平面位置 j0 + f·nl；座標 idx·Delta；通量 D | 保留 | 保留 |
| tests/level2_fresnel.py:37–58, 65–69, 105–128 | layout(nl)、通量 D；D7 已改判定基準 | 保留 | 保留 |
| tests/level_cpml_stability.py:46, 67 | meta Delta 只用於紀錄與標題 | 保留 | 保留 |
| tests/level3_meep.py:32 | `--nl 20` | 保留 | 保留 |
| tests/make_report.py:63, 83–85 | 報告字串 Δ=λ0/nl | Part II 另寫產生器 | 保留 |
| fdtd_io_2.py、fdtd3d_oblique_2.txt、meep_ref_2.txt | 舊副本（無任何檔案 import 或引用；fdtd3d_oblique_2.txt 與 fdtd3d_oblique.c 內容相同） | 不動，不在任何路徑上 | 不需改 |

## 3. prompt 指定清單對照

| prompt 項目 | 位置 |
|---|---|
| 全域常數 D、dt/D、c = dt/D 的使用處 | C:63, 243, 259, 272–273, 317–320, 330, 341, 464, 475, 512–513, 666–667 |
| ikxD、ikzD | C:224, 235–236, 245–248, 261–264 |
| analytic_val、dirichlet_y、inc_amp | C:370–374, 392–398, 214–219 |
| aux_update_E/H 與 aux 的 1D TF/SF | C:239–268 |
| 主網格 TF/SF 修正係數 | C:538–552 |
| CPML 剖面（現況：以物理座標 (j+½h)·D 定義，但 d、σ_max 用 D） | C:329–350 |
| Courant | C:273, 275 |
| DFT 切面 metadata | C:565–572, 597–613, 742–748 |
| Python 座標、相位擬合、散度、Poynting 內插、通量積分 | fdtd_io.py:88–109；analyze.py:148–164, 205；fdtd_theory.py:137–223；compare.py 全檔 |
| meep_ref.py | 全檔（保留給舊關卡 3） |
| aux 中 y 差分與相量項共用 c | C:243–248, 259–264（最重要的必改項） |
| dft_accumulate 半格分量 j = Ny | C:609, 612 |

## 4. 統計

- C：37 列，其中「必改」34、「保留」3。
- Python：25 列，其中「必改」3（fdtd_io 座標與快取、analyze 通量），「保留 + 新增」7，「保留」14，「不需改」1。
- 原則：舊的 uniform 測試與工具繼續使用舊函式，確保舊回歸數值不變；nonuniform 分析一律走新函式，座標一律讀 `grid_used.json`。`tests/nu_lint_coords.py` 只掃描新檔案（nu_*、grid_gen、tmm、rcwa、meep_ref_uniform/transform、analyze 新函式區段）。
