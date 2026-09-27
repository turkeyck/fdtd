# 使用說明：輸入物理參數，直接執行 3D 斜入射 FDTD（均勻與非均勻網格）

這份說明給「想算一個結構的反射、透射或繞射效率」的使用者。你只需要寫一個 JSON 設定檔（波長、入射角、各層厚度與折射率、網格密度），`simulate.py` 會自動產生網格、選時間步、執行 C 求解器、量測 R/T，並和解析參考值比對。

- 網格類型用一個欄位切換：`"mesh": {"type": "uniform"}` 或 `"nonuniform"`。兩者走同一個 C 程式碼路徑（均勻是非均勻的特例，已驗證逐位元相同）。
- 驗證狀態見 [`validation_report.md`](../validation_report.md)；規格見 [`SPEC_nonuniform.md`](../SPEC_nonuniform.md)。

---

## 1. 環境與編譯

求解器是 C（double、OpenMP），分析是 Python（只用 numpy + matplotlib）。在 WSL Ubuntu 或 Linux 執行：

```bash
cd /mnt/c/Users/Sasha/onedrive_copy/coding/claude/FDTD/3D_oblique
make
```

`simulate.py` 若找不到執行檔也會自動 `make`。Meep 只在驗證（第 3 關）才需要，一般模擬不需要。

## 2. 三步驟快速開始

```bash
python3 simulate.py --example film > my_film.json     # 1. 取得範本（film | uniform | stack | spectrum | grating）
python3 simulate.py my_film.json --dry-run            # 2. 只建網格：看實際角度、格數、Δt、預估時間
python3 simulate.py my_film.json                      # 3. 執行
```

執行後的輸出（以 `examples/film_nonuniform.json`，633 nm、80 nm TiO2（n = 2.0）在玻璃（n = 1.46）上、θ = 30° 為例）：

```
λ = 633 nm  mesh = nonuniform  injection = aux line (exact)
  angle: requested θ = 30°, φ = 0°; realized θ = 30.0000°, φ = 0.000° (m = 1, n = 0, Lx = 1266 nm, Lz = 63.3 nm)
  cells 40 x 398 x 2 = 3.19e+04; Δy 7.273–15.83 nm (max ratio 1.090); Δx min 31.65 nm, Δz 31.65 nm
  Δt = 0.01992 fs (106 steps/period), 9540 steps; estimated 0.1 min per run
  s: R = 0.268586  T = 0.731416  R+T−1 = +1.5e-06 | TMM R = 0.269024  T = 0.730976 | exact discrete R = 0.268585
     (FDTD − discrete +1.2e-06, discretization -4.4e-04)
  p: R = 0.166054  T = 0.833947  R+T−1 = +1.3e-06 | TMM R = 0.166248  ...
```

每一欄的意義見第 6 節。所有範例都在 [`examples/`](../examples/)。

## 3. 設定檔欄位

長度一律可寫成帶單位的字串（`"633nm"`、`"0.08um"`、`"1.2mm"`）；只寫數字時單位是 nm。沒寫的欄位用預設值。

| 欄位 | 意義 | 預設 |
|---|---|---|
| `name` | 這次模擬的名字（輸出資料夾名稱） | 檔名 |
| `wavelength` | 真空波長；給一個字串，或一個清單（做光譜，每個波長一次執行） | 必填 |
| `angle_deg` | 入射角 θ（度），從傳播軸 +y 量起；0 < θ < 90 | 必填 |
| `azimuth_deg` | 方位角 φ = atan2(kz, kx)（度） | 0 |
| `polarizations` | `["s"]`、`["p"]` 或 `["s", "p"]` | `["s","p"]` |
| `layers` | 從入射側往下的各層：`{"name", "thickness", "n"}`；光柵層見第 7 節 | `[]` |
| `substrate_n` | 基板（出射側半無限介質）折射率 | 1.0 |
| `mesh.type` | `"uniform"` 或 `"nonuniform"` | `nonuniform` |
| `mesh.ppw` | 每個「介質內波長」的格數（points per wavelength） | 40 |
| `mesh.transverse` | x/z 間距：`"lambda/20"` 或長度 | 非均勻 λ/20；均勻 = y 的 Δ |
| `mesh.r_max` | 非均勻網格相鄰格比例上限 | 1.1 |
| `run.periods` | 模擬的光學週期數 | 90 |
| `run.dft_periods` | 最後多少週期做 DFT（取穩態相量） | 30 |
| `run.courant` | Courant 數 S（Δt = S·√3/√(Σ1/Δ²_min) 後再取每週期整數步） | 0.5 |
| `geometry.sf`, `tf`, `substrate` | 散射場區、真空全場區、基板區長度（單位 λ0） | 1, 3, 3 |
| `geometry.pml_cells` | 兩端 CPML 格數 | 20 |
| `angle.tol_deg` | 角度量化容許誤差（度） | 0.05 |
| `angle.max_period_cells` | x 或 z 週期的最大格數 | 400 |
| `output.dir` | 輸出資料夾 | `runs/user/<name>` |
| `output.field_map` | 是否輸出 x–y 截面場圖 | true |

命令列選項：`--dry-run`（只規劃）、`--pol s`（覆寫偏振）、`--example NAME`（印範本）。

## 4. 座標、角度與偏振約定

- y 是傳播軸（層的法線方向），光從 y 小的一側（真空）入射，往 +y 進入各層與基板；x、z 是週期方向。
- θ 從 +y 量起；φ = atan2(kz, kx)。s 偏振：E 垂直入射面；p 偏振：E 在入射面內。
- **角度會被量化。** x/z 是一般週期邊界（沒有 Bloch 相位），所以 kx = 2πm/Lx、kz = 2πn/Lz，Lx、Lz 必須是整數格。程式自動選最小的 (m, Lx)、(n, Lz)，讓實際角度誤差在 `angle.tol_deg` 以內，並印出實際 θ、φ 與週期長度。φ = 0 時 z 只有 2 格（z 方向沒有變化），計算量最小。
- 比較或做光譜時，請用印出的「realized」角度解讀結果；TMM／RCWA 參考值已經用實際角度計算。

## 5. 均勻 vs 非均勻網格：怎麼選

| | `uniform` | `nonuniform` |
|---|---|---|
| y 間距 | 全域同一個 Δ = λ/(ppw·n_max)，由最高折射率決定 | 每個介質各自 λ/(n·ppw)；介面兩側以等比例（≤ r_max）漸變 |
| 介面位置 | 必須落在格點上：厚度被**四捨五入到整數格**（程式會印出實際厚度） | 精確：每層厚度照輸入值，介面一定在主節點上 |
| x、z | 預設與 y 同一個 Δ（真正的立方格） | 預設 λ/20；光柵層的 x 在脊與槽內各自加密並漸變 |
| 適用 | 簡單結構、要與舊版結果逐位元比對 | 高折射率薄膜、多層膜、要求精確厚度、要省格數 |

同一片 80 nm TiO2 膜（633 nm，θ = 30°），相對於連續 TMM 的離散誤差：

| 網格（範例檔） | 間距 | 總格數 | 厚度 | R 相對 TMM 的誤差（s） |
|---|---|---|---|---|
| uniform（`film_uniform.json`，ppw 20） | 全域 15.8 nm（x、y、z） | 80 × 325 × 2 = 5.2e4 | 80 → 79.12 nm | −5.7e-3（相對同樣 79.12 nm 的 TMM） |
| nonuniform（`film_nonuniform.json`，ppw 40） | 膜內 7.3 nm、真空 15.8 nm、基板 10.8 nm；x/z 31.7 nm | 40 × 398 × 2 = 3.2e4 | 80 nm（精確） | −4.4e-4 |

非均勻網格把細格只放在需要的地方：總格數較少，誤差小 13 倍，而且沒有厚度四捨五入。兩者的 FDTD 與各自的精確離散解都一致到約 1e-6，差別完全來自網格本身。

**時間步。** Δt 由最小格決定（CFL），再取成「每週期整數步」，讓 DFT 在單一頻率精確。非均勻網格的最小格越小，步數越多：`--dry-run` 會印出步數與預估時間。

**加密收斂的做法。** 把 `mesh.ppw` 加倍重跑；非均勻網格的誤差是 O(Δ²)（驗證第 2 關實測階數 2.00–2.03），所以 R_∞ ≈ R_2 + (R_2 − R_1)/3 是很好的外插值。

## 6. 輸出怎麼讀

### 平面多層膜

| 印出的量 | 意義 | 正常範圍 |
|---|---|---|
| `R`, `T` | FDTD 量到的反射率（真空區前行／後行波擬合）與透射率（基板區守恆通量 / 入射通量） | — |
| `R+T−1` | 能量守恆（R、T 是獨立量測的） | ~1e-6（無損介質） |
| `TMM R, T` | 連續（真實物理）轉移矩陣結果，用實際角度 | — |
| `exact discrete R` | 同一個 Yee 網格的精確離散解（不含時間步進誤差以外的任何近似） | — |
| `FDTD − discrete` | 求解器與它自己的離散模型的差：檢查執行是否到穩態、注入是否正確 | ≲ 1e-5 |
| `discretization` | 離散解 − 連續 TMM = 網格造成的物理誤差；加密 ppw 來降低 | 依網格 |

均勻網格若有厚度四捨五入，會另外印出「在要求厚度下的 TMM」以及 FDTD 與它的差。

### 檔案

`runs/user/<name>/`：

- `summary.json`：每個波長、偏振的所有數字與網格資訊。
- `<λ>nm_<pol>/field_xy.png`：z = 0 的 x–y 截面，各電場分量在 t = 0 的實部（DFT 相量）；黑線是介面、綠虛線是 TF/SF 面、灰色是 PML。座標用實際（非均勻）格點位置。
- `spectrum.png`：`wavelength` 是清單時，R(λ)、T(λ) 與 TMM 曲線。
- 每個執行資料夾：`grid.json`（輸入網格）、`grid_used.json`（**所有後處理座標的唯一來源**）、`meta.json`（Δt、k、入射振幅…）、`dft_*.bin`（complex128 相量）、`log.csv`。

相同參數再跑一次會直接重用結果（參數、網格雜湊與執行檔都沒變時）。

## 7. 光柵（x 方向 lamellar，階段 B）

在 `layers` 裡放一層：

```json
{"name": "grating", "thickness": "300nm",
 "grating": {"period": "2000nm", "duty": 0.5, "n_ridge": 2.0, "n_groove": 1.0}}
```

- 光柵週期就是 x 週期 Lx，因此 **kx = 2πm/Λ**：可實現的入射只有 sin θ cos φ = mλ/Λ（m 整數，可為 0）。程式會選最接近的 m 並印出實際角度；kz 由 z 週期自由量化。
- 每週期一個脊（x ∈ [0, duty·Λ)）；只支援一層光柵，上下可以有其他平面層。
- 非均勻網格：x 節點落在脊的兩個邊上，脊內 λ/(n_ridge·ppw)、槽內 λ/(n_groove·ppw)，以 r_max 漸變；此時使用解析平面波注入（`inc=p`），並多跑一次同網格的真空正規化執行。
- 均勻網格：x 間距 Λ/N_x，占空比四捨五入到整數格（程式印出實際值），使用精確的輔助線注入。
- 輸出：各傳播階的反射／透射效率、總和、以離散守恆通量計算的能量檢查；若結構只有「一層光柵 + 基板」，另外列出 RCWA（161 階）參考值與最大差。
- 實測準確度（`examples/lamellar_grating.json` 的幾何：Λ = 2λ、脊 n = 2、厚 0.3λ、基板 1.46、A2 錐形入射，s 偏振）：

  | 網格 | 週期數（DFT） | 各階效率與 RCWA 最大差 | 能量（守恆通量）− 1 |
  |---|---|---|---|
  | uniform，ppw 16（厚度 300 → 312.5 nm，RCWA 也用 312.5 nm） | 200（60） | 5.5e-3 | −1.6e-5 |
  | nonuniform，ppw 20，x/z λ/20 | 200（60） | 4.0e-3 | −7.0e-5 |
  | nonuniform，ppw 40，z λ/40（驗證 B3） | 450（150） | 9.6e-4 | −1.2e-5 |

  非均勻網格從 ppw 20 到 40 的誤差比約 4，符合二階收斂。範例檔用 ppw 20、200 週期（每個偏振約 9 分鐘，含真空正規化執行）；要 1e-3 等級請用 ppw 40、`transverse` λ/40、450 週期（每個偏振約 2.5 小時）。
- **穩態。** 光柵可能有導模共振、衰減很慢。驗證時（`examples/lamellar_grating.json` 的幾何）90 週期的能量誤差約 1e-3，450 週期（DFT 取最後 150）降到 1e-5 等級。若能量檢查超過 1e-3，程式會提醒；請增加 `run.periods` 與 `run.dft_periods`，看效率是否還在變。

## 8. 限制（程式會檢查並給出錯誤訊息）

- 入射介質是真空；介質無損、無色散、非磁性，n ≥ 1（每個波長一個固定 n；做光譜時色散要自己為每個波長寫一份設定）。
- 0 < θ < 90°，實際 sin θ < 0.95；不支援正入射（s/p 基底在正入射退化）。
- 角度被 x/z 週期量化（第 4 節）；光柵只能在 sin θ cos φ = mλ/Λ。
- 單頻（連續波）：每個波長一次執行。
- 均勻網格上比 Δ/2 薄的層無法表示（會報錯並建議改用非均勻網格）。

## 9. 常見問題

| 訊息 | 原因與處理 |
|---|---|
| `angle must be in (0, 90)` | 不支援正入射；用很小的角度（例如 1°）時週期會很長，注意格數 |
| `no period fits within --max-cells` / 實際角度誤差大 | 放寬 `angle.tol_deg`、提高 `angle.max_period_cells`，或換 `mesh.transverse` |
| `layer ... thinner than Δ/2` | 均勻網格太粗：提高 `mesh.ppw` 或改 `nonuniform` |
| `layer L=... too thin for neighbour spacing` | 非均勻漸變放不進薄層：提高 `mesh.ppw` 或 `mesh.r_max`（≤ 1.2 建議） |
| `warning: energy off by > 1e-3` | 還沒到穩態（共振結構）：增加 `run.periods` |
| 執行時間比預估長很多 | 預估假設機器閒置。兩個 OpenMP 模擬同時跑會互相拖慢數十倍（實測約 30 倍）：請一次只跑一個，或用 `OMP_NUM_THREADS=4` 分配核心 |

## 10. 進階：直接使用各元件

- `grid_gen.py`：y 方向以 block 描述（`uniform` / `layer` / `grade`），x/z 用 `uniform_axis` 或 `periodic_axis_layers`；`python3 grid_gen.py --check FILE` 檢查網格檔。
- C 求解器：`./fdtd3d_oblique key=value ...`。舊版均勻介面仍可直接用：`nl=20 eps2=2.25 y1=60 pol=s ...`（半無限介質）；網格檔：`grid=FILE`。所有參數見 [`README.md`](../README.md)。
- `oblique.py`：真空中斜入射平面波的場圖（物理單位輸入）。
- 參考解：`tmm.py`（連續 TMM、精確離散 Yee TMM）、`rcwa.py`（conical RCWA）。
- 驗證：`python3 tests/run_nu_regression.py`（凍結雜湊、舊回歸、第 0–3 關、階段 B）。
