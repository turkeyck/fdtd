# 規格整理 v 1.2.0

> 主題：在既有 3D 斜入射 FDTD（C、double；x/z PBC、y CPML、相量化複數 aux line 的 TF/SF）上加入 tensor-product 非均勻網格。
> 模式：互動模式。研究 checkpoint 已交付，使用者已於 2026-09-25 確認下列決策。本文件不取代既有 `SPEC.md`（均勻版），兩者並存；本文件的決策 D7 會修改舊版關卡 2 的判定基準。

```md
[關鍵概念定義]（查證日期 2026-09-25）
- Tensor-product 非均勻 Yee 網格：每軸一條 1D 主節點陣列。對偶節點取相鄰主節點中點。
  y 方向：主間距 h_{j+½} = y_{j+1} − y_j，對偶間距 d_j = y_{j+½} − y_{j−½} = (h_{j−½} + h_{j+½})/2。
  規則：y 差分寫進 H（Hx、Hz）時除以 h_{j+½}；寫進 E（Ex、Ez）時除以 d_j。
- 權重內積：W_E = 主邊長 × 對偶面積，W_H = 對偶邊長 × 主面積。在這組權重下 C_E 與 C_H 互為伴隨，
  蛙跳格式守恆 ½Σ εW_E|E^n|² + ½Σ μW_H H^{n+½}·H^{n−½}。現有程式的 ½Σ(εE^n·E^{n+1} + |H^{n+½}|²)Δ³ 是它的均勻特例。
- 相量化 aux line 的精確約化：x、z 均勻且為 PBC 時，e^{i(kx x + kz z)} 是離散 curl 的特徵函數，
  所以 3D 方程可精確約化成 y 向 1D 複數方程；y 是否均勻不影響這一點。只要 x 或 z 非均勻，約化就失效。
- 均勻子區離散色散：(2/h)² sin²(ky h/2) = ω̃²εμ − K̃x² − K̃z²。K̃y ≡ (2/h) sin(ky h/2) 在各子區相同，改變的只有 ky。
- 介面 ε 平均：介面落在主節點 j，切向 E 用 ε_j = (h_{j−½}ε₁ + h_{j+½}ε₂)/(h_{j−½} + h_{j+½})。
  現有 ifmode=a（算術平均）是它的均勻特例。
- Meep：官方 FAQ 明寫 "Meep only supports uniform Cartesian grids"。所以比對走 (A) 收斂值、(B) 變換光學等效材料兩條路。
```

```md
[競品 / 類似服務比較]
| 對象 | 非均勻網格做法 | 優點 | 缺點 | 本案借鏡 |
|---|---|---|---|---|
| Meep | 不支援 | Bloch k_point、epsilon_diag/mu_diag、material_function | 只能做收斂比對或變換材料 | 第 3 關 (A)(B) |
| Tidy3D AutoGrid | min_steps_per_wvl 預設 10（官方建議 ≥20），max_scale 預設 1.4 | 依材料折射率分格，並對齊層狀介面 | 閉源，TF/SF 細節不可見 | grid_gen 參數化（每材料 PPW + 相鄰比上限） |
| openEMS | SmoothMeshLines，ratio 預設 1.5 | 在固定線之間平滑插格 | 無斜入射 Bloch TF/SF | 「固定線 + 幾何級數插格」演算法骨架 |
| Lumerical FDTD | 自動非均勻 + 加密區 | 成熟 | 閉源 | 報告中 Δt 代價的呈現方式 |
```

```md
[GitHub / 開源 repo 比較]
| Repo | 技術路線 | 亮點 | 侷限 | 本案借鏡 |
|---|---|---|---|---|
| NanoComp/meep | C++ 均勻 Yee | k_point 不含 2π、Courant 可調、eps_averaging 可關 | 無非均勻網格 | meep_ref_uniform.py、meep_ref_transform.py |
| thliebig/openEMS + CSXCAD | C++ 非均勻 EC-FDTD | 每軸 1D 係數陣列 | 無相量 aux line | cE_y[j]/cH_y[j] 佈局 |
| gprMax | 均勻 Yee + GPU | GPU kernel 結構 | 不支援非均勻 | 日後移植 GPU 的記憶體佈局 |
```

```md
[建議方案與已確認決策]
- 範圍：階段 A（只有 y 非均勻）為主；階段 B（x、z 非均勻）寫成完整規格（D4）。
- D1（Q1=b）突變 r=4：理論值改用「精確離散接面反射」R_disc（tmm.py 的離散 Yee 約化解），
  門檻 |R_meas − R_disc|/|R_disc| < 1e-3。相速 Fresnel 估計與比值（預估 ≈ 3.0）只列 INFO。
- D2（Q2=c）漸變 r=1.1，底格 λ0/20：門檻 |R_meas − R_disc|/|R_disc| < 1e-3，
  且 |R| 在 λ0/20 → λ0/40 → λ0/80 細化下嚴格下降。−60 dB 只列 INFO（預估 −51.4/−56.6/−77.6 dB）。
- D3（Q3=ok）漸變區角度：|θ_meas − θ_pred(h_loc)| ≤ 0.1·|θ_pred − θ_cont|，且 |θ_meas − θ_cont| 的收斂階數擬合 ∈ [1.8, 2.2]。
- D5 交付檔名：SPEC_nonuniform.md（專案根目錄）。
- D6 第 2 關預設網格：y 方向每材料 PPW = 40（Δ ≤ λ0/(40n)），x/z = λ0/20。預估 |ΔR_s| ≈ 8e-5；PPW = 20 只用於收斂圖。
- D7 舊關卡 2：採舊報告選項 B，1% 門檻改對「精確離散 Fresnel」比較（目前 < 2e-7）。
  連續 Fresnel 的差距降為 INFO。舊驗證因此可以全部 PASS，並作為回歸基線。
- 規格補充（非門檻變更，執行前寫死）：aux line 複製主網格 y 節點 [0, Ny]，兩端以端點間距延伸；
  另設 aux_ref 線（同節點、同 ε、同 CPML），作為 3D 場的逐點精確參考。
- 待使用者事後確認（D8，本規格先採用，可否決）：graded 網格會產生真實的數值反射並進入 SF 區，
  所以第 1-1 關的「洩漏」定義為「SF 場 − aux_ref 預測散射場」。原字面量「SF 最大場」列 INFO。理由見 §7.3。
```

### 規格變更紀錄（實作期間）

| 編號 | 日期 | 內容 | 依據 | 狀態 |
|---|---|---|---|---|
| D9 | 2026-09-25 | 1-3b 收斂階數的細化系列改成 x/z 與 y 一起細化（`L1_taper_r1.1_base{10,20,40}_xz`；原訂 20/40/80 因 base80_xz 需約 2650 萬格 × 3.5 萬步（約 10 小時/次），在任何執行前改為 10/20/40，預測階數 2.003）；D3 比值判準仍同時在 x/z 固定系列上判定 | Stage 1 預檢：x/z 固定時橫向色散形成誤差底限，預測階數 0.26；x/z 一起細化（prompt 原文「Δ 加倍細化」）預測 2.000 | 待使用者事後確認 |
| D10 | 2026-09-25 | 2-3 的 y 細化系列改為 PPW=20 網格的巢狀二分（`L2_{F1,F5}_bis{1,2,4,8}`），Δt 同比例縮小，x/z 固定，階數對 Richardson R_∞(Δx) | 預檢：PPW 系列薄膜格數 4/7/13/26 不是 2 倍關係，誤差不單調（預測 2.9/2.4）；二分系列預測 2.00–2.01。修訂（執行前）：k 改為 1、2、4（k=8 每次約 2.5 小時），階數用三層估計 p = log₂((R₁−R₂)/(R₂−R₄))，預測 2.006–2.025 | 待使用者事後確認 |
| D11 | 2026-09-25 | 2-4 對照組改為均勻 Δ = 0.08/(M+½)，兩介面固定落在格子 0.25 與 0.75 位置，ε 點取樣（階梯），誤差對 D10 的 R_∞ | 預檢：「偏移 0.3Δ」在各解析度的偏移比例不同，預測階數 < 0；兩介面同比例時只是平移薄膜（2 階）；0.25/0.75 預測 0.97。修訂（執行前）：mh 取 4.5、9.5、18.5（37.5 每次約 4.4 小時），誤差對 D10 三層 R_∞，預測 0.96–0.97 | 待使用者事後確認 |
| D13 | 2026-09-26 | 3A-1 薄 cell 與完整 cell 的比對改在解析度 40（原 80），門檻 \|ΔR\| < 1e-6 不變 | 實測 Meep 單執行緒約 190 ns/格·步：完整 2×3 cell 在 80 約 12 小時/次，在 40 約 45 分鐘；於任何 Meep 執行前決定 | 待使用者事後確認 |
| D14 | 2026-09-26 | 3A-3 本程式外插值改用 `L2_{F1,F5}_bisxz{1,2}`（y 二分且 x/z 間距同時減半、Δt ∝ 1/k）的二層 Richardson（2 階，階數由 2-3 與理論確立），門檻 1e-4 不變 | PPW 20–160 且 x/z 同步細化需上億格；三層組合系列約 3500 萬格 × 2.4 萬步（約 6 小時）；預測 \|R_∞ − R_TMM\| ≤ 6e-6 | 待使用者事後確認 |
| D15 | 2026-09-26 | B3-3 RCWA 階數改為 41、81、161、321（門檻「最後兩階差 < 1e-5」不變），參考解用 321 階；光柵：週期 Λ = Lx = 2λ0、脊（n = 2）位於 x ∈ [0, 1)λ0、槽為真空、厚 0.3λ0、真空上覆、n = 1.46 基板（原規格未定基板） | 任何光柵 FDTD 前實測：41→81 階變化 3.7e-5、81→161 為 8.4e-6、161→321 為 2.7e-6 | 待使用者事後確認 |
| D16 | 2026-09-26 | 階段 B 測試配置：B0-2 用 x、z 皆漸變的 `B_ops_xz_graded_pec`；B1/B2 系列 `B_{vac,film}_k{1,2,4}` 採古典入射（m=1、n=0、θ=30°）、z 為薄的均勻格、x 漸變、y 均勻且對齊薄膜，s 偏振；B3 光柵維持 A2 錐形入射 | A2 角度（Lz=3）下三層 x/z 細化在 k=4 每次約 2500 萬格、每層 8 次執行；注入誤差的來源是 x 非均勻。於任何階段 B 執行前決定 | 待使用者事後確認 |
| D17 | 2026-09-26 | 1-4 判定用的基底 40/80 漸變執行，遠端 PML 維持固定物理厚度 3λ0（120／240 格；基底 20 的 60 格即 3λ0），門檻不變 | 失敗流程：60 格 PML 在各基底的回波約 1e-7，是基底 80（\|r\| = 1.5e-5）的 1%；改 3λ0 後回波 1.5e-9、相對誤差 1.7e-4（results/nu_failures.json） | 已採用（失敗流程） |
| D18 | 2026-09-26 | 3A Meep 解析度改為 100、200、400（所有介面落在像素邊界），門檻不變 | 失敗流程：Meep 1.34 在非對齊解析度下次像素平均未生效，40–320 只有一階收斂（F5 在 160 誤差 −1e-2）；對齊解析度為 2 階（F5 s：−5.8e-4、−1.4e-4） | 已採用（失敗流程） |
| D19 | 2026-09-26 | B3 判定用的光柵執行（注入 (i)）延長為 450 週期，DFT 取最後 150 週期（真空正規化仍 90／30；(ii) 以 90／30 列 INFO），門檻不變 | 失敗流程：90／30 時兩種注入的 Σ−1 都約 −1e-3、T0 差 −1.4e-3；RCWA 顯示 f ≈ 1.039 有導模共振，慢衰減的振鈴洩漏進 30 週期的 DFT 窗；180／30 已降到 −3.2e-4 | 已採用（失敗流程） |
| D12（實作註） | 2026-09-25 | 所有非均勻執行使用 Δt = T0/⌈T0/Δt_C⌉（每週期整數步，≤ §7.5 的 Courant 公式值，S_eff ≤ 0.5），與舊碼每週期 40 步的作法相同；理論值以此 Δt 重算（`v2/` 前綴），在任何量測前凍結。另加 `proj=1`：每一 y 列對橫向 Floquet 相位投影後的 DFT | 非整數週期的 DFT 窗會讓實數場的負頻映像以約 3e-3 洩漏進相量，遠大於 1-2（1e-8 rad）與 1-3a（1e-6）的門檻 | 已採用 |
| 實作註 | 2026-09-26 | 0-5 的「錯誤均勻算子」對照改為依偏振取法向分量不為 0 的場（s：div H；p：div E），並兩種偏振都判定；主判準（兩種算子外的正確算子 < 1e-10）不變 | s 偏振 Ey ≡ 0，div E 對 y 度規完全不敏感（實測兩種算子同為 6.6e-14），不能用來證明敏感度 | 已採用 |
| 實作註 | 2026-09-26 | 2-6（INFO）的「相同誤差」改以 y 離散誤差 \|R − R_∞(Δx)\| 比較（兩種網格共用 x/z = λ0/20，橫向誤差相同） | F1 預設網格的總誤差 7.4e-5 來自 y 誤差（+7.5e-5）與橫向誤差（−1.5e-4）相消；任何均勻網格的橫向誤差底限都高於此值，比較失去意義 | 已採用 |
| 實作註 | 2026-09-26 | 1-1 依凍結定義的時窗判定；另列整段執行（入射振幅已達 1）的 D8 洩漏與字面量為 INFO | L1 網格的時窗只有約 7 週期，入射仍在 erf 起步段（振幅約 1e-6），時窗內的數值沒有資訊量；1-0 已在整段執行上界住差值 | 已採用 |
| 實作註 | 2026-09-25 | 1-3b 前行波分解作用在主節點切向 E（s、p 皆同），與 FDTD 取樣位置一致 | 預檢：用對偶節點 H 作 p 偏振輸入時比值 ≈ 9；用主節點 E_t 時與 s 相同（≤ 4.5e-3） | 已採用 |
| 實作註 | 2026-09-25 | grid 的 sha256 由 Python（grid_gen.check、fdtd_io）驗證；C 只做結構檢查（節點數、均勻區、比例） | 避免在 C 內實作 sha256 | 已採用 |

以上都沒有改任何門檻數值，也沒有改理論值。D9–D16、D12 與 1-3b、sha256 兩條實作註是在任何非均勻 FDTD 執行之前，由 `tests/nu_predict.py` 的精確離散預檢或耗時估計發現的量測設計修正；D17–D19 與 0-5、2-6、1-1 三條實作註是在執行後經失敗流程（假設→最小實驗→修正→重跑）決定的，只改量測配置（PML 厚度、Meep 解析度、執行長度、對照場與 INFO 列的比較基準），證據記在上表「依據」欄，失敗流程的完整紀錄在 `results/nu_failures.json`。

---

## 技術規格文件

### 1. 假設與前提

| 編號 | 內容 | 類型 | 若不成立 |
|---|---|---|---|
| A1 | 單位沿用現有程式：c = ε0 = μ0 = 1、λ0 = 1、T0 = 1、ω0 = 2π。長度以 λ0、時間以 T0 報告 | 已確認（程式碼） | — |
| A2 | 預設角度沿用舊設定：m = 1、n = 1、Lx = 2、Lz = 3 → θ = 36.94°，φ = atan2(kz, kx) = 33.69° | 已確認（舊報告） | 所有理論值要重算 |
| A3 | Courant 數 S 的意義改為 Δt = S·√3 / sqrt(1/Δx_min² + 1/Δy_min² + 1/Δz_min²)。均勻時退化成 Δt = S·Δ，所以舊參數 S = 0.5 的結果逐位一致 | 規格決策 | 回歸 parity 失敗 |
| A4 | 5 層堆疊定義為 H L H L H：n_H = 2.0、d_H = 0.08λ0；n_L = 1.46、d_L = 0.12λ0；基板 n = 1.46 延伸進 PML | 假設（prompt 未定義） | 改 grids/stack5.json 即可 |
| A5 | 建置與執行環境沿用舊報告：WSL2 Ubuntu、gcc 13.3、`-O3 -march=native -std=c99 -fopenmp`；Meep 1.34.0（conda-forge、nompi）經 MEEP_PYTHON 呼叫 | 已確認（舊報告） | Meep 缺失時第 3 關為 BLOCKED |
| A6 | Python 只用 numpy + matplotlib；Meep 腳本例外，需要 `meep`。RCWA、TMM、線性求解都自行以 numpy 實作（不用 scipy） | 已確認（prompt） | — |
| A7 | 「全部舊驗證 PASS」依 D7 生效後的舊測試判定 | 已確認（D7） | — |
| A8 | 所有新門檻在第一次執行前寫死於 `tests/thresholds_nu.json`；執行後不得修改（只能由使用者決策改） | 規則 | 違反即為 FAIL |

### 2. 背景

現有程式已完成均勻網格版的全部驗證：注入精確到 4e-16、PML −110 dB、ky 相對誤差 2.5e-13。關卡 2 過去以連續 Fresnel 判定為 FAIL，D7 決策後改判。
現有程式把「均勻」寫死在很多地方：單一全域 D 與 ch = dt/D；aux line 的 c = dt/D 同時乘 y 差分與 ikxD；CPML 的 σ_max 用 D；能量乘 D³；散度除以 D；Python 端約 141 處用 Delta 或 idx*D。
薄膜堆疊需要在高折射率薄層內加密 y 網格，但不想讓整個計算域都變細。所以要加入 y 方向非均勻網格，同時保證均勻路徑只是非均勻路徑的特例，舊驗證全數維持。

### 3. 目標

| 編號 | 目標 | 可量化成功標準 |
|---|---|---|
| G1 | 單一程式路徑同時支援 uniform 與 nonuniform | 等間距 grid 檔與舊版凍結 binary 的所有場最大相對差 < 1e-12；舊驗證全 PASS |
| G2 | 非均勻 Yee 離散的數學性質被證明並以數值驗證 | 伴隨性 < 1e-13；1e5 步能量漂移 < 1e-10；Δt_max 預測準確（0.99 穩定、1.02 發散） |
| G3 | 相量 aux line 在 y 非均勻時仍是精確約化 | SF 洩漏（扣除預測散射場）< 1e-10；x–z 相位殘差 < 1e-8 rad |
| G4 | 數值色散與反射可被預測 | 子區 ky 相對誤差 < 1e-6；突變與漸變反射對 R_disc 相對差 < 1e-3 |
| G5 | 薄膜物理正確，且非均勻網格有效率優勢 | 預設網格 \|R−R_TMM\|、\|T−T_TMM\| < 1e-3；\|R+T−1\| < 1e-5；對齊介面收斂階數 ≈ 2 |
| G6 | 與 Meep 獨立交叉驗證 | Meep Richardson 外插值與 TMM 差 < 1e-4；變換光學 DFT 切向場相對 L2 < 1e-3 |
| G7（階段 B） | x、z 非均勻 + 光柵 | 各繞射階效率誤差 < 1e-3（對 RCWA）；\|ΣR_p + ΣT_p − 1\| < 1e-5 |

### 4. 範圍

**本次要做：**
- `grid_gen.py`：兩種模式。(a) 介面驅動 + 幾何級數漸變；(b) 映射模式 y = f(u)，給第 3 關 (B) 用。輸出 grid JSON。
- C 求解器改成每軸 1D 係數陣列，以 grid 檔驅動；未給 grid 檔時，在程式內產生均勻節點（feature flag `mesh=uniform|file`）。
- aux line 與 aux_ref 線改為非均勻 y 節點；主網格與 aux 的 TF/SF 係數用局部間距。
- CPML 以物理距離定義；Courant 以最小間距公式計算，並報告代價。
- 修正 `dft_accumulate` 在 j = Ny 讀半格分量的錯誤。
- 後處理：座標一律讀檔；非均勻散度、距離加權 Poynting、三點 ky 估計、局部前行波分解。
- `tmm.py`：連續 TMM（精確解）與離散 Yee 約化解 R_disc。
- 第 0～3 關全部測試、圖表與報告；階段 B 完整實作（注入 i、ii 必做，iii 選做）、Floquet 分解、`rcwa.py` 與光柵驗證。

### 5. 不做什麼

- 不做非 tensor-product 網格（局部 subgridding、八叉樹、conformal mesh）。
- 不做寬頻脈衝源的角度校正。固定 k_t 下只取 f0 分量，其餘頻率只在 Meep normalization run 中出現。
- 不做 GPU 移植；只保證記憶體佈局（每軸 1D 係數陣列、場陣列索引不變）方便日後移植。
- 不放寬任何門檻、不修改理論值。唯一例外是本文件開頭 D1–D8 已記錄的使用者決策。
- 不改 PBC ghost-layer 機制與 OpenMP 平行化策略。

### 6. Persona

| 角色 | 關心的事 | 痛點 | 本規格給的東西 |
|---|---|---|---|
| 計算電磁工程師（主要使用者，也是實作者） | 薄膜結構能用較少格點算準 | 均勻細網格太慢；非均勻網格常偷偷壞掉 TF/SF 與能量 | 精確約化的 aux line、可預測的色散與反射、逐關門檻 |
| 審查者（PI、指導教授、合作者） | 結果可信、可重現 | 看不出數字是不是挑出來的 | 事先寫死的門檻檔、驗證表、INFO 與 PASS 分開 |
| 下一位開發者 / coding agent | 改東西不會壞掉舊功能 | 均勻假設散落各處 | 盤點表、回歸 parity 測試、lint 禁止 idx*D |

### 7. 系統說明與核心設計

#### 7.1 網格與係數（階段 A）

- 主節點 y_j（j = 0..Ny）放 Ex、Ez、Hy；對偶節點 y_{j+½}（j = 0..Ny−1）放 Ey、Hx、Hz。
- 係數陣列（長度 Ny+1，SY 佈局不變）：

| 陣列 | 定義 | 使用處 |
|---|---|---|
| `cH_y[j]` | Δt/(μ h_{j+½}) | Hx 的 δyEz 項、Hz 的 δyEx 項；CPML 的 ψH 同樣乘它 |
| `cE_y[j]` | Δt/(ε_t[j] d_j) | Ex 的 δyHz 項、Ez 的 δyHx 項；CPML 的 ψE 同樣乘它 |
| `cE_tx[j]`, `cE_tz[j]` | Δt/(ε_t[j] Δx)、Δt/(ε_t[j] Δz) | Ez 的 δxHy 項、Ex 的 δzHy 項 |
| `cE_nx[j]`, `cE_nz[j]` | Δt/(ε_n[j+½] Δx)、Δt/(ε_n[j+½] Δz) | Ey 的 δxHz、δzHx 項 |
| 純量 `cH_x`, `cH_z` | Δt/(μΔx)、Δt/(μΔz) | Hy 的 δxEz（cH_x）與 δzEx（cH_z）；Hx 的 δzEy（cH_z）；Hz 的 δxEy（cH_x） |
| `epsT[j]`, `epsN[j]` | 主節點 ε（介面用 h 加權平均）、對偶節點 ε（單一材料） | 係數與能量權重 |
| `WE_y[j]`, `WH_y[j]` | d_j、h_{j+½}（y 方向權重因子；x、z 權重為 Δx、Δz） | 能量、伴隨性、散度 |

- 均勻特例：grid 檔直接給 h 陣列（等於 D，不以節點相減求得），所以 d_j = (D + D)/2 = D 逐位精確。
- 更新式由「共用係數乘差分之差」改成「兩個係數各乘一個差分」，例如 `ex += cE_y[j]*dHz − cE_tz[j]*dHy`。浮點運算順序改變造成的差距 ≤ 數個 ulp/步，回歸門檻 1e-12 可以吸收。

#### 7.2 aux line（相量化複數 1D 線）

- 節點：global j ∈ [0, Ny] 完全複製主網格節點；j < 0 以 h_{½} 均勻延伸，j > Ny 以 h_{Ny−½} 均勻延伸。長度沿用 `half = nsteps/2 + 20` 格。Yee 模板每步最多傳遞一格，這個索引因果界與間距無關。
- 更新式：y 差分乘局部 `Δt/h` 或 `Δt/d`；x、z 相量項乘 `Δt·iK̃x`、`Δt·iK̃z`（K̃ 由 Δx、Δz 決定），兩者拆開，不再共用 c。
- aux 源 j_a = j0 − na：要求 j_a ± 5 格均勻（間距 Δ_a）。inc_amp 的 ky 用 ky(Δ_a)，y 座標取 aux 節點座標。
- aux 本身的 1D TF/SF 係數使用 j_a 處局部 Δt/Δ_a。
- aux 材料為真空，只提供入射場。
- **aux_ref 線（新增，旗標 `auxref=1`）**：節點、ε_t/ε_n、CPML（σ、κ、α、ψ）都與主網格相同。在 j0 接收與主網格相同的 TF/SF 修正（入射值取自 aux line）。
  由於 CPML 只作用於 y、x/z 為 PBC、方程為線性，aux_ref 是整個主網格對 (kx, kz) 的精確約化，主網格場應等於 aux_ref × e^{i(kx x + kz z)}（只差捨入）。

#### 7.3 洩漏量測定義（D8）

漸變區會產生真實的數值反射，約 2.7e-3，通過 j0 進入 SF 區。這是正確的散射場，不是 TF/SF 洩漏。
所以第 1-1 關的量測是 max_SF |F_main − Re[aux_ref_SF · e^{i(k_t·x − ω0 t)}]| / E0，時窗照原文「入射波抵達遠端 PML 之前」（以 c 與 10% 餘裕計算）。
原字面量（SF 最大場 / E0）照樣記錄並列 INFO；另跑一組 TF 區全均勻的對照網格，以原字面量判定 PASS/FAIL。

#### 7.4 CPML

- PML 區與相鄰內部 ≥ 5 格都要均勻，間距分別為 Δ_pml,lo、Δ_pml,hi。
- ρ = 距 PML 內界面的物理距離，以 y_j、y_{j+½} 計算；d = N_pml·Δ_pml。
- σ_max = sig_fac·(m+1)/(η_loc·Δ_pml)，η_loc = 1/n_loc；遠端 n_loc 為基板折射率。κ、α 剖面同樣以 ρ/d 定義。
- ψ 仍以「差分」形式儲存，乘上局部 `cE_y` 或 `cH_y`。PML 區均勻，所以與舊碼的形式一致。

#### 7.5 Courant

Δt = S·√3/sqrt(1/Δx_min² + 1/Δy_min² + 1/Δz_min²)，其中 S < 1/√3。
報告 Δt 代價 = Δt_base/Δt，Δt_base 是以最大間距計算的 Δt。例：第 1 關 λ0/20 → λ0/80 漸變，Δt = 0.01021 T0，代價 2.45 倍步數。

#### 7.6 階段 B（x、z 非均勻）

- 係數陣列擴充為 `cE_x[i]`、`cH_x[i]`、`cE_z[k]`、`cH_z[k]`（x、z 方向週期）。週期接縫的對偶間距 d_0 = (h_{Nx−½} + h_{½})/2。
- 偵測規則：若 x 或 z 非均勻（max|h − mean h|/mean h > 1e-14），而 `inc=a` 或 `auxref=1`，程式以 exit code 2 中止並印出：
  `ERROR: phasor aux line requires uniform x and z (x nonuniform=%d, z nonuniform=%d); use inc=p|j|m`。
- 注入方式：
  - `inc=p`：在 TF/SF 面以實際分量座標取樣連續解析平面波（乘 ramp）。
  - `inc=j`：電流片 + normalization run。
  - `inc=m`（選做）：以 x–z 橫向離散 Bloch 本徵模態作為模態注入，再走 1D aux line，洩漏應回到機器精度。

### 8. 核心流程設計

1. **盤點（Stage 0）**：產出 `docs/inventory_nonuniform.md`，列出每一處「檔案:行號 | 假設內容 | 修改方式」。凍結舊碼 → `ref/fdtd3d_oblique_uniform_ref.c`。套用 D7 後舊回歸全 PASS。
2. **推導（Stage 1）**：`docs/derivation_nonuniform.md` 完成 §12 的 5 項推導。`tmm.py`、`fdtd_theory.py` 產出所有理論值。事先預測寫入 `results/nu_predictions.json`。
3. **網格（Stage 2）**：`grid_gen.py` 產生 `grids/*.json`，附檢查器與 lint。
4. **實作（Stage 3–5）**：C 改係數陣列 → 回歸 parity → 第 0 關單元測試 → aux/TF/SF/CPML 非均勻化 → 第 1 關。
5. **物理（Stage 6）**：薄膜與 5 層堆疊、TMM、收斂、誤差預算、效率 → 第 2 關。
6. **Meep（Stage 7）**：(A) 收斂值、(B) 變換光學 → 第 3 關。
7. **階段 B（Stage 8–9）**：只有在 A 全部 PASS 後才開始。
8. **回歸（Stage 10）+ 文件（Stage 11）**。
9. **失敗分支**：任一項 FAIL → 寫假設 → 設計最小實驗 → 修正 → 重跑全部測試（含舊測試）。每次紀錄寫入報告的「未通過紀錄」表；門檻不得改。

### 9. 開發應注意重點與應避開誤區

| 誤區 | 後果 | 正確做法 |
|---|---|---|
| y 差分與相量項共用 c = dt/D | 非均勻時 aux 不再是精確約化，洩漏約 1e-3 | y 用 Δt/h 或 Δt/d，相量用 Δt·iK̃ |
| 從節點相減得到均勻 h | 回歸 parity 在 1e-16 就開始漂移 | grid 檔直接存 h 陣列；節點由 h 累加（math.fsum） |
| 後處理用 idx*D 算座標 | 非均勻時圖與擬合全錯 | 一律讀 `out/grid_used.json`；lint 測試擋下 `*D`、`*Delta` |
| 介面 ε 仍用算術平均 | 非均勻網格只剩 1 階收斂 | 用 h 加權平均 |
| ky 用單點相位差量測 | 反射污染約 5e-3，遠大於 1e-6 門檻 | 均勻子區用三點恆等式 cos(ky h) = (E_{j+1} + E_{j−1})/(2E_j)；漸變區用局部前行波分解 |
| 以「½ 平均」內插 H 算 S_y | prompt 禁止 | 距離加權內插到 E 節點。理論上兩者在 1D 皆守恆（離散 Wronskian），½ 平均只列 INFO |
| 看完結果再挑量測區或 N_pml | 違反工作規則 | 量測區、N_pml、時窗全部寫在 `thresholds_nu.json`，先 commit 再跑 |
| 在 x/z 非均勻時仍用 aux line | 靜默錯誤 | 偵測並中止（§7.6） |
| 改動 log.csv 欄位順序 | 舊測試讀錯欄 | 只能在最後追加欄位 |

### 10. 任務模型與資訊優先級

本專案沒有互動式畫面。使用者的「介面」是命令列、驗證表與圖。

#### 任務模型表
| 層級 | 內容 | 為何屬於這一層 | 是否必須首屏支援 |
|---|---|---|---|
| 唯一主目標 | 讓 y 非均勻網格在不破壞舊功能的前提下，通過第 0–3 關 | 這是 prompt 的主要目標 | 是：報告第一個表就是各關總表 |
| 次目標 | 產出可預測的色散與反射理論，並被量測驗證 | 決定網格設計規則（r_max、PPW） | 是：第 1 關表 |
| 低頻目標 | 誤差預算分離、效率比較、Meep 變換光學 | 研究價值，不影響放行以外的使用 | 否：報告後段 |
| 罕見目標 | 階段 B 光柵、橫向 Bloch 模態注入 | 選做延伸 | 否：報告附錄 |

#### 資訊分類表
| 資訊項目 | 分類 | 使用頻率 | 是否首屏必須 | 不顯示的風險 |
|---|---|---|---|---|
| 各關 PASS/FAIL/BLOCKED 總表 | status-feedback | 高 | 是 | 無法判斷能否進下一關 |
| 驗證表（項目/理論/量測/門檻/判定 + 網格資訊） | action-critical | 高 | 是 | 無法追溯 |
| 圖（θ(y)、收斂 log-log、瞬時場） | decision-supporting | 中 | 否 | 看不出趨勢 |
| 盤點表、推導 | reference | 低 | 否 | 審查無依據 |
| 未通過紀錄 | exception-handling | 只在 FAIL 時 | 否（FAIL 時置頂） | FAIL 被掩蓋 |
| 執行環境、指令、耗時 | audit-history | 低 | 否 | 無法重現 |

#### 資訊架構表（validation_report.md）
| 資訊項目 | 使用頻率 | 是否首屏必須 | 所屬任務階段 | 顯示條件 | 建議容器 | 是否可收合 |
|---|---|---|---|---|---|---|
| 總結論 + 各關總表 | 高 | 是 | resolved / blocked | 永遠 | 頁首表格 | 否 |
| 未通過紀錄 | 只在 FAIL | 是（FAIL 時） | blocked | 有 FAIL | 頁首下方表格 | 否 |
| 各關驗證表 | 高 | 否 | validating / resolved | 永遠 | 各關章節 | 否 |
| 圖 | 中 | 否 | resolved | 永遠 | 表下方 | 否 |
| 盤點、推導、diff 說明 | 低 | 否 | reference | 永遠 | 附錄章節（`<details>`） | 是 |
| 環境與重現 | 低 | 否 | audit | 永遠 | 附錄 | 是 |

### 11. 狀態模型與揭露策略

每一個驗證項目與每一關都有下列狀態，寫在 `results/nu_*.json` 的 `status` 欄位。

#### 狀態矩陣
| State | 進入條件 | 使用者此刻目標 | 必顯資訊 | 隱藏資訊 | 主 CTA | 離開條件 |
|---|---|---|---|---|---|---|
| empty（not_run） | 尚未執行 | 知道還缺什麼 | 項目名、門檻 | 量測值 | `python3 tests/run_nu_regression.py` | 開始執行 |
| drafting（predicted） | 已在 nu_predictions.json 寫入理論值 | 確認理論值在跑前就固定 | 理論值、來源公式 | — | 執行該關 | 開始模擬 |
| validating（running） | 模擬中 | 等待 | 進度（步數/總步數） | 結果 | — | 模擬結束或中止 |
| resolved（pass） | 所有門檻達成 | 進下一關 | 理論/量測/門檻/網格資訊 | — | 下一關 | — |
| blocked（fail / blocked） | 任一門檻未達，或外部環境缺失（Meep） | 找原因 | 差距、未通過紀錄 | 後續關卡的放行 | 寫假設 → 最小實驗 | 修正後重跑全部測試並 PASS |
| submitted（reported） | make_report 產出報告 | 交付 | 總表 + 全部表與圖 | — | — | 程式或門檻變動 → 回到 empty |

#### Progressive disclosure 規則
- 報告頁首只放：總結論、各關總表、有 FAIL 時的未通過紀錄。
- INFO 列與 PASS/FAIL 列分開。INFO 不得出現在判定欄。
- 盤點表、推導、diff 說明放附錄，用 `<details>` 收合。
- 前一關未 PASS 時，後面各關若為了診斷而執行，一律標「僅供歸因，不構成放行」。

### 12. 推導規格（Stage 1 產出，寫進 `docs/derivation_nonuniform.md`，最後併入報告）

| # | 推導 | 必須包含的結論（驗收時逐條核對） |
|---|---|---|
| 1 | 非均勻 Yee 更新式 | 6 分量逐一寫出：每個差分項用 h 或 d 或 Δx、Δz，以及 ε_t 或 ε_n；CPML ψ 的係數；TF/SF 修正位置 |
| 2 | 散度、權重內積、伴隨、能量 | 散度：div E 在主節點 = δxEx/Δx + (Ey_{j+½} − Ey_{j−½})/d_j + δzEz/Δz；div H 在胞心，y 項除以 h_{j+½}。證明 ⟨E, C_H H⟩_{W_E} = ⟨C_E E, H⟩_{W_H}（逐項求和交換，PBC 與 PEC 邊界項為 0）。證明 𝓔^n = ½ΣεW_E\|E^n\|² + ½ΣμW_H H^{n+½}·H^{n−½} 守恆，並寫出 Δt 條件使 𝓔 為正定 |
| 3 | aux line 精確約化 | 代入 F(x,y,z) = f(y)e^{i(kx x + kz z)}，x、z 的 δ 化為 iK̃（因為 Δx、Δz 均勻、PBC）；y 的 δ 保持原樣，所以 1D 方程對任意 y 節點都與 3D 等價。說明 x、z 非均勻時 e^{ikx x} 不是 δx 的特徵函數 |
| 4 | TF/SF 係數 | 主網格：Hx[j0−1] += cH_y[j0−1]·EzI、Hz[j0−1] −= cH_y[j0−1]·ExI、Ex[j0] −= cE_y[j0]·HzI、Ez[j0] += cE_y[j0]·HxI；aux 源 inc_amp 使用 ky(Δ_a) |
| 5 | 色散預測 | ky(h) = (2/h) asin((h/2)·sqrt(ω̃²εμ − K̃t²))，θ(y) = atan(\|K̃t\|/ky_local)；數值 Fresnel 估計 n_num = sqrt(kt² + ky²)/k0（INFO）；R_disc 的定義與算法（§13.3） |

### 13. 核心資料模型

#### 13.1 grid JSON（`grids/<name>.json`，由 grid_gen.py 產生；C 與 Python 都讀）

```json
{
  "format": "fdtd-grid", "version": 1, "units": "lambda0",
  "name": "L1_taper_r1.1_base20",
  "generator": {"mode": "interfaces", "ppw": 20, "r_max": 1.1, "args": "..."},
  "x": {"uniform": true, "h": [0.05, "..."], "L": 2.0},
  "z": {"uniform": true, "h": [0.05, "..."], "L": 3.0},
  "y": {"h": [0.05, "..."], "y0": 0.0},
  "eps_y": [1.0, "..."],
  "interfaces_j": [412, 419],
  "zones": {"npml_lo": 20, "npml_hi": 60, "j0": 40, "ja": 0, "uniform_min_cells": 5},
  "hash": "sha256 of canonical h/eps arrays"
}
```

- `eps_y[m]`：第 m 個主格（y_m 到 y_{m+1}）的相對介電常數。ε_t 與 ε_n 由 C 依 §7.1 規則導出。
- C 讀檔：只支援本格式的「平面 key → 數字或數字陣列」子集（最多一層巢狀物件）。自寫約 150 行的讀取器，遇到未知 key 直接忽略，缺必要 key 則中止。
- 驗證器（`grid_gen.py --check`）：
  - Σh 與 L 相差 < 1e-13；
  - 相鄰比 ≤ r_max(1+1e-12)；
  - 介面 j 必須是主節點；
  - j0 ± 5、ja ± 5、兩端 PML 與相鄰 5 格內必須均勻（相對差 < 1e-14）；
  - 每一格 h ≤ λ0/(n·PPW)。

#### 13.2 執行輸出（每個 run 目錄）

| 檔案 | 內容 | 變更 |
|---|---|---|
| `meta.json` | 既有欄位 + `mesh`, `grid_file`, `grid_hash`, `dt_courant`, `dt_penalty`, 每軸 `Delta_min/Delta_max/r_max`, `auxref` | 追加欄位；舊欄位 `Delta` 在 file 模式下填 Δ_ref 並加 `Delta_is_reference: true` |
| `grid_used.json` | 各軸主/對偶節點座標、h、d、epsT、epsN、各分量有效 j 範圍 | 新增；後處理唯一座標來源 |
| `dft_<slice>.bin` | 形狀不變。半格分量在 j = Ny 寫 NaN | 修正 dft_accumulate |
| `dft_auxref.bin` | aux_ref 線 6 分量相量 | 新增 |
| `log.csv` | 既有欄位 + `W_mod`（prompt 的能量式）, `divE_TF_nu`, `divH_TF_nu` | 只在最後追加 |
| `fields_n<step>.bin` | 選定步數的全場（`dump_at=` 清單） | 新增，給回歸與能量測試用 |

#### 13.3 理論值資料 `results/nu_predictions.json`

每筆記錄：`{id, quantity, value, formula_ref, grid_hash, created_utc}`。
在第一次模擬前由 `tests/nu_predict.py` 產生，之後唯讀。若 grid 改動，必須換新 id，舊 id 保留。
R_disc 的算法：以與 C 相同的 y 節點、ε_t、ε_n、Δt、K̃x、K̃z 組成 1D 離散約化方程（s：E_s 在主節點；p：H_s 在對偶節點）。兩端以均勻區的精確離散出射模態封閉，用自寫 Thomas 三對角解法求 R、T。

#### 13.4 門檻檔 `tests/thresholds_nu.json`

所有 §20 的門檻、量測區定義（索引範圍規則，不是數值）、時窗規則、N_pml 設定。Stage 0 建立，commit 後唯讀。
任何修改都必須在報告的「規格變更紀錄」註明使用者決策編號。

### 14. State 管理與持久化

| 項目 | 暫存 / 持久 | 位置 | 恢復方式 |
|---|---|---|---|
| 模擬結果 | 持久 | `runs/<key>/`，key = hash(參數 + grid_hash + binary mtime) | fdtd_io 發現同 key 且 binary 未更新就重用；`--fresh` 強制重跑 |
| grid 檔 | 持久（版本化） | `grids/` | 修改即產生新 hash → 舊 run 自動失效；刪除 grid 檔時，引用它的 run 在下次回歸時標 stale 並重跑 |
| 理論值 | 持久、唯讀 | `results/nu_predictions.json` | 新增只能 append |
| 驗證結果 | 持久 | `results/nu_gate*.json` | make_report 重新組表 |
| 重新開始 | — | `python3 tests/run_nu_regression.py --fresh` | 清除 runs 快取並全部重跑；不刪 predictions 與 thresholds |

新增、修改、刪除對照：
- grid：新增用 grid_gen；修改是重新產生，會得到新 hash；刪除要一併刪除引用它的 runs。
- run：新增由 fdtd_io 觸發；參數改變時會自動重跑；以 `--fresh` 或刪目錄移除。
- prediction：只能 append。修改必須另起新 id 並記錄使用者決策。

### 15. 專案目錄規劃

```text
3D_oblique/
├── fdtd3d_oblique.c              # 求解器（單一路徑：uniform 為 nonuniform 特例）
├── Makefile                      # 新增 target: ref（建置凍結參考 binary）
├── ref/
│   ├── fdtd3d_oblique_uniform_ref.c   # Stage 0 凍結的舊碼（唯讀）
│   └── nonuniform.diff                # Stage 11 產出：ref → 新碼的完整 diff
├── grid_gen.py                   # 網格產生器 + 驗證器（--check）
├── grids/                        # *.json，檔名 = 用途_參數（例：L1_taper_r1.1_base20.json）
├── fdtd_theory.py                # 既有 + 非均勻色散、θ(y)、數值 Fresnel 估計
├── tmm.py                        # 連續 TMM、離散 Yee 約化 R_disc（s/p）
├── rcwa.py                       # 階段 B：conical RCWA（Li 正確分解法則）
├── fdtd_io.py                    # 既有 + 讀 grid_used.json、快取 key 含 grid_hash
├── analyze.py                    # 既有 + 非均勻散度、三點 ky、前行波分解、距離加權 S_y、Floquet 分解
├── meep_ref.py                   # 舊關卡 3 用（保留不動）
├── meep_ref_uniform.py           # 第 3 關 (A)
├── meep_ref_transform.py         # 第 3 關 (B)
├── compare.py                    # 既有 + 第 3 關 (A)(B) 比對
├── docs/
│   ├── inventory_nonuniform.md   # 盤點表
│   ├── derivation_nonuniform.md  # 5 項推導
│   └── DIFF_NOTES.md             # C 修改逐段說明
├── tests/
│   ├── (既有 level*.py 保留)
│   ├── thresholds_nu.json        # 門檻（唯讀）
│   ├── nu_predict.py             # 產生 nu_predictions.json
│   ├── nu_lint_coords.py         # 禁止 idx*D
│   ├── nu_gate0_regression.py / nu_gate0_adjoint.py / nu_gate0_energy.py / nu_gate0_stability.py / nu_gate0_divergence.py
│   ├── nu_gate1_vacuum.py        # 1-0..1-7
│   ├── nu_gate2_film.py          # 薄膜、5 層、收斂、誤差預算、效率
│   ├── nu_gate3_meep.py          # (A)(B)
│   ├── nuB_gate*.py              # 階段 B
│   ├── run_nu_regression.py      # 先跑舊 run_all_regression，再跑全部 nu_*
│   └── make_report.py            # 產生合併版 validation_report.md
├── runs/  results/  figures/nu/  # 產出物（可重建）
├── AGENTS.md / CLAUDE.md         # 長期規則（Stage 0 建立）
└── validation_report.md          # Part I 舊均勻版（D7 後）+ Part II 非均勻
```

命名原則：
- 新測試以 `nu_` 開頭（階段 B 用 `nuB_`），圖放在 `figures/nu/<gate>/`。
- 結果 JSON 命名為 `results/nu_gate<k>.json`。
- 舊檔名不改，避免舊測試失效。

### 16. 模組與資料流

| 模組 | 責任 | 輸入 | 輸出 |
|---|---|---|---|
| grid_gen.py | 產生並驗證網格 | 介面、材料、PPW、r_max，或映射 f(u) | grids/*.json |
| fdtd3d_oblique.c | 時域求解 | key=value 參數 + grid 檔 | runs/<key>/ |
| fdtd_theory.py / tmm.py / rcwa.py | 理論值 | grid 檔、角度、偏振 | nu_predictions.json |
| fdtd_io.py | 執行與快取、讀檔 | 參數 | numpy 陣列 + 座標 |
| analyze.py | 量測 | run 輸出 | 量測值 |
| meep_ref_*.py | Meep 參考 | 參數、grid 檔（映射） | runs/meep_nu/ |
| compare.py | 第 3 關比對 | 兩邊輸出 | 量測值 |
| tests/nu_*.py | 判定 | 量測、理論、門檻 | results/nu_gate*.json、圖 |
| make_report.py | 報告 | results、圖、docs | validation_report.md |

```svg
<svg viewBox="0 0 1200 720" xmlns="http://www.w3.org/2000/svg" font-family="Noto Sans TC, sans-serif">
  <rect width="1200" height="720" fill="#f8fafc"/>
  <text x="600" y="40" text-anchor="middle" font-size="24" fill="#0f172a" font-weight="bold">非均勻網格 FDTD：模組與資料流</text>
  <!-- 輸入層 -->
  <rect x="40" y="80" width="240" height="120" rx="10" fill="#e0f2fe" stroke="#0369a1"/>
  <text x="160" y="115" text-anchor="middle" font-size="18" fill="#0c4a6e" font-weight="bold">grid_gen.py</text>
  <text x="160" y="145" text-anchor="middle" font-size="14" fill="#0c4a6e">介面、PPW、r_max / f(u)</text>
  <text x="160" y="170" text-anchor="middle" font-size="14" fill="#0c4a6e">→ grids/*.json（含 hash）</text>
  <rect x="40" y="240" width="240" height="100" rx="10" fill="#e0f2fe" stroke="#0369a1"/>
  <text x="160" y="275" text-anchor="middle" font-size="18" fill="#0c4a6e" font-weight="bold">thresholds_nu.json</text>
  <text x="160" y="305" text-anchor="middle" font-size="14" fill="#0c4a6e">門檻、量測區、時窗（唯讀）</text>
  <!-- 理論層 -->
  <rect x="340" y="80" width="260" height="120" rx="10" fill="#ede9fe" stroke="#6d28d9"/>
  <text x="470" y="115" text-anchor="middle" font-size="18" fill="#4c1d95" font-weight="bold">理論：fdtd_theory / tmm / rcwa</text>
  <text x="470" y="145" text-anchor="middle" font-size="14" fill="#4c1d95">ky(h)、θ(y)、R_disc、R_TMM</text>
  <text x="470" y="170" text-anchor="middle" font-size="14" fill="#4c1d95">→ nu_predictions.json（先於模擬）</text>
  <!-- 求解層 -->
  <rect x="340" y="240" width="520" height="220" rx="12" fill="#fef3c7" stroke="#b45309"/>
  <text x="600" y="272" text-anchor="middle" font-size="18" fill="#78350f" font-weight="bold">fdtd3d_oblique.c（mesh=uniform|file，同一路徑）</text>
  <rect x="360" y="290" width="150" height="70" rx="8" fill="#fffbeb" stroke="#b45309"/>
  <text x="435" y="320" text-anchor="middle" font-size="14" fill="#78350f">係數陣列</text>
  <text x="435" y="342" text-anchor="middle" font-size="12" fill="#78350f">cE_y[j] cH_y[j] ε_t ε_n</text>
  <rect x="525" y="290" width="150" height="70" rx="8" fill="#fffbeb" stroke="#b45309"/>
  <text x="600" y="320" text-anchor="middle" font-size="14" fill="#78350f">3D Yee + CPML</text>
  <text x="600" y="342" text-anchor="middle" font-size="12" fill="#78350f">物理距離剖面</text>
  <rect x="690" y="290" width="150" height="70" rx="8" fill="#fffbeb" stroke="#b45309"/>
  <text x="765" y="320" text-anchor="middle" font-size="14" fill="#78350f">aux 線 + aux_ref 線</text>
  <text x="765" y="342" text-anchor="middle" font-size="12" fill="#78350f">dt/h 與 dt·iK̃ 分開</text>
  <rect x="360" y="375" width="480" height="65" rx="8" fill="#fffbeb" stroke="#b45309"/>
  <text x="600" y="402" text-anchor="middle" font-size="14" fill="#78350f">TF/SF（j0 局部間距）· DFT（j=Ny 修正）· 能量 W_mod · 非均勻散度</text>
  <text x="600" y="425" text-anchor="middle" font-size="12" fill="#78350f">輸出 meta.json、grid_used.json、dft_*.bin、log.csv</text>
  <!-- 參考層 -->
  <rect x="920" y="80" width="240" height="120" rx="10" fill="#dcfce7" stroke="#15803d"/>
  <text x="1040" y="115" text-anchor="middle" font-size="18" fill="#14532d" font-weight="bold">Meep（外部）</text>
  <text x="1040" y="145" text-anchor="middle" font-size="14" fill="#14532d">meep_ref_uniform.py</text>
  <text x="1040" y="170" text-anchor="middle" font-size="14" fill="#14532d">meep_ref_transform.py</text>
  <rect x="920" y="240" width="240" height="100" rx="10" fill="#dcfce7" stroke="#15803d"/>
  <text x="1040" y="275" text-anchor="middle" font-size="18" fill="#14532d" font-weight="bold">凍結參考 binary</text>
  <text x="1040" y="305" text-anchor="middle" font-size="14" fill="#14532d">ref/…uniform_ref.c</text>
  <!-- 分析層 -->
  <rect x="340" y="520" width="520" height="80" rx="10" fill="#f1f5f9" stroke="#334155"/>
  <text x="600" y="552" text-anchor="middle" font-size="18" fill="#0f172a" font-weight="bold">fdtd_io → analyze / compare → tests/nu_*.py</text>
  <text x="600" y="580" text-anchor="middle" font-size="14" fill="#0f172a">量測值 vs 理論值 vs 門檻 → results/nu_gate*.json + figures/nu/</text>
  <rect x="340" y="630" width="520" height="60" rx="10" fill="#fee2e2" stroke="#b91c1c"/>
  <text x="600" y="667" text-anchor="middle" font-size="18" fill="#7f1d1d" font-weight="bold">make_report.py → validation_report.md</text>
  <!-- 箭頭 -->
  <defs><marker id="a" markerWidth="10" markerHeight="10" refX="9" refY="5" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#334155"/></marker></defs>
  <line x1="280" y1="140" x2="340" y2="140" stroke="#334155" stroke-width="2" marker-end="url(#a)"/>
  <line x1="280" y1="170" x2="360" y2="290" stroke="#334155" stroke-width="2" marker-end="url(#a)"/>
  <line x1="280" y1="300" x2="340" y2="540" stroke="#334155" stroke-width="2" marker-end="url(#a)"/>
  <line x1="470" y1="200" x2="470" y2="520" stroke="#6d28d9" stroke-width="2" stroke-dasharray="6 4" marker-end="url(#a)"/>
  <line x1="600" y1="460" x2="600" y2="520" stroke="#334155" stroke-width="2" marker-end="url(#a)"/>
  <line x1="1040" y1="200" x2="820" y2="520" stroke="#15803d" stroke-width="2" marker-end="url(#a)"/>
  <line x1="1040" y1="340" x2="860" y2="545" stroke="#15803d" stroke-width="2" marker-end="url(#a)"/>
  <line x1="600" y1="600" x2="600" y2="630" stroke="#334155" stroke-width="2" marker-end="url(#a)"/>
</svg>
```

### 17. 介面設計（命令列參數與 Python 函式）

本專案沒有網路服務，「API」指 C 的命令列介面與 Python 模組的公開函式。

#### 17.1 C 新增/變更參數（`./fdtd3d_oblique key=value ...`）

| key | 預設 | 說明 | 相容性 |
|---|---|---|---|
| `mesh` | `uniform` | `uniform`：由 nl、Lx、Lz、npml、sf、tf 內部產生均勻節點；`file`：讀 `grid=` | 舊指令不變 |
| `grid` | — | grid JSON 路徑（mesh=file 必填） | 新 |
| `dtfac` | 1.0 | Δt = dtfac × Courant 公式值；只供穩定邊界測試 | 新 |
| `dt` | — | 直接指定 Δt（優先於 dtfac），給第 3 關 (B) 對齊 Meep | 新 |
| `init` | `0` | 新增 `r`：隨機初始場（`seed=`），給能量守恆與穩定測試 | 擴充 |
| `auxref` | 0 | 1：啟用 aux_ref 線並輸出 dft_auxref.bin | 新 |
| `dump_at` | — | 逗號清單，在這些步數輸出全場 | 新 |
| `inc` | `a` | 階段 B 新增 `p`（解析取樣）、`j`（電流片）、`m`（Bloch 模態，選做） | 擴充 |
| `divop` | `nu` | `nu`：非均勻散度；`uni`：錯誤的均勻算子（只供對照組） | 新 |

- `eps2`、`y1`、`ifmode` 在 mesh=file 時不得同時給定（材料改由 grid 的 eps_y 決定）；同時給定則中止。
- `ifmode=s`（階梯）保留，用於「介面不對齊」對照組。

#### 17.2 Python 公開函式（穩定介面）

| 模組.函式 | 簽名 | 回傳 |
|---|---|---|
| `grid_gen.make_interface_grid` | `(layers, ppw, r_max=1.1, zones, x_spec, z_spec)` | dict（grid JSON） |
| `grid_gen.make_mapped_grid` | `(f, du, u_range, ...)` | dict |
| `grid_gen.check` | `(grid) → list[str]` | 違規清單，空清單才可用 |
| `tmm.tmm_continuous` | `(layers, n_in, n_out, kt, k0, pol)` | (r, t, R, T) |
| `tmm.tmm_discrete` | `(grid, dt, Kx, Kz, pol) → (R, T, fields)` | 與 C 相同離散的精確解，附節點場 |
| `fdtd_theory.ky_local` | `(h, eps, dt, Kt)` | ky |
| `analyze.ky_three_point` | `(F, h) → (ky, imag_residual)` | 最小平方擬合 cos(ky h) |
| `analyze.forward_phase` | `(E_t, H_t, grid, Z_tilde) → phase` | 局部前行波相位 |
| `analyze.div_nu` / `analyze.div_uniform_wrong` | `(fields, grid)` | 散度陣列 |
| `analyze.sy_weighted` | `(dft, grid)` | 各 y 平面 S_y |
| `analyze.floquet` | `(F, xnodes, znodes, weights, kx, kz, orders)` | 各階相量 |
| `fdtd_io.run` / `fdtd_io.load` | 既有簽名 + `grid=` | 同既有 + `coords` |

### 18. 錯誤處理、回退策略與可觀測性

| 情境 | 偵測 | 行為 | 訊息範例 |
|---|---|---|---|
| grid 檔不存在或格式錯 | C 讀檔 | exit 2 | `ERROR: grid file: missing key "y.h"` |
| grid 違反均勻區規則 | C 讀檔後重檢（與 grid_gen.check 同規則） | exit 2 | `ERROR: j0±5 not uniform (max rel diff 3.2e-02)` |
| x/z 非均勻卻用 aux | setup | exit 2 | §7.6 訊息 |
| 無傳播 ky（局部 h 太大） | ky_local | exit 2，並印出是哪一段 h | `ERROR: no propagating ky in cells 200..215 (h=0.1)` |
| Δt 超過 Courant | setup | 若 dtfac ≤ 1 則 exit 2；dtfac > 1 只警告（穩定測試用） | `WARNING: dt exceeds Courant bound by 2.0%` |
| 數值發散 | 每步 isfinite | exit 2 + 記錄步數（穩定測試把它當預期結果） | `instability detected at step 187` |
| Meep 不存在 | meep 腳本 | exit 3 → BLOCKED | `BLOCKED: MEEP_PYTHON not found` |
| 門檻檔被修改 | run_nu_regression 比對 hash | FAIL | `FAIL: thresholds_nu.json hash changed without decision record` |

可觀測性：
- log.csv 每步記錄 SF 最大場、兩種能量、散度。
- meta.json 記錄全部參數、網格摘要、Δt 代價、執行時間、執行緒數。
- 回歸摘要寫入 `results/nu_regression.log`。
- 每個結果 JSON 附 `git_rev` 或原始檔 sha256，以及執行指令。

### 19. 通知與背景執行

- 長時間模擬（第 1-7 關 20000 步、第 2 關收斂系列、第 3 關 Meep 320 解析度）由 `run_nu_regression.py` 依序執行。每一步寫一行摘要，含 PASS/FAIL/BLOCKED 與耗時。
- 預設第一個 FAIL 就停止。`--continue` 只供診斷，結果標「不構成放行」。
- 中斷處理：Ctrl-C 後未完成的 run 目錄缺 meta.json，下次視為無效並重跑；已完成的 run 保留。
- 重試：同一 run 不自動重試（決定性程式重跑結果相同）。FAIL 一律走「假設 → 最小實驗」流程。
- 完成通知：回歸結束時印出總表，exit code 0/1/3（PASS/FAIL/BLOCKED）。

### 20. 驗證關卡（完整門檻與量測定義）

所有表格欄位：項目 | 理論值 | 量測值 | 門檻 | PASS/FAIL，另附網格資訊（Δ_min、Δ_max、r_max、Δt、N_pml）。所有數字附單位（λ0、T0、rad、dB 或無因次）。

#### 20.1 第 0 關（不跑物理模擬）

| # | 項目 | 設定 | 門檻 | 對照 / 備註 |
|---|---|---|---|---|
| 0-1a | 回歸 parity | 等間距 grid 檔（mesh=file）與 mesh=uniform，分別對凍結 ref binary，比較第 1、10、100、1000、5000 步的全場。s/p 各一；含 CPML、TF/SF、介質 | max\|ΔF\|/max\|F_ref\| < 1e-12（每分量每檢查步） | NaN 位置（無效半格）排除 |
| 0-1b | 舊驗證 | `tests/run_all_regression.py --fresh`（D7 生效後） | 全 PASS | Meep 缺失則 BLOCKED |
| 0-2 | 伴隨性 | 隨機 E、H，Level 1 漸變網格與 r=4 網格，PBC x/z + PEC y | \|⟨E,C_H H⟩_{W_E} − ⟨C_E E,H⟩_{W_H}\| / (‖E‖‖C_H H‖) < 1e-13 | 以錯誤權重（全 1）計算應 ≫ 1e-13，列 INFO |
| 0-3 | 能量守恆 | npml=0（PEC）、inc=0、init=r、1e5 步、Level 1 漸變網格 | 𝓔 的相對漂移 max\|𝓔^n − 𝓔^0\|/𝓔^0 < 1e-10 | 舊能量式同時記錄（INFO） |
| 0-4 | 穩定邊界 | Python power iteration：M = ε⁻¹W_E⁻¹C_E^T W_H μ⁻¹ C_E（實際以兩次 curl 作用實作），Rayleigh 商相對變化 < 1e-10 才停 → Δt_max = 2/sqrt(λ_max)。C 以 dt = 0.99Δt_max、1.02Δt_max 各跑 20000 步 | 0.99：全程有限且 𝓔 漂移 < 1e-10；1.02：20000 步內出現非有限值，或 max\|F\| 超過初值 1e6 倍；Courant 公式 Δt ≤ Δt_max | 報告 Δt_courant/Δt_max |
| 0-5 | 散度 | 第 1 關漸變網格（真空）實際注入 run，量測區：TF 內 j ∈ [j0+2, Ny−N_pml,hi−2] | max\|div E\|/(\|K̃\|·max\|E\|) < 1e-10，H 同 | 對照：`divop=uni` 結果 ≥ 1e-6（事先定義的「明顯非零」），證明敏感度 |

#### 20.2 第 1 關（階段 A，真空）

網格（寫入 `grids/`）：
- `L1_taper_r1.1_base{20,40,80}`：底格 λ0/b；x/z = λ0/20（細化系列中 x/z 固定）。
  y 配置：近端 PML 20 | SF 20 | j0 | 均勻 2λ0 | r = 1.1 縮到 λ0/(4b) | 保持 0.5λ0 | r = 1.1 放回 | 均勻 3λ0 | 遠端 PML 60。
- `L1_abrupt_r4_base20`：同上配置，但縮放改成單格突變（λ0/20 ↔ λ0/80）。
- `L1_uniform_ctrl`：全均勻 λ0/20，給 1-1 字面量對照。
- 反射量測（1-4）與 PML（1-6）一律用遠端 N_pml = 60 與 20 兩組，事先寫死：判定用 60，20 只列 INFO。

| # | 項目 | 量測方法（事先寫死） | 門檻 |
|---|---|---|---|
| 1-0 | 精確約化（新增） | 所有 y 平面 max\|F_main − Re[aux_ref·e^{iφ}]\|/E0 | < 1e-12 |
| 1-1 | 洩漏（D8） | max_SF \|F_main − aux_ref 預測散射場\|/E0，時窗 = 入射波前抵達遠端 PML 前 | < 1e-10；字面量 SF max 列 INFO；`L1_uniform_ctrl` 以字面量判定 < 1e-10 |
| 1-2 | x–z 平面性 | 每個 y 平面，對 unwrap 相位做 kx x + kz z + c 最小平方擬合，全部分量 | 殘差 RMS < 1e-8 rad |
| 1-3a | 子區 ky | 三點恆等式；量測點為三點都落在同一均勻子區、且不含 j0±1 與 PML 的全部節點；E_t 分量 | \|ky_meas − ky(h)\|/ky(h) < 1e-6（每個子區） |
| 1-3b | 漸變區 θ(y) | 局部前行波分解相位，h_loc = 兩相鄰對偶節點的距離 d_j（D3） | \|θ_meas − θ_pred\| ≤ 0.1\|θ_pred − θ_cont\|；\|θ_meas − θ_cont\| 對 b = 20/40/80 的收斂階數 ∈ [1.8, 2.2] |
| 1-3c | 圖 | θ(y) vs 預測曲線；x–y 切面瞬時場（equal aspect、真實座標）疊上各子區預測等相位線 | 圖存在且由程式產生 |
| 1-4a | r=4 突變反射（D1） | 在 j0 與漸變起點之間的均勻區做前行/反射兩波擬合 → \|R\| | \|R_meas − R_disc\|/\|R_disc\| < 1e-3；Fresnel 估計與比值列 INFO（預估 3.0） |
| 1-4b | r=1.1 漸變反射（D2） | 同上 | \|R_meas − R_disc\|/\|R_disc\| < 1e-3（b = 20、40、80 各自）；且 \|R\|_{20} > \|R\|_{40} > \|R\|_{80}；−60 dB 列 INFO |
| 1-4c | 反射 vs r_max 與解析度 | r ∈ {1.05, 1.1, 1.2, 1.5, 4} × b ∈ {20, 40, 80} | 報告表 + 圖（INFO） |
| 1-5 | S_y 守恆 | 時間平均 S_y，H 以真實座標距離加權內插到 E 主節點；TF 區所有 y 平面 | (max − min)/mean < 1e-4；½ 平均列 INFO |
| 1-6 | PML | 在漸變區後方均勻區擬合反向波 → 遠端 PML 反射，N_pml ∈ {10, 20, 30, 60} | N_pml=20 時與舊均勻版相差 ≤ 3 dB，且 < −40 dB（舊門檻）；其他 N_pml 列 INFO |
| 1-7 | 穩定性 | ≥ 20000 步，off_t 關源 | 從能量峰值單調衰減到 < 1e-8×峰值；第二半段最大值 ≤ 起點值（沿用舊版已接受的解讀） |

#### 20.3 第 2 關（階段 A，薄膜）

- 結構 F1：真空 | n = 2.0、0.08λ0 薄膜 | n = 1.46 基板一路延伸進遠端 PML（σ_max 用 η0/1.46）。結構 F5：H L H L H（A4）+ 同基板。
- 偏振 s、p；角度 A2。預設網格（D6）：y 每材料 PPW = 40，x/z = λ0/20。
- R、T 以離散守恆通量（E 在主節點、H 在對偶節點，距離加權）在反射側與透射側量測，並經 normalization（真空 run）正規化。

| # | 項目 | 門檻 |
|---|---|---|
| 2-1 | 預設網格 \|R − R_TMM\|、\|T − T_TMM\|（F1、F5 × s、p） | < 1e-3 |
| 2-2 | \|R + T − 1\| | < 1e-5 |
| 2-3 | 介面對齊非均勻網格，R 對 Δy 的收斂階數（PPW = 20、40、80、160，x/z 固定；對 R_∞ = Richardson 外插值計算） | ∈ [1.8, 2.2]（以 prompt「≈ 2」事先定義） |
| 2-4 | 對照：介面不對齊的均勻網格（ifmode=s，介面偏移 0.3Δ） | 階數 ∈ [0.8, 1.2]（「≈ 1」）；log-log 圖 |
| 2-5 | 誤差預算：只細化 Δy、只細化 Δx 與 Δz | 報告 y 離散誤差與 K̃t vs kt 橫向誤差兩項，並與 tmm_discrete 預測對比（INFO 表 + 圖） |
| 2-6 | 效率：達到相同 TMM 誤差時的格點數與執行時間 | 報告表（INFO），同執行緒數 |
| 2-7（新增，診斷） | FDTD vs R_disc（同網格） | 列 INFO；若 2-1 FAIL，用它區分程式錯誤與離散化本質 |

#### 20.4 第 3 關（Meep）

(A) 收斂值比對：
- 薄 3D Bloch cell（x、z 各 4 個格點）；k_point = (kx/2π, 0, kz/2π)。
- 先以解析度 80 比對薄 cell 與完整 Lx × Lz cell，事先定義一致標準為 \|ΔR\| < 1e-6。
- 解析度 40、80、160、320，GaussianSource(f0 = 1, fwidth = 0.1) + normalization run，DFT 只取 f0；Richardson 外插。

| # | 項目 | 門檻 |
|---|---|---|
| 3A-1 | 薄 cell vs 完整 cell（解析度 80） | \|ΔR\| < 1e-6 |
| 3A-2 | Meep 外插值 vs TMM | < 1e-4 |
| 3A-3 | 本程式外插值（PPW 20/40/80/160、x/z 同步細化）vs Meep 外插值 | < 1e-4（事先定義「指向同一外插值」） |
| 3A-4 | 本程式預設網格 vs Meep 外插值 | < 1e-3 |

(B) 變換光學：
- y = f(u)，Δu = λ0/20；s(u) = 1 − (3/4)·B(u)，B 為 C∞ bump（在薄膜附近 = 1，兩端 PML 區 = 0）。
- material_function 回傳 mp.Medium(epsilon_diag = (εs, ε/s, εs), mu_diag = (s, 1/s, s))，eps_averaging = False；Meep Courant = Δt_own/Δu。
- 位置探針測試：先確認 Meep 各分量的 Yee 位置與本程式的 u 節點對應正確。

| # | 項目 | 門檻 |
|---|---|---|
| 3B-1 | 真空：u 映射回 y 後，DFT 切向場（Ex、Ez、Hx、Hz）相對 L2 差 | < 1e-3 |
| 3B-2 | 薄膜：同上 | < 1e-3 |
| 3B-3 | 殘差歸因：比較 s(u_j) 與 d_j/Δu、s(u_{j+½}) 與 h/Δu 的差，預測殘差量級；另跑離散度規變體（整數節點 s = d_j/Δu，半格 s = h/Δu） | INFO：預期殘差大幅下降；量化並說明 |
| 3B-4 | 兩邊穩定 | 兩邊都跑完且為有限值 |

#### 20.5 階段 B（x、z 非均勻，完整規格）

| # | 項目 | 門檻 |
|---|---|---|
| B0-1 | x/z 非均勻 + inc=a 或 auxref=1 → 中止 | exit 2 + 指定訊息（測試逐字比對） |
| B0-2 | 第 0 關 0-2、0-3、0-4、0-5 在 x、z 非均勻網格重跑 | 同第 0 關門檻 |
| B1-1 | 注入 (i)：SF 洩漏對 Δ（x/z 同步細化 λ0/20、40、80） | 收斂階數 ∈ [1.8, 2.2]；絕對值列 INFO |
| B1-2 | 注入 (ii)：電流片 + normalization；R、T 與 (i) 的差 | 差值收斂階數 ∈ [1.8, 2.2] |
| B1-3（選做） | 注入 (iii) 模態注入洩漏 | < 1e-10 |
| B2-1 | x–z 相位殘差 RMS | 收斂階數 ∈ [1.8, 2.2] |
| B2-2 | Floquet 純度：以對偶間距為求積權重分解，非入射階功率比例 | 報告值與收斂階數（INFO） |
| B3-1 | 光柵：x 向 lamellar，占空比 0.5，n = 2，厚 0.3λ0，x 節點對齊邊緣；conical 入射（A2 角度） | 各傳播階效率 \|η_FDTD − η_RCWA\| < 1e-3 |
| B3-2 | 能量 | \|ΣR_p + ΣT_p − 1\| < 1e-5（「= 1」依第 2 關精度事先定義） |
| B3-3 | RCWA 階數收斂：N = 11、21、41、81 | 最後兩階的效率差 < 1e-5，作為「足夠階數」的事先定義 |

- 階段 B 預設網格：y 同 D6；x 在脊內 λ0/(2·40)、槽內 λ0/40，漸變 r_max = 1.1；z = λ0/40。
- 若 B3-1 在預設網格 FAIL，照工作規則處理，不得另選網格。

### 21. 盤點表（Stage 0 驗證並補齊；以下 C 部分已依現有原始碼預填）

| 檔案:行號 | 假設內容 | 修改方式 |
|---|---|---|
| fdtd3d_oblique.c:63 | 全域 D（=1/nl） | 只保留為 Δ_ref（uniform 產生器、meta 相容） |
| fdtd3d_oblique.c:162–168 | ky_disc 用 D 與 P.S | ky_local(h, n)；ω̃ 以 Δt 計算 |
| fdtd3d_oblique.c:182–184 | amplitudes 的 K̃ 全用 D | K̃x 用 Δx、K̃z 用 Δz、K̃y 用 Δ_a |
| fdtd3d_oblique.c:215 | inc_amp 的 y = (j+OFF)·D | 取 aux 節點座標 |
| fdtd3d_oblique.c:228–236 | aux 長度與 ikxD、ikzD（含 D） | 長度規則保留（索引因果）；改存 iK̃x、iK̃z |
| fdtd3d_oblique.c:243–248, 259–264 | c = dt/D 同乘 y 差分與相量項 | 拆成 cy_aux[u] 與 Δt·iK̃ |
| fdtd3d_oblique.c:251–252, 266–267 | aux 1D TF/SF 係數 c | 用 j_a 處 Δt/Δ_a |
| fdtd3d_oblique.c:272–279 | D = 1/nl、dt = S·D、Nx = Lx·nl | Courant 公式（A3）；Nx 由 grid |
| fdtd3d_oblique.c:314 | 介面算術平均 | h 加權平均 |
| fdtd3d_oblique.c:317–320 | ceE* = dt/(εD)、ch = dt/D | §7.1 係數陣列 |
| fdtd3d_oblique.c:330–341 | CPML d = npml·D、y = (j+½h)·D、σ_max 含 D | 物理座標、Δ_pml |
| fdtd3d_oblique.c:360 | 相位表 x = (i+½)D | x、z 節點陣列 |
| fdtd3d_oblique.c:371, 392–398 | analytic_val、dirichlet_y 用 D | 節點陣列 |
| fdtd3d_oblique.c:425–442, 464–494 | update_H/E 單一係數 | 分項係數 |
| fdtd3d_oblique.c:464, 475, 512–513 | 能量權重 D³ | W_E、W_H；新增 W_mod |
| fdtd3d_oblique.c:541–550 | TF/SF 係數 ch、ceEx[j0] | cH_y[j0−1]、cE_y[j0] |
| fdtd3d_oblique.c:597–613 | dft_accumulate type 1/2 對 Ey/Hx/Hz 取 j = Ny | 半格分量 j = Ny 寫 NaN；meta 記有效範圍 |
| fdtd3d_oblique.c:648–667 | div_max 除以 D | 非均勻算子（divop） |
| fdtd3d_oblique.c:720 | meta "Delta" | 追加 mesh 欄位 + grid_used.json |
| fdtd3d_oblique.c:774–780 | init=b 的 y = j·D | 節點陣列 |
| Python（analyze、compare、fdtd_io、fdtd_theory、meep_ref、oblique、tests/*，約 141 處） | Delta、idx·D、resolution=nl | Stage 0 逐行列出；改為讀 grid_used.json；meep_ref.py 保留給舊關卡 |

### 22. 報告與圖表呈現（UI 設計）

- **風格定調**：學術工具感。清楚、可列印、沒有裝飾。
- **色彩策略**：
  - 主色深藍 #1d4ed8：本程式量測；
  - 輔色灰 #475569：理論與預測曲線（虛線）；
  - 輔色綠 #15803d：Meep；
  - 強調色紅 #b91c1c：只用在 FAIL 標記與門檻線。
  - 同一張圖最多 4 條曲線；色盲友善時以線型區分。
- **比例規則**：x–y 切面與任何場圖一律 `set_aspect('equal')`，以真實座標繪製，縮放時維持原始寬高比。log-log 圖的兩軸採相同十進位跨度。
- **動畫**：`matplotlib.animation.PillowWriter`，GIF，≤ 200 幀。
- **分步導覽**：報告以關卡分章（Part II：第 0 → 1 → 2 → 3 → 階段 B），同頁分段，不用 tabs。每章順序：驗證表 → 圖 → INFO 表。

#### UI 元件清單（報告構件）
| 元件 | 用途 | 關鍵屬性 | 產生者 |
|---|---|---|---|
| 總表 | 各關結論 | 關卡、PASS/FAIL/BLOCKED | make_report |
| 驗證表 | 判定 | 項目、理論、量測、門檻、判定、網格資訊 | tests/nu_* |
| INFO 表 | 非門檻數據 | 同上但無判定 | tests/nu_* |
| 未通過紀錄 | 失敗處理 | 關卡、現象、假設、實驗、修正、結果 | 手寫 + JSON |
| θ(y) 圖 | 1-3b | 量測點、預測線、連續線 | nu_gate1 |
| 瞬時場圖 | 1-3c | equal aspect、等相位線 | nu_gate1 |
| 收斂圖 | 1-3b、2-3、2-4、3A、B | log-log、參考斜率線 | nu_gate* |
| 反射圖 | 1-4c | R vs r_max、解析度 | nu_gate1 |

#### UI 事件回報（log 事件）
| 事件 | 觸發 | 欄位 | 用途 |
|---|---|---|---|
| `run_start` / `run_end` | 每次 C 或 Meep 執行 | key、參數、grid_hash、耗時、exit code | 重現與快取 |
| `gate_item` | 每個判定項 | id、理論、量測、門檻、status | 驗證表 |
| `gate_summary` | 每關結束 | 關卡、status | 總表 |
| `threshold_check` | 回歸開始 | thresholds hash | 防竄改 |

#### UI ↔ 介面對照
| 報告區塊 | 呼叫 | 資料來源 |
|---|---|---|
| 第 0 關表 | tests/nu_gate0_*.py | fields_n*.bin、log.csv、Python 算子 |
| 第 1 關表與圖 | nu_gate1_vacuum.py → analyze.* | dft_*.bin、dft_auxref.bin、grid_used.json、nu_predictions.json |
| 第 2 關 | nu_gate2_film.py → tmm.* | 同上 + TMM |
| 第 3 關 | nu_gate3_meep.py → meep_ref_* → compare.* | Meep 輸出 |
| 附錄 | make_report 讀 docs/*.md | 盤點、推導、DIFF_NOTES |

#### 報告版面示意（SVG）
```svg
<svg viewBox="0 0 1440 960" xmlns="http://www.w3.org/2000/svg" font-family="Noto Sans TC, sans-serif">
  <rect width="1440" height="960" fill="#ffffff"/>
  <rect x="80" y="40" width="1280" height="70" fill="#f1f5f9" stroke="#334155"/>
  <text x="110" y="85" font-size="26" fill="#0f172a" font-weight="bold">驗證報告 Part II：非均勻網格　總結論：PASS / FAIL / BLOCKED</text>
  <rect x="80" y="130" width="620" height="220" fill="#eff6ff" stroke="#1d4ed8"/>
  <text x="100" y="165" font-size="20" fill="#1e3a8a" font-weight="bold">各關總表（唯一視覺重點）</text>
  <text x="100" y="200" font-size="16" fill="#1e3a8a">第 0 關 ● PASS　第 1 關 ● PASS　第 2 關 ● …</text>
  <text x="100" y="230" font-size="16" fill="#1e3a8a">第 3 關 ● …　階段 B ● …</text>
  <text x="100" y="270" font-size="14" fill="#475569">網格資訊：Δ_min、Δ_max、r_max、Δt、Δt 代價</text>
  <rect x="740" y="130" width="620" height="220" fill="#fef2f2" stroke="#b91c1c" stroke-dasharray="8 4"/>
  <text x="760" y="165" font-size="20" fill="#7f1d1d" font-weight="bold">未通過紀錄（只在有 FAIL 時顯示）</text>
  <text x="760" y="200" font-size="15" fill="#7f1d1d">現象 → 假設 → 最小實驗 → 修正 → 重跑結果</text>
  <rect x="80" y="380" width="1280" height="200" fill="#ffffff" stroke="#334155"/>
  <text x="100" y="415" font-size="20" fill="#0f172a" font-weight="bold">第 1 關驗證表：項目 | 理論值 | 量測值 | 門檻 | 判定</text>
  <line x1="100" y1="430" x2="1340" y2="430" stroke="#cbd5e1"/>
  <text x="100" y="460" font-size="15" fill="#334155">1-4a 突變 r=4 反射 | R_disc | R_meas | 相對差 &lt; 1e-3 | PASS</text>
  <text x="100" y="490" font-size="15" fill="#334155">1-3b θ(y) 偏移 | θ_pred | θ_meas | ≤ 0.1|θ_pred−θ_cont| | PASS</text>
  <rect x="80" y="600" width="620" height="300" fill="#f8fafc" stroke="#334155"/>
  <text x="100" y="635" font-size="18" fill="#0f172a" font-weight="bold">圖：θ(y) 量測 vs 預測</text>
  <polyline points="120,860 200,850 300,800 400,760 500,800 600,850 680,860" fill="none" stroke="#1d4ed8" stroke-width="3"/>
  <polyline points="120,858 200,848 300,798 400,762 500,798 600,848 680,858" fill="none" stroke="#475569" stroke-width="2" stroke-dasharray="6 4"/>
  <rect x="740" y="600" width="620" height="300" fill="#f8fafc" stroke="#334155"/>
  <text x="760" y="635" font-size="18" fill="#0f172a" font-weight="bold">圖：x–y 瞬時場（等比例）+ 等相位線</text>
  <rect x="800" y="660" width="500" height="220" fill="#dbeafe" stroke="#1d4ed8"/>
  <line x1="800" y1="880" x2="1000" y2="660" stroke="#475569" stroke-dasharray="6 4"/>
  <line x1="950" y1="880" x2="1150" y2="660" stroke="#475569" stroke-dasharray="6 4"/>
  <text x="80" y="940" font-size="14" fill="#64748b">附錄（收合）：盤點表、推導、DIFF_NOTES、環境與重現</text>
</svg>
```

#### UI 狀態保存與重新開始
- 模擬結果與判定結果都保存在 `runs/` 與 `results/`，中斷後重跑只補未完成的 run。
- 「重新開始」＝ `run_nu_regression.py --fresh`：清除 runs 快取，保留 grids、predictions、thresholds。
- 若要重做理論值，必須另起新 prediction id 並記錄原因。

### 23. 非功能需求

| 類別 | 要求 |
|---|---|
| 效能 | 均勻路徑的每步時間不得比 ref 慢超過 10%（Level 1 λ0/20 網格，同執行緒數）；第 0–2 關完整回歸在 8 執行緒工作站 ≤ 12 小時（估計，列 INFO） |
| 記憶體 | 係數陣列為 O(Ny + Nx + Nz)，場陣列佈局與 ID 巨集不變 |
| 可靠性 | 決定性：同參數兩次執行，場逐位相同（OpenMP 只影響能量等 reduction，差距 < 1e-14） |
| 可攜性 | C99 + libm；Windows MinGW 建置說明保留；Python 3.10+、numpy、matplotlib |
| 可維護性 | 均勻不另寫一條分支；係數集中在 setup_coeffs()；lint 禁止 idx*D |
| 安全性 | 不讀寫專案外路徑；grid 讀取器對陣列長度做上限檢查（≤ 1e7），避免惡意或錯誤檔案耗盡記憶體 |

### 24. 測試案例

| ID | 測試 | 輸入 | 預期 |
|---|---|---|---|
| T-G1 | grid 介面對齊 | 薄膜 0.08λ0、PPW 40 | 介面節點精確（\|y − y_if\| < 1e-15），薄膜 7 格 |
| T-G2 | grid 相鄰比 | r_max 1.1 | 所有 h_{j+1}/h_j ∈ [1/1.1, 1.1] |
| T-G3 | grid 均勻區 | 故意把 j0 放在漸變區 | check 回傳違規，C 拒絕執行 |
| T-L1 | lint | 新增含 `j*D` 的 Python | nu_lint_coords FAIL |
| T-C1 | 參數衝突 | mesh=file + eps2=2 | exit 2 |
| T-C2 | 階段 B 中止 | x 非均勻 + inc=a | exit 2，訊息逐字相符 |
| T-D1 | dft 修正 | xy 切面 Ey 在 j = Ny | NaN |
| T-TMM1 | tmm_discrete 均勻特例 | 舊 n=1.5 半空間 | 與舊「精確離散 Fresnel」相差 < 1e-12 |
| T-TMM2 | tmm_discrete vs 1D 時域 | 漸變網格 | 與 aux_ref 的 DFT 相差 < 1e-9 |
| T-TMM3 | tmm_continuous | 單介面 | 與 Fresnel 公式相差 < 1e-14 |
| T-A1 | 三點估計 | 合成 a e^{iky} + b e^{−iky} | ky 誤差 < 1e-13 |
| T-A2 | 前行波分解 | 均勻區合成場 | 與三點估計相差 < 1e-10 |
| T-A3 | 前行波分解預檢 | tmm_discrete 給出的漸變區精確場 | 滿足 D3 門檻；若不滿足，停下回報使用者（理論/估計器問題，不是程式錯誤） |
| T-R1 | RCWA | 均勻層（無光柵） | 與 TMM 相差 < 1e-12 |

### 25. Gherkin / BDD

```gherkin
Feature: 均勻路徑是非均勻路徑的特例
  Scenario: 等間距 grid 檔與凍結版逐步一致
    Given 以 grid_gen 產生 Δ=λ0/20 的等間距 grid 檔
    And 凍結參考 binary 已由 ref/fdtd3d_oblique_uniform_ref.c 建置
    When 兩者以相同參數各跑 1、10、100、1000、5000 步並輸出全場
    Then 每個分量在每個檢查步的 max|ΔF|/max|F_ref| < 1e-12

Feature: 相量 aux line 的適用範圍
  Scenario: x 非均勻時拒絕使用 aux line
    Given grid 檔的 x 間距不全相等
    When 以 inc=a 執行
    Then 程式以 exit code 2 中止
    And stderr 含 "phasor aux line requires uniform x and z"

Feature: 漸變區數值反射可預測
  Scenario: r=1.1 漸變
    Given L1_taper_r1.1_base20 網格與遠端 N_pml=60
    And nu_predictions.json 已含 R_disc
    When 執行真空注入並做兩波擬合
    Then |R_meas − R_disc|/|R_disc| < 1e-3
```

### 26. Edge / Abuse cases

| 情境 | 處理 |
|---|---|
| r_max < 1 或 PPW ≤ 0 | grid_gen 拒絕 |
| 介面距離小於一個最小格 | grid_gen 拒絕並說明需要的 PPW |
| 薄膜厚度不是間距整數倍 | 取 ceil(厚度/Δ_max) 格，均分厚度（例：0.08/7） |
| j0 或 ja 落在漸變區 | check 違規 → 拒絕 |
| 局部 h 使 ky 成為虛數（掠射或太粗） | C 中止並列出格段 |
| 使用者手改 grid 檔但沒改 hash | C 重算 hash，不符時中止 |
| 手改 thresholds_nu.json | 回歸比對 hash → FAIL |
| 事後挑選量測區 | 量測區定義寫在門檻檔、以索引規則表示；程式不接受手動範圍參數 |
| Meep 版本不同 | 記錄版本；結果差異超過門檻時走 FAIL 流程，不以版本為由放行 |
| 極端比例（r=4）時 CPML 相鄰區不均勻 | check 違規 |
| x/z 非均勻 + auxref=1 | 同 aux 中止 |
| 20000 步後能量在捨入底限小幅上升 | 依 1-7 的事先解讀判定，逐樣本版列 INFO |

### 27. 建議補充的功能（非本次範圍）

| 功能 | 價值 | 優先度 |
|---|---|---|
| grid 自動最佳化：給目標誤差，搜尋最小格點數 | 效率 | 中 |
| aux_ref 的頻域版本（直接用 tmm_discrete 取代時域 1D） | 更快的參考解 | 低 |
| GPU kernel（利用每軸 1D 係數佈局） | 速度 | 低 |
| 色散補償的 x/z 係數（讓 K̃t = kt） | 消除橫向誤差 | 中 |

### 28. 驗收條件（總結）

1. §20 所有判定項 PASS；BLOCKED 只允許出現在 Meep 環境缺失，且總結論因此為 BLOCKED。
2. 舊回歸（D7 後）全 PASS，0-1a parity < 1e-12。
3. 交付物齊全：grid_gen.py、fdtd3d_oblique.c + ref/nonuniform.diff + docs/DIFF_NOTES.md、analyze.py、tmm.py、meep_ref_uniform.py、meep_ref_transform.py、compare.py、validation_report.md（含盤點表、全部推導、驗證表與圖）；階段 B 另含 rcwa.py。
4. 報告每張表附網格資訊與單位；INFO 與判定分開；前一關未 PASS 的結果標「不構成放行」。
5. thresholds_nu.json 與 nu_predictions.json 的 hash 自 Stage 0/1 起未被無紀錄地修改。

### 29. 風險與未決事項

| 編號 | 風險 | 證據 / 預估 | 處理 |
|---|---|---|---|
| R1 | D3 的前行波分解在漸變區可能仍受反射污染 | 反射 2.7e-3 對應相位污染約 5e-3（相對），與效應 1.5e-3 同量級；前行波分解預計把污染壓低一個數量級以上 | T-A3 在跑 FDTD 前先用精確離散場預檢；不過則停下回報 |
| R2 | 2-3 收斂階數在 PPW 20–160 間可能不單調 | 初步 1D 估算中，PPW 80 的誤差大於 PPW 40（y 誤差與 x/z 橫向誤差相消） | 階數對 R_∞(Δx) 計算（x/z 固定），而不是對 TMM；事先寫死 |
| R3 | 1-4 的 1e-3 相對門檻 ≈ 絕對 3e-6，接近 N_pml=20 的 PML 回波 | 舊報告 −110.5 dB = 3e-6 | 判定用 N_pml = 60（事先寫死） |
| R4 | Meep 在高對比的 ε'、μ' 下需要更小的 Courant | s ∈ [1/4, 1] | 以 Δt_own/Δu 設定並確認穩定（3B-4） |
| R5 | 3B 位置對應錯誤造成假殘差 | Meep 以 cell 中心為原點 | 位置探針測試先行 |
| R6 | 階段 B 的 RCWA conical + Li 分解實作難度高 | — | T-R1 + 階數收斂 B3-3 |
| R7 | D8 的洩漏定義是規格解讀 | 見 §7.3 | 使用者可否決；否決時改回字面量，預期 FAIL |
| R8 | 第 2 關預設網格 PPW 40 的判定依賴 s 偏振 1D 預估 | 預估 8e-5 | p 偏振與 F5 在 Stage 1 用 tmm_discrete 預測；若預測 > 1e-3，先回報使用者再開跑 |

---

## 非技術規格文件

### 這份規格是寫給誰看的

寫給要判斷「這次升級可不可信」的人：指導教授、計畫審查人、合作的研究者。你懂光學與電磁學，但不需要會寫程式。
讀完之後，你應該知道：這次改了什麼、怎麼證明沒有改壞、每一關看什麼、什麼情況算失敗，以及最後會拿到哪些文件與圖。

### 這個工具能做什麼

原本的模擬器把空間切成一樣大小的小方格，來計算一道斜射的光。要算很薄的膜（例如只有十分之一波長厚的高折射率鍍膜），就必須把整個空間都切得很細，計算時間會大幅增加。
升級之後，只有在需要的地方（薄膜裡面、材料交界附近）才把格子切細，其他地方維持原本的大小。
升級也保證：如果你把格子設成一樣大，結果和舊版完全相同，差距小於一兆分之一。

### 你會怎麼使用它

1. 描述結構：告訴格子產生工具每一層的材料、厚度，以及每個波長至少要幾格。工具自動把格子從粗慢慢變細，每一步最多變化 10%，並讓材料交界剛好落在格線上。它會存成一份座標檔，之後所有計算與畫圖都讀這份檔案，不再自己用格子編號乘寬度去猜位置。
2. 執行模擬：給模擬器這份座標檔，照原本的方式指定入射角、偏振、計算長度。
3. 看結果：執行驗證流程，它會逐關檢查，最後產生一份驗證報告，裡面有表格與圖。
4. 想從頭再來時，下達「重新全部計算」的指令。舊的計算紀錄會清掉重算，但事先寫下的理論值與及格標準不會被動到。

### 你會看到哪些主要畫面

- **報告第一頁**：一張總表，列出每一關是「通過」、「未通過」或「無法執行」。如果有未通過，下方緊接著一張「未通過紀錄」：看到什麼現象、猜測原因、做了什麼小實驗、怎麼修、重跑結果。
- **每一關的驗證表**：每一行是一個檢查項目，依序是理論上應該是多少、實際量到多少、及格標準、結果。每張表都附上當時格子最小多大、最大多大、相鄰格子最多差幾倍。
- **圖**：
  - 光在不同格子大小區域的行進角度，實測點對上預測曲線；
  - 某一瞬間的光場切面，疊上預測的波峰線；
  - 誤差隨格子變細而下降的對數圖。
- **附錄**（預設收起來）：舊程式哪些地方假設了「格子一樣大」的清單、全部數學推導、程式修改說明。

### 畫面風格與色彩

- 風格：學術報告，乾淨、可列印。
- 深藍色：我們的模擬結果。灰色虛線：理論預測。綠色：對照軟體 Meep 的結果。紅色只用來標示「未通過」與及格線。
- 所有光場圖都照真實長寬比例畫，放大縮小時不會被拉扁或拉長。

### 畫面示意（SVG）

```svg
<svg viewBox="0 0 1440 960" xmlns="http://www.w3.org/2000/svg" font-family="Noto Sans TC, sans-serif">
  <rect width="1440" height="960" fill="#ffffff"/>
  <text x="720" y="60" text-anchor="middle" font-size="30" fill="#0f172a" font-weight="bold">一關一關檢查，前一關全過才進下一關</text>
  <rect x="60" y="120" width="250" height="300" rx="16" fill="#eff6ff" stroke="#1d4ed8" stroke-width="2"/>
  <text x="185" y="165" text-anchor="middle" font-size="22" fill="#1e3a8a" font-weight="bold">第 0 關</text>
  <text x="185" y="205" text-anchor="middle" font-size="17" fill="#1e3a8a">不跑真正的光</text>
  <text x="185" y="240" text-anchor="middle" font-size="16" fill="#334155">格子一樣大時</text>
  <text x="185" y="265" text-anchor="middle" font-size="16" fill="#334155">和舊版完全相同</text>
  <text x="185" y="300" text-anchor="middle" font-size="16" fill="#334155">能量不會憑空增減</text>
  <text x="185" y="335" text-anchor="middle" font-size="16" fill="#334155">找出最大安全時間步</text>
  <rect x="340" y="120" width="250" height="300" rx="16" fill="#eff6ff" stroke="#1d4ed8" stroke-width="2"/>
  <text x="465" y="165" text-anchor="middle" font-size="22" fill="#1e3a8a" font-weight="bold">第 1 關</text>
  <text x="465" y="205" text-anchor="middle" font-size="17" fill="#1e3a8a">真空中的光</text>
  <text x="465" y="240" text-anchor="middle" font-size="16" fill="#334155">光沒有漏到不該去的區域</text>
  <text x="465" y="265" text-anchor="middle" font-size="16" fill="#334155">角度偏差符合預測</text>
  <text x="465" y="300" text-anchor="middle" font-size="16" fill="#334155">格子變大小處的</text>
  <text x="465" y="325" text-anchor="middle" font-size="16" fill="#334155">微小反射可預測</text>
  <rect x="620" y="120" width="250" height="300" rx="16" fill="#eff6ff" stroke="#1d4ed8" stroke-width="2"/>
  <text x="745" y="165" text-anchor="middle" font-size="22" fill="#1e3a8a" font-weight="bold">第 2 關</text>
  <text x="745" y="205" text-anchor="middle" font-size="17" fill="#1e3a8a">薄膜與多層膜</text>
  <text x="745" y="240" text-anchor="middle" font-size="16" fill="#334155">反射與穿透</text>
  <text x="745" y="265" text-anchor="middle" font-size="16" fill="#334155">對上教科書精確解</text>
  <text x="745" y="300" text-anchor="middle" font-size="16" fill="#334155">反射 + 穿透 = 1</text>
  <text x="745" y="335" text-anchor="middle" font-size="16" fill="#334155">比均勻格子更省</text>
  <rect x="900" y="120" width="250" height="300" rx="16" fill="#ecfdf5" stroke="#15803d" stroke-width="2"/>
  <text x="1025" y="165" text-anchor="middle" font-size="22" fill="#14532d" font-weight="bold">第 3 關</text>
  <text x="1025" y="205" text-anchor="middle" font-size="17" fill="#14532d">和 Meep 對照</text>
  <text x="1025" y="240" text-anchor="middle" font-size="16" fill="#334155">兩套程式越算越細</text>
  <text x="1025" y="265" text-anchor="middle" font-size="16" fill="#334155">收斂到同一個答案</text>
  <text x="1025" y="300" text-anchor="middle" font-size="16" fill="#334155">「拉伸空間」的</text>
  <text x="1025" y="325" text-anchor="middle" font-size="16" fill="#334155">等效寫法光場一致</text>
  <rect x="1180" y="120" width="220" height="300" rx="16" fill="#f8fafc" stroke="#64748b" stroke-width="2" stroke-dasharray="8 5"/>
  <text x="1290" y="165" text-anchor="middle" font-size="22" fill="#334155" font-weight="bold">延伸階段</text>
  <text x="1290" y="205" text-anchor="middle" font-size="17" fill="#334155">左右方向也不等寬</text>
  <text x="1290" y="240" text-anchor="middle" font-size="16" fill="#334155">光柵的各繞射方向</text>
  <text x="1290" y="265" text-anchor="middle" font-size="16" fill="#334155">亮度對上精確解</text>
  <text x="720" y="490" text-anchor="middle" font-size="22" fill="#0f172a" font-weight="bold">格子只在需要的地方變細</text>
  <rect x="160" y="520" width="1120" height="200" fill="#f8fafc" stroke="#334155"/>
  <g stroke="#94a3b8">
    <line x1="200" y1="520" x2="200" y2="720"/><line x1="260" y1="520" x2="260" y2="720"/><line x1="320" y1="520" x2="320" y2="720"/>
    <line x1="380" y1="520" x2="380" y2="720"/><line x1="440" y1="520" x2="440" y2="720"/><line x1="490" y1="520" x2="490" y2="720"/>
    <line x1="535" y1="520" x2="535" y2="720"/><line x1="575" y1="520" x2="575" y2="720"/><line x1="610" y1="520" x2="610" y2="720"/>
    <line x1="640" y1="520" x2="640" y2="720"/><line x1="665" y1="520" x2="665" y2="720"/><line x1="690" y1="520" x2="690" y2="720"/>
    <line x1="715" y1="520" x2="715" y2="720"/><line x1="740" y1="520" x2="740" y2="720"/><line x1="765" y1="520" x2="765" y2="720"/>
    <line x1="790" y1="520" x2="790" y2="720"/><line x1="820" y1="520" x2="820" y2="720"/><line x1="855" y1="520" x2="855" y2="720"/>
    <line x1="895" y1="520" x2="895" y2="720"/><line x1="940" y1="520" x2="940" y2="720"/><line x1="990" y1="520" x2="990" y2="720"/>
    <line x1="1050" y1="520" x2="1050" y2="720"/><line x1="1110" y1="520" x2="1110" y2="720"/><line x1="1170" y1="520" x2="1170" y2="720"/>
  </g>
  <rect x="640" y="520" width="150" height="200" fill="#fde68a" fill-opacity="0.55"/>
  <text x="715" y="625" text-anchor="middle" font-size="18" fill="#78350f" font-weight="bold">薄膜</text>
  <text x="330" y="760" text-anchor="middle" font-size="17" fill="#334155">空氣：一般大小的格子</text>
  <text x="715" y="760" text-anchor="middle" font-size="17" fill="#334155">薄膜：細格子，交界落在格線上</text>
  <text x="1080" y="760" text-anchor="middle" font-size="17" fill="#334155">基板：依折射率決定格子大小</text>
  <text x="720" y="830" text-anchor="middle" font-size="18" fill="#b91c1c">相鄰格子大小最多只差 10%，避免在格子變化處產生多餘反射</text>
  <text x="720" y="880" text-anchor="middle" font-size="16" fill="#64748b">每張表都寫明當時的最小格、最大格與相鄰變化上限</text>
</svg>
```

### 操作流程

1. 先把舊版本原封不動保存一份，當作比較用的「標準答案」，並確認舊版所有檢查都通過。
2. 列出舊程式裡所有「假設格子一樣大」的地方，逐一寫明要怎麼改。
3. 在寫程式之前，先把全部理論推導與預測值寫好並封存，之後不能再改。
4. 改寫模擬器，再逐關檢查：第 0 關 → 第 1 關 → 第 2 關 → 第 3 關。前一關有任何一項沒過，就不進下一關。
5. 如果某項沒過：先寫下猜測的原因，設計一個最小的小實驗驗證，修好之後把所有檢查（含舊版的）從頭再跑一次。
6. 全部通過後，產出驗證報告與修改說明。
7. 主要目標完成後，才開始延伸階段（左右方向也用不等寬格子，並驗證光柵）。

### 你會看到的提示語

| 情況 | 提示語範例 |
|---|---|
| 成功 | 「第 1 關：14 項全部通過（最小格 λ0/80、最大格 λ0/20、相鄰變化 ≤ 1.1 倍）」 |
| 成功 | 「格子設為一樣大時，與舊版最大差距 3×10⁻¹⁵，通過（標準 10⁻¹²）」 |
| 等待中 | 「正在計算：薄膜，s 偏振，第 8,000 / 20,000 步（已用 6 分鐘）」 |
| 等待中 | 「對照軟體正在計算最細的一組（解析度 320），預計較久」 |
| 失敗 | 「第 1-4 項未通過：格子突變處的反射量測值與預測相差 2.3×10⁻³，標準 10⁻³。請看未通過紀錄。」 |
| 失敗 | 「座標檔錯誤：光的入射面附近 5 格內的格子不一樣大，請重新產生座標檔。」 |
| 無法執行 | 「找不到對照軟體 Meep，第 3 關標示為無法執行，總結論不能算通過。」 |
| 拒絕執行 | 「左右方向的格子不等寬時，不能使用目前的入射光產生方式，已停止。請改用延伸階段的入射方式。」 |

### 限制與注意事項

- 主要目標只讓「上下方向」（光穿過膜層的方向）的格子不等寬；左右兩個方向維持等寬，因為目前的入射光產生方式只有在左右等寬時才完全精確。
- 格子從粗變細的地方，一定會有極微小的反射。規格事先把這個反射量算好，實測只要跟預測一致就算正確，不要求它是零。
- 用來產生正確入射光的小區域（入射面附近、光源附近、吸收邊界）必須維持等寬格子。座標檔工具會自動檢查。
- 格子最細的地方決定了每一步能前進的時間。格子切得越細，需要的步數越多。報告會寫出這個代價，例如第 1 關約多 2.45 倍步數。
- 對照軟體 Meep 本身只能用一樣大的格子，所以對照方式是「兩邊越算越細，看是否收斂到同一個答案」，以及「把空間拉伸改寫成等效材料」。
- 所有及格標準在開始計算前就寫死封存。事後不能為了通過而放寬，也不能改理論值或挑量測位置。本規格裡少數幾項經過你確認的調整，都已寫明原因。

### 成功完成後會得到什麼

- 一個升級後的模擬器：可以只在需要處加細格子；格子設成一樣大時與舊版一致。
- 一個格子產生工具：給定膜層與精度，自動產生座標檔。
- 一份自寫的多層膜精確解工具，用來對答案。
- 與對照軟體 Meep 的兩種對照程式。
- 一份驗證報告：總表、每關驗證表與圖、全部推導、舊程式假設清單、修改說明。
- 延伸階段完成時，另有光柵的精確解工具與對應驗證。

### 常見問題與錯誤提示

| 問題 | 回答 |
|---|---|
| 為什麼不直接把整個空間切細？ | 可以，但計算量會暴增。第 2 關會實際比較：要達到同樣精度，兩種做法各需要多少格子、多少時間。 |
| 格子變化處的反射會不會毀了結果？ | 相鄰格子只差 10% 時，反射非常小，而且可以事先算準。第 1 關就是在驗證「算得準」。 |
| 舊版原本有一項沒過，現在怎麼處理？ | 依你的決定，改成與「在同樣格子下的精確答案」比較（目前差距小於一千萬分之二）。和教科書連續答案的差距仍然會記錄，只是不當作及格與否的依據。 |
| 看到「無法執行」怎麼辦？ | 通常是對照軟體沒安裝。安裝後重跑即可；這段期間總結論不能算通過。 |
| 看到「座標檔錯誤」怎麼辦？ | 代表入射面、光源或吸收邊界附近的格子不等寬。用格子產生工具重新產生，不要手動修改座標檔。 |

---

## Codex / Claude Code 分階段開發計畫

**共用不可變規則**（Stage 0 寫入 `AGENTS.md` 給 Codex、`CLAUDE.md` 給 Claude Code，之後每個 stage 都必須遵守）：

1. 規格來源：`SPEC_nonuniform.md`。門檻只在 `tests/thresholds_nu.json`，理論值只在 `results/nu_predictions.json`，兩者都不得無使用者決策地修改。
2. 均勻是非均勻的特例：不得出現 `if (uniform) {舊碼} else {新碼}` 的分叉更新式。
3. Python 只用 numpy + matplotlib（Meep 腳本例外）；動畫用 `matplotlib.animation.PillowWriter`，不用 imageio。座標一律讀 `grid_used.json`，禁止 idx*D。
4. 失敗時：寫假設 → 最小實驗 → 修正 → 重跑全部測試（含舊測試），並把紀錄寫入報告的「未通過紀錄」。不得放寬門檻、不得改理論值、不得挑量測區。
5. 每一關產出驗證表（項目 | 理論值 | 量測值 | 門檻 | PASS/FAIL）與圖，數字附單位與網格資訊（Δ_min、Δ_max、r_max）。
6. 前一關未全部 PASS，不得宣稱下一關通過；診斷用的執行標「不構成放行」。
7. 本專案不涉及生成式 AI 輸出，沒有 Streaming 需求。

### Stage 0：基線凍結、舊關卡 2 改判（D7）、盤點表

- 目標：凍結舊碼；讓舊回歸在 D7 規則下全 PASS；產出完整盤點表；建立門檻檔與長期規則檔。這一階段不改任何求解器程式碼。
- 前置條件：現有專案可 `make`；WSL 或 Linux；Meep 環境可用（否則第 3 關記 BLOCKED）。

Codex Instructions
```text
[建議貼用方式] 在專案根目錄開新 Codex 任務，整段貼上。先把「共用不可變規則」寫入 AGENTS.md。
[任務範圍] 基線凍結 + D7 改判 + 盤點表 + 門檻檔。禁止修改 fdtd3d_oblique.c。
[需修改/新增的檔案清單]
  新增：AGENTS.md、ref/fdtd3d_oblique_uniform_ref.c（逐位複製現有 fdtd3d_oblique.c）、docs/inventory_nonuniform.md、tests/thresholds_nu.json
  修改：Makefile（新增 target ref → 建置 fdtd3d_oblique_ref，旗標與主程式相同）、tests/level2_fresnel.py（D7）、SPEC.md（在關卡 2 加「2026-09-25 使用者決策 D7」修訂註記）
[具體步驟]
  1. 若不是 git repo，先 git init 並 commit 現況作為 baseline tag `pre-nonuniform`。
  2. 複製舊碼到 ref/，Makefile 加 ref target，確認兩個 binary 對同參數輸出逐位相同。
  3. D7：level2_fresnel.py 的 1% 判定改對「精確離散 Fresnel」（fdtd_theory 已有），連續 Fresnel 相對誤差改列 INFO。不動其他測試。
  4. 跑 python3 tests/run_all_regression.py --fresh，必須全 PASS（Meep 缺失時為 BLOCKED，要回報）。
  5. 盤點：逐檔（C 與所有 .py，含 tests/）列出「檔案:行號 | 假設內容 | 修改方式」。C 部分以 SPEC_nonuniform.md §21 為起點逐行核對；Python 用 grep（Delta、*D、D*、/D、resolution、nl）逐筆分類：必改 / 不需改（說明理由）。
  6. 依 SPEC §20 寫 tests/thresholds_nu.json：每個項目 id、門檻值、量測區規則、時窗規則、N_pml；記錄 sha256。
[輸出格式要求] 盤點表為 Markdown 表格，每列一個位置；最後附統計（必改 N 處、不需改 M 處）。
[測試要求] 舊回歸 --fresh 全 PASS；ref binary 與現 binary 對 3 組參數輸出 diff 為空。
[驗收標準 DoD] 盤點表涵蓋 SPEC §21 所列項目與 prompt 指定清單（D、dt/D、ikxD/ikzD、analytic_val、dirichlet_y、inc_amp、aux、TF/SF、CPML、Courant、DFT metadata、全部 Python 後處理、meep_ref.py）；thresholds_nu.json 已 commit；validation_report.md Part I 關卡 2 顯示 PASS。
```

Claude Code Instructions
```text
[建議貼用方式] 在 Claude Code 專案根目錄貼上。先把「共用不可變規則」寫入 CLAUDE.md（專案層級），再開始。
[任務範圍] 同 Codex 版：基線凍結、D7 改判、盤點、門檻檔；不得修改 fdtd3d_oblique.c。
[需修改/新增的檔案清單] CLAUDE.md、ref/fdtd3d_oblique_uniform_ref.c、Makefile、tests/level2_fresnel.py、SPEC.md（修訂註記）、docs/inventory_nonuniform.md、tests/thresholds_nu.json
[具體步驟]
  1. 用 Grep 找出 C 與 Python 中所有 D、Delta、dt/D、ikxD、ikzD、resolution、nl 的使用處；用 Read 核對上下文後填入盤點表。
  2. 凍結舊碼、Makefile 加 ref target；以 Bash 比對兩個 binary 的輸出逐位相同。
  3. 套用 D7（只改 level2 判定基準與 INFO 列），重跑舊回歸 --fresh（長時間執行用 run_in_background，完成後讀結果）。
  4. 寫 thresholds_nu.json（SPEC §20 全部項目），記錄 sha256 到 CLAUDE.md 的「凍結雜湊」段落。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版；另外 CLAUDE.md 含共用規則與凍結雜湊。
```

- 風險與回滾：D7 只改判定、不改數值。若舊回歸不是全 PASS，先回報，不進 Stage 1。回滾 = `git checkout pre-nonuniform`。

### Stage 1：推導與理論值（寫程式前）

- 目標：完成 SPEC §12 的 5 項推導；實作 tmm.py（連續 + 離散約化）與 fdtd_theory 的非均勻函式；在任何模擬前產生 nu_predictions.json。
- 前置條件：Stage 0 DoD 達成。

Codex Instructions
```text
[建議貼用方式] 新任務貼上；確認 AGENTS.md 已含共用規則。
[任務範圍] 推導文件 + 理論工具 + 預測值。不改 C。
[需修改/新增的檔案清單] docs/derivation_nonuniform.md、tmm.py、fdtd_theory.py（新增函式，不改既有簽名）、tests/nu_predict.py、results/nu_predictions.json、tests/test_tmm.py
[具體步驟]
  1. 寫 5 項推導（SPEC §12 表格每一列的「必須包含的結論」逐條出現）；伴隨性以求和分部推導，能量守恆寫出正定條件。
  2. tmm.tmm_continuous(layers, n_in, n_out, kt, k0, pol)：s/p 精確解。
  3. tmm.tmm_discrete(grid, dt, Kx, Kz, pol)：以 grid 的 h、d、ε_t（介面 h 加權）、ε_n 建 1D 離散約化方程，兩端以均勻區精確離散出射模態封閉，自寫 Thomas 三對角解法（numpy only），回傳 R、T 與節點場。
  4. fdtd_theory：ky_local(h, eps, dt, Kt)、theta_local、fresnel_numeric_estimate（INFO 用）。
  5. nu_predict.py：依 SPEC §20 列出的網格（先用 grid_gen 尚未存在的話，以函式內建最小產生器建立相同 h 陣列，Stage 2 以後改讀 grids/）計算全部理論值：子區 ky、θ_pred、θ_cont、R_disc（1-4a、1-4b，b=20/40/80）、R_TMM（F1、F5、s/p）、第 2 關預設網格 tmm_discrete 預估。
  6. 若任何預估顯示門檻必然 FAIL（例如 2-1 的 p 偏振預估 > 1e-3），停下並回報使用者，不要自行換網格。
[輸出格式要求] nu_predictions.json 每筆 {id, quantity, value, formula_ref, grid_hash, created_utc}。
[測試要求] T-TMM1（均勻特例 vs 舊精確離散 Fresnel < 1e-12）、T-TMM3、T-A1、T-A3（前行波分解預檢）。
[驗收標準 DoD] 推導 5 項齊全；測試 PASS；nu_predictions.json 已 commit 並記錄 sha256；預檢結果寫入報告草稿。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上；Claude Code 會讀 CLAUDE.md。建議先進 Plan mode 列出推導章節與函式簽名，確認後實作。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟]
  1. 先寫 derivation_nonuniform.md，每一節結尾用「結論：」條列可被測試的式子。
  2. 實作 tmm.py 與 fdtd_theory 新函式；每個函式寫 docstring 註明對應推導小節。
  3. 寫 nu_predict.py 產生預測；遇到必然 FAIL 的預估時，用 AskUserQuestion 回報，並附數字與可選方案。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：預測值一旦 commit 即唯讀。錯誤只能新增更正 id 並記錄原因。回滾 = revert 本 stage commit（尚未有 C 變更）。

### Stage 2：網格產生器與座標規則

- 目標：grid_gen.py 兩種模式 + 驗證器；產生 SPEC §20 所需全部 grids；Python 座標 lint。
- 前置條件：Stage 1 DoD。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] 只做網格與 lint，不改 C。
[需修改/新增的檔案清單] grid_gen.py、grids/*.json、tests/test_grid_gen.py、tests/nu_lint_coords.py
[具體步驟]
  1. make_interface_grid：輸入介面位置、各材料 n、PPW、r_max（預設 1.1）、zones（npml_lo/hi、sf、j0 位置、ja、uniform_min_cells=5）、x/z 規格。介面精確落在主節點；各段取 ceil(厚度/Δ_max) 格等分；幾何級數漸變不超過 r_max；h 陣列直接輸出，節點以 math.fsum 累加。
  2. make_mapped_grid(f, du, u_range)：第 3 關 (B) 用，y_j = f(u_j)，s(u) = 1 − 0.75·bump(u)。
  3. check(grid)：SPEC §13.1 全部規則，回傳違規清單。CLI：python3 grid_gen.py --check file.json。
  4. 產生：L1_taper_r1.1_base{20,40,80}、L1_abrupt_r4_base20、L1_uniform_ctrl、L1 反射掃描 r∈{1.05,1.1,1.2,1.5,4}、F1/F5 × PPW{20,40,80,160}、F1 不對齊均勻對照、3B 映射網格、等間距 parity 網格（h 全等於 1/20，與舊 nl=20 版面相同）。
  5. nu_lint_coords.py：掃描新增或修改的 .py，找出 idx*D、*Delta、Delta*、/nl 等樣式；允許清單只有 grid_gen 的均勻產生函式。
  6. 重跑 nu_predict.py，改讀 grids/ 的檔案；若 grid_hash 與 Stage 1 預測不同，新增 id（不覆寫）。
[輸出格式要求] grid JSON 依 SPEC §13.1。
[測試要求] T-G1、T-G2、T-G3、T-L1。
[驗收標準 DoD] 所有 grids 通過 check；lint 對現有新檔 PASS、對故意違規樣本 FAIL。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上；在 CLAUDE.md 追加「座標只從 grid_used.json 取得」。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–6；完成後以 matplotlib 畫每個 grid 的 h(y) 圖，存到 figures/nu/grids/，供人工檢視。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版 + h(y) 圖存在。
```

- 風險與回滾：grid 變動會讓預測 id 失效，必須新增 id。回滾 = 刪除 grids 並 revert。

### Stage 3：C 係數陣列化、讀 grid、回歸 parity、DFT 修正（第 0 關 0-1）

- 目標：C 改成每軸係數陣列，單一路徑支援 mesh=uniform|file；修正 dft_accumulate；通過 0-1a 與 0-1b。這一階段的 aux line 仍是均勻 y 版本，只做係數拆分（y 與相量項分開）。
- 前置條件：Stage 2 DoD。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。這是風險最高的階段，要小步 commit。
[任務範圍] fdtd3d_oblique.c 重構 + grid 讀取 + dft 修正 + 回歸比對。
[需修改/新增的檔案清單] fdtd3d_oblique.c、fdtd_io.py（grid 參數、grid_used.json、快取 key 含 grid_hash）、tests/nu_gate0_regression.py、docs/DIFF_NOTES.md（持續累積）
[具體步驟]
  1. 新增 setup_coeffs()：由 h_x、h_y、h_z 與 eps_y 建 SPEC §7.1 全部陣列；mesh=uniform 時在程式內產生 h = D（直接賦值，不相減）。
  2. update_H/E 改成分項係數（例：ex += cE_y[j]*dHz − cE_tz[j]*dHy）；CPML ψ 仍存差分，乘 cE_y 或 cH_y。
  3. 新增 JSON 子集讀取器（僅本格式；陣列長度上限 1e7；缺 key 則 exit 2），讀後重跑 check 規則與 hash。
  4. 參數：mesh、grid、dtfac、dt、dump_at、divop；mesh=file 與 eps2/y1 同給時 exit 2。
  5. Δt 依 A3 公式；meta.json 追加欄位；輸出 grid_used.json。
  6. dft_accumulate：半格分量在 j = Ny 寫 NaN；meta 記各分量有效 j 範圍；檢查舊 Python 是否讀過該位置（盤點表已列），有就修正並記錄。
  7. aux：c 拆成 cy = Δt/Δ（本 stage 仍均勻）與 Δt·iK̃；確認數值不變。
  8. nu_gate0_regression.py：以 dump_at=1,10,100,1000,5000 分別跑新 binary（mesh=uniform 與 mesh=file 等間距 grid）與 ref binary（ref 以 nsteps=N、dump=1 各跑一次），比較全場；s/p、真空/介質、含 CPML 與 TF/SF。
  9. 重跑舊回歸 --fresh。
[輸出格式要求] results/nu_gate0.json 的 0-1a、0-1b 項；DIFF_NOTES 每段說明「舊 → 新 → 理由 → 對應推導小節」。
[測試要求] 0-1a < 1e-12；0-1b 全 PASS；T-C1、T-D1；均勻路徑每步時間不超過 ref 的 1.10 倍。
[驗收標準 DoD] 上述全過；程式中沒有 uniform/nonuniform 分叉更新式（code review 檢查）。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上；先用 Plan mode 列出要改的函式與順序，每改完一個函式就跑一次 parity（小網格、100 步）。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–9。長時間回歸用 run_in_background。parity 失敗時先用二分法找出第一個出現差異的步數與分量，再修。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：浮點運算順序造成的差距若 > 1e-12，先檢查是否誤用非對稱公式，不得調門檻。回滾 = revert 到 Stage 2 commit；ref binary 永遠可用。

### Stage 4：第 0 關其餘項目（伴隨、能量、穩定邊界、散度）

- 目標：0-2 到 0-5 全 PASS。
- 前置條件：Stage 3 DoD。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] Python 離散算子 + C 的 init=r、W_mod、divop + 4 個測試。
[需修改/新增的檔案清單] analyze.py（curl_E、curl_H、weights、div_nu、div_uniform_wrong）、fdtd3d_oblique.c（init=r、seed、W_mod 欄位、divE/H_TF_nu 欄位）、tests/nu_gate0_adjoint.py、nu_gate0_energy.py、nu_gate0_stability.py、nu_gate0_divergence.py
[具體步驟]
  1. Python 以 grid_used.json 建 C_E、C_H（PBC x/z、PEC y），權重依 SPEC §12-2；隨機場檢查伴隨性（0-2），另以全 1 權重列 INFO。
  2. C：W_mod = ½ΣεW_E|E^n|² + ½ΣμW_H H^{n+½}·H^{n−½}，log.csv 只在最後追加欄位。
  3. 0-3：npml=0、inc=0、init=r、1e5 步、漸變網格；讀 W_mod 算相對漂移。
  4. 0-4：power iteration（Rayleigh 商相對變化 < 1e-10 才停）→ Δt_max；C 用 dt= 跑 0.99 與 1.02 各 20000 步；判定依 SPEC §20.1；報告 Δt_courant/Δt_max。
  5. 0-5：第 1 關漸變網格真空注入 run，div_nu 在量測區 < 1e-10（正規化依 SPEC）；divop=uni 對照 ≥ 1e-6。
[輸出格式要求] results/nu_gate0.json 追加 0-2..0-5；圖：能量漂移 vs 步數、power iteration 收斂。
[測試要求] 0-2..0-5 全 PASS；0-1 重跑仍 PASS。
[驗收標準 DoD] 第 0 關表全 PASS。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–5；0-4 的 1.02 run 預期會以 exit 2 結束，測試要把它判為「預期發散 → PASS」，並記錄發散步數。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：散度量測需要 Stage 5 的非均勻注入，所以 0-5 可以在 Stage 5 完成後補跑；第 0 關在 0-5 PASS 前不算通過。回滾 = revert 本 stage。

### Stage 5：非均勻 aux line、aux_ref、TF/SF、CPML → 第 1 關

- 目標：aux 與 aux_ref 依 SPEC §7.2；TF/SF 局部係數；CPML 物理距離；第 1 關 1-0 到 1-7 全 PASS。
- 前置條件：Stage 4 的 0-1..0-4 PASS。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] C 的 aux/aux_ref/TF-SF/CPML 非均勻化 + 第 1 關全部量測與圖。
[需修改/新增的檔案清單] fdtd3d_oblique.c、analyze.py（ky_three_point、forward_phase、sy_weighted、two_wave_fit、phase_planarity）、tests/nu_gate1_vacuum.py、docs/DIFF_NOTES.md
[具體步驟]
  1. aux 節點 = 主網格 [0, Ny] + 兩端均勻延伸；y 差分乘 Δt/h 或 Δt/d，相量項乘 Δt·iK̃；ja 處局部 1D TF/SF；inc_amp 用 ky(Δ_a) 與實際座標。
  2. aux_ref（auxref=1）：同節點、同 ε、同 CPML；在 j0 接收與主網格相同的修正；輸出 dft_auxref.bin。
  3. 主網格 TF/SF 係數用 cH_y[j0−1]、cE_y[j0]。
  4. CPML 依 SPEC §7.4（物理距離、Δ_pml、遠端 n_loc）。
  5. 先跑第 0 關全部（含補跑 0-5），再跑第 1 關；量測區、時窗、N_pml 一律讀 thresholds_nu.json。
  6. 圖：θ(y) vs 預測與連續曲線；x–y 瞬時場（真實座標、equal aspect）疊子區預測等相位線；反射 vs r_max/解析度；PML vs N_pml；能量衰減。可選 GIF 用 PillowWriter。
[輸出格式要求] results/nu_gate1.json；每列附 Δ_min、Δ_max、r_max、Δt、Δt 代價、N_pml。
[測試要求] 1-0..1-7 全 PASS；第 0 關與舊回歸重跑 PASS。
[驗收標準 DoD] 第 1 關表全 PASS；INFO 表含 Fresnel 估計比值、−60 dB 對照、字面量洩漏。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–6。先用 1-0（主網格 = aux_ref × 相位）當除錯主工具：不過時，逐分量、逐 y 找第一個偏離的節點。任何 FAIL 依 CLAUDE.md 的失敗流程處理，並寫入 results/nu_failures.json。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：R1（D3 估計器）、R3（PML 回波）。若 1-3b FAIL 且 T-A3 預檢已通過，就是程式問題；若預檢也不過，停下回報。回滾 = revert 到 Stage 4。

### Stage 6：第 2 關薄膜、5 層堆疊、收斂、誤差預算、效率

- 目標：2-1 到 2-7。
- 前置條件：第 0、1 關全 PASS。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] 介質 run、R/T 量測、收斂與效率分析。
[需修改/新增的檔案清單] analyze.py（守恆通量 R/T、Richardson、階數擬合）、tests/nu_gate2_film.py、figures/nu/gate2/*
[具體步驟]
  1. F1、F5 × s/p，預設網格（D6：y 每材料 PPW=40、x/z=λ0/20），含真空 normalization run；遠端 PML 匹配基板。
  2. R、T 以離散守恆通量量測；對 R_TMM、T_TMM 判定 2-1、2-2；2-7 對 R_disc 列 INFO。
  3. 2-3：PPW 20/40/80/160（x/z 固定），對 Richardson R_∞ 擬合階數。2-4：ifmode=s、介面偏移 0.3Δ 的均勻網格，同樣擬合；畫 log-log。
  4. 2-5：只細化 Δy / 只細化 Δx、Δz；以 tmm_discrete 預測兩個分量，列表對比。
  5. 2-6：以 2-1 的誤差為目標，由均勻網格收斂曲線內插所需 Δ，實跑，比較格點數與牆鐘時間（同執行緒數）。
[輸出格式要求] results/nu_gate2.json；收斂圖含參考斜率 1 與 2。
[測試要求] 2-1..2-4 PASS；前面各關重跑 PASS。
[驗收標準 DoD] 第 2 關表全 PASS；INFO 表齊全。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上；收斂系列很耗時，用 run_in_background 平行排程，但同一時間只跑一個 OpenMP 程序，避免計時失真。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–5。若 2-1 FAIL：先看 2-7。FDTD 與 R_disc 相符，代表離散化本質（R8），回報使用者；不相符代表程式錯誤，照失敗流程處理。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：R2（不單調）、R8（p 偏振預估）。回滾 = revert 本 stage 的分析程式；求解器不應在本 stage 修改，若要修改，視為回到 Stage 5。

### Stage 7：第 3 關 Meep（A）收斂值、（B）變換光學

- 目標：3A-1..3A-4、3B-1..3B-4。
- 前置條件：第 0–2 關全 PASS；MEEP_PYTHON 可用（否則 BLOCKED）。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。Codex sandbox 若無 Meep，寫完腳本後回報 BLOCKED，不得跳過比對邏輯。
[任務範圍] meep_ref_uniform.py、meep_ref_transform.py、compare.py 擴充、tests/nu_gate3_meep.py。
[需修改/新增的檔案清單] 上述四檔 + figures/nu/gate3/*
[具體步驟]
  1. (A) 薄 Bloch cell（x、z 各 4 格點），k_point = mp.Vector3(kx/2π, 0, kz/2π)；GaussianSource(f0=1, fwidth=0.1)；normalization run（真空）+ load_minus_flux_data；DFT/flux 只取 f0。
  2. 解析度 80 做薄 cell vs 完整 cell（3A-1）；40/80/160/320 → Richardson → 3A-2；本程式 PPW 20/40/80/160（x/z 同步細化）→ 外插 → 3A-3；預設網格 → 3A-4。
  3. (B) material_function 回傳 mp.Medium(epsilon_diag=(εs, ε/s, εs), mu_diag=(s, 1/s, s))，eps_averaging=False，Courant=Δt_own/Δu；先做位置探針（已知函數取樣）確認 Yee 位置與 u 節點對應；本程式用 make_mapped_grid 的 grid 且 dt=Δt_own。
  4. 比較 DFT 切向場（法向分量先除以 s；u 映射回 y），相對 L2；真空與薄膜。
  5. 3B-3：量化 s(u_j) − d_j/Δu、s(u_{j+½}) − h/Δu，並跑離散度規變體；列 INFO 並寫歸因說明。
[輸出格式要求] results/nu_gate3.json；Meep 版本、解析度、耗時。
[測試要求] 3A、3B 判定項全 PASS 或 BLOCKED（僅限環境缺失）。
[驗收標準 DoD] 第 3 關表完成；殘差歸因段落完成。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上；Meep 以 MEEP_PYTHON 指定的直譯器在背景執行（run_in_background），完成後讀輸出。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–5。Meep 在某解析度不穩定時，記錄並回報，不自行改 Courant 以外的設定。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：R4、R5。Meep 腳本與求解器獨立，回滾只影響本 stage 檔案。

### Stage 8：階段 B 基礎（x/z 係數陣列、中止偵測、注入 i/ii/iii）

- 目標：B0-1、B0-2、B1-1..B1-3、B2-1、B2-2。
- 前置條件：階段 A（第 0–3 關）全部 PASS（第 3 關若為 BLOCKED，須使用者明確同意才開始）。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] C 擴充 x/z 係數陣列與週期接縫、inc=p/j/m、中止偵測；Floquet 分解。
[需修改/新增的檔案清單] fdtd3d_oblique.c、grid_gen.py（x/z 非均勻）、analyze.py（floquet）、tests/nuB_gate0.py、tests/nuB_gate1.py、docs/DIFF_NOTES.md
[具體步驟]
  1. cE_x[i]、cH_x[i]、cE_z[k]、cH_z[k]；x/z 均勻時陣列值全相同 → 階段 A 全部測試必須不變（重跑）。
  2. 中止偵測（SPEC §7.6 訊息逐字）。
  3. inc=p：在 TF/SF 面以分量實際座標取樣解析平面波 × ramp。inc=j：電流片 + normalization run。inc=m（選做）：x–z 橫向離散 Bloch 本徵問題（numpy 以 power/inverse iteration 求最接近平面波的模態）→ 模態 + 1D aux。
  4. B0-2：第 0 關 0-2..0-5 在 x/z 非均勻網格重跑。
  5. B1-1、B1-2、B2-1：x/z 同步細化 λ0/20、40、80，擬合階數；B2-2 Floquet 以對偶間距為求積權重。
[輸出格式要求] results/nuB_gate0.json、nuB_gate1.json。
[測試要求] B0、B1、B2 判定項 PASS；階段 A 全部回歸 PASS；T-C2。
[驗收標準 DoD] 上述全過。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–5；步驟 1 完成後立刻跑完整階段 A 回歸，再開始 2–5。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：x/z 陣列化可能拖慢均勻路徑；非功能需求（≤ 10%）必須重測。回滾 = revert 到 Stage 7。

### Stage 9：階段 B 光柵與 RCWA

- 目標：B3-1..B3-3。
- 前置條件：Stage 8 DoD。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] rcwa.py（conical、Li 正確 Fourier 分解法則、numpy only）+ 光柵 FDTD + 比對。
[需修改/新增的檔案清單] rcwa.py、tests/test_rcwa.py、grids/grating_*.json、tests/nuB_gate3_grating.py
[具體步驟]
  1. RCWA：均勻層特例與 TMM 相差 < 1e-12（T-R1）；階數 N = 11、21、41、81 收斂（B3-3）。
  2. 光柵網格：x 節點對齊脊邊緣；脊內 λ0/(2·40)、槽內 λ0/40、r_max = 1.1；y 同 D6；z = λ0/40。
  3. FDTD 以 inc=p（或 inc=m）注入，Floquet 分解反射與透射各傳播階效率。
  4. B3-1、B3-2 判定。
[輸出格式要求] results/nuB_gate3.json；每階效率表。
[測試要求] B3 全 PASS。
[驗收標準 DoD] 階段 B 表全 PASS。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–4；RCWA 先以文獻中的 lamellar 基準例（若有）自我核對，並把來源記在 docs。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：R6。若 B3-1 在預設網格 FAIL，不得換網格，照失敗流程處理。rcwa.py 獨立，回滾不影響求解器。

### Stage 10（倒數第二）：完整回歸、邊界測試、決定性

- 目標：一次跑完舊回歸 + 全部 nu 測試 + 邊界案例；確認決定性與門檻未被竄改。
- 前置條件：Stage 0–9 各自 DoD。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] tests/run_nu_regression.py 與邊界測試，不新增功能。
[需修改/新增的檔案清單] tests/run_nu_regression.py、tests/test_edge_cases.py
[具體步驟]
  1. run_nu_regression.py：檢查 thresholds_nu.json、nu_predictions.json 的 sha256 → 舊 run_all_regression → nu_lint → nu_gate0..3 → nuB_*；第一個 FAIL 停止；--fresh、--continue 語意同舊版。
  2. 邊界測試：SPEC §26 全部情境（非法 grid、介面太近、j0 在漸變區、ky 虛數、hash 被改、x/z 非均勻 + aux）。
  3. 決定性：同參數跑兩次，場逐位相同；OpenMP 1 與 8 執行緒的場逐位相同。
  4. --fresh 完整跑一次，記錄總耗時。
[輸出格式要求] results/nu_regression.log；最後一行 REGRESSION: PASS|FAIL|BLOCKED。
[測試要求] 全部 PASS（或僅因 Meep 缺失而 BLOCKED）。
[驗收標準 DoD] --fresh 全綠；邊界測試全部依規格行為。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上；完整回歸用 run_in_background，完成通知後讀 log。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–4。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：只新增測試，回滾無副作用。

### Stage 11（最終）：文件化與交付

- 目標：validation_report.md（Part I 舊均勻版 D7 後 + Part II 非均勻，含盤點表、全部推導、驗證表與圖）、ref/nonuniform.diff、DIFF_NOTES、README 更新。
- 前置條件：Stage 10 REGRESSION: PASS（或僅 Meep BLOCKED，並在總結論如實標示）。

Codex Instructions
```text
[建議貼用方式] 新任務貼上。
[任務範圍] 只寫文件與報告產生器，不改求解器與門檻。
[需修改/新增的檔案清單] tests/make_report.py、validation_report.md、ref/nonuniform.diff、docs/DIFF_NOTES.md、README.md、SPEC_nonuniform.md（只追加「實作結果」連結）
[具體步驟]
  1. make_report.py 從 results/*.json 產生全部表格：頁首總表、未通過紀錄（有才顯示）、各關驗證表 + INFO 表 + 圖、附錄（盤點表、推導、DIFF_NOTES、環境與重現、規格變更紀錄 D1–D8）。
  2. git diff pre-nonuniform -- fdtd3d_oblique.c > ref/nonuniform.diff；DIFF_NOTES 每段對應 diff hunk。
  3. README：新增 mesh=file 用法、grid_gen 範例、新測試指令、Δt 代價說明。
  4. 以 grep 檢查報告中每個判定列都有單位與網格資訊。
[輸出格式要求] 報告數字全部由程式產生，不得手填。
[測試要求] make_report 重跑兩次輸出相同；報告內圖片連結全部存在。
[驗收標準 DoD] SPEC_nonuniform.md §28 驗收條件 1–5 全部成立。
```

Claude Code Instructions
```text
[建議貼用方式] 貼上。
[任務範圍] 同 Codex 版。
[需修改/新增的檔案清單] 同 Codex 版。
[具體步驟] 同 Codex 版 1–4；完成後把 CLAUDE.md 的「目前進度」更新為完成，並列出後續建議（SPEC §27）。
[輸出格式要求] 同 Codex 版。
[測試要求] 同 Codex 版。
[驗收標準 DoD] 同 Codex 版。
```

- 風險與回滾：文件錯誤不影響程式；以 make_report 重建。
