# fdtd3d_oblique.c：非均勻網格修改說明

基準：`ref/fdtd3d_oblique_uniform_ref.c`（git tag `pre-nonuniform`）。完整 diff：`ref/nonuniform.diff`（Stage 11 產生）。
原則：只有一條程式路徑。`mesh=uniform` 在程式內產生舊版版面的均勻節點；`mesh=file` 讀 `grid_gen.py` 的 JSON。
兩者走同一套更新式。已驗證 mesh=uniform 與凍結版在所有檢查步逐位相同（gate 0-1a）。

| 區段 | 舊 → 新 | 理由 | 推導 |
|---|---|---|---|
| 參數 | 新增 `mesh`、`grid`、`dtfac`、`dt`、`auxref`、`ref_every`、`seed`、`divop`、`dump_at`、`mode=e`、`eig_maxit/tol`、`proj`、`planar`、`jsrc`；`init=r`；`inc=p/j` | SPEC §17.1；舊參數與預設值不變 | — |
| JSON 讀取器 | 新增約 150 行的子集讀取器（物件、數字、數字陣列；未知 key 忽略；陣列長度上限） | 讀 grid 檔；sha256 由 Python 驗證（實作註） | — |
| 幾何 | 全域 `D` → 每軸 `h`（主間距）、`d`（對偶間距）、主/對偶座標陣列；x、z 週期接縫的 d 用兩端 h；y 牆邊 d = 相鄰 h | 非均勻網格 | §0 |
| 係數 | 單一 `ch`、`ceE*` → `cHy[j]=Δt/h`、`ay[j]=Δt/(ε_t d)`、`cEyn[j]`、`cHz[k]`，以及精確的間距比表（`rExz`、`rHxz`、`rEzx`、`rHzx`、`rEyk`、`rHyk`、`zfac`）。y 係數提出來，另一項乘間距比；均勻時比值恰為 1.0，所以舊的浮點運算逐位重現 | 每個差分除以正確的 h 或 d | §1 |
| 材料 | 介面主節點的 ε_t 由算術平均改為 h 加權平均（相等間距時就是算術平均）；新增每 (i, j) 的 `eTx/eTz/eNy` 與對應係數 `ayx/ayz/cEy2`，給階段 B 的光柵（層狀介質時是 1D 值的逐位複本） | 二階收斂；x–y 相依材料 | §1 |
| 更新式 | 6 分量都改成「y 係數 × (y 差分 − 比值 × 另一差分)」；CPML ψ 仍以差分儲存，乘同一個 y 係數 | 同上 | §1 |
| Courant | `dt = S·D` → Δt = S√3/sqrt(Σ 1/Δ_min²)（三軸最小間距相同時恰為 S·Δ）；另可用 `dt=` 或 `dtfac=` 指定 | SPEC §7.5；D12 的每週期整數步由 Python 傳入 `dt=` | §2 |
| ky、振幅 | `ky_disc` → `ky_local(h, n)`：y 用 aux 源的局部間距 Δ_a，x、z 用各自的間距 | 數值色散 | §5 |
| aux line | 節點 = 主網格 [0, Ny]，兩端以端點間距延伸；`c = dt/D` 拆成 `acH[u] = Δt/h`、`acE[u] = Δt/d` 與 `ikxH = iK̃x h` 等因子，使 y 差分用局部間距、相量項等於 Δt·iK̃；1D TF/SF 係數用源處間距 | 精確約化 | §3、§4 |
| aux_ref | 新增：與主網格同節點、同 ε、同 CPML 的複數 1D 線，在 j0 接收同樣的 TF/SF 修正；log 記錄主網格與它在 SF 區、全域的最大差 | gate 1-0、1-1（D8） | §3 |
| TF/SF | 係數改用 `cHy[j0−1]`、`ayx/ayz[·, j0]` | 局部間距 | §4 |
| CPML | 剖面改由物理座標計算；σ_max 用各端 PML 間距與該端折射率；兩端層數可不同 | SPEC §7.4 | §1 |
| 能量 | 舊式能量加入 W 權重；新增 `W_mod`（½ΣεW_E\|E^n\|² + ½ΣμW_H H^{n+½}·H^{n−½}），在 H 更新時累計 | gate 0-3 | §2 |
| 散度 | `div_max` 改成非均勻算子；`divop=u` 保留錯誤的均勻算子作對照 | gate 0-5 | §2 |
| DFT | 半格分量在 j = Ny 不再累計，輸出 NaN（修正 bug）；新增 `proj`（每列對橫向 Floquet 相位的加權投影）與 `planar`（每個 y 平面的相位殘差） | 盤點表；D12 | — |
| 輸出 | `meta.json` 追加欄位；新增 `grid_used.json`（所有座標、間距、ε），後處理唯一的座標來源；`fields_n<step>.bin`、`auxref_n<step>.bin` | SPEC §13.2 | — |
| 階段 B | x 或 z 非均勻時，`inc=a` 或 `auxref=1` 以固定訊息中止；`inc=p`（在分量實際座標取樣連續平面波）、`inc=j`（電流片）；光柵材料 | SPEC §7.6 | §3 |
| 模態注入（D22） | 新增 `inc=m`：每個非均勻週期軸解 A φ = κ² W φ（shift-invert 子空間迭代 + Rayleigh–Ritz，主網格的 h、d 權重），每個乘積模態一條 aux line（共用 y 幾何 `acH/acE`，自己的 K̃、ky、E0、H0），TF/SF 入射值 = Σ Re[aux_t · w_t · 剖面]；輸出 `modes.json`。`inc=a` 的程式碼與資料完全未動 | 失敗流程 B3-2 (s)；SPEC B1-3 | §3b |
| 本徵模式 | `mode=e`：以求解器自己的更新核心做 power iteration，求 λ_max 與 Δt_max = 2/sqrt(λ_max)（CPML 關閉、PEC 牆） | gate 0-4 | §2 |
| 無源執行 | mesh=file、inc=0、init≠a 時不建立平面波（kx = kz = 0），可用薄 x–z 格 | gate 0-3/0-4 | — |
