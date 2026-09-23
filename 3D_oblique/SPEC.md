# 規格整理 v 1.2.0

> 研究前置（單回合模式，依使用者選擇：Claude Code 專屬分階段計畫、單回合直接產出完整規格）。以下四個 code block 為研究摘要與假設列表，供事後追溯；使用者已確認方向，不再等待回覆。

```text
[研究1] 關鍵概念定義
- Yee FDTD：Yee(1966) 提出的交錯網格（staggered grid）有限差分時域法，E、H 分量在空間與時間上各偏移半格/半步，
  是本專案數值核心，來源：K.S. Yee, IEEE Trans. Antennas Propag., 1966；Taflove & Hagness,
  "Computational Electrodynamics: The Finite-Difference Time-Domain Method"（標準教科書，公式與本 spec 的
  離散色散關係、CPML 係數形式一致）。
- CPML（Convolutional PML）：Roden & Gedney (2000) 提出的 PML 變體，對消散/色散介質與低頻分量吸收效果更好，
  以 σ、κ、α 三參數與遞迴卷積實作，避免 Berenger 原始 PML 在斜向入射時的數值不穩定。本 spec 採用的
  σ_max ≈ 0.8(m+1)/(η0Δ) 多項式 grading 是 CPML 文獻中常見的經驗設定（m=3 時最常見）。
- TF/SF（Total-Field/Scattered-Field）：在計算域內劃出一個邊界面（本例為 y=y0 單一平面），內側為「入射+散射」
  總場，外側為純散射場；透過在邊界面切向分量加減入射場修正項來注入平面波，避免電流片源造成的雙向輻射與
  DC 電荷累積問題。
- 離散平面波精確解：連續色散關係 ky = sqrt(k0²−kx²−kz²) 不是 Yee 格點方程的精確解；必須用
  [sin(ωΔt/2)/(cΔt)]² = Σ[sin(k_iΔ/2)/Δ]² 反解 ky，否則 TF/SF 注入會在 SF 區漏出非零殘留場
  （leakage），這是本 spec 「關卡0」與「關卡1-1」要量化證明的核心風險點。
- 週期邊界 + 固定切向波數：因為 x、z 用 PBC 且 kx=2πm/Lx、kz=2πn/Lz 已被網格週期性鎖定為離散值，
  只要用 CW（連續波）+ ramp 而非寬頻脈衝，就能維持角度嚴格單一，這是 spec 明確禁止脈衝源（除非先處理
  角度色散問題）的原因。
```

```text
[研究2] 相似工具 / 業界比較（供理解本專案定位，非採購比較）
| 工具 | 授權 | 語言 | 斜向入射作法 | PML 型式 | 備註 |
|---|---|---|---|---|---|
| MEEP (NanoComp/meep) | 開源 GPL | C++/Python | k_point Bloch 週期 + 電流片源（可加 ramp/DFT） | 標準 PML（非 CPML 選項可切換） | 本專案 Level 3 的比對基準，本規格採用其 k_point 以 2π 為單位、Bloch 相位 exp(2πi k·r) 的慣例（來源：MEEP 官方文件 Python_User_Interface / Exploiting_Symmetry） |
| openEMS (thliebig/openEMS) | 開源 GPL | C++ (Matlab/Python 介面) | EC-FDTD，多用埠激發非 TF/SF | CPML | 工業級開源工具，證實 CPML 多項式 grading 手法業界常見 |
| gprMax | 開源 GPL | Python+Cython(+CUDA) | 平面波/點源，GPR 應用為主 | CPML | 驗證流程（含解析解比對）與本 spec 的分關卡驗證精神類似 |
| Lumerical FDTD | 商用 | GUI+腳本 | 斜向平面波源（Bloch/週期 BC） | PML | 商用參考，不在交付範圍內，僅供概念比對 |
- 本專案與上述工具的差異：本專案要求「離散精確解注入」並用「代入 Yee 更新式殘差 < 1e-12」證明，
  這比一般教學型 FDTD 專案的驗收基準嚴格很多（多數開源專案僅以視覺化或粗略能量守恆驗收）。
```

```text
[研究3] 相似 GitHub repo（供實作慣例參考，非依賴）
- NanoComp/meep — https://github.com/NanoComp/meep（本專案 Level 3 直接依賴其 Python 介面）
- gprMax/gprMax — https://github.com/gprMax/gprMax（CPML 係數命名、驗證報告格式可參考）
- thliebig/openEMS — https://github.com/thliebig/openEMS（EC-FDTD 與 CPML C++ 實作慣例參考）
- flaport/fdtd — https://github.com/flaport/fdtd（輕量 Python FDTD，僅供交錯網格索引寫法對照，
  不作為效能或正確性基準，功能遠不如本專案要求的驗證深度）
- 結論：本專案不「抄」任何 repo 的驗證方法；上述 repo 僅用來確認「CPML 多項式 grading」「TF/SF 面注入」
  「Yee 交錯索引」是業界標準寫法，不是本專案自創手法。
```

```text
[研究4] 建議方案與已確認事項（AskUserQuestion 已收斂，以下為最終決策記錄）
- 平台範圍：只產出 Claude Code 專屬分階段 instructions，不產出 Codex 版本。
  → 影響：「## Claude Code 分階段開發計畫」取代 skill 預設的雙平台標題。
- 研究/確認流程：單回合直接產出完整規格，不再等待第二輪確認。
- 以下為單回合模式下必須明示的假設（未逐一詢問使用者，皆標記為假設，執行中若與現實衝突以現實為準並回報）：
  A1. 白話規格讀者 = 具電磁/物理背景但非本專案程式撰寫者的共同研究者或計畫審查人（例如 PI、指導教授），
      非完全外行；因此「反射率」「相位」等物理詞彙保留，但避免 API/schema/DFT 內部實作用語。
  A2. C 語言標準 = C99/C11 + 標準函式庫（math.h、stdio.h）為主，不依賴 POSIX-only API，以便在 Windows
      （MinGW/gcc 或 WSL）與 Linux 都能編譯；Makefile 需同時給出兩種建置說明。
  A3. Python 依賴嚴格限制為 numpy + matplotlib（使用者已明講），不用 scipy；因此最小平方擬合（波前平面性）
      須用 numpy 正規方程（normal equation）手刻，不用 scipy.optimize；ky 求解須用使用者給的 arcsin 閉式解，
      不用數值求根器。
  A4. meep_ref.py 需要 mpi4py/meep 套件環境；本規格假設該環境由使用者另行提供或安裝，Claude Code 階段計畫
      只負責撰寫腳本與比對邏輯，若環境缺失以「環境缺失」回報，不得跳過比對邏輯本身的正確性。
  A5. 效能不是本專案的驗收指標（use case 是驗證/教學/研究用途，非量產級效能），因此不要求 OpenMP/SIMD，
      但程式碼不得刻意寫成無法在合理時間內跑完（Ly≥10λ0+2×N_pml、20000 步）的量級，20x/λ 網格下應可在
      一般工作站數分鐘內完成單次真空模擬。
  A6. 若第 1、2、3 關某一項在合理範圍內找不到 bug 而反覆 FAIL（例如 Meep 版本差異造成的已知限制），
      依 spec 的「工作規則」寫出假設 → 最小實驗 → 修正 → regression，最多重試 3 輪，第 3 輪仍 FAIL
      則在 validation_report.md 中列為 BLOCKED 並說明根因，不得放寬門檻或悄悄跳過。
```

---

## 技術規格文件

### 0. 文件範圍與交付物總覽

本規格描述一個「3D 均勻網格 Yee FDTD + CPML + TF/SF 斜向平面波注入」的 C 模擬引擎，搭配 Python
後處理／自我驗證，並與 Meep 進行交叉比對。規格本身**不包含**實際的物理推導過程（六分量公式、TF/SF
修正式的符號推導等）——那是「Stage 0」的必要產出物，本規格只定義**這些推導必須滿足的形式、公式骨架與
驗收門檻**，確保不同時間執行的 Claude Code session 都能做出可比對、可驗收的結果。

**最終交付物（固定路徑，位於專案根目錄）：**

| 檔案 | 內容 | 產出關卡 |
|---|---|---|
| `derivation.md` | 六分量 Yee 座標表、離散平面波六分量公式、TF/SF 修正式與正負號推導、PBC 索引規則、各驗證項目理論值公式 | Stage 0 |
| `fdtd3d_oblique.c` | 3D Yee FDTD 主程式（PBC + CPML + TF/SF + on-the-fly DFT + 介質層） | Stage 1–5 |
| `Makefile` | Windows(MinGW)/Linux 雙平台建置指令 | Stage 1 |
| `analyze.py` | Python 後處理：讀取 DFT 二進位 + JSON metadata、跑關卡 0/1/2 全部驗證與繪圖 | Stage 0–5 |
| `meep_ref.py` | Meep 參考模擬（真空、介質、DFT 取場、通量） | Stage 6 |
| `compare.py` | 讀取本程式與 Meep 的輸出，做 Level 3 (a)–(e) 全部比對 | Stage 6 |
| `validation_report.md` | 全部推導、每關驗證表（項目/理論值/量測值/門檻/PASS-FAIL）與對應圖 | Stage 7 |
| `tests/level{0,1,2,3}_*.py` | 各關卡對應的獨立驗證腳本（`analyze.py`/`compare.py` 呼叫這些模組） | 各 Stage |

### 1. 目標與非目標

**目標**
- G1：對真空中沿任意 (m, n) 週期鎖定角度入射的平面波，FDTD 注入必須是「Yee 離散方程的精確解」，
  以殘差 < 1e-12 證明，而非以連續解近似後「看起來還行」。
- G2：真空自我驗證（洩漏、波前平面性、色散收斂、橫波性、阻抗能量、PML 反射、長時間穩定性、對稱性）
  8 大項全數量化通過。
- G3：加入 n=1.5 半無限介質後，R/T 與 Fresnel/Snell 解析解比對誤差 <1%，並呈 O(Δ²) 收斂。
- G4：與 Meep 在真空 ky、波前角度、Fresnel R/T、DFT 場逐點比對、PML 反射記錄五項上交叉驗證，
  差異都要歸因到具體來源。

**非目標（Out of scope，明確排除）**
- NG1：不支援寬頻脈衝源（除非未來需求明確要求先解決角度隨頻散的問題，本規格版本不處理）。
- NG2：不支援非均勻網格、非直角座標、GPU/多執行緒平行化、任意形狀散射體。
- NG3：不做效能優化與量產化打包（無需安裝腳本、無需跨語言 binding）。
- NG4：不提供 GUI；輸出僅二進位 + JSON metadata + matplotlib 圖。
- NG5：Brewster 角驗證為選做（G3 的加分項），未通過不影響整體關卡放行。

### 2. 名詞與符號表

| 符號 | 意義 | 單位 |
|---|---|---|
| Δ | 均勻網格間距 Δx=Δy=Δz | m（以 λ0 為單位時無因次） |
| Δt | 時間步長，由 Courant 數 S=cΔt/Δ=0.5 決定 | s |
| λ0, ω0, k0 | 自由空間波長、角頻率、波數，k0=2π/λ0=ω0/c | — |
| m, n | x、z 方向週期模數（整數），kx=2πm/Lx、kz=2πn/Lz | — |
| θ, φ | 入射極角（自 +y 軸量起）、方位角 φ=atan2(kz,kx) | rad |
| K̃_i | 離散波數 (2/Δ)·sin(k_iΔ/2)，i=x,y,z | 1/m |
| ω̃ | 離散角頻率 (2/Δt)·sin(ω0Δt/2) | rad/s |
| ŝ, p̂ | s/p 偏振單位向量，ŝ∝K̃×ŷ，p̂∝ŝ×K̃ | — |
| N_pml | CPML 層數（每端） | 格 |
| j0 | TF/SF 注入平面的 y 格 index | — |
| j1 | 介質介面的 y 格 index（Level 2） | — |

### 3. 幾何、網格與預設參數

| 項目 | 預設值 | 規則 |
|---|---|---|
| Δ | λ0/20 | 色散收斂測試另需 λ0/10、λ0/40 兩組 |
| Lx, Lz | 2λ0, 3λ0 | 必須是 Δ 的整數倍（否則 PBC 與 kx/kz 週期鎖定失效），啟動時需檢查並在不整除時直接報錯終止，不得自動取整靜默修改 |
| m, n | 1, 1 | 對稱性測試需另跑 (−1, 1) |
| N_pml | 20 | CPML grading m=3, σ_max≈0.8(m+1)/(η0Δ), 含 κ、α |
| TF 區長度 | ≥10λ0 | y 方向扣除兩端 PML 後的自由傳播區間 |
| S（Courant 數）| 0.5 | 必須 <1/√3，啟動時檢查 |
| 傳播條件 | \|k_t\|<ω0/c，需留裕度 | 啟動時檢查 kx²+kz²<（ω0/c）² 且留至少 5% 裕度，否則報錯終止（避免消逝波/根號為負） |

**邊界**：x、z 為純 PBC；y 兩端 CPML。

### 4. Yee 分量座標表（Stage 0 必須產出，此處定義驗收格式）

`derivation.md` 與程式碼須包含以下表格（(i,j,k) 為格點整數 index，n 為時間步）：

| 分量 | 空間位置 | 時間取樣 |
|---|---|---|
| Ex | (i+½, j, k) | nΔt |
| Ey | (i, j+½, k) | nΔt |
| Ez | (i, j, k+½) | nΔt |
| Hx | (i, j+½, k+½) | (n+½)Δt |
| Hy | (i+½, j, k+½) | (n+½)Δt |
| Hz | (i+½, j+½, k) | (n+½)Δt |

驗收方式：Stage 0 的關卡0殘差測試（見 §11）即是對此表的間接驗證——若座標表錯誤，代入殘差不可能 <1e-12。

### 5. 激發訊號與離散平面波規格

**時間訊號**：CW + raised-cosine ramp（≥10 個週期），在 ω0 做 on-the-fly DFT（累加 `field·exp(-iω0 nΔt)`，
不得先存整段時間序列再離線 DFT，因為 20000 步 ×多平面×六分量的全時間序列會佔用不必要記憶體，且
「on-the-fly」是使用者明確要求）。

**ky 求解（禁止連續色散近似）**：

```
[sin(ω0Δt/2)/(cΔt)]² = [sin(kxΔ/2)/Δ]² + [sin(kyΔ/2)/Δ]² + [sin(kzΔ/2)/Δ]²
```

給定 kx、kz，解出 sin(kyΔ/2)，須檢查其平方落在 [0,1]（否則表示該 (m,n,Δ,Δt) 組合在此頻率下無傳播解，
必須在啟動時報錯而非產生 NaN）。ky 用 arcsin 閉式解：`ky = (2/Δ)·arcsin(sqrt(...))`，正根取傳播方向為 +y。

**離散量定義**：K̃_i=(2/Δ)sin(k_iΔ/2)（i=x,y,z）、ω̃=(2/Δt)sin(ω0Δt/2)。

**偏振**：ŝ∝K̃×ŷ；p̂∝ŝ×K̃；E0 依所測偏振（s 或 p）取 ŝ 或 p̂ 方向，使 K̃·E0=0（離散高斯定律）；
H0=K̃×E0/(μ0ω̃)。

**六分量公式**：`derivation.md` 須依 §4 座標表，對每個分量寫出形如
`E_x(i+½,j,k,n) = E0_x · cos(K̃x·(i+½)Δ + K̃y·jΔ + K̃z·kΔ − ω̃·nΔt)`（H 分量同理，時間用 (n+½)Δt、
相位可能有 ±Δ/2 或 ±Δt/2 的偏移，正負號需推導並附殘差測試佐證）的明確表達式，六個分量缺一不可。

### 6. TF/SF 修正

在 y=y0（E 節點）與 y=y0−Δ/2（H 節點，即 j0−½）上，對切向分量做修正：

| 受影響分量 | 位置 | 修正方向 |
|---|---|---|
| Ex, Ez | j0 | 用入射場修正相鄰 H 更新式中缺的鄰居項 |
| Hx, Hz | j0−½ | 用入射場修正相鄰 E 更新式中缺的鄰居項 |

`derivation.md` 須明確寫出每一項修正的正負號（TF 側加、SF 側減，或相反，依更新式方向而定），並在
關卡0/關卡1-1測試中以數值方式佐證正負號正確（正負號錯誤會導致洩漏或能量不守恆，測不出來就是符號錯）。

**明確禁止**：不得用電流片源（J）替代 TF/SF 面修正法。若之後有人想改用電流片，必須先在 `derivation.md`
補充「面內 ∇·J≠0 電荷累積」與「時間訊號 DC 分量造成靜電殘留場」兩個問題的分析，否則視為不合規變更。

### 7. PBC 索引規則

`derivation.md` 與程式碼註解須明確寫出交錯分量在邊界的 wrap 規則，例如：

- Nx 格點於 x 方向週期：Ex 在 i=0..Nx−1（cell 為主，非節點為主的量），Hz、Hy 需要的 `E(i=Nx)` 一律
  wrap 為 `E(i=0)`；`H(i=−1)` wrap 為 `H(i=Nx−1)`。
- z 方向同理以 Nz 為週期。
- y 方向兩端為 CPML，不做週期 wrap。

驗收方式：關卡0 的精確代入殘差測試涵蓋邊界格點（i=0、i=Nx−1、k=0、k=Nz−1），若 wrap 錯誤，邊界格點
殘差會明顯大於內部格點，門檻仍是 <1e-12，但需額外在報告中列出「邊界格點 vs 內部格點」殘差對照。

### 8. CPML 規格

- 僅施加於 y 方向兩端，各 N_pml=20 層。
- 多項式 grading，m=3：`σ(x)=σ_max·(x/N_pml Δ)^m`，`σ_max≈0.8(m+1)/(η0Δ)`。
- 含 κ（一般由 1 grading 到 κ_max，建議 κ_max∈[5,15]，需在 `derivation.md` 註明實際採用值與理由）、
  α（用於低頻穩定性，建議 α_max 一般在 0~π f0 ε0 量級，需註明採用值）。
- x、z 方向不施加 PML（純 PBC）。
- 介質延伸（Level 2）：ε 需延伸進遠端 PML 內部，PML 係數以介質中的 η、c 對應調整（`derivation.md`
  需說明介質 PML 的等效阻抗如何處理）。

### 9. DFT 輸出與 Metadata 格式

在指定的一個 y 平面、一個 x–y 切面、一個 x–z 切面上，對六分量做 ω0 的 on-the-fly DFT，輸出：

- 二進位檔（建議 `.bin`，double precision，行主序，檔名含平面類型與座標）。
- 對應 JSON metadata，至少包含：`{Delta, dt, Nx, Ny, Nz, N_pml, plane_type, plane_index, component_offsets:
  {Ex:[...], ..., Hz:[...]}, omega0, courant, m, n, kx, kz, ky_discrete, timestamp}`。
- `analyze.py` 讀檔時必須用 metadata 中的 `component_offsets` 對齊座標，不得寫死假設。

### 10. 驗證關卡總覽

四個關卡（0/1/2/3）為嚴格遞進關係：前一關全數 PASS 才能進下一關；任何一項 FAIL 依「工作規則」流程
（假設→最小實驗→修正→全量 regression）處理，不得放寬門檻、修改理論值或挑對自己有利的量測區域（§0 前言
的不可變規則同樣適用於本規格本身）。

#### 關卡 0：精確注入證明（不跑完整模擬）

| 項目 | 理論值 | 門檻 | 對照組 |
|---|---|---|---|
| max\|E_update−E_analytic\|/max\|E\| | 0 | <1e-12 | 連續 ky 版本，殘差應為 O((k0Δ)²) |
| H 同上 | 0 | <1e-12 | 同上 |
| 解析場離散 ∇·E | 0 | <1e-12 | — |
| 解析場離散 ∇·H | 0 | <1e-12 | — |

#### 關卡 1：真空自我驗證（8 項）

| # | 項目 | 理論值/預期 | 門檻 |
|---|---|---|---|
| 1 | SF 區洩漏（因果時窗內） | 0 | max\|E_SF\|/\|E_inc\| <1e-10；連續 ky 對照組洩漏隨 Δ² 縮放（兩解析度驗證） |
| 2 | 波前平面性 | 擬合 (kx,ky,kz) = 離散理論值 | 相對誤差 <1e-6；殘差 RMS <1e-3 rad；PML 回波需量化並扣除或加長 Ly |
| 3 | 色散收斂 | 與連續 ky 誤差 ~O(Δ²) | log-log 斜率 2±0.1；與離散 ky 誤差 <1e-6（三種 Δ：λ/10, λ/20, λ/40） |
| 4 | 橫波性/高斯定律 | E·K̃=H·K̃=E·H≈0，∇·E=∇·H≈0 | 相對值 <1e-10，且不隨時間增長 |
| 5a | 阻抗 | \|H\|/\|E\|=\|K̃\|/(μ0ω̃) | 一致到 1e-8；與 1/η0 差距為 O((k0Δ)²)（需報告，非門檻） |
| 5b | 能量守恆 | S_y 沿 TF 區多平面不變 | 相對變化 <1e-4；S 與 k 夾角需報告 |
| 6 | PML 反射 | 隨 N_pml 增加而降低 | 目標 <−40 dB |
| 7 | 長時間穩定性 | 關源後能量單調衰減 | ≥20000 步，衰減至 <1e-8×峰值，後期不增長 |
| 8 | 對稱性 | (m,n)→(−m,n) 為 x 鏡像 | s、p 偏振都需驗證，鏡像誤差需量化（沿用 §11 通用數值誤差門檻 <1e-6 相對值，除非另有物理理由） |

#### 關卡 2：介質物理定律（Fresnel/Snell）

| 項目 | 理論值 | 門檻 |
|---|---|---|
| R（s、p 偏振） | Fresnel 解析解 | 誤差 <1%（Δ=λ0/20），且 O(Δ²) 收斂 |
| T（s、p 偏振） | Fresnel 解析解 | 誤差 <1%，O(Δ²) 收斂 |
| R+T | 1 | 誤差 <1e-4 |
| Snell（介質中 ky） | 介質離散色散關係預測值 | 沿用 <1e-6 相對誤差（與關卡1-3一致基準） |
| Brewster 低谷（選做） | p 偏振 R 在 Brewster 角附近出現低谷 | 選做，不影響放行 |

介面 ε 處理（切向算術平均 / 法向調和平均，或錯開介面位置避開 Ey）須在 `derivation.md` 明確記錄選用方式，
因為會直接影響與 Meep 的比對結果。

#### 關卡 3：與 Meep 比對

| 項目 | 門檻 | 備註 |
|---|---|---|
| (a) 真空 ky | 相對誤差 <1e-6 | 同 Δ、Δt 下比對；不符先查單位/Courant，不先懷疑物理 |
| (b) 波前角度與相位殘差 | 沿用關卡1-2 門檻 | — |
| (c) Fresnel R/T | 與 Meep 差 <0.5%，各自與解析解比對 | Meep 需設 `eps_averaging=False` 或本程式採相同介面平均策略 |
| (d) DFT 場圖逐點比對 | 真空中相對 L2 差 <1e-3 | 需先對齊 Yee 座標、歸一化振幅與相位 |
| (e) PML 反射 | 不設門檻，只記錄 | 兩者 PML 實作不同 |

Meep 設定重點：`k_point=mp.Vector3(kx/2π,0,kz/2π)`（Meep 以 2π 為單位，Bloch 相位 `exp(2πi k·r)`，
來源：Meep 官方文件 Python User Interface / Exploiting Symmetry）；`mp.PML(d, direction=mp.Y)`；
`force_complex_fields=True`；來源需在 `meep_ref.py` 註解中附官方文件連結；`amp_func=exp(2πi k·r)` 的
振幅依 s/p 偏振設定；用 `add_dft_fields`/`get_dft_array`/`get_array_metadata` 取場與座標。

### 11. 通用規則

- 所有數字結果須附單位與所用網格解析度（Δ、Δt）。
- 不得為了通過測試而放寬門檻、修改理論值，或挑選有利的量測區域/時間窗。
- FAIL 處理流程固定：寫出假設 → 設計最小實驗 → 修正 → 重跑全部既有測試（不得只重跑失敗項）。
- 每一關結束須產出：一張驗證表（項目｜理論值｜量測值｜門檻｜PASS/FAIL）＋對應圖，寫入
  `validation_report.md` 對應章節（不得在該關全數 PASS 前先寫「已完成」）。

### 12. 非功能需求

| 類別 | 要求 |
|---|---|
| 數值精度 | C 端一律 double precision；不得用 float 累積 DFT 或做殘差比較 |
| 可重現性 | 不使用隨機數；若未來新增隨機測試（無此需求）須固定 seed 並記錄 |
| 可攜性 | C99/C11 + 標準函式庫；Windows(MinGW)/Linux 皆可編譯（假設 A2） |
| Python 依賴 | 僅 numpy + matplotlib；動畫僅用 `matplotlib.animation.PillowWriter` |
| 記憶體 | on-the-fly DFT，不得整段時間序列常駐記憶體 |
| 錯誤處理 | 啟動時的參數檢查（Lx/Lz 整除、Courant<1/√3、傳播條件裕度、ky 根號域）失敗須直接報錯終止，不得靜默修正 |

### 13. 風險與依賴

| 風險 | 影響 | 緩解 |
|---|---|---|
| ky 閉式解在極端 (m,n,Δ) 組合下根號為負 | 程式無法啟動 | 啟動時檢查並報錯，非執行期 NaN |
| PML 回波污染波前擬合 | 關卡1-2 誤判 FAIL | 先量化回波量級，必要時加長 Ly，不得只換一段「乾淨」區域 |
| Meep 環境缺失（假設 A4） | 關卡3 無法執行 | 明確回報環境缺失，不得略過比對邏輯撰寫 |
| 介面 ε 平均策略選擇影響 Meep 比對 | 關卡3(c)(d) 誤差超標 | `derivation.md` 記錄策略，`meep_ref.py` 對應設定 `eps_averaging` |
| Lx、Lz 非 Δ 整數倍 | kx/kz 週期鎖定失效，PBC 出錯 | 啟動即檢查並報錯終止 |

### 14. 邊界/異常情境（Edge & Abuse Cases）

| 情境 | 預期行為 |
|---|---|
| \|kx\|²+\|kz\|²過於接近 (ω0/c)²（掠射角附近） | 啟動即報錯（傳播條件裕度不足），不得產生 NaN 或消逝波 |
| N_pml 設太小（如 <5） | 允許執行但 PML 反射門檻（<−40dB）大機率 FAIL，報告需誠實記錄，不得調整門檻 |
| TF 區長度不足 10λ0 | 啟動時警告並可選擇報錯（建議：報錯終止，因會影響洩漏因果時窗判定） |
| (m,n) 使 kx 或 kz 為 0（正入射的退化情形） | ŝ、p̂ 定義中 K̃×ŷ 可能退化，`derivation.md` 需說明此邊界情形是否在支援範圍內（若不支援，於文件與程式檢查中明確排除並報錯，不得靜默給錯誤偏振） |
| Meep 版本不同造成 k_point 慣例差異 | `meep_ref.py` 需以官方文件連結佐證慣例，若實測不符先懷疑單位換算，並在 `validation_report.md` 記錄查證過程 |

---

## 非技術規格文件

> 讀者設定：計畫主管、共同研究者、審查人——懂物理與「這個模擬對不對」的意義，但不需要知道程式怎麼寫。
> 以下用白話說明「我們在做什麼、怎麼知道它是對的、什麼時候算完成」。

### 這個東西是做什麼的

我們要做一套會自己「證明自己算對了」的電磁波模擬程式：模擬一道光（電磁波）從某個角度斜斜地射入真空，
之後再射進一塊玻璃般的介質（折射率 1.5），看它怎麼反射、怎麼折射。這套程式是用一種叫 FDTD
（把空間和時間切成很多小格子、一步一步推進）的方法自己刻出來的，之後拿業界常用的模擬軟體 Meep 來對答案，
確認我們刻的東西跟別人做的一致。

### 為什麼要一關一關驗證，不能寫完就說完成

FDTD 模擬最怕的狀況是「畫面看起來很合理，但其實算錯了」——例如角度差一點、能量偷偷變多或變少、
邊界處理錯誤讓假訊號跑進來。這些錯誤肉眼很難看出來，所以我們不會用「看起來對」當作完成標準，
而是設計一整套數字化的檢查項目（我們稱為「關卡」），每一項都有明確的及格標準（門檻），全部通過才能
進到下一步。如果沒通過，不能放寬標準去「湊過關」，只能真的把問題找出來修掉。

### 四個關卡在檢查什麼（白話版）

1. **關卡 0：注入是不是精確的**
   在正式開始跑模擬之前，先用紙筆（其實是程式）驗證「我們準備要送進模擬裡的那道波」，跟這套模擬程式
   自己的計算規則是完全吻合的，誤差要小到幾乎等於零（1 兆分之一等級）。如果這一步就沒過，後面全部
   都不用跑，因為地基是歪的。

2. **關卡 1：在什麼都沒有的空間裡，波的行為對不對**
   包含 8 個子項目，白話理解：
   - 波有沒有「漏」到不該出現的地方（洩漏）。
   - 波前（波峰所在的面）是不是真的平的、角度對不對。
   - 換更細的網格，答案會不會如預期地變得更準（收斂）。
   - 波是不是真的是橫波（電場、磁場、傳播方向互相垂直，這是電磁波的基本性質）。
   - 電場和磁場的比例對不對（阻抗），能量沿路有沒有守恆。
   - 邊界吸收層（模擬邊界的「消音棉」）夠不夠安靜，不會反彈假訊號回來。
   - 長時間跑（兩萬步以上）不會爆炸或發散。
   - 把入射角左右鏡射，結果也應該跟著鏡射，這是基本的物理對稱性檢查。

3. **關卡 2：加入介質之後，反射與折射的比例對不對**
   放進一塊折射率 1.5 的介質，量測反射了多少、穿透了多少，跟課本上的 Fresnel 公式（描述光在介面反射/
   折射比例的經典公式）比對，誤差要在 1% 以內，而且反射加穿透的能量要精確等於 100%。

4. **關卡 3：跟業界工具 Meep 對答案**
   把同樣的設定丟進 Meep 這套廣泛使用的模擬軟體，比較兩者在真空中的波數、波前角度、反射穿透率、
   還有逐點的場圖，確認差異都在合理誤差內，而且每一個差異都要講得出「為什麼」（不能只說「有點差但沒關係」）。

### 完成後會拿到什麼

- 一份會自己驗證的模擬程式（C 語言）。
- 一套 Python 分析工具，可以讀模擬結果、畫圖、自動跑完上面四關的檢查。
- 一份跟 Meep 比對用的腳本。
- 一份完整報告 `validation_report.md`，裡面有每一關的驗證表格（該項目、理論上應該是多少、實際量到多少、
  及格標準、有沒有過）和對應的圖，任何人都可以照著報告重新檢查一次結論。

### 完成標準

每一關的每一項都要「有數字、有圖、有沒有過」三件事齊全，且全部標記為通過，才算這一關完成；四關全部
完成，加上最終報告整理好，才算整個專案完成。中途若某一項卡住，會先寫下「我們認為問題出在哪裡」，
再設計一個小實驗去驗證這個猜測，而不是直接調整及格標準讓它「看起來過了」。

---

## Claude Code 分階段開發計畫

> 依使用者選擇，僅提供 Claude Code 專屬版本（不含 Codex 版本）。每個 Stage 都是可獨立驗收、可回滾的
> vertical slice；前一 Stage 的 DoD 全部通過才進下一 Stage。Stage 0 為推導與骨架，Stage 6 為倒數第二的
> 完整測試／迴歸關，Stage 7 為最終文件交付關。

### Stage 0 — 專案骨架 + 推導文件 + 關卡 0

**目標**：建立檔案骨架、產出 `derivation.md`（含 §4–§7 要求的全部推導與表格）、實作最小 Python
（`tests/level0_exact_injection.py`）驗證離散平面波精確代入 Yee 更新式的殘差 <1e-12。此階段**不需要**
完整 C 引擎，只需要能代入單步更新式的最小驗證腳本（可先用 Python 實作單步驗證，C 引擎於 Stage 1 建立）。

**涉及檔案**：`derivation.md`、`tests/level0_exact_injection.py`、`README.md`（專案總覽骨架）、
`Makefile`（先放空殼，Stage 1 補內容）。

**Claude Code instructions（可直接貼用）**：

```
你是資深計算電磁工程師。這是 SPEC.md 定義的 3D Yee FDTD 斜向入射專案的第 0 階段。

任務：
1. 讀取專案根目錄的 SPEC.md，理解 §2–§7（符號表、幾何、Yee 座標表、激發與 ky 求解、TF/SF 修正、PBC 索引）。
2. 撰寫 derivation.md，內容必須包含且僅包含以下五節，每節都要有完整數學推導（不是抄 SPEC 的骨架，
   是把骨架填成真正的推導）：
   (1) 六個分量的 Yee 座標表（含推導為什麼是這個偏移，不是只抄表格）
   (2) 離散平面波的六分量明確表達式（含 K̃_i、ω̃、ŝ/p̂ 偏振定義、H0=K̃×E0/(μ0ω̃) 的推導）
   (3) TF/SF 修正式，含正負號推導過程（要說明為什麼是加不是減）
   (4) PBC 索引方式（i=0、i=Nx-1 邊界的具體 wrap 規則，x、z 各自寫一次）
   (5) 關卡0、關卡1 每一項驗證的理論值計算公式（照抄 SPEC.md §10 的表格逐項展開成可執行的公式）
3. 用 Python（僅 numpy）實作 tests/level0_exact_injection.py：
   - 建一個小型 3D 網格（例如 Nx=Ny=Nz=8 即可，重點是測代入邏輯不是效能）。
   - 用 derivation.md 推出的六分量公式產生解析場（含 PBC）。
   - 用同一組解析場代入一步 Yee 更新式（照 SPEC.md §7 的 PBC wrap 規則），算出更新後的場。
   - 計算 max|E_update - E_analytic|/max|E|，H 同樣計算，斷言 <1e-12。
   - 計算解析場的離散散度 ∇·E、∇·H，斷言 <1e-12。
   - 對照組：把 ky 換成連續色散 sqrt(k0²-kx²-kz²) 重跑一次，印出殘差量級，斷言它明顯大於離散版本
     （量級應接近 O((k0Δ)²)，不需要精確斷言數值，但要印出來供人工確認）。
4. 在 README.md 寫最小骨架：專案一句話介紹、檔案清單、如何跑 tests/level0_exact_injection.py。

嚴禁事項：不得為了讓殘差變小而偷偷改用連續 ky 當作「離散」版本；不得省略對照組；不得跳過散度檢查。
```

**測試指令**：
```bash
python tests/level0_exact_injection.py
```

**DoD**：
- `derivation.md` 五節齊全，且六分量公式、TF/SF 修正式、PBC 規則皆為具體數學式（非佔位文字）。
- `tests/level0_exact_injection.py` 執行通過：E、H 殘差 <1e-12，∇·E、∇·H <1e-12，對照組殘差明顯較大
  且量級與 O((k0Δ)²) 大致相符（同一數量級即可，非精確等式）。

**風險與回滾**：若殘差怎麼調都過不了 1e-12，優先懷疑 (1) Yee 座標偏移打錯 (2) H 的半步時間相位漏加。
不得放寬門檻到 1e-9 之類；找不到原因就先在 derivation.md 記錄假設與已排除的可能性，重新走一次
「假設→最小實驗→修正」流程，不進入 Stage 1。

---

### Stage 1 — C 語言 3D Yee 核心引擎（純 PBC，無 CPML/TF-SF）

**目標**：把 Stage 0 驗證過的公式與索引規則搬進 C，建立可編譯執行的 3D Yee 更新迴圈（x、z 為 PBC，
y 暫時用簡單截斷邊界，CPML 留到 Stage 2），並把關卡 0 的驗證延伸成多步驟版本確認引擎本身無誤。

**涉及檔案**：`fdtd3d_oblique.c`（初版：grid 配置、Yee 更新式、PBC index wrap）、`Makefile`（補齊
Windows/Linux 雙平台建置目標）、`tests/level0b_multistep.py`（呼叫編譯後的執行檔，比對多步後的場）。

**Claude Code instructions**：

```
延續 Stage 0 的 derivation.md，開始寫 fdtd3d_oblique.c 的第一版核心引擎。

任務：
1. 用 C99，double precision，實作：
   - 3D 陣列配置（Ex, Ey, Ez, Hx, Hy, Hz），依 derivation.md 的 Yee 座標表決定各分量陣列大小。
   - 標準 Yee 更新式（先不含 CPML，y 方向暫時用最簡單的邊界，例如場在最外層固定為 0，因為這階段
     只測 PBC 與核心迴圈，不測邊界吸收）。
   - x、z 方向的 PBC index wrap，嚴格照 derivation.md §4 的規則寫（用註解引用 derivation.md 的章節）。
   - 啟動時的參數檢查（SPEC.md §3、§13）：Lx/Lz 是否為 Δ 整數倍、Courant<1/√3、傳播條件裕度，
     任何一項不合格直接印錯誤訊息並以非 0 狀態碼結束，不做靜默修正。
2. 加一個「debug 模式」：可以從命令列參數指定用 derivation.md 的解析公式把整個網格初始化為解析場
   （不是從 TF/SF 注入，是直接灌整個網格，純粹測試更新式本身），推進 N 步後把場輸出成二進位檔。
3. 寫 Makefile：`make` 目標同時給出 Linux (gcc) 與 Windows (MinGW gcc) 的建置說明（可用同一份
   Makefile，或是文件註明兩種呼叫方式）。
4. 寫 tests/level0b_multistep.py：呼叫上面的 debug 模式跑 50 步，讀輸出二進位檔，跟 Python 端算出的
   解析場（用同一組公式推進 50 步）比較，殘差門檻沿用 <1e-12（純 PBC、無 PML、無 TF/SF 的封閉系統，
   精確解在此假設下仍應嚴格成立）。額外印出「邊界格點 vs 內部格點」殘差分開統計（SPEC.md §7 要求）。

嚴禁事項：不得引入任何隨機數；不得用 float；不得依賴 POSIX-only 標頭。
```

**測試指令**：
```bash
make
python tests/level0b_multistep.py
```

**DoD**：編譯成功（Linux 與 Windows/MinGW 至少各驗證一種，若環境只有一種可用，另一種以「文件審查」
方式確認 Makefile 語法正確並記錄未實測）；50 步後 E/H 殘差 <1e-12；邊界格點殘差與內部格點同量級
（不得系統性偏大，若偏大即 PBC wrap 有誤）。

**風險與回滾**：若邊界格點殘差明顯大於內部格點，優先檢查 index wrap 的 off-by-one；不得先懷疑「PML
還沒加所以邊界本來就會怪」——這階段 y 方向邊界固定為 0 只影響最外層一兩格，x/z 的 PBC 邊界應該和內部
一樣精確。

---

### Stage 2 — CPML（y 方向兩端）

**目標**：在 y 方向兩端加入 CPML（§8），先確認數值穩定（長時間跑不發散），完整的 -40dB 反射量測留到
Stage 4（需要 TF/SF 源才能做有意義的反射量測）。

**涉及檔案**：`fdtd3d_oblique.c`（新增 CPML 輔助陣列與更新式）、`tests/level_cpml_stability.py`。

**Claude Code instructions**：

```
在 fdtd3d_oblique.c 中加入 CPML，範圍限定 y 方向兩端各 N_pml=20 層，依 SPEC.md §8：
- 多項式 grading m=3，σ_max≈0.8(m+1)/(η0Δ)，實作 σ、κ、α 三參數的分層 profile，並在程式碼註解寫出
  實際採用的 κ_max、α_max 數值與選擇理由（回填進 derivation.md 一個新小節「CPML 參數選擇」）。
- 用標準 CPML 遞迴卷積（auxiliary variable ψ）更新法，只作用在 y 方向的空間微分項。
- x、z 方向維持純 PBC，不受 CPML 影響。

寫 tests/level_cpml_stability.py：不加任何源，只把中心區域一小塊初始化非零場（例如一個高斯脈衝點），
跑 5000 步，確認：
1. 全域場不發散（每 500 步印一次最大場值，應呈下降或持平趨勢，不應成長）。
2. 最終場能量遠小於初始（因為波會跑到 PML 被吸收）。

這階段先不要求 -40dB 的精確反射係數（那需要 Stage 4 的 TF/SF 源與 DFT 才能精確量測），只確認 CPML
本身寫得穩定、不發散。
```

**測試指令**：
```bash
make
python tests/level_cpml_stability.py
```

**DoD**：5000 步內場值單調下降或持平，無發散（NaN/Inf 檢查需內建在 C 程式中，一旦偵測到 NaN 立即報錯
退出並印出發生的步數與位置）。

**風險與回滾**：若初期就發散，優先檢查 α、κ 的數值範圍（過大/過小都會不穩定）；不得靠調小 Δt 掩蓋
CPML 本身寫錯的問題（Courant 數已固定為 0.5，不得為了穩定而私自更改）。

---

### Stage 3 — TF/SF 注入（CW+ramp 源、on-the-fly DFT）

**目標**：在 y=y0 實作 TF/SF 面注入（§5、§6），CW+raised-cosine ramp 訊號，並在指定平面做 ω0 on-the-fly
DFT，輸出二進位+JSON metadata（§9）。

**涉及檔案**：`fdtd3d_oblique.c`（TF/SF 修正項、DFT 累加、輸出模組）、`tests/level1_leakage.py`。

**Claude Code instructions**：

```
在 fdtd3d_oblique.c 中實作：
1. y=y0 平面的 TF/SF 注入，依 derivation.md §3（TF/SF 修正式）對 Ex、Ez（j0）與 Hx、Hz（j0-½）做修正，
   正負號必須與 derivation.md 完全一致（程式碼註解引用推導章節）。
2. CW + 10 週期 raised-cosine ramp 的時間訊號生成函式。
3. On-the-fly DFT：對六分量在 (1) 一個 y=const 平面 (2) 一個 x-y 切面 (3) 一個 x-z 切面，
   每步累加 field·exp(-iω0·nΔt)（實部虛部分開累加，double precision），模擬結束時寫出二進位檔
   +對應 JSON metadata（SPEC.md §9 的欄位需齊全，尤其 component_offsets 一定要照 Yee 座標表填）。
4. 命令列參數需可指定 m, n, 偏振（s/p）、輸出平面座標。

寫 tests/level1_leakage.py：
- 跑一次模擬，時間長度只到「入射波抵達遠端 PML 之前」的因果時窗（用 Ly、注入平面位置、群速度算出
  這個時窗，寫成明確公式，不要用經驗值瞎猜）。
- 讀 SF 區（y<y0 的一側）的 DFT 結果，算 max|E_SF|/|E_inc| ，斷言 <1e-10。
- 對照組：把 C 程式加一個編譯開關或執行期參數，改用連續 ky 算 TF/SF 注入相位，重跑，記錄洩漏量。
- 用兩種解析度（Δ=λ0/20 與 λ0/40）各跑一次對照組，確認洩漏量隨 Δ² 縮放（log-log 兩點斜率約 2）。
```

**測試指令**：
```bash
make
python tests/level1_leakage.py
```

**DoD**：離散 ky 版本洩漏 <1e-10；連續 ky 對照組洩漏隨 Δ² 縮放（兩點估計斜率落在合理範圍，例如 1.5–2.5，
非嚴格 2±0.1，因只有兩個資料點；精確的 2±0.1 斜率驗證留到 Stage 4 的色散收斂項，用三種解析度）。

**風險與回滾**：若離散版本洩漏也超標，優先重查 TF/SF 正負號（Stage 0 的符號推導）與 DFT 累加是否漏了
ramp 期間的暫態（因果時窗的起點需晚於 ramp 結束，否則暫態本身就會被誤判為洩漏）。

---

### Stage 4 — 完整關卡 1（8 項全數）

**目標**：完成 SPEC.md §10 關卡 1 的全部 8 項驗證，含波前擬合、色散收斂（三種 Δ）、橫波性/高斯定律、
阻抗能量、PML 反射 dB、20000 步穩定性、(m,n)→(−m,n) 對稱性（s、p 偏振）。

**涉及檔案**：`analyze.py`（新增波前擬合、色散收斂、橫波性、能量、PML 反射、對稱性各模組）、
`tests/level1_*.py`（拆成多個子腳本或一個主腳本呼叫多個函式，需與 `analyze.py` 共用邏輯，不得重寫兩份）、
圖檔輸出到 `figures/level1/`。

**Claude Code instructions**：

```
基於 Stage 3 產出的 DFT 二進位+JSON，實作 analyze.py 的核心分析函式（僅用 numpy+matplotlib），並跑完
SPEC.md §10 關卡1 的 8 個項目：

1. 波前平面性：讀 TF 區某個切面的 DFT 相位，unwrap 後，用 numpy 手刻正規方程（A^T A x = A^T b，
   不用 scipy）做 3D 平面 φ=k·r+φ0 最小平方擬合。輸出擬合 (kx,ky,kz) 與理論值的相對誤差（斷言<1e-6）、
   殘差 RMS（斷言<1e-3 rad）。先量化 PML 回波（比較有無 PML 情形下的相位殘差差異，或用時間窗排除
   PML 建立前的暫態），選擇扣除受污染區域或加長 Ly，並在報告中記錄選擇理由（不得為了通過默默丟資料）。
   畫兩張圖：x-y 切面瞬時場疊理論等相位線（應為平行斜直線）、x-z 平面相位殘差圖。

2. 色散收斂：用 Δ=λ0/10, λ0/20, λ0/40 各跑一次 Stage 1-3 的完整流程（呼叫 fdtd3d_oblique 執行檔），
   量測波前擬合得到的 ky，分別跟連續 ky 與離散理論 ky 比較誤差，畫 log-log 圖，用 numpy.polyfit
   對 log(Δ) vs log(誤差) 做線性擬合，取斜率，斷言在 2±0.1；離散 ky 誤差斷言全部 <1e-6。

3. 橫波性與高斯定律：在 TF 區內（排除 TF/SF 平面兩側各一格、排除 PML）算離散 ∇·E、∇·H 相對值
   （斷言<1e-10，且需驗證多個時間點確認不隨時間增長），以及 E·K̃、H·K̃、E·H 是否≈0。

4. 阻抗與能量：算 |H|/|E| 跟 |K̃|/(μ0ω̃) 比較（斷言一致到1e-8），並報告與 1/η0 的差距（非門檻，僅報告）；
   實作 Yee 分量內插到同一點的方法（明確寫出用什麼內插規則，例如相鄰兩點平均），算時間平均 Poynting
   S=½Re(E×H*)，在 TF 區多個 y 平面上比較 S_y 相對變化（斷言<1e-4），報告 S 方向與 k 方向夾角。

5. PML 反射：用 Stage 2/3 的 SF 區穩態場當作 PML 反射訊號，換算反射係數 dB，跑不同 N_pml（例如
   10、20、30）觀察趨勢，目標 N_pml=20 時 <−40dB。

6. 20000 步穩定性：跑 ≥20000 步，中途關掉源（ramp 結束後可選擇直接停止注入或反向 ramp down，
   需在 derivation.md 註明選擇），確認總能量單調衰減到 <1e-8×峰值，且後期不增長（用 C 程式輸出
   逐步能量到一個 log 檔，Python 端檢查單調性，允許小幅數值雜訊但不得有系統性上升趨勢）。

7. 對稱性：另跑 (m,n)=(-1,1)，比較與 (1,1) 的場圖是否為 x 方向鏡像（s、p 偏振各跑一次），用相對誤差
   量化鏡像程度並斷言 <1e-6（或說明若因數值路徑不同而略寬鬆的理由，需在報告記錄）。

每一項都要在 validation_report.md 對應章節寫一張表（項目|理論值|量測值|門檻|PASS/FAIL）+對應圖，
存到 figures/level1/。全部 8 項 PASS 才算 Stage 4 完成。
```

**測試指令**：
```bash
python analyze.py --level 1 --all
```

**DoD**：SPEC.md §10 關卡1 表格全部 8 項 PASS，圖與表格已寫入 `validation_report.md` 對應段落（草稿
即可，Stage 7 會做最終整理）。

**風險與回滾**：任一項 FAIL 時，先確認不是關卡0/Stage1-3 遺留的 bug 復發（跑一次 regression：
`tests/level0_exact_injection.py`、`tests/level0b_multistep.py`、`tests/level1_leakage.py`），
確認地基仍然穩固後才在關卡1本身找問題。

---

### Stage 5 — 介質層與 Fresnel/Snell 驗證（關卡 2）

**目標**：在 TF 區 y=y1 放置 n=1.5 半無限介質（延伸進遠端 PML），實作 §10 關卡2 的 R/T/Snell 驗證。

**涉及檔案**：`fdtd3d_oblique.c`（介質 ε 陣列、介面平均策略）、`analyze.py`（新增通量/R/T 計算）、
`tests/level2_fresnel.py`。

**Claude Code instructions**：

```
在 fdtd3d_oblique.c 中加入介質支援：
1. 在 y=y1（TF 區內某平面）之後（含遠端 PML 內部）把 ε 設為 1.5²，其餘保持真空 ε0。
2. 依 SPEC.md §10 關卡2 選定介面 ε 處理策略（切向 Ex/Ez 算術平均、法向 Ey 調和平均，或錯開介面位置
   避免 Ey 落在介面上），把選擇與理由寫進 derivation.md 新章節「介面 ε 處理」，這個選擇會影響 Stage 6
   跟 Meep 的比對，必須跟 meep_ref.py 的 eps_averaging 設定對應。
3. 命令列參數加入介質開關與 n 值。

在 analyze.py 加入：
- SF 區通量（反射）與介質內通量（穿透）計算，需先定義清楚「通量」怎麼從 DFT 場算（用 Poynting
  在垂直於 y 的平面積分）。
- R = SF通量/入射通量，T = 介質內通量/入射通量。

寫 tests/level2_fresnel.py：
- s、p 偏振各跑一次，Δ=λ0/20，算 R、T，跟 Fresnel 解析公式（需在腳本內明確寫出解析公式，附上角度
  換算）比較，斷言誤差<1%。
- 額外用 Δ=λ0/10 跑一次，確認誤差呈 O(Δ²) 收斂（跟 λ0/20 的結果比較，誤差應明顯縮小，量級對應）。
- 斷言 R+T 在誤差 1e-4 內等於 1。
- Snell 檢查：量測介質內 DFT 相位得到的 ky，跟介質中離散色散關係（同 SPEC.md §5 公式但 c 換成
  介質中相速度）比較。
- （選做，不影響放行）掃描接近 Brewster 角的 (m, Lx) 組合，畫 p 偏振 R vs 角度圖，確認低谷存在。
```

**測試指令**：
```bash
python analyze.py --level 2 --all
```

**DoD**：SPEC.md §10 關卡2 表格全部必要項（R、T 誤差<1%且 O(Δ²) 收斂、R+T誤差<1e-4、Snell）PASS；
Brewster 選做項執行與否都需在報告中註明（做了就附圖，沒做就寫「選做未執行」，不得假裝做了）。

**風險與回滾**：若 R+T≠1 誤差過大，優先檢查介面 ε 平均策略是否與通量計算的取樣點一致（例如 Ey 沒有
正確平均會導致法向能流計算系統性偏差）。

---

### Stage 6 — Meep 比對 + 全量 Regression（倒數第二關：完整測試關）

**目標**：完成 SPEC.md §10 關卡3 全部 (a)-(e) 比對項目，並把 Stage 0–5 已建立的**全部**測試重新跑一次，
確認新增的介質/比對程式碼沒有破壞先前已 PASS 的項目。

**涉及檔案**：`meep_ref.py`、`compare.py`、`tests/run_all_regression.py`（依序呼叫 Stage 0–6 全部測試腳本）。

**Claude Code instructions**：

```
1. 撰寫 meep_ref.py：
   - 用跟本程式相同的 Lx, Ly, Lz、resolution=1/Δ、Courant=0.5。
   - mp.PML(d, direction=mp.Y)，只在 y 方向。
   - k_point = mp.Vector3(kx/(2*pi), 0, kz/(2*pi))（Meep 以 2π 為單位，Bloch 相位 exp(2πi k·r)，
     來源見 SPEC.md §10 關卡3 備註引用的 Meep 官方文件；用量到的相位梯度反驗證這個慣例，不要只信文件）。
   - force_complex_fields=True。
   - 源用 mp.Source + amp_func=lambda pos: cmath.exp(2j*pi*(kx/(2*pi)*pos.x + kz/(2*pi)*pos.z))，
     振幅依 s/p 偏振設定各分量；用 ContinuousSource(含 ramp) 或窄頻 GaussianSource+DFT。
   - 注意 Meep 平面波電流片源會同時往 ±y 發射，需在分析時只取關注方向的資料，或設計對應的
     TF/SF 等效比較窗（在腳本註解說明怎麼處理這個差異）。
   - 用 add_dft_fields/get_dft_array 取 DFT 場，get_array_metadata 取實際座標；用 add_flux 做
     normalization run 量 R/T；設 eps_averaging=False（或依 derivation.md 的介面平均策略決定）。
   - 若執行環境缺少 meep/mpi4py，腳本需在開頭做 import 檢查，缺少時印清楚的錯誤訊息並以非 0
     狀態碼結束，不得靜默跳過。

2. 撰寫 compare.py，讀本程式與 meep_ref.py 的輸出，做：
   (a) 真空 ky 比對，斷言相對誤差<1e-6，不符先檢查單位與 Courant 設定。
   (b) 波前角度與相位殘差比對（沿用 analyze.py 的擬合函式，不要重寫一份）。
   (c) Fresnel R/T 比對，斷言與 Meep 差<0.5%，且各自與解析解比較。
   (d) DFT 場逐點比對：對齊 Yee/Meep 座標（Meep 原點在 cell 中心，需要座標平移），用參考點歸一化
       振幅相位後，算真空中相對 L2 差，斷言<1e-3。
   (e) PML 反射：只記錄兩者數值，不設門檻，寫入報告。
   所有差異都要在報告輸出中標註可能來源（單位、網格偏移、介面平均、源型式、PML 差異）。

3. 撰寫 tests/run_all_regression.py，依序呼叫：
   level0_exact_injection.py → level0b_multistep.py → level1_leakage.py →
   analyze.py --level 1 --all → analyze.py --level 2 --all → compare.py --all
   任何一個失敗就停止並印出是哪一關、哪一項失敗，不得繼續跑後面的當作沒事。
```

**測試指令**：
```bash
python tests/run_all_regression.py
```

**DoD**：SPEC.md §10 關卡3 (a)-(d) 全部 PASS，(e) 已記錄（不要求 PASS/FAIL）；`run_all_regression.py`
全綠，代表 Stage 0–5 的既有結論在加入 Stage 6 程式碼後仍然成立。

**風險與回滾**：Meep 環境缺失時，(a)-(e) 標記為 BLOCKED 並寫明原因，但 regression 中 Stage 0–5 部分
仍必須全數 PASS 才能進 Stage 7；BLOCKED 的比對項目需在 Stage 7 的 validation_report.md 中誠實列出，
不得假裝通過。

---

### Stage 7 — 最終文件交付關

**目標**：彙整所有推導、驗證表與圖，產出最終 `validation_report.md`，完成 README 建置/執行說明，
確認全部交付物齊全可重現。

**涉及檔案**：`validation_report.md`（最終版）、`README.md`（最終版）、`figures/`（整理歸檔）。

**Claude Code instructions**：

```
1. 彙整 validation_report.md，結構：
   - 摘要（一段話說明整體結論：全部關卡是否通過，若有 BLOCKED 項目在最上方列出）。
   - derivation.md 全部內容原文收錄或引用（五節推導）。
   - 關卡0、關卡1（8項）、關卡2、關卡3（5項）各自的驗證表 + 對應圖，格式統一為
     「項目｜理論值(附單位與解析度)｜量測值(附單位與解析度)｜門檻｜PASS/FAIL/BLOCKED」。
   - 「已知限制與後續工作」小節：列出 Brewster 選做是否執行、Meep 環境限制、任何 BLOCKED 項目
     的根因分析。
2. 更新 README.md：專案簡介、目錄結構、建置指令（Windows/Linux）、如何重跑
   tests/run_all_regression.py、如何單獨重跑某一關。
3. 把 figures/ 底下的圖依關卡分類整理，確認 validation_report.md 中的圖片連結都有效
   （相對路徑，不得用絕對路徑）。
4. 最後跑一次 tests/run_all_regression.py 確認交付前狀態仍然全綠（或已記錄的 BLOCKED 狀態未變壞）。

嚴禁事項：不得在有任何一項標記 FAIL（非 BLOCKED）的情況下產出「結論：全部通過」的摘要。
```

**測試指令**：
```bash
python tests/run_all_regression.py
```

**DoD**：`validation_report.md` 涵蓋全部四關、所有表格與圖齊全，摘要結論與實際 PASS/FAIL/BLOCKED
狀態一致；`README.md` 可讓一個沒參與開發的人照著指令重新建置並重跑驗證。

**風險與回滾**：若彙整時發現先前某關的圖或數字對不上目前程式碼版本（例如中途改過 CPML 參數但沒
重新產圖），必須重新執行對應關卡重新產出，不得沿用舊圖硬套新結論。

