# 非均勻（tensor-product）Yee 網格推導

單位：c = ε0 = μ0 = 1、λ0 = 1、ω0 = 2π；ε、μ 為相對值。相量約定 f(t) = Re[F e^{−iωt}]。
本文件是 SPEC_nonuniform.md §12 的 5 項推導。每一節最後的「結論」條列可被測試的式子，測試檔對應寫在括號中。

## 0. 記號

- 每軸主節點：x_i、y_j、z_k。對偶節點取中點：x_{i+½} = (x_i + x_{i+1})/2，其餘兩軸同理。
- 主間距 h^a_{m+½} = a_{m+1} − a_m；對偶間距 d^a_m = a_{m+½} − a_{m−½} = (h^a_{m−½} + h^a_{m+½})/2。
- 階段 A：h^x ≡ d^x ≡ Δx，h^z ≡ d^z ≡ Δz，只有 y 非均勻。以下 y 軸簡寫 h_{j+½} ≡ h^y_{j+½}、d_j ≡ d^y_j。
- 分量位置（與現有程式 OFF 表一致）：

| 分量 | 位置 | 時間 |
|---|---|---|
| Ex | (x_{i+½}, y_j, z_k) | n |
| Ey | (x_i, y_{j+½}, z_k) | n |
| Ez | (x_i, y_j, z_{k+½}) | n |
| Hx | (x_i, y_{j+½}, z_{k+½}) | n+½ |
| Hy | (x_{i+½}, y_j, z_{k+½}) | n+½ |
| Hz | (x_{i+½}, y_{j+½}, z_k) | n+½ |

- 差分：δ_a F 表示沿 a 軸相鄰兩個樣本相減（前差或後差依位置而定，無歧義）。

## 1. 非均勻 Yee 更新式

**規則**：差分的結果落在某軸的對偶位置時，除以該軸的主間距 h；落在主位置時，除以對偶間距 d。
理由：兩個主節點樣本的差除以它們的距離 h，得到中點（對偶點）的二階精確導數；兩個對偶樣本的差除以它們的距離 d，得到主節點的導數。對偶節點取中點時，這是一階（主節點端）與二階（對偶節點端）局部截斷誤差；對應整體 supraconvergence（Monk & Süli 1994）。

μ ∂\_t H = −∇×E：

| 分量 | 更新式（Δt 乘右邊） | y 係數 | x/z 係數 |
|---|---|---|---|
| Hx | Hx −= (Δt/μ)[(Ez_{j+1} − Ez_j)/h_{j+½} − (Ey_{k+1} − Ey_k)/h^z_{k+½}] | h | h^z |
| Hy | Hy −= (Δt/μ)[(Ex_{k+1} − Ex_k)/h^z_{k+½} − (Ez_{i+1} − Ez_i)/h^x_{i+½}] | — | h^z、h^x |
| Hz | Hz −= (Δt/μ)[(Ey_{i+1} − Ey_i)/h^x_{i+½} − (Ex_{j+1} − Ex_j)/h_{j+½}] | h | h^x |

ε ∂\_t E = ∇×H：

| 分量 | 更新式 | 材料 | y 係數 | x/z 係數 |
|---|---|---|---|---|
| Ex | Ex += (Δt/ε_t[j])[(Hz_{j+½} − Hz_{j−½})/d_j − (Hy_{k+½} − Hy_{k−½})/d^z_k] | ε_t（主節點） | d | d^z |
| Ey | Ey += (Δt/ε_n[j+½])[(Hx_{k+½} − Hx_{k−½})/d^z_k − (Hz_{i+½} − Hz_{i−½})/d^x_i] | ε_n（對偶節點） | — | d^z、d^x |
| Ez | Ez += (Δt/ε_t[j])[(Hy_{i+½} − Hy_{i−½})/d^x_i − (Hx_{j+½} − Hx_{j−½})/d_j] | ε_t | d | d^x |

程式係數（SPEC §7.1）：cH_y[j] = Δt/h_{j+½}，cE_y[j] = Δt/(ε_t[j] d_j)，cE_tx[j] = Δt/(ε_t[j]Δx)，cE_tz[j] = Δt/(ε_t[j]Δz)，
cE_nx[j] = Δt/(ε_n[j]Δx)，cE_nz[j] = Δt/(ε_n[j]Δz)，cH_x = Δt/Δx，cH_z = Δt/Δz。

**材料平均**：介面只放在主節點 y_I。對 Ampère 式在 Ex 的對偶胞 [y_{I−½}, y_{I+½}] 上積分（切向 E 在介面連續），
∫ε dy / d_I = (h_{I−½}ε₁ + h_{I+½}ε₂)/(2d_I)，所以 ε_t[I] = (h_{I−½}ε₁ + h_{I+½}ε₂)/(h_{I−½} + h_{I+½})。
均勻時就是算術平均（舊 ifmode=a）。ε_n 位於對偶節點，永遠落在單一材料內部。

**CPML**：只在 y。∂\_y → (1/κ)∂\_y + ψ。ψ 的遞迴作用在差分（未除以間距），再與 y 差分一起乘 cE_y 或 cH_y：
ψ^{n} = b ψ^{n−1} + c·δ_y F，更新式中 y 項為 (1/κ)δ_yF + ψ。PML 區間距均勻（Δ_pml），剖面以物理距離 ρ 定義：
σ(ρ) = σ_max(ρ/d)^m、κ(ρ) = 1 + (κ_max − 1)(ρ/d)^m、α(ρ) = α_max(1 − ρ/d)，d = N_pml·Δ_pml，σ_max = sig_fac(m+1)/(η_loc Δ_pml)，η_loc = 1/n_loc。

**結論 1**（tests/nu_gate0_regression.py、nu_gate0_adjoint.py）
- 上兩表的 12 個差分項中，除以 h 的是：Hx 的 δyEz、Hz 的 δyEx，以及 x/z 版本的 δxEy、δxEz、δzEx、δzEy；除以 d 的是 Ex、Ez 的 δyH，以及 x/z 版本的 δxH、δzH。
- ε_t[I] = (h_{I−½}ε₁ + h_{I+½}ε₂)/(h_{I−½} + h_{I+½})。
- 均勻特例的更新式與舊碼只差浮點運算順序。

## 2. 散度、權重內積、伴隨性、能量守恆

**權重**（E：主邊長 × 對偶面積；H：對偶邊長 × 主面積）：

| 分量 | 權重 W |
|---|---|
| Ex | h^x_{i+½} · d_j · d^z_k |
| Ey | d^x_i · h_{j+½} · d^z_k |
| Ez | d^x_i · d_j · h^z_{k+½} |
| Hx | d^x_i · h_{j+½} · h^z_{k+½} |
| Hy | h^x_{i+½} · d_j · h^z_{k+½} |
| Hz | h^x_{i+½} · h_{j+½} · d^z_k |

內積 ⟨E, E'⟩\_{W_E} = Σ_分量 Σ_點 W·E·E'，H 同理。ε、μ 另外乘在內積中：‖E‖²_ε = ⟨εE, E⟩\_{W_E}。

**算子**：C_E = ∇\_h× 作用在 E、結果在 H 位置（§1 的 H 更新式中方括號，含正負號）；C_H = ∇\_h× 作用在 H、結果在 E 位置。

**伴隨性**：⟨E, C_H H⟩\_{W_E} = ⟨C_E E, H⟩\_{W_H}。
證明（以 Ex–Hz 的 y 項為例，其餘 11 項相同）：
- 左邊該項 = Σ Ex_{j}·(Hz_{j+½} − Hz_{j−½})/d_j · h^x d_j d^z = Σ h^x d^z Ex_j(Hz_{j+½} − Hz_{j−½})。
- 對 j 做分部求和 = −Σ h^x d^z Hz_{j+½}(Ex_{j+1} − Ex_j) + 邊界項。
- 右邊對應項：(C_E E)\_{Hz} 中的 −(Ex_{j+1} − Ex_j)/h_{j+½}，權重 h^x h_{j+½} d^z，得 −Σ h^x d^z Hz_{j+½}(Ex_{j+1} − Ex_j)。兩邊相等。
- 關鍵：E 權重中的 d_j 與差分的 1/d_j 相消；H 權重中的 h_{j+½} 與 1/h_{j+½} 相消，兩者剩下相同的 h^x d^z。
- 邊界項：y 方向 PEC（Ex、Ez 在 j = 0、Ny 為 0）使邊界項為 0；x、z 方向 PBC 使求和循環、無邊界項。

**能量守恆**：蛙跳格式 E^{n+1} = E^n + Δt ε^{−1}C_H H^{n+½}，H^{n+½} = H^{n−½} − Δt μ^{−1}C_E E^n。定義
𝓔^n = ½‖E^n‖²_ε + ½⟨μH^{n+½}, H^{n−½}⟩\_{W_H}。則
𝓔^{n+1} − 𝓔^n = ½⟨ε(E^{n+1} − E^n), E^{n+1} + E^n⟩ + ½⟨μH^{n+½}, H^{n+3/2} − H^{n−½}⟩
= ½Δt⟨C_H H^{n+½}, E^{n+1} + E^n⟩\_{W_E} − ½Δt⟨H^{n+½}, C_E(E^{n+1} + E^n)⟩\_{W_H} = 0（由伴隨性）。

**正定條件**：H^{n−½} = H^{n+½} + Δt μ^{−1}C_E E^n，所以 𝓔^n = ½‖E^n‖²_ε + ½‖H^{n+½}‖²_μ + ½Δt⟨H^{n+½}, C_E E^n⟩。
令 A = μ^{−½}C_E ε^{−½}（在加權空間中），‖A‖² = λ_max(M)，M = ε^{−1}C_H μ^{−1}C_E。
由 Cauchy–Schwarz，𝓔^n ≥ ½(1 − Δt‖A‖/2)(‖E‖² + ‖H‖²)。所以 Δt < 2/sqrt(λ_max) 時 𝓔 正定；這也是蛙跳格式的精確穩定界：
Δt_max = 2/sqrt(λ_max(M))。

**最小間距公式是充分條件**：
- tensor-product 網格上，不同軸的差分算子作用在不同索引，彼此可交換。
- 因此離散恆等式 C_H C_E = −L_h + G_h D_h 成立。L_h 是逐分量的非均勻 Laplacian；D_h 是散度，G_h = −D_h^† 是梯度。
- 所以 C_H C_E = −L_h − D_h^†D_h ≤ −L_h（在加權內積下）。
- −L_h 是三個 1D 算子 −(1/w)δ(δ/g) 之和。對稱化後，每一列的對角元素 ≤ 2/Δ_min²、兩個非對角元素各 ≤ 1/Δ_min²（因 h, d ≥ Δ_min）。Gershgorin 給出 λ_max(−L_h) ≤ Σ_a 4/Δ_min,a²。
- 真空或 ε, μ ≥ 1 時，λ_max(M) ≤ 4Σ_a 1/Δ_min,a²。所以 Δt ≤ 1/sqrt(Σ_a 1/Δ_min,a²) ⇒ Δt ≤ Δt_max。

**離散散度**：
- div E 在主節點 (i, j, k)：(Ex_{i+½} − Ex_{i−½})/d^x_i + (Ey_{j+½} − Ey_{j−½})/d_j + (Ez_{k+½} − Ez_{k−½})/d^z_k。介質中改用 ε 乘各分量（div D）。
- div H 在胞心 (i+½, j+½, k+½)：(Hx_{i+1} − Hx_i)/h^x_{i+½} + (Hy_{j+1} − Hy_j)/h_{j+½} + (Hz_{k+1} − Hz_k)/h^z_{k+½}。
- 由 D_h C_H = 0、D_h C_E = 0（混合差分可交換），無源區保持散度不變。

**結論 2**（nu_gate0_adjoint.py、nu_gate0_energy.py、nu_gate0_stability.py、nu_gate0_divergence.py）
- ⟨E, C_H H⟩\_{W_E} = ⟨C_E E, H⟩\_{W_H}，以上表權重；全 1 權重時不成立。
- 𝓔^n = ½Σ εW_E|E^n|² + ½Σ μW_H H^{n+½}·H^{n−½} 在 PEC/PBC 無源腔中精確守恆。
- Δt_max = 2/sqrt(λ_max(M))；Δt_Courant = 1/sqrt(Σ 1/Δ_min²) ≤ Δt_max。
- 散度算子如上；以 d_j 換成 Δ 的「均勻算子」在漸變區不為 0。

## 3. aux line 是精確約化

令 3D 場 F_c(x, y, z, t) = Re[f_c(y, t) e^{i(kx x_c + kz z_c)}]，x_c、z_c 為分量 c 自身的 x、z 位置，kx Lx、kz Lz ∈ 2πℤ（滿足 PBC）。

- **x 方向**：Δx 均勻，所以 (e^{ikx x_{i+1}} − e^{ikx x_i})/Δx = e^{ikx x_{i+½}}·(2i/Δx) sin(kxΔx/2) = iK̃x e^{ikx x_{i+½}}。對偶差分也得到同一個 iK̃x。
  所以對這個函數族，δ_x/Δx 等於「乘以常數 iK̃x」。z 方向同理（iK̃z）。
- **y 方向**：所有係數（h、d、ε_t、ε_n、CPML 的 b、c、κ）只依賴 y，不依賴 x、z，所以 y 差分完全不受影響。
- **線性、實係數**：若複數場滿足格式，它的實部也滿足。
- **結論**：3D 格式對這族場精確約化成 y 向 1D 複數格式，對任意 y 節點陣列成立。x、z 導數換成 iK̃x、iK̃z 並乘 Δt；y 差分使用局部 Δt/h、Δt/d。

**x 非均勻時失效**：(e^{ikx x_{i+1}} − e^{ikx x_i})/h^x_{i+½} = iK̃(h^x_{i+½}) e^{ikx x_{i+½}}，其中 K̃(h) = (2/h) sin(kx h/2) 隨 i 改變，
所以 e^{ikx x} 不再是 δ_x 的特徵函數，無法抽出共同因子；相量化 aux line 不成立（階段 B 必須中止）。

**aux 節點範圍**：aux 在 global j ∈ [0, Ny] 使用主網格節點，j < 0、j > Ny 以端點間距均勻延伸。
TF/SF 修正只讀 aux 在 j0 與 j0 − ½ 的值，所以 aux 至少要在 [j_a, j0 + 1] 與主網格相同（prompt 的最低要求）。
在真空中複製到 Ny 則讓 aux 對主網格的 TF 區也精確。aux_ref 線再加上主網格的 ε 與 CPML，對整個主網格精確。

**結論 3**（nu_gate1_vacuum.py 1-0、1-1、1-2）
- 真空、x/z 均勻：F_main = Re[aux_ref · e^{i(kx x_c + kz z_c)}]，除捨入外相等（< 1e-12）。
- 每個 y 平面上相位對 kx x + kz z 精確線性（殘差只有捨入）。
- aux 更新式：y 差分乘 Δt/h（H）或 Δt/d（E），相量項乘 Δt·iK̃x、Δt·iK̃z，不可共用 c。

## 4. TF/SF 修正係數

主網格 j0（TF 側第一個主節點）：
- H 側：Hx、Hz 在 j0 − ½ 屬於 SF。它們的 y 差分讀到 j0 的總場 E，必須減去入射場：
  Hx[j0−1] += cH_y[j0−1]·Ez_inc(j0)，Hz[j0−1] −= cH_y[j0−1]·Ex_inc(j0)。係數 Δt/h_{j0−½}。
- E 側：Ex、Ez 在 j0 屬於 TF。它們的 y 差分讀到 j0 − ½ 的散射場 H，必須加上入射場：
  Ex[j0] −= cE_y[j0]·Hz_inc(j0−½)，Ez[j0] += cE_y[j0]·Hx_inc(j0−½)。係數 Δt/(ε_t[j0] d_{j0})。
- 入射值取自 aux line（與主網格同節點、同 Δt、同 K̃），所以修正對任意間距都精確。SPEC 仍要求 j0 ± 5 均勻（保守規則，也方便量測）。

aux 源（j_a，1D TF/SF）：
- 入射值是解析離散平面波 inc_amp，只在均勻區是離散方程的精確解。所以 j_a ± 5 必須均勻（間距 Δ_a）。
- ky 用 ky(Δ_a)，y 座標用 aux 節點座標，修正係數用 Δt/Δ_a。

**結論 4**（nu_gate1_vacuum.py 1-1）：上式四個修正係數；aux 源的 ky = ky_local(Δ_a)。

## 5. 數值色散、角度、反射

**均勻子區**（間距 h、介電常數 ε）：
(2/h)² sin²(ky h/2) = ω̃² εμ − K̃t²，其中 ω̃ = (2/Δt) sin(ωΔt/2)，K̃t² = K̃x² + K̃z²，K̃x = (2/Δx) sin(kxΔx/2)。
所以 ky(h) = (2/h) asin((h/2) sqrt(ω̃²εμ − K̃t²))。K̃y ≡ (2/h) sin(ky h/2) = sqrt(ω̃²εμ − K̃t²) 與 h 無關。

**局部角度**：θ(y) = atan(|K̃t|/ky_local)；連續角 θ_cont = atan(kt/ky_cont)，ky_cont = sqrt(k0²ε − kt²)。
θ_pred(h) − θ_cont = O(h²) + O(Δx², Δz², Δt²)。細化序列固定 x/z 時，對 h 的部分是 O(h²)。

**局部前行波分解**（1-3b 的量測方法）：
- 任一在同一節點族上取樣、滿足均勻子區離散方程的場 u_j = a e^{iκy_j} + b e^{−iκy_j}，在對偶點 j+½ 定義
  A_{j+½} = ½[(u_j + u_{j+1})/(2cos(κh/2)) + (u_{j+1} − u_j)/(iK̃y h)]。
- 前行波兩項都等於 a e^{iκy_{j+½}}；反向波兩項相消。所以 A 只含前行波，相位不受反射影響。
- κ 只出現在實數因子 cos(κh/2)，K̃y 在各子區相同。
- 量測 ky_meas(j) = arg(A_{j+½}/A_{j−½})/d_j，h_loc = d_j；θ_meas = atan(|K̃t|/ky_meas)。

**離散約化散射問題**（R_disc，tmm.tmm_discrete）：
- 由 §3，給定 (kx, kz) 後 3D 格式等價於 1D。在 x–z 平面把座標轉到 t̂ = K̃t/|K̃t| 方向後分成兩組：s（e_s、h_t、h_y）與 p（h_s、e_t、e_y）。
- 消去代數分量後，兩組都寫成三點形式 (1/w_m)[(u_{m+1} − u_m)/g_{m+½} − (u_m − u_{m−1})/g_{m−½}] + q_m u_m = 0：

| 偏振 | u（位置） | w_m | g_{m+½} | q_m |
|---|---|---|---|---|
| s | e_s（主節點 j） | d_j | h_{j+½} | ω̃² ε_t[j] − K̃t² |
| p | h_s（對偶節點 j+½） | h_{j+½} | ε_t[j+1] d_{j+1} | ω̃² − K̃t²/ε_n[j+½] |

- 均勻區兩者都回到上面的色散式。
- 守恆通量（離散 Wronskian）：Φ = Im(ū_m(u_{m+1} − u_m))/g_{m+½}，對所有 m 相同。所以
  R = |r|²，T = |t|²·[sin(k_R h_R)/g_R]/[sin(k_L h_L)/g_L]，R + T = 1（無損）。
- 兩端為均勻區，以精確離散模態封閉：入射端 u_{−1} = u_0 e^{iκh} − 2i sin(κh)（入射振幅 1），出射端 u_{N+1} = u_N e^{iκh}。三對角系統以 Thomas 演算法求解。

**數值 Fresnel 估計（INFO）**：
- n_num,i = sqrt(kt² + ky_i²)/k0，代入 Fresnel：s 為 r = (ky₁ − ky₂)/(ky₁ + ky₂)，p 為 r = (ky₁/n₁² − ky₂/n₂²)/(ky₁/n₁² + ky₂/n₂²)。
- 預研究（SPEC 研究 checkpoint）顯示，突變時 |R_disc|/|R_估計| ≈ 3。多出的部分來自接面處一階截斷誤差（d_j ≠ h），不是相速差。

**結論 5**（nu_predict.py、nu_gate1_vacuum.py 1-3、1-4；nu_gate2_film.py）
- ky(h) 與 θ(h) 的公式如上；子區量測以三點恆等式 cos(ky h) = (u_{j+1} + u_{j−1})/(2u_j)。
- R_disc、T_disc 由上表三點問題精確求得，均勻特例等於舊 discrete_fresnel（T-TMM1）。
- 漸變區 ky 以局部前行波分解量測。
