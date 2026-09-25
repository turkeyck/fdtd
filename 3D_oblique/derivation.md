# 推導文件（derivation.md）

本文件是 SPEC.md「寫碼之前先輸出推導」的正式產出。§1–§5 對應 SPEC 要求的五項推導；§6 起為實作
過程中補上的設計推導（入射場時間包絡、CPML、介面 ε、通量定義），各節標明新增於哪個 Stage。

**單位**：全專案採正規化單位 $`c=\varepsilon_0=\mu_0=1`$（因此 $`\eta_0=1`$），$`\lambda_0=1`$，
故 $`\omega_0=k_0=2\pi`$、週期 $`T_0=1`$。預設 $`\Delta=\lambda_0/20`$，$`S=c\Delta t/\Delta=0.5`$，
$`\Delta t=\Delta/2=1/40`$，每週期 $`T_0/\Delta t=N_\lambda/S=40`$ 步（$`N_\lambda=\lambda_0/\Delta`$；$`S=0.5`$ 時對 $`N_\lambda=10,20,40`$ 皆為整數，便於整數週期 DFT）。

---

## 1. Yee 分量座標表

格點整數索引 $`(i,j,k)`$，實體座標 $`(x,y,z)=(i\Delta, j\Delta, k\Delta)`$；時間索引 $`n`$。

| 分量 | 空間位置 | 時間 | 陣列索引範圍 |
|---|---|---|---|
| $`E_x`$ | $`(i+\tfrac12,\ j,\ k)`$ | $`n\Delta t`$ | $`i\in[0,N_x),\ j\in[0,N_y],\ k\in[0,N_z)`$ |
| $`E_y`$ | $`(i,\ j+\tfrac12,\ k)`$ | $`n\Delta t`$ | $`j\in[0,N_y)`$ |
| $`E_z`$ | $`(i,\ j,\ k+\tfrac12)`$ | $`n\Delta t`$ | $`j\in[0,N_y]`$ |
| $`H_x`$ | $`(i,\ j+\tfrac12,\ k+\tfrac12)`$ | $`(n+\tfrac12)\Delta t`$ | $`j\in[0,N_y)`$ |
| $`H_y`$ | $`(i+\tfrac12,\ j,\ k+\tfrac12)`$ | $`(n+\tfrac12)\Delta t`$ | $`j\in[0,N_y]`$ |
| $`H_z`$ | $`(i+\tfrac12,\ j+\tfrac12,\ k)`$ | $`(n+\tfrac12)\Delta t`$ | $`j\in[0,N_y)`$ |

**為什麼是這個偏移**：Ampère 定律 $`\varepsilon\,\partial_t E_x = \partial_y H_z-\partial_z H_y`$ 要以
二階中心差分近似，$`E_x`$ 必須位於 $`\partial_y H_z`$ 與 $`\partial_z H_y`$ 兩個差分的中點：
$`H_z`$ 需在 $`E_x`$ 的 $`y\pm\tfrac12`$、$`H_y`$ 需在 $`E_x`$ 的 $`z\pm\tfrac12`$。對六個分量同時要求「每個
旋度分量都是其兩個鄰居的中心差分」，唯一的自洽解（差一個整體平移）就是：$`E_\alpha`$ 只在自己的
方向 $`\alpha`$ 偏移半格；$`H_\alpha`$ 在另外兩個方向各偏移半格。時間上 Faraday/Ampère 交錯（蛙跳），
$`E`$ 在整數步、$`H`$ 在半整數步，使時間導數也是中心差分。整個格式對空間與時間都是二階精度。

**y 方向邊界**：$`j=0`$ 與 $`j=N_y`$ 兩個外表面放 PEC（切向 $`E_x,E_z\equiv0`$，不更新），其內側是
CPML；$`H_y`$ 位於整數 $`j`$ 但只含 $`x,z`$ 差分，所以在 $`j=0,N_y`$ 仍可正常更新。

**Stage 0 驗證**：若本表任一偏移錯誤，§2 的解析場代入一步更新後殘差不可能降到 $`10^{-15}`$ 量級。
關卡 0 實測殘差 $`\approx 2\times10^{-15}`$（見 `results/level0.json`），間接證明本表正確。

---

## 2. 離散平面波（六分量精確解）

### 2.1 差分算子作用在指數上的本徵值

對 $`f=e^{i k_x x}`$，Yee 中心差分（間隔 $`\Delta`$、中心在 $`x`$）：

```math
\frac{f(x+\tfrac\Delta2)-f(x-\tfrac\Delta2)}{\Delta}= i\,\underbrace{\frac{2}{\Delta}\sin\frac{k_x\Delta}{2}}_{\tilde K_x}\,e^{ik_x x}.
```

時間上對 $`e^{-i\omega t}`$：$`\dfrac{g(t+\tfrac{\Delta t}{2})-g(t-\tfrac{\Delta t}{2})}{\Delta t}=-i\tilde\omega\,e^{-i\omega t}`$，
$`\tilde\omega=\dfrac{2}{\Delta t}\sin\dfrac{\omega\Delta t}{2}`$。

因此把 $`\mathbf F=\mathrm{Re}\{\mathbf F_0\,e^{i(\mathbf k\cdot\mathbf r-\omega t)}\}`$（每個分量在自己的
Yee 點與時間取樣）代入離散 Maxwell，等同於連續 Maxwell 中把 $`\nabla\to i\tilde{\mathbf K}`$、
$`\partial_t\to -i\tilde\omega`$：

```math
\tilde\omega\,\mu_0\mathbf H_0=\tilde{\mathbf K}\times\mathbf E_0,\qquad
-\tilde\omega\,\varepsilon\,\mathbf E_0=\tilde{\mathbf K}\times\mathbf H_0 .
```

**注意**：相位裡用的是真實波數 $`\mathbf k=(k_x,k_y,k_z)`$ 與 $`\omega_0`$；$`\tilde{\mathbf K}`$、$`\tilde\omega`$
只出現在振幅／偏振關係中。

### 2.2 色散關係與 $`k_y`$ 閉式解

兩式相消並用 $`\tilde{\mathbf K}\cdot\mathbf E_0=0`$，得 $`|\tilde{\mathbf K}|^2=\mu_0\varepsilon\,\tilde\omega^2`$，即

```math
\Big[\frac{\sin(\omega_0\Delta t/2)}{c_m\Delta t}\Big]^2=\sum_{i=x,y,z}\Big[\frac{\sin(k_i\Delta/2)}{\Delta}\Big]^2,\qquad c_m=c/n .
```

已知 $`k_x=2\pi m/L_x`$、$`k_z=2\pi n/L_z`$：

```math
Q\equiv\sin^2\frac{k_y\Delta}{2}=\Big(\frac{n}{S}\Big)^2\sin^2\frac{\omega_0\Delta t}{2}-\sin^2\frac{k_x\Delta}{2}-\sin^2\frac{k_z\Delta}{2},\qquad
k_y=\frac{2}{\Delta}\arcsin\sqrt{Q}.
```

程式在啟動時檢查 $`0\lt Q\le1`$（$`Q\le0`$：該模態在此網格上為消逝波；$`Q\gt 1`$：無實數解），不合格直接報錯
結束，不產生 NaN。取正根，代表往 $`+y`$ 傳播。另外檢查連續傳播條件 $`|k_t|\lt (1-0.05)\,n\,k_0`$（保留 5% 裕度，
避免掠射角附近群速度趨近 0）。

**預設參數數值**（$`m=n=1`$，$`L_x=2`$，$`L_z=3`$，$`\Delta=\lambda_0/20`$，$`S=0.5`$）：

| 量 | 值 |
|---|---|
| $`k_x/k_0,\ k_z/k_0`$ | 0.5, 1/3 |
| $`\lvert k_t\rvert/k_0`$ | 0.60093（裕度 40%） |
| $`\theta=\mathrm{atan2}(\lvert k_t\rvert,k_y)`$ | 36.895°；$`\varphi=\mathrm{atan2}(k_z,k_x)=33.690°`$ |
| $`k_y`$（離散） | 5.0297668864 /λ0 |
| $`k_y`$（連續 $`\sqrt{k_0^2-k_t^2}`$） | 5.0221830272 /λ0 |
| $`\tilde{\mathbf K}`$ | (3.13836383, 5.01652259, 2.09343825) |
| $`\lvert\tilde{\mathbf K}\rvert=\tilde\omega`$ | 6.2767276582 |

### 2.3 偏振基底

```math
\hat s=\frac{\tilde{\mathbf K}\times\hat y}{|\tilde{\mathbf K}\times\hat y|}=\frac{(-\tilde K_z,\,0,\,\tilde K_x)}{\tilde K_t},\qquad
\hat p=\frac{\hat s\times\tilde{\mathbf K}}{|\hat s\times\tilde{\mathbf K}|}=\frac{(-\tilde K_x\tilde K_y,\ \tilde K_t^2,\ -\tilde K_z\tilde K_y)}{\tilde K_t\,|\tilde{\mathbf K}|},
```

其中 $`\tilde K_t=\sqrt{\tilde K_x^2+\tilde K_z^2}`$。兩者皆與 $`\tilde{\mathbf K}`$ 正交，故 $`\tilde{\mathbf K}\cdot\mathbf E_0=0`$
（離散 Gauss 定律）。**退化情形**：$`m=n=0`$（正入射）時 $`\tilde K_t=0`$，基底無定義；程式直接報錯拒絕，
不支援此情形。

```math
\mathbf E_0=E_0\,\hat e,\ \ \hat e\in\{\hat s,\hat p\},\qquad \mathbf H_0=\frac{\tilde{\mathbf K}\times\mathbf E_0}{\mu_0\tilde\omega}.
```

預設（$`E_0=1`$）：s：$`\mathbf E_0=(-0.55492,0,0.83190)`$，$`\mathbf H_0=(0.66488,-0.60103,0.44351)`$；
p：$`\mathbf E_0=(-0.66488,0.60103,-0.44351)`$，$`\mathbf H_0=(-0.55492,0,0.83190)`$。

### 2.4 六分量明確表達式

令 $`\phi(x,y,z,t)=k_x x+k_y y+k_z z-\omega_0 t+\phi_0`$（本專案 $`\phi_0=0`$），則

```math
\begin{aligned}
E_x\big|_{i+\frac12,j,k}^{n}&=E_{0x}\cos\phi\big((i+\tfrac12)\Delta,\ j\Delta,\ k\Delta,\ n\Delta t\big)\\
E_y\big|_{i,j+\frac12,k}^{n}&=E_{0y}\cos\phi\big(i\Delta,\ (j+\tfrac12)\Delta,\ k\Delta,\ n\Delta t\big)\\
E_z\big|_{i,j,k+\frac12}^{n}&=E_{0z}\cos\phi\big(i\Delta,\ j\Delta,\ (k+\tfrac12)\Delta,\ n\Delta t\big)\\
H_x\big|_{i,j+\frac12,k+\frac12}^{n+\frac12}&=H_{0x}\cos\phi\big(i\Delta,\ (j+\tfrac12)\Delta,\ (k+\tfrac12)\Delta,\ (n+\tfrac12)\Delta t\big)\\
H_y\big|_{i+\frac12,j,k+\frac12}^{n+\frac12}&=H_{0y}\cos\phi\big((i+\tfrac12)\Delta,\ j\Delta,\ (k+\tfrac12)\Delta,\ (n+\tfrac12)\Delta t\big)\\
H_z\big|_{i+\frac12,j+\frac12,k}^{n+\frac12}&=H_{0z}\cos\phi\big((i+\tfrac12)\Delta,\ (j+\tfrac12)\Delta,\ k\Delta,\ (n+\tfrac12)\Delta t\big)
\end{aligned}
```

（$`\mathbf E_0,\mathbf H_0`$ 為實向量，六分量同相；程式實作見 `fdtd_theory.analytic()`。）

### 2.5 由此推出的恆等式（驗證用）

1. **離散 Gauss**：節點 $`(i,j,k)`$ 的 $`\nabla_h\!\cdot\mathbf E = i\tilde{\mathbf K}\cdot\mathbf E_0\,e^{i\phi}=0`$；
   胞心的 $`\nabla\cdot\mathbf H=i\tilde{\mathbf K}\cdot\mathbf H_0 e^{i\phi}=0`$（因 $`\mathbf H_0\propto\tilde{\mathbf K}\times\mathbf E_0`$）。
2. **$`\mathbf E\cdot\mathbf H=0`$**：$`\mathbf E_0\cdot(\tilde{\mathbf K}\times\mathbf E_0)=0`$。
3. **阻抗（對 SPEC 的更正）**：$`|\mathbf H_0|/|\mathbf E_0|=|\tilde{\mathbf K}|/(\mu_0\tilde\omega)`$，而色散關係正是
   $`|\tilde{\mathbf K}|=\tilde\omega/c_m`$，所以

   $`\displaystyle \frac{|\mathbf H_0|}{|\mathbf E_0|}=\frac{1}{\mu_0 c_m}=\frac{1}{\eta}\quad\text{（恆等，非近似）}.`$

   SPEC 寫「與 $`1/\eta_0`$ 的差距應為 $`O((k_0\Delta)^2)`$」只對「內插到同一點後」的振幅成立（內插引入
   $`\cos(k_d\Delta/2)`$ 因子，見 §5 L1-5）；對各分量在自身 Yee 點的原始振幅，差距理論值是 0。
   關卡 1 兩者都量，理論值各自依此設定，不修改門檻。
4. **能流方向**：原始振幅的 $`\tfrac12\mathbf E_0\times\mathbf H_0=\tilde{\mathbf K}|\mathbf E_0|^2/(2\mu_0\tilde\omega)`$ 平行
   $`\tilde{\mathbf K}`$；Yee 群速度 $`v_{g,i}=c_m^2\Delta t\sin(k_i\Delta)/(\Delta\sin\omega_0\Delta t)`$ 平行
   $`(\sin k_x\Delta,\sin k_y\Delta,\sin k_z\Delta)`$。預設參數下 $`\angle(\tilde{\mathbf K},\mathbf k)=0.0499^\circ`$、
   $`\angle(\mathbf v_g,\mathbf k)=0.2003^\circ`$，兩者都是 $`O((k\Delta)^2)`$。

---

## 3. TF/SF 修正（$`y=y_0=j_0\Delta`$ 單一平面）

**區域劃分**：總場（TF）＝ $`y\ge y_0`$ 的所有節點；散射場（SF）＝ $`y\lt y_0`$。因此
$`E_x,E_z,H_y`$（整數 $`j`$）：$`j\ge j_0`$ 為 TF；$`E_y,H_x,H_z`$（$`j+\tfrac12`$）：$`j\ge j_0`$ 為 TF、
$`j\le j_0-1`$（即 $`y_0-\tfrac\Delta2`$）為 SF。

只有「更新式跨越 $`y_0`$ 的 $`y`$ 差分」需要修正。$`E_y`$、$`H_y`$ 的更新式只含 $`x,z`$ 差分（同一個 $`y`$），
不跨界，**不需修正**。其餘四個：

**(a) $`H_x`$ at $`(i,\,j_0-\tfrac12,\,k+\tfrac12)`$（SF 節點）**：
$`H_x^{n+\frac12}=H_x^{n-\frac12}-\frac{\Delta t}{\mu_0}\big[\frac{E_z|_{j_0}-E_z|_{j_0-1}}{\Delta}-\partial_zE_y\big]`$。
SF 節點需要的是散射場 $`E_z^{\rm scat}|_{j_0}=E_z|_{j_0}-E_z^{\rm inc}|_{j_0}`$，陣列存的是總場，所以

```math
H_x|_{j_0-\frac12}^{n+\frac12}\ \mathrel{+}=\ \frac{\Delta t}{\mu_0\Delta}\,E_z^{\rm inc}\big|_{i,\,j_0,\,k+\frac12}^{\,n}.
```

**(b) $`H_z`$ at $`(i+\tfrac12,\,j_0-\tfrac12,\,k)`$（SF）**：
$`H_z^{n+\frac12}=H_z^{n-\frac12}-\frac{\Delta t}{\mu_0}\big[\partial_xE_y-\frac{E_x|_{j_0}-E_x|_{j_0-1}}{\Delta}\big]`$，
把 $`E_x|_{j_0}`$ 換成 $`E_x|_{j_0}-E_x^{\rm inc}`$：

```math
H_z|_{j_0-\frac12}^{n+\frac12}\ \mathrel{-}=\ \frac{\Delta t}{\mu_0\Delta}\,E_x^{\rm inc}\big|_{i+\frac12,\,j_0,\,k}^{\,n}.
```

**(c) $`E_x`$ at $`(i+\tfrac12,\,j_0,\,k)`$（TF）**：
$`E_x^{n+1}=E_x^{n}+\frac{\Delta t}{\varepsilon}\big[\frac{H_z|_{j_0+\frac12}-H_z|_{j_0-\frac12}}{\Delta}-\partial_zH_y\big]`$，
TF 節點需要總場 $`H_z|_{j_0-\frac12}+H_z^{\rm inc}`$：

```math
E_x|_{j_0}^{n+1}\ \mathrel{-}=\ \frac{\Delta t}{\varepsilon\Delta}\,H_z^{\rm inc}\big|_{i+\frac12,\,j_0-\frac12,\,k}^{\,n+\frac12}.
```

**(d) $`E_z`$ at $`(i,\,j_0,\,k+\tfrac12)`$（TF）**：
$`E_z^{n+1}=E_z^{n}+\frac{\Delta t}{\varepsilon}\big[\partial_xH_y-\frac{H_x|_{j_0+\frac12}-H_x|_{j_0-\frac12}}{\Delta}\big]`$：

```math
E_z|_{j_0}^{n+1}\ \mathrel{+}=\ \frac{\Delta t}{\varepsilon\Delta}\,H_x^{\rm inc}\big|_{i,\,j_0-\frac12,\,k+\frac12}^{\,n+\frac12}.
```

**正負號的規律**：修正項等於「該差分中跨界鄰居的係數 × 入射場」並帶上把總場↔散射場換算所需的符號：
SF 節點用到 TF 鄰居時減去入射場，TF 節點用到 SF 鄰居時加上入射場；再乘上原更新式中該鄰居前面的
符號（$`\partial_y`$ 前的正負與 $`\pm1/\Delta`$），就得到上面 $`+,-,-,+`$。

**時間取樣**：(a)(b) 在 H 更新時使用 $`E^{\rm inc}`$ 在 $`n\Delta t`$；(c)(d) 在 E 更新時使用 $`H^{\rm inc}`$ 在
$`(n+\tfrac12)\Delta t`$——與被修正的更新式所使用的場同一時刻。

**充要條件**：若 $`\mathbf F^{\rm inc}`$ 在這四個跨界節點的更新模板上**逐步**滿足離散齊次 Maxwell
方程，則散射場的源項恆為 0，SF 區場精確為 0（洩漏 = 捨入誤差）。這是 §6 設計的出發點。

---

## 4. PBC 索引（x、z 週期）

陣列 $`i\in[0,N_x)`$、$`k\in[0,N_z)`$。差分有兩種方向：

- **向後差分**（目標點在整數位置、來源在 $`+\tfrac12`$ 位置）：$`(F[i]-F[i-1])/\Delta`$，wrap $`F[-1]\equiv F[N_x-1]`$。
- **向前差分**（目標點在 $`+\tfrac12`$ 位置、來源在整數位置）：$`(F[i+1]-F[i])/\Delta`$，wrap $`F[N_x]\equiv F[0]`$。

逐項列出（z 方向同理，以 $`N_z`$ 為週期）：

| 更新 | x 差分 | 在 $`i`$ 邊界取值 | z 差分 | 在 $`k`$ 邊界取值 |
|---|---|---|---|---|
| $`E_x(i+\frac12,j,k)`$ | — | — | $`H_y[k]-H_y[k-1]`$ | $`k=0`$ 用 $`H_y[N_z-1]`$ |
| $`E_y(i,j+\frac12,k)`$ | $`H_z[i]-H_z[i-1]`$ | $`i=0`$ 用 $`H_z[N_x-1]`$ | $`H_x[k]-H_x[k-1]`$ | $`k=0`$ 用 $`H_x[N_z-1]`$ |
| $`E_z(i,j,k+\frac12)`$ | $`H_y[i]-H_y[i-1]`$ | $`i=0`$ 用 $`H_y[N_x-1]`$ | — | — |
| $`H_x(i,j+\frac12,k+\frac12)`$ | — | — | $`E_y[k+1]-E_y[k]`$ | $`k=N_z-1`$ 用 $`E_y[0]`$ |
| $`H_y(i+\frac12,j,k+\frac12)`$ | $`E_z[i+1]-E_z[i]`$ | $`i=N_x-1`$ 用 $`E_z[0]`$ | $`E_x[k+1]-E_x[k]`$ | $`k=N_z-1`$ 用 $`E_x[0]`$ |
| $`H_z(i+\frac12,j+\frac12,k)`$ | $`E_y[i+1]-E_y[i]`$ | $`i=N_x-1`$ 用 $`E_y[0]`$ | — | — |

以 SPEC 舉的例子：$`i=0`$ 處 $`E_y`$ 更新中的 $`\partial H_z/\partial x=(H_z[0]-H_z[N_x-1])/\Delta`$，
$`H_z[N_x-1]`$ 位於 $`x=(N_x-\tfrac12)\Delta\equiv-\tfrac12\Delta\pmod{L_x}`$，正是 $`x=0`$ 左側半格。

**精確性條件**：$`k_xL_x=2\pi m`$、$`k_zL_z=2\pi n`$ ⇒ $`e^{ik_x(x+L_x)}=e^{ik_xx}`$，解析場在網格上嚴格週期，
wrap 不引入任何誤差；因此場可用實數（不需要 Bloch 相位）。C 程式以 ghost 層實作：每次 H 更新前把
$`E[N_x]\leftarrow E[0]`$、$`E[\cdot][\cdot][N_z]\leftarrow E[\cdot][\cdot][0]`$，E 更新前把
$`H[-1]\leftarrow H[N_x-1]`$、$`H[\cdot][\cdot][-1]\leftarrow H[\cdot][\cdot][N_z-1]`$，與上表等價。

**Stage 0 驗證**：關卡 0 分開統計邊界格（$`i\in\{0,N_x-1\}`$ 或 $`k\in\{0,N_z-1\}`$）與內部格殘差，
兩者同為 $`\sim2\times10^{-15}`$，無系統性偏大。

---

## 5. 各驗證項目的理論值公式

記號：$`\mathbf k_d=(k_x,k_y^{\rm disc},k_z)`$，$`k_y^{\rm cont}=\sqrt{k_0^2-k_t^2}`$。所有「相對」量的分母在各項中註明。

### 關卡 0
| 項目 | 理論值 | 說明 |
|---|---|---|
| $`\max\lvert E^{\rm upd}-E^{\rm an}\rvert/\max\lvert E\rvert`$ | 0（捨入 $`\sim10^{-15}`$） | §2.1 本徵值關係 |
| 同上（H） | 0 | 同上 |
| $`\max\lvert\nabla_h\!\cdot\mathbf E\rvert/(\lvert\tilde{\mathbf K}\rvert\max\lvert E\rvert)`$ | 0 | §2.5-1；分母是單一差分項的量級 |
| 對照 C1（只把 $`k_y`$ 換成連續值，$`\tilde{\mathbf K}`$ 由連續 $`k_y`$ 算） | E 殘差 $`=\Delta t\,\lvert\tilde\omega^2-\lvert\tilde{\mathbf K}_c\rvert^2\rvert/\tilde\omega`$（相對）；H 殘差 $`=0`$（$`\mathbf H_0`$ 由同一個 $`\tilde{\mathbf K}_c`$ 定義，Faraday 自動滿足） | 閉式解，程式逐一比對 |
| 對照 C2（完全連續平面波 $`\mathbf E_0\perp\mathbf k`$、$`\mathbf H_0=\mathbf k\times\mathbf E_0/\mu_0\omega`$） | E、H 殘差皆 $`\ne0`$，另 $`\nabla\cdot\mathbf E\ne0`$ | 閉式解見 `predicted_control()` |
| 對照組縮放 | 每步殘差 $`\propto(k_0\Delta)^2\cdot(\omega_0\Delta t)`$，即每弧度相位推進 $`O((k_0\Delta)^2)`$ | SPEC 的「$`O((k_0\Delta)^2)`$」指每弧度；每步因 $`\Delta t\propto\Delta`$ 是 $`O(\Delta^3)`$ |

### 關卡 1
| # | 項目 | 理論值公式 |
|---|---|---|
| 1 | SF 洩漏 | 0（使用 §6 的精確入射場）；因果時窗終點 $`t_{\rm end}=t_{\rm front}+(y_{\rm far}-y_0)/c+(y_{\rm far}-y_{\rm obs})/c`$，以 $`c`$ 為所有物理訊號群速度上界；洩漏量測在 SF 區（排除 PML）。對照組（解析場 × ramp，連續 $`k_y`$）穩態洩漏幅度 $`\propto\lvert k_y^{\rm cont}-k_y^{\rm disc}\rvert\propto\Delta^2`$ |
| 2 | 波前擬合 | $`\hat{\mathbf k}=\mathbf k_d`$；殘差 0 |
| 3 | 色散收斂 | $`k_y^{\rm disc}-k_y^{\rm cont}\approx-\dfrac{\Delta^2}{24k_y}\Big(S^2\omega^4/c^4-\sum_ik_i^4\Big)`$（由 $`\sin^2a\approx a^2-a^4/3`$ 展開）；預設 $`=3.013\,\Delta^2`$，$`\Delta=\lambda/20`$ 時 $`7.61\times10^{-3}`$（實際 $`7.58\times10^{-3}`$）；斜率 2 |
| 4 | Gauss／橫波 | $`\nabla\cdot\mathbf E=\nabla\cdot\mathbf H=0`$；$`\hat{\mathbf E}\cdot\tilde{\mathbf K}=\hat{\mathbf H}\cdot\tilde{\mathbf K}=\hat{\mathbf E}\cdot\hat{\mathbf H}^*=0`$（$`\hat{\ }`$＝去掉 $`e^{i\mathbf k\cdot\mathbf r_c}`$ 後的複振幅） |
| 5a | 阻抗 | 原始振幅：$`\lvert\mathbf H\rvert/\lvert\mathbf E\rvert=\lvert\tilde{\mathbf K}\rvert/(\mu_0\tilde\omega)\equiv1/\eta_0`$；胞心內插後：$`\lvert\mathbf H'\rvert/\lvert\mathbf E'\rvert=\sqrt{\sum_c(a_cH_{0c})^2}/\sqrt{\sum_c(a_cE_{0c})^2}`$，$`a_c`$ 為內插因子（$`E_x`$：$`\cos\frac{k_y\Delta}2\cos\frac{k_z\Delta}2`$，$`E_y`$：$`\cos\frac{k_x\Delta}2\cos\frac{k_z\Delta}2`$，$`E_z`$：$`\cos\frac{k_x\Delta}2\cos\frac{k_y\Delta}2`$，$`H_x`$：$`\cos\frac{k_x\Delta}2`$，$`H_y`$：$`\cos\frac{k_y\Delta}2`$，$`H_z`$：$`\cos\frac{k_z\Delta}2`$），與 $`1/\eta_0`$ 差 $`O((k\Delta)^2)`$ |
| 5b | 通量守恆 | 離散守恆通量（§9）$`\Phi(y)`$ 對 $`y`$ 為常數；理論值 $`\Phi_{\rm inc}=\tfrac12\cos(k_y\Delta/2)(E_{0z}H_{0x}-E_{0x}H_{0z})L_xL_z`$ |
| 5c | $`\mathbf S`$ 方向 | 胞心內插 $`\mathbf S'=\tfrac12\mathbf E'\times\mathbf H'`$，理論方向可由 $`a_c`$ 精確算出；另報 $`\angle(\mathbf S,\mathbf k)`$、$`\angle(\tilde{\mathbf K},\mathbf k)`$、$`\angle(\mathbf v_g,\mathbf k)`$ |
| 6 | PML 反射 | 無閉式；$`R_{\rm dB}=20\log_{10}(\lvert B\rvert/\lvert A\rvert)`$，$`B`$＝SF 區穩態反向平面波振幅，目標 $`\lt -40`$ dB |
| 7 | 穩定性 | 關源後總能量不增、衰減至 $`\lt 10^{-8}\times`$ 峰值 |
| 8 | 對稱性 | 鏡射 $`x\to-x`$：$`E_x'(x)=-\sigma E_x(-x)`$、$`E_{y,z}'(x)=\sigma E_{y,z}(-x)`$、$`H_x'=\sigma H_x(-x)`$、$`H_{y,z}'=-\sigma H_{y,z}(-x)`$，$`\sigma=-1`$（s，因 $`\hat s(-m)=-M\hat s(m)`$）、$`+1`$（p）；格點映射：半整數 $`x`$ 分量 $`i\to N_x-1-i`$，整數 $`x`$ 分量 $`i\to(N_x-i)\bmod N_x`$。Yee 格對此鏡射是等變的 ⇒ 理論差 0 |

### 關卡 2
| 項目 | 理論值 |
|---|---|
| $`R_s,R_p,T_s,T_p`$ | 連續 Fresnel：$`\cos\theta_1=k_y^{\rm cont}/k_0`$，$`n\sin\theta_2=\sin\theta_1`$；$`r_s=\frac{\cos\theta_1-n\cos\theta_2}{\cos\theta_1+n\cos\theta_2}`$，$`r_p=\frac{n\cos\theta_1-\cos\theta_2}{n\cos\theta_1+\cos\theta_2}`$；預設 $`R_s=0.0700`$、$`R_p=0.0179`$ |
| $`R+T`$ | 1（以 §9 守恆通量量測時為恆等式） |
| 介質中 $`k_y`$ | §2.2 公式取 $`n=1.5`$ |

### 關卡 3
Meep 與本程式同為標準 Yee、同 $`\Delta`$、同 $`\Delta t`$ ⇒ 真空 $`k_y`$ 理論上逐位相同（兩者都滿足 §2.2）。

---

## 6. 入射場的時間包絡：一維模態輔助線（Stage 0 推導，Stage 3 實作）

### 6.1 問題：解析平面波 × ramp 不是精確解

令 $`\mathbf F^{\rm inc}=g(t)\,\mathbf P`$，$`\mathbf P`$ 為 §2 的精確 CW 解。E 更新的殘差

```math
\mathbf R_E=g_{n+1}\mathbf P^{n+1}-g_n\mathbf P^{n}-\tfrac{\Delta t}{\varepsilon}\nabla\times(g_{n+\frac12}\mathbf P_H^{n+\frac12})
=(g_{n+1}-g_{n+\frac12})\mathbf P^{n+1}+(g_{n+\frac12}-g_n)\mathbf P^{n}\approx\tfrac{\Delta t}{2}g'(t)(\mathbf P^{n+1}+\mathbf P^n).
```

ramp 期間 $`\mathbf R_E\ne0`$，TF/SF 面就是一個大小 $`\sim g'/\omega_0\sim1/(\omega_0T_{\rm ramp})\approx1.6\times10^{-2}`$
（10 週期 ramp）的等效源，洩漏遠大於 $`10^{-10}`$。這不是實作錯誤，是「解析式 × 包絡」本身不滿足離散方程。

### 6.2 解法：同一 $`(k_x,k_z)`$ 模態的精確一維化

因 $`x,z`$ 週期且只有單一模態，入射場可寫成

```math
F_c^{\rm inc}(i,j,k,n)=\mathrm{Re}\big\{a_c[j](t_c)\,e^{i(k_xx_c+k_zz_c)}\big\},
```

$`x_c,z_c`$ 為分量 $`c`$ 自己的 Yee 座標。§2.1 的本徵關係對 $`x,z`$ 差分**逐點精確**成立，所以 3D Yee
方程等價於下列複數一維方程（$`y`$ 交錯與 3D 相同：$`a_{E_x},a_{E_z},a_{H_y}`$ 在整數 $`j`$，其餘在 $`j+\frac12`$）：

```math
\begin{aligned}
a_{H_x}[j+\tfrac12]&\mathrel{-}=\tfrac{\Delta t}{\mu_0}\Big[\tfrac{a_{E_z}[j+1]-a_{E_z}[j]}{\Delta}-i\tilde K_z\,a_{E_y}[j+\tfrac12]\Big]\\
a_{H_y}[j]&\mathrel{-}=\tfrac{\Delta t}{\mu_0}\Big[i\tilde K_z\,a_{E_x}[j]-i\tilde K_x\,a_{E_z}[j]\Big]\\
a_{H_z}[j+\tfrac12]&\mathrel{-}=\tfrac{\Delta t}{\mu_0}\Big[i\tilde K_x\,a_{E_y}[j+\tfrac12]-\tfrac{a_{E_x}[j+1]-a_{E_x}[j]}{\Delta}\Big]\\
a_{E_x}[j]&\mathrel{+}=\tfrac{\Delta t}{\varepsilon}\Big[\tfrac{a_{H_z}[j+\frac12]-a_{H_z}[j-\frac12]}{\Delta}-i\tilde K_z\,a_{H_y}[j]\Big]\\
a_{E_y}[j+\tfrac12]&\mathrel{+}=\tfrac{\Delta t}{\varepsilon}\Big[i\tilde K_z\,a_{H_x}[j+\tfrac12]-i\tilde K_x\,a_{H_z}[j+\tfrac12]\Big]\\
a_{E_z}[j]&\mathrel{+}=\tfrac{\Delta t}{\varepsilon}\Big[i\tilde K_x\,a_{H_y}[j]-\tfrac{a_{H_x}[j+\frac12]-a_{H_x}[j-\frac12]}{\Delta}\Big]
\end{aligned}
```

（例：3D 的 $`(E_y[k+1]-E_y[k])/\Delta`$ 作用在 $`e^{ik_zk\Delta}`$ 上等於 $`i\tilde K_z e^{ik_z(k+\frac12)\Delta}`$，正好是 $`H_x`$
的 $`z`$ 座標，所以係數裡不會出現額外相位。）

**做法**：在這條複數線上，於 $`j_a=j_0-N_a`$ 處用一維 TF/SF 注入「§2 解析平面波 × $`g(t)`$」。ramp 造成的
不精確只在 $`j_a`$ 產生一個一維散射場，但此後整條線上的 $`a_c`$ 在 $`j_a`$ 以外的每一個更新模板上都**精確**滿足
齊次方程。3D TF/SF（§3）使用 $`a_c[j_0]`$、$`a_c[j_0-\tfrac12]`$ 重建的 $`\mathbf F^{\rm inc}`$，因此 §3 的充要條件
逐步成立 ⇒ 3D 洩漏只剩捨入誤差，**與時間訊號形狀無關**。穩態時線上的場就是 §2 的解析解（因為一維注入
使用離散 $`k_y`$），所以 SPEC「入射場＝§2 六分量公式」在穩態逐字成立，ramp 期間則由輔助線保證精確。

**參數規則**：
- $`N_a\ge2`$ 即足以精確；取 $`N_a=40`$（$`2\lambda_0`$），讓一維注入點的準靜態殘留（p 偏振時 $`\nabla\cdot\mathbf J\ne0`$
  的電荷累積，場沿 $`y`$ 以 $`e^{-|k_t||y-y_a|}`$ 衰減）在 $`j_0`$ 處再衰減 $`e^{-7.5}\approx5\times10^{-4}`$ 倍。
- 輔助線兩端為 PEC，長度使反射來不及回到 $`j_0`$：Yee 的數值影響錐每步最多 1 格，端點距注入點至少
  $`N_{\rm steps}/2+N_a+10`$ 格；每步只更新 $`|j-j_a|\le n+2`$ 的光錐範圍。成本 $`O(N_{\rm steps}^2)`$ 個複數運算，遠小於 3D。
- 這不是電流片源：3D 網格內仍是單一平面 TF/SF。輔助線內的一維 TF/SF 等效源會使 p 偏振出現面電荷，
  但它位於 3D 計算域之外，其準靜態場本身也是齊次方程的精確解，不造成 3D 洩漏。

**保留解析模式作對照**：`inc=analytic` 直接用 §2 公式 × $`g(t)`$（可選 `ky=continuous`），用於關卡 1-1
的對照組與證明 §6.1 的推論。

---

---

## 7. CPML（y 兩端；Stage 2 實作）

座標伸縮 $`\partial_y\to\frac{1}{s_y}\partial_y`$，$`s_y=\kappa+\dfrac{\sigma}{\alpha+i\omega\varepsilon_0}`$（CFS）。遞迴卷積形式
（Roden & Gedney 2000）：對每個含 $`\partial_y`$ 的更新項

```math
\frac{1}{s_y}\partial_yF\ \to\ \frac{1}{\kappa}\partial_yF+\psi,\qquad
\psi^{n}=b\,\psi^{n-1}+c\,\partial_yF,\quad b=e^{-(\sigma/\kappa+\alpha)\Delta t/\varepsilon_0},\quad
c=\frac{\sigma\,(b-1)}{\kappa(\sigma+\kappa\alpha)} .
```

四個輔助量：$`\psi_{E_xy}`$（$`\partial_yH_z`$）、$`\psi_{E_zy}`$（$`\partial_yH_x`$）、$`\psi_{H_xy}`$（$`\partial_yE_z`$）、$`\psi_{H_zy}`$
（$`\partial_yE_x`$）；更新式中的符號與原 $`\partial_y`$ 項相同，例如
$`E_x\mathrel{+}=\frac{\Delta t}{\varepsilon}\big[\frac1\kappa\partial_yH_z+\psi_{E_xy}-\partial_zH_y\big]`$、
$`E_z\mathrel{+}=\frac{\Delta t}{\varepsilon}\big[\partial_xH_y-(\frac1\kappa\partial_yH_x+\psi_{E_zy})\big]`$。
程式中 $`\psi`$ 以「差分單位」（$`\psi\Delta`$）儲存。

**剖面**（$`\rho`$ = 節點深入 PML 的距離，$`d=N_{\rm pml}\Delta`$，每個節點在自己的 $`y`$ 位置取值，E 節點整數、H 節點半整數）：
$`\sigma=\sigma_{\max}(\rho/d)^{m}`$，$`\kappa=1+(\kappa_{\max}-1)(\rho/d)^{m}`$，$`\alpha=\alpha_{\max}(1-\rho/d)`$，$`m=3`$。

| 參數 | 值 | 理由 |
|---|---|---|
| $`\sigma_{\max}`$ | $`0.8(m+1)/(\eta\Delta)`$，$`\eta=\eta_0/n`$ | SPEC 指定的最佳值；介質端除以 $`n`$，使每格衰減與真空相同（衰減率 $`\propto n\,\sigma\eta_0`$，Taflove 式 $`\sigma_{\rm opt}\propto1/\sqrt{\varepsilon_r}`$） |
| $`\kappa_{\max}`$ | 5 | 讓近掠射/消逝成分也被拉伸；關卡 1-6 另報 $`\kappa_{\max}=1,10`$ 的敏感度 |
| $`\alpha_{\max}`$ | $`0.05\,\omega_0\varepsilon_0=0.1\pi`$ | CFS 極點，改善低頻與晚期穩定性；遠低於 $`\omega_0`$ 以免削弱工作頻率吸收；關卡 1-6 另報 $`\alpha_{\max}=0`$ |
| 外邊界 | PEC | $`j=0,N_y`$ 的 $`E_x,E_z\equiv0`$ |

**Stage 2 驗證**：(a) 無 PML 的 PEC/PBC 腔體中，Yee 守恆能量
$`W^{n+\frac12}=\tfrac12\sum(\varepsilon\,\mathbf E^n\!\cdot\mathbf E^{n+1}+\mu|\mathbf H^{n+\frac12}|^2)\Delta^3`$
1000 步漂移 $`8.9\times10^{-16}`$（證明能量診斷正確；若用 $`\tfrac12\sum(E^2+H^2)`$，會因 E、H 半步錯開而出現 $`2\omega`$ 振盪，
第一次 Stage 2 就是因此誤判「成長」）；(b) 有 CPML 時總能量逐步單調不增（最大 $`+6\times10^{-16}`$ = 捨入），
5000 步衰減到 $`4\times10^{-9}`$。

---

## 9. 守恆通量與 Poynting 向量的內插（Stage 4）

### 9.1 為什麼 $`S_y`$ 不能用胞心內插
把六分量都平均到胞心再算 $`\tfrac12\mathrm{Re}(\mathbf E\times\mathbf H^*)`$，每個分量乘上不同的
$`\cos(k_d\Delta/2)`$ 因子（§5 L1-5a）。這對方向與阻抗是 $`O((k\Delta)^2)`$ 的偏差；但在真空與介質中 $`k_y`$ 不同，
會讓 $`R+T`$ 產生 $`\sim1.5\times10^{-2}`$ 的假誤差（$`\cos(k_y^{\rm vac}\Delta/2)=0.9921`$ vs
$`\cos(k_y^{\rm med}\Delta/2)=0.9768`$），遠超 $`10^{-4}`$ 門檻。所以通量另用下面的精確定義。

### 9.2 離散守恆通量
時諧、無損、無源區域內，離散方程給 $`\nabla_h\times\mathbf H=-i\tilde\omega\varepsilon\mathbf E`$、
$`\nabla_e\times\mathbf E=i\tilde\omega\mu\mathbf H`$，所以
$`\mathrm{Re}\sum_V[\mathbf E^*\!\cdot(\nabla_h\times\mathbf H)-\mathbf H\cdot(\nabla_e\times\mathbf E)^*]=0`$。
$`x,z`$ 週期求和使橫向差分項相消；$`y`$ 方向分部求和（summation by parts）只留下兩個端面項，因此

```math
\Phi(j)=\tfrac12\mathrm{Re}\sum_{i,k}\Big[E_z|_{j}\,H_x^*|_{j+\frac12}-E_x|_{j}\,H_z^*|_{j+\frac12}\Big]\Delta^2
```

對 $`j`$ **精確為常數**（對任何滿足離散方程的場，包括正反向波疊加與介面上平均過的 $`\varepsilon`$）。
配對的分量在 $`x,z`$ 上同點（$`E_z`$ 與 $`H_x`$ 都在 $`(i,\,\cdot\,,k+\frac12)`$；$`E_x`$ 與 $`H_z`$ 都在 $`(i+\frac12,\,\cdot\,,k)`$），
只差 $`y`$ 半格——不需要任何內插。這就是 C 程式 y-平面 DFT 的配對方式（平面 $`j`$ 存整數分量於 $`j`$、半整數分量於 $`j+\frac12`$）。
單一平面波的理論值：$`\Phi_{\rm inc}=\tfrac12\cos(k_y\Delta/2)(E_{0z}H_{0x}-E_{0x}H_{0z})L_xL_z`$。

正反向波的交叉項對 $`j`$ 以 $`e^{\pm2ik_yj\Delta}`$ 振盪，而 $`\Phi`$ 恆定 ⇒ 交叉項必為 0；所以
$`\Phi_{\rm TF}=\Phi_{\rm inc}+\Phi_{\rm refl}`$、$`\Phi_{\rm SF}=\Phi_{\rm refl}`$，$`R+T=1`$ 在此定義下是恆等式。

### 9.3 Poynting 向量方向（需要三個分量時）
把六分量各自以兩點平均移到胞心 $`(i+\frac12,j+\frac12,k+\frac12)`$：$`E_x`$ 平均 $`(j,j+1)\times(k,k+1)`$ 四點；
$`E_y`$ 平均 $`(i,i+1)\times(k,k+1)`$；$`E_z`$ 平均 $`(i,i+1)\times(j,j+1)`$；$`H_x`$ 平均 $`(i,i+1)`$；$`H_y`$ 平均 $`(j,j+1)`$；
$`H_z`$ 平均 $`(k,k+1)`$。再算 $`\mathbf S=\tfrac12\mathrm{Re}(\mathbf E\times\mathbf H^*)`$。對平面波這個內插的結果可精確預測
（分量乘上 §5 的 $`a_c`$），所以方向同時與「內插理論」與 $`\mathbf k`$、$`\tilde{\mathbf K}`$、$`\mathbf v_g`$ 比較。

---

---

## 8. 介面上的 ε（Stage 5）

**位置選擇**：介面放在整數 $`y`$ 平面 $`y_1=j_1\Delta`$（$`j_1=j_0+3N_\lambda`$），也就是切向 $`E_x,E_z`$ 所在的平面；
法向 $`E_y`$ 在 $`j\pm\frac12`$，**永遠不落在介面上**。

| 分量 | 位置 | 取值 | 理由 |
|---|---|---|---|
| $`E_x,E_z`$（切向） | $`j=j_1`$ | $`\varepsilon=\tfrac12(1+n^2)`$（算術平均） | 對偶胞 $`[y_1-\frac\Delta2,y_1+\frac\Delta2]`$ 一半真空一半介質；切向 $`E`$ 連續 ⇒ 胞內平均 $`D_t=\langle\varepsilon\rangle E_t`$（並聯電容） |
| $`E_x,E_z`$ | $`j\lt j_1`$ / $`j\gt j_1`$ | 1 / $`n^2`$ | 整個對偶胞在單一介質 |
| $`E_y`$（法向） | $`j+\frac12\lt j_1`$ / $`\gt j_1`$ | 1 / $`n^2`$ | 對偶胞 $`[j,j+1]\Delta`$ 整個在單一介質，無需調和平均（若介面切過 $`E_y`$ 才需 $`\langle\varepsilon^{-1}\rangle^{-1}`$，因法向 $`D`$ 連續、串聯電容） |
| PML 內 | $`j\gt N_y-N_{\rm pml}`$ | $`n^2`$，$`\sigma_{\max}`$ 除以 $`n`$ | 介質延伸進遠端 PML（SPEC 要求），見 §7 |

對平面介面，切向算術平均 + 法向不落介面，給出 $`O(\Delta^2)`$ 的反射係數誤差；對照組 `ifmode=s`（介面平面上切向
直接取 $`n^2`$，即階梯近似）只有 $`O(\Delta)`$。關卡 2 兩者都跑以示差異。

**與 Meep 的對應**：Meep 的 subpixel averaging（Kottke 各向異性平均）對「與網格對齊的平面介面」，在每個 Yee
分量自己的體素上做切向 $`\langle\varepsilon\rangle`$、法向 $`\langle\varepsilon^{-1}\rangle^{-1}`$。本專案的介面落在 Meep 的整數
網格面上（Meep 原點在 cell 中心，網格點在 $`\Delta`$ 的整數倍，已由 `get_array_metadata`／陣列形狀確認），
所以 Meep `eps_averaging=True` 給出與本程式**相同**的 $`\varepsilon`$ 分佈；關卡 3 以此為主要比對，並另跑
`eps_averaging=False` 作為差異歸因。

### 8.1 Yee 晶格的精確（離散）Fresnel 係數

因 $`k_x,k_z`$ 固定，s、p 在晶格上也互不耦合（$`\hat s\propto\tilde{\mathbf K}\times\hat y=(-\tilde K_z,0,\tilde K_x)`$ 與 $`\tilde K_y`$ 無關，
兩側相同）。取介面 $`y=0`$：
- 真空側（$`y\le0`$ 的 E 節點與 $`y\le-\frac\Delta2`$ 的 H 節點）：入射 $`(\mathbf E_i,\mathbf H_i)e^{ik_yy}`$ ＋ 反射 $`r(\mathbf E_r,\mathbf H_r)e^{-ik_yy}`$；
- 介質側（$`y\ge0`$、$`y\ge\frac\Delta2`$）：$`t(\mathbf E_t,\mathbf H_t)e^{ik_y^{(m)}y}`$，$`k_y^{(m)}`$ 由 §2.2 取 $`n=1.5`$；
- 各波皆為 §2 的離散平面波（$`\mathbf H=\tilde{\mathbf K}\times\mathbf E/\mu_0\tilde\omega`$，偏振由各自 $`\tilde{\mathbf K}`$ 定義）。

介面平面上兩個條件決定 $`r,t`$：(i) 切向 $`E`$ 連續（$`E_x,E_z`$ 節點只有一個值）；(ii) 該平面的 Ampère 更新式
（$`\varepsilon=\varepsilon_{\rm if}`$）：$`-i\tilde\omega\varepsilon_{\rm if}E_x=\frac{H_z(\frac\Delta2)-H_z(-\frac\Delta2)}{\Delta}-i\tilde K_zH_y(0)`$、
$`-i\tilde\omega\varepsilon_{\rm if}E_z=i\tilde K_xH_y(0)-\frac{H_x(\frac\Delta2)-H_x(-\frac\Delta2)}{\Delta}`$，其中 $`H_y(0)`$ 只由切向 $`E(0)`$ 決定。
投影到切向方向 $`\hat u`$（s：$`\hat s`$；p：$`\hat k_t`$）得到 $`2\times2`$ 線性方程。$`R=-\Phi_r/\Phi_i`$、$`T=\Phi_t/\Phi_i`$ 以 §9 的守恆通量計算，
$`R+T=1`$ 到捨入誤差（已數值確認）。實作：`fdtd_theory.discrete_fresnel()`。

**結果（$`\theta_1=36.9^\circ`$，$`n=1.5`$，算術平均介面）**：

| Δ | $`R_s`$（離散） | 相對連續 Fresnel | $`R_p`$（離散） | 相對連續 Fresnel |
|---|---|---|---|---|
| λ0/10 | 0.053767 | −23.2% | 0.009809 | −45.1% |
| λ0/20 | 0.066145 | −5.50% | 0.015830 | −11.4% |
| λ0/40 | 0.069042 | −1.36% | 0.017356 | −2.85% |
| λ0/160 | 0.069932 | −0.084% | 0.017832 | −0.18% |

誤差精確呈 $`O(\Delta^2)`$（每次減半約 ÷4），但係數大：$`R`$ 是兩個相近導納之差的平方，介質內每格相位
$`nk_0\Delta=0.47`$ rad，晶格色散對 $`r`$ 的相對影響被放大。**在 Δ=λ0/20 時，標準 Yee＋任何非人為調參的介面處理都達不到
「$`R`$ 與連續 Fresnel 相對誤差 < 1%」**：算術平均 −5.5%/−11.4%，調和平均 −3.8%/−7.5%，階梯（任一側）+5.7%/+14.8%；
要讓離散 $`R`$ 等於連續值需 $`\varepsilon_{\rm if}\approx1.19`$（s）／$`1.21`$（p）——偏振相依的湊數，不是方法。
$`T`$ 則在 λ0/20 時只差 +0.41%（s）、+0.20%（p）。Meep（獨立實作）在 λ0/10 得 $`R_s=0.053738`$，與本節理論差 0.05%，
佐證這是 Yee 離散化本身的性質，而非本程式的錯誤。
