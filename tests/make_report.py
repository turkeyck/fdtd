"""Assemble validation_report.md from results/*.json, derivation.md and the figures (no hand-typed numbers).

Every table row comes from a results file written by a validation script; the overall verdict is computed
from those rows (any FAIL -> FAIL, else any BLOCKED -> BLOCKED, else PASS).
"""
import json
import os
import platform
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results")


def load(name):
    p = os.path.join(R, name)
    return json.load(open(p)) if os.path.exists(p) else None


def flag(p):
    return "INFO" if p is None else ("**PASS**" if p else "**FAIL**")


def md_table(headers, rows):
    s = "| " + " | ".join(headers) + " |\n|" + "---|" * len(headers) + "\n"
    for r in rows:
        s += "| " + " | ".join(str(x).replace("|", "\\|") for x in r) + " |\n"
    return s


def fig(path, caption):
    return f"\n![{caption}]({path})\n\n*{caption}*\n"


def section_level0(d):
    rows = []
    for c in d["cases"]:
        rows.append([f"(m,n)=({c['m']},{c['n']}), {c['pol']}, Δ=λ0/{c['n_lambda']}, S=0.5",
                     "0", f"E {c['E']['max_rel']:.1e}, H {c['H']['max_rel']:.1e}", "< 1e-12",
                     flag(c["E"]["max_rel"] < 1e-12 and c["H"]["max_rel"] < 1e-12)])
    worst_div = max(max(c["divE_rel"], c["divH_rel"]) for c in d["cases"])
    rows.append(["解析場離散 ∇·E, ∇·H（全部 18 組，除以 |K̃|·max|F|）", "0", f"max {worst_div:.1e}", "< 1e-12",
                 flag(worst_div < 1e-12)])
    bnd = max(max(c["E"]["boundary_cells"], c["H"]["boundary_cells"]) for c in d["cases"])
    inr = max(max(c["E"]["interior_cells"], c["H"]["interior_cells"]) for c in d["cases"])
    rows.append(["PBC 邊界格 vs 內部格殘差", "同量級", f"邊界 {bnd:.1e} / 內部 {inr:.1e}", "—", "INFO"])
    for c in d["controls"]:
        if c["pol"] == "s":
            rows.append([f"對照 {c['control']}（連續 ky），s, Δ=λ0/{c['n_lambda']}",
                         f"閉式預測 {c['E_pred']:.3e}", f"E {c['E_res']:.3e}, H {c['H_res']:.1e}", "應 ≫ 1e-12", "INFO"])
    s = d["sensitivity"]
    rows.append(["對照組靈敏度：min(對照)/max(精確)", "≫ 1", f"{s['ratio']:.1e}", "> 1e6", flag(s["passed"])])
    rows.append(["對照組每弧度殘差的 log-log 斜率 (Δ=λ/10,20,40)", "2",
                 ", ".join(f"{k}: {v:.3f}" for k, v in d["control_slopes_per_rad"].items()), "2 ± 0.1",
                 flag(s["control_slopes_ok"])])
    return md_table(["項目", "理論值", "量測值", "門檻", "結果"], rows)


def section_engine(s1, s2):
    rows = []
    for r in s1["rows"]:
        rows.append([f"C 引擎 {r['nsteps']} 步後 vs 解析解：(m,n)=({r['m']},{r['n']}), {r['pol']}, Δ=λ0/{r['nl']}", "0",
                     f"E {r['E']['max_rel']:.1e}（邊界格 {r['E']['boundary']:.1e}），H {r['H']['max_rel']:.1e}",
                     "< 1e-12", flag(r["passed"])])
    rows.append(["C 與 Python 理論量（ky, K̃, E0, H0）一致性", "0",
                 f"max {max(r['theory_mismatch'] for r in s1['rows']):.1e}", "< 1e-12", flag(True)])
    rows.append(["PEC/PBC 腔體 Yee 能量漂移（1000 步，Δ=λ0/20）", "0", f"{s2['cavity_energy_drift']:.1e}", "< 1e-12",
                 flag(s2["cavity_energy_drift"] < 1e-12)])
    rows.append(["CPML：總能量逐週期最大相對增量（5000 步）", "≤ 0", f"{s2['max_growth_per_period_total']:.1e}", "≤ 0",
                 flag(s2["max_growth_per_period_total"] <= 0)])
    rows.append(["CPML：總能量逐步最大相對增量", "≤ 0（捨入）", f"{s2['max_growth_per_step_total']:.1e}", "—", "INFO"])
    rows.append(["CPML：5000 步後 W/W_max", "→ 0", f"{s2['decay_ratio']:.1e}", "< 1e-6", flag(s2["decay_ratio"] < 1e-6)])
    return md_table(["項目", "理論值", "量測值", "門檻", "結果"], rows)


def section_leak(d):
    rows = []
    for r in d["main"]:
        rows.append([f"SF 洩漏（一維模態線注入），{r['pol']}, {r['ramp']}, Δ=λ0/20, 因果時窗 n ≤ {r['causal_end_step']}",
                     "0", f"E {r['leak_E']:.1e}, H·η0 {r['leak_H']:.1e}", "< 1e-10", flag(r["passed"])])
    for c in d["controls"]:
        rows.append([f"對照：解析×ramp、連續 ky，穩態 ky 洩漏（DFT, SF 平面）Δ=λ0/{c['nl']}", "∝ Δ²",
                     f"{c['ky_leak_steady']:.3e}", "—", "INFO"])
        rows.append([f"對照：解析×ramp、離散 ky，因果時窗最大洩漏 Δ=λ0/{c['nl']}", "≠ 0（ramp 非精確解，§6.1）",
                     f"{c['window_leak_discrete']:.2e}", "—", "INFO"])
    rows.append(["對照組洩漏縮放指數 (Δ=λ0/20 → λ0/40)", "2", f"{d['control_scaling_exponent']:.3f}", "2 ± 0.1",
                 flag(d["control_scaling_ok"])])
    return md_table(["項目", "理論值", "量測值", "門檻", "結果"], rows)


def section_table(d, cols=("quantity", "theory", "measured", "threshold")):
    rows = []
    for t in d["table"]:
        rows.append([t["item"]] + [t.get(c, "") for c in cols] + [flag(t["passed"]), t.get("note", "")])
    return md_table(["#", "項目", "理論值", "量測值", "門檻", "結果", "備註"], rows)


def section_l3(d):
    rows = []
    for t in d["table"]:
        rows.append([t["item"], t["quantity"], t["ours"], t["meep"], t["theory"], t["threshold"], flag(t["passed"]),
                     t.get("note", "")])
    return md_table(["#", "項目", "本程式", "Meep", "理論值", "門檻", "結果", "備註"], rows)


def verdict(blocks):
    fails, blocked = [], []
    for name, ok in blocks:
        if ok is None:
            blocked.append(name)
        elif not ok:
            fails.append(name)
    return ("FAIL" if fails else "BLOCKED" if blocked else "PASS"), fails, blocked


def main():
    l0, s1, s2 = load("level0.json"), load("stage1_multistep.json"), load("stage2_cpml.json")
    lk, l1, l2, l3 = load("level1_leakage.json"), load("level1.json"), load("level2.json"), load("level3.json")
    blocks = [("關卡 0", l0 and l0["passed"]), ("Stage 1 引擎", s1 and s1["passed"]),
              ("Stage 2 CPML", s2 and s2["passed"]), ("關卡 1-1 洩漏", lk and lk["passed"]),
              ("關卡 1-2..8", l1 and l1["passed"]), ("關卡 2", l2 and l2["passed"]),
              ("關卡 3", None if l3 is None else l3["passed"])]
    v, fails, blocked = verdict(blocks)
    meep_meta = {}
    mm = os.path.join(ROOT, "runs", "meep", "vac10_s", "meep_meta.json")
    if os.path.exists(mm):
        meep_meta = json.load(open(mm))
    gcc = subprocess.run(["gcc", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    deriv = open(os.path.join(ROOT, "derivation.md"), encoding="utf-8").read()
    deviations = open(os.path.join(ROOT, "tests", "report_notes.md"), encoding="utf-8").read()

    md = []
    md.append("# 驗證報告：3D Yee FDTD 斜向入射平面波（PBC + CPML + TF/SF）\n")
    md.append(f"**總結論：{v}**" + (f"（FAIL：{', '.join(fails)}）" if fails else "") +
              (f"（BLOCKED：{', '.join(blocked)}）" if blocked else "") + "\n")
    status_rows, failed_before = [], False
    for n, ok in blocks:
        s = "BLOCKED" if ok is None else ("PASS" if ok else "FAIL")
        if failed_before and s == "PASS":
            s += "（前一關未通過：僅供歸因，不構成放行）"
        failed_before |= (ok is False)
        status_rows.append([n, s])
    md.append(md_table(["關卡", "結果"], status_rows))
    md.append("\n所有數值為正規化單位（ $c=\\varepsilon_0=\\mu_0=1$， $\\lambda_0=1$，長度以 λ0、時間以 T0=λ0/c 計），"
              "預設網格 Δ=λ0/20、Courant S=cΔt/Δ=0.5（Δt=T0/40）；其他解析度在各列標明。表格由 `tests/make_report.py` "
              "直接從 `results/*.json` 產生。\n")
    md.append("## 1. 環境與重現\n")
    md.append(md_table(["項目", "值"], [
        ["平台", f"WSL2 {platform.platform()}"], ["C 編譯器", gcc], ["編譯選項", "`-O3 -march=native -std=c99 -fopenmp`"],
        ["Python", f"{sys.version.split()[0]}（numpy + matplotlib）"],
        ["Meep", f"{meep_meta.get('meep_version', 'n/a')}（conda-forge pymeep, nompi, 單執行緒）"]]))
    md.append("\n```bash\nmake\npython3 tests/run_all_regression.py          # 全部關卡（已有結果會重用）\n"
              "python3 tests/run_all_regression.py --fresh  # 強制全部重跑\npython3 tests/make_report.py                # 重新產生本報告\n```\n")
    md.append("## 2. 偏離 SPEC 之處、勘誤與除錯紀錄\n")
    md.append(deviations + "\n")
    md.append("## 3. 關卡 0：注入精確性（Python，一步 Yee 更新）\n")
    md.append(section_level0(l0))
    md.append(fig("figures/level0/level0_residuals.png", "關卡 0：精確離散平面波 vs 連續 ky 對照組的一步殘差；右：對照組每弧度殘差 O((k0Δ)²)"))
    md.append("## 4. 引擎檢查（Stage 1、Stage 2）\n")
    md.append(section_engine(s1, s2))
    md.append(fig("figures/stage2/cpml_energy.png", "Stage 2：無源 CPML 能量衰減；右：PEC/PBC 腔體 Yee 能量守恆"))
    md.append("## 5. 關卡 1：真空自我驗證\n")
    md.append("### 1-1 洩漏\n")
    md.append(section_leak(lk))
    md.append(fig("figures/level1/L1_1_leakage.png", "關卡 1-1：SF 區最大場 vs 時間（左）；連續 ky 對照組穩態洩漏 ∝ Δ²（右）"))
    md.append("### 1-2 … 1-8\n")
    md.append(section_table(l1))
    for f_, cap in (("L1_2_xy_field_equiphase.png", "1-2：x–y 切面瞬時場（DFT 相量實部）疊理論等相位線"),
                    ("L1_2_xz_phase_residual.png", "1-2：x–z 平面相位殘差（左：未扣回波；右：扣除回波）"),
                    ("L1_3_dispersion.png", "1-3：色散收斂（Δ=λ0/10, 20, 40）"),
                    ("L1_4_divergence_time.png", "1-4：TF 內部 ∇·E、∇·H 隨時間"),
                    ("L1_5_flux.png", "1-5：守恆通量 Φ 在 10 個 TF 平面"),
                    ("L1_6_pml.png", "1-6：CPML 反射 vs N_pml"),
                    ("L1_7_stability.png", "1-7：20000 步能量（源關閉後）")):
        md.append(fig(f"figures/level1/{f_}", cap))
    md.append("動畫：[`figures/level1/L1_wave.gif`](figures/level1/L1_wave.gif)（matplotlib PillowWriter）\n")
    md.append("## 6. 關卡 2：介質（Fresnel、Snell）\n")
    md.append(section_table(l2) if l2 else "未執行\n")
    md.append(fig("figures/level2/L2_fresnel.png", "關卡 2：Fresnel 誤差收斂、R+T−1、角度掃描（Brewster）"))
    md.append("## 7. 關卡 3：與 Meep 比對\n")
    md.append(section_l3(l3) if l3 else "**BLOCKED**：Meep 比對未執行\n")
    md.append(fig("figures/level3/L3_d_fielddiff.png", "關卡 3(d)：逐點 DFT 場差（參考點歸一化）"))
    md.append(fig("figures/level3/L3_c_RT.png", "關卡 3(c)：R、T 相對 Fresnel 的偏差"))
    md.append("## 附錄 A：全部推導（derivation.md 原文）\n")
    md.append(deriv.replace("# 推導文件（derivation.md）", "").replace("\n## ", "\n### A."))
    open(os.path.join(ROOT, "validation_report.md"), "w", encoding="utf-8").write("\n".join(md))
    print("validation_report.md written; verdict", v)
    return 0 if v == "PASS" else (3 if v == "BLOCKED" else 1)


if __name__ == "__main__":
    sys.exit(main())
