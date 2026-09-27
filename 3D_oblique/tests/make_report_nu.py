"""Build validation_report.md = Part I (legacy uniform report, tests/make_report.py, after decision D7) +
Part II (nonuniform grid, SPEC_nonuniform.md). Every number comes from results/*.json; nothing is typed by hand.

    python3 tests/make_report_nu.py
"""
import json
import os
import platform
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
GATES = [
    ("gate0_regression", "第 0 關 0-1：回歸 parity 與舊驗證"),
    ("gate0_adjoint", "第 0 關 0-2：伴隨性"),
    ("gate0_energy_stability", "第 0 關 0-3、0-4：能量守恆與穩定邊界"),
    ("gate0_divergence", "第 0 關 0-5：散度"),
    ("gate1", "第 1 關（階段 A）：真空自我驗證"),
    ("gate2", "第 2 關（階段 A）：薄膜與 5 層堆疊 vs TMM"),
    ("gate3", "第 3 關：與 Meep 比對"),
    ("B_gate0", "階段 B：中止偵測與第 0 關"),
    ("B_gate1", "階段 B：注入誤差、波前與 Floquet 純度"),
    ("B_gate3", "階段 B：光柵 vs RCWA"),
]
FIGS = [
    ("figures/nu/gate0/energy_eig.png", "0-3/0-4：離散能量漂移與 power iteration 收斂"),
    ("figures/nu/gate1/theta_y.png", "1-3b：θ(y) 量測（點）與預測（虛線）"),
    ("figures/nu/gate1/snapshot_xy.png", "1-3c：x–y 瞬時場（真實座標、等比例）與預測等相位線"),
    ("figures/nu/gate1/reflection_scan.png", "1-4c：漸變反射 vs r_max 與解析度"),
    ("figures/nu/gate1/pml.png", "1-6：遠端 PML 反射 vs N_pml"),
    ("figures/nu/gate1/energy_decay.png", "1-7：關源後能量衰減"),
    ("figures/nu/gate2/convergence.png", "2-3/2-4：R 對 Δy 的收斂（對齊非均勻 vs 不對齊均勻對照）"),
    ("figures/nu/gate3/meep_convergence.png", "3A：Meep 收斂"),
    ("figures/nu/grids", "網格間距 h(y)"),
    ("figures/stageB/grating_orders.png", "B3：光柵各繞射階效率 vs RCWA"),
]


def load(name):
    p = os.path.join(RES, f"nu_{name}.json")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def esc(s):
    return str(s).replace("|", "\\|").replace("\n", " ")


def table(rows):
    out = ["| 項目 | 量 | 理論值 | 量測值 | 門檻 | 判定 | 網格 | 備註 |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        v = "INFO" if r["passed"] is None else ("**PASS**" if r["passed"] else "**FAIL**")
        out.append(f"| {esc(r['item'])} | {esc(r['quantity'])} | {esc(r['theory'])} | {esc(r['measured'])} | "
                   f"{esc(r['threshold'])} | {v} | {esc(r.get('grid', ''))} | {esc(r.get('note', ''))} |")
    return "\n".join(out)


def status(d):
    if d is None:
        return "NOT RUN"
    if d.get("blocked"):
        return "BLOCKED"
    return "PASS" if d.get("passed") else "FAIL"


def demote(md):
    return re.sub(r"^(#+) ", lambda m: "#" + m.group(1) + " ", md, flags=re.M)


def main():
    subprocess.run([sys.executable, os.path.join(ROOT, "tests", "make_report.py")], cwd=ROOT, check=True)
    part1 = open(os.path.join(ROOT, "validation_report.md"), encoding="utf-8").read()
    res = {k: load(k) for k, _ in GATES}
    st = {k: status(v) for k, v in res.items()}
    stageA = [k for k, _ in GATES if not k.startswith("B_")]
    overall = "PASS" if all(st[k] == "PASS" for k, _ in GATES) else (
        "BLOCKED" if all(st[k] in ("PASS", "BLOCKED") for k, _ in GATES) else "FAIL")
    md = ["# 驗證報告：3D 斜入射 FDTD（均勻網格 Part I ＋ 非均勻網格 Part II）", ""]
    md.append(f"**Part II（非均勻網格）總結論：{overall}**　階段 A："
              f"{'PASS' if all(st[k] == 'PASS' for k in stageA) else 'NOT PASS'}；階段 B："
              f"{'PASS' if all(st[k] == 'PASS' for k, _ in GATES if k.startswith('B_')) else 'NOT PASS'}")
    md.append("")
    md.append("| 關卡 | 結果 |\n|---|---|")
    for k, title in GATES:
        md.append(f"| {title} | {st[k]} |")
    fails = [(k, r) for k, _ in GATES if res[k] for r in res[k]["table"] if r["passed"] is False]
    fl = os.path.join(RES, "nu_failures.json")
    if fails or os.path.exists(fl):
        md.append("\n## 未通過紀錄\n")
        if fails:
            md.append(table([r for _, r in fails]))
        if os.path.exists(fl):
            with open(fl, encoding="utf-8") as fh:
                F = json.load(fh)
            md.append("\n| 關卡 | 現象 | 假設 | 最小實驗 | 修正 | 結果 |\n|---|---|---|---|---|---|")
            for f in F:
                md.append(f"| {esc(f['gate'])} | {esc(f['symptom'])} | {esc(f['hypothesis'])} | {esc(f['experiment'])} | "
                          f"{esc(f['fix'])} | {esc(f['result'])} |")
    md.append("\n## 規格變更紀錄（使用者決策 D1–D8，實作期間修正 D9–D19，使用者決策 D20–D21）\n")
    spec = open(os.path.join(ROOT, "SPEC_nonuniform.md"), encoding="utf-8").read()
    m = re.search(r"### 規格變更紀錄（實作期間）\n\n(.*?)\n\n以上都沒有改任何門檻數值", spec, re.S)
    if m:
        md.append(m.group(1))
    md.append("\n所有門檻值在執行前寫入 `tests/thresholds_nu.json`；理論值在執行前寫入 `results/nu_predictions.json`（`v2/` 前綴）。"
              "兩檔的 sha256 記錄在 `CLAUDE.md`，`tests/run_nu_regression.py` 每次先比對。")
    md.append("\n## 環境\n")
    try:
        gcc = subprocess.run(["gcc", "--version"], capture_output=True, text=True).stdout.splitlines()[0]
    except Exception:
        gcc = "?"
    md.append(f"| 項目 | 值 |\n|---|---|\n| 平台 | {platform.platform()} |\n| C 編譯器 | {gcc} |\n"
              f"| 編譯選項 | `-O3 -march=native -std=c99 -fopenmp` |\n| Python | {platform.python_version()} (numpy, matplotlib) |")
    for k, title in GATES:
        md.append(f"\n## {title}\n")
        if res[k] is None:
            md.append("（尚未執行）")
            continue
        md.append(f"結果：**{st[k]}**（{res[k]['created_utc']}）\n")
        md.append(table(res[k]["table"]))
    md.append("\n## 圖\n")
    for path, cap in FIGS:
        full = os.path.join(ROOT, path)
        if os.path.isdir(full):
            for f in sorted(os.listdir(full)):
                if f.endswith(".png"):
                    md.append(f"![{cap}: {f}]({path}/{f})")
        elif os.path.exists(full):
            md.append(f"**{cap}**\n\n![{cap}]({path})\n")
    for title, doc in (("附錄 A：均勻假設盤點表", "docs/inventory_nonuniform.md"),
                       ("附錄 B：推導", "docs/derivation_nonuniform.md"),
                       ("附錄 C：C 程式修改說明", "docs/DIFF_NOTES.md")):
        txt = open(os.path.join(ROOT, doc), encoding="utf-8").read()
        md.append(f"\n<details><summary>{title}（{doc}）</summary>\n\n{demote(demote(txt))}\n\n</details>\n")
    md.append("\n---\n\n# Part I：均勻網格驗證（舊版，D7 生效後重新產生）\n")
    md.append(demote(part1))
    with open(os.path.join(ROOT, "validation_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(md) + "\n")
    print("validation_report.md written; Part II overall:", overall)


if __name__ == "__main__":
    main()
