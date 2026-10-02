# CLAUDE.md（Claude Code 長期規則）

## 共用不可變規則（SPEC_nonuniform.md）

1. 規格來源：`SPEC_nonuniform.md`。門檻只在 `tests/thresholds_nu.json`，理論值只在 `results/nu_predictions.json`；兩者不得在沒有使用者決策紀錄的情況下修改（`tests/run_nu_regression.py` 會比對 sha256）。
2. 均勻是非均勻的特例：不得出現 `if (uniform) {舊碼} else {新碼}` 的分叉更新式。舊碼凍結在 `ref/fdtd3d_oblique_uniform_ref.c`（`make ref`）。
3. Python 只用 numpy + matplotlib（Meep 腳本例外）；動畫用 `matplotlib.animation.PillowWriter`，不用 imageio。nonuniform 分析的座標一律讀 run 目錄的 `grid_used.json`，禁止 idx*D（`tests/nu_lint_coords.py`）。
4. 失敗時：寫假設 → 最小實驗 → 修正 → 重跑全部測試（含舊測試），紀錄寫入 `results/nu_failures.json` 與報告的「未通過紀錄」。不得放寬門檻、不得改理論值、不得挑量測區。
5. 每一關產出驗證表（項目 | 理論值 | 量測值 | 門檻 | PASS/FAIL）與圖，數字附單位與網格資訊（Δ_min、Δ_max、r_max）。
6. 前一關未全部 PASS，不得宣稱下一關通過；診斷用的執行標「不構成放行」。
7. 本專案不涉及生成式 AI 輸出，沒有 Streaming 需求。

## 執行環境

- 求解器與 Python 測試在 WSL Ubuntu 執行（gcc 13.3、numpy、matplotlib）；專案路徑 `/mnt/c/Users/Sasha/onedrive_copy/coding/claude/FDTD/3D_oblique`。
- Meep 1.34：`MEEP_PYTHON`（預設 `~/micromamba/envs/mp/bin/python`）。
- 舊回歸：`python3 tests/run_all_regression.py --fresh`；新回歸：`python3 tests/run_nu_regression.py`。

## 凍結雜湊

- tests/thresholds_nu.json: b49f0e3e4fa15d67a9bbb4066c96311525e4da979265bb5de6a4a4bafb1cd839（含 amendments D9–D22；D17–D19 經失敗流程；D20–D22 為使用者決策，D22 於 2026-10-01）
- results/nu_predictions.json: 7cbfa0900a2c66619fc9d3887fef3ad8a693ed55c5219a00b64960182dd70b87（v2/ = D12 的 Δt 規則）
