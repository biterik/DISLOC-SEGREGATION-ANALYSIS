"""Plotting module (house style): fonts >= 13 pt, colour-blind-safe fixed colours, symmetric scales for signed maps."""

import numpy as np
import pytest

pytest.importorskip("matplotlib")

from disegreg import plotting  # noqa: E402
from disegreg.cli import main  # noqa: E402
from disegreg.results import GAMMA_CURVE_COLUMNS, STAGE3_COLUMNS, write_table  # noqa: E402
from test_across import write_run  # noqa: E402
from test_manifest_stage2 import demo  # noqa: E402


def test_symmetric_limits_and_signed_map(tmp_path):
    assert plotting.symmetric_limits([-0.2, 0.5, np.nan]) == (-0.5, 0.5)
    assert plotting.symmetric_limits([0.0, 0.0]) == (-1.0, 1.0)
    assert plotting.symmetric_limits([-3.0, 1.0]) == (-3.0, 3.0)
    v = np.array([[0.1, -0.4], [0.2, 0.05]])
    fig = plotting.signed_map([0, 1, 2], [0, 1, 2], v, "Y (Å)", "Z (Å)", "excess")
    mesh = fig.axes[0].collections[0]
    assert mesh.norm.vmin == -0.4 and mesh.norm.vmax == 0.4          # symmetric about zero
    cmap = mesh.cmap
    mid = cmap(0.5)[:3]
    assert np.allclose(mid, np.array([0xf0, 0xef, 0xec]) / 255, atol=0.02)  # grey, not a hue, at zero
    assert plotting.check_fonts(fig) == []
    assert plotting.save(fig, tmp_path / "map", ("png", "pdf"))[1].is_file()


def test_fixed_colour_per_character_and_font_rule(tmp_path):
    rows = [{"character": ch, "potential": "P", "delta": d, "value": v, "error": 0.01}
            for ch, v0 in (("d60", 0.3), ("d0", 0.1)) for d, v in ((-0.5, v0), (-0.25, 2 * v0))]
    fig = plotting.vs_delta(rows, "c_C")
    handles, labels = fig.axes[0].get_legend_handles_labels()
    colours = {lab: h[0].get_color() for h, lab in zip(handles, labels)}
    assert colours == {plotting.CHARACTER_NAME["d0"]: plotting.CHARACTER_COLOUR["d0"],
                       plotting.CHARACTER_NAME["d60"]: plotting.CHARACTER_COLOUR["d60"]}  # colour follows the entity
    assert plotting.check_fonts(fig) == []
    fig.axes[0].set_xlabel("tiny", fontsize=8)
    with pytest.raises(ValueError, match="fonts below"):
        plotting.save(fig, tmp_path / "x")


def test_palette_is_the_validated_okabe_ito_order():
    assert plotting.PALETTE[:4] == ("#0072B2", "#D55E00", "#009E73", "#E69F00")
    assert len(set(plotting.MARKERS)) == len(plotting.MARKERS)


def test_plot_cli_all_kinds(tmp_path):
    # gamma curve
    curve = [{"run_id": "RUN", "R": R, "basis": b, "ref": ref, "role": "", "mean": 0.1 * R, "sem": 0.05, "n_samples": 3,
              "flat": True, "definition_set": "D1", "tool_version": "t"}
             for R in (5.0, 10.0, 20.0, 40.0) for b in ("r", "r_dxa") for ref in ("c_B", "x_iso")]
    write_table(tmp_path / "RUN_gammaR_D1.csv", curve, {}, columns=GAMMA_CURVE_COLUMNS)
    assert main(["plot", "gamma", str(tmp_path / "RUN_gammaR_D1.csv"), "--format", "png"]) == 0
    assert (tmp_path / "RUN_gammaR_D1.png").is_file()
    # time series with its Stage-3 window
    m = demo()
    m["run_id"] = "RUN"
    s2 = [{"run_id": "RUN", "snapshot": f"s{k}", "timestep": 100000 * k, "quantity": "c", "region": "C", "basis": "",
           "ref": "", "value": 0.2 + 0.01 * np.sin(k), "definition_set": "D1", "tool_version": "t"} for k in range(1, 11)]
    write_table(tmp_path / "RUN_stage2_D1.csv", s2, {"manifest_data": m})
    s3 = [{"run_id": "RUN", "quantity": "c", "region": "C", "basis": "", "ref": "", "mean": 0.2, "sem": 0.003,
           "n_samples": 5, "sem_plateau": True, "change": 0.001, "change_se": 0.004, "flat": True, "flat_rev3_sem": True,
           "t_first": 600000.0, "t_last": 1000000.0, "definition_set": "D1", "tool_version": "t"}]
    write_table(tmp_path / "RUN_stage3_D1.csv", s3, {"manifest_data": m, "stage3": {
        "time_unit": "timestep", "window_t_start": 600000.0, "window_t_end": 1000000.0}}, columns=STAGE3_COLUMNS)
    assert main(["plot", "timeseries", str(tmp_path / "RUN_stage2_D1.csv"), "--quantity", "c:C"]) == 0
    assert (tmp_path / "RUN_stage2_D1_c-C.png").is_file()
    # results table -> delta and regions; Stage 4 -> uptake
    files = [write_run(tmp_path / "s3", f"r{k}", d, 2e4 * (1e-3 * np.exp(d)) ** 3, 1e-3 * np.exp(d), 1e-3 * np.exp(d))
             for k, d in enumerate((-1.0, -0.5, -0.25))] if (tmp_path / "s3").mkdir() is None else []
    assert main(["table", str(tmp_path / "s3"), "--out", str(tmp_path / "all.csv")]) == 0
    assert main(["plot", "delta", str(tmp_path / "all.csv"), "--quantity", "c:C", "--log"]) == 0
    assert main(["plot", "regions", str(tmp_path / "all.csv"), "--run", "r1"]) == 0
    assert main(["across", str(tmp_path / "s3"), "--out", str(tmp_path / "s4.csv")]) == 0
    assert main(["plot", "uptake", str(tmp_path / "s4.csv"), "--format", "png", "pdf"]) == 0
    assert (tmp_path / "s4_uptake.pdf").is_file() and len(files) == 3
