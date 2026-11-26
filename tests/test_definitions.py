import pytest

from disegreg.cli import main
from disegreg.definitions import available, load_definitions


def test_d1_loads_with_decided_values():
    d = load_definitions("D1")
    assert d.label == "D1" and d.status.startswith("rev. 4")
    assert len(d.sha256) == 64
    reg = d["regions"]
    assert reg["far_field_radius"] == 60.0       # decided 22 Nov
    assert reg["surface_exclusion"] == 15.0
    assert reg["core_radius"] == 5.0 and reg["shell_radius"] == 10.0 and reg["fault_half_width"] == 4.0
    assert d["distance"]["polyline_spacing"] == 0.25
    fb = d["distance"]["fallback_axis"]
    assert fb["min_slab_fraction"] == 0.5 and fb["active_below_dxa_coverage"] == 0.5
    assert d["order"]["shell2_cutoff_a"] == 1.112
    comp = d["composition"]                      # rev. 4 (Erik 24 Nov)
    assert comp["headline_radius"] == 20.0 and comp["check_radius"] == 40.0 and comp["partial_merge_distance"] == 2.0
    assert comp["gamma_curve_radii"][0] == 5 and comp["gamma_curve_radii"][-1] == 60
    ref = d["reference"]
    assert ref["without_bulk_control_report_against"] == ["c_B", "x_iso"]
    assert ref["paper_hold_without_bulk_control"] == [{"potential_contains": "ACE", "host": "Ni", "solute": "Cu"}]
    assert d["timeseries"]["drift_criterion"] == "change_se"
    assert "D1" in available()


def test_definitions_are_read_only_copies():
    d = load_definitions()
    d["regions"]["far_field_radius"] = 1.0
    assert d["regions"]["far_field_radius"] == 60.0


def test_unknown_label():
    with pytest.raises(FileNotFoundError):
        load_definitions("D99")


def test_cli_version(capsys):
    assert main(["version"]) == 0
    assert "definition set D2" in capsys.readouterr().out


def test_d2_is_d1_plus_new_sections():
    """D2 (Erik 25 Nov 2026) = D1 rev. 4 unchanged + the sections decided from P1-P4 and E1; P5 is not in D2."""
    import yaml

    d1, d2 = load_definitions("D1"), load_definitions("D2")
    assert d2.label == "D2" and d2["extends"] == "D1" and available() == ["D1", "D2"]
    for k, v in d1.params.items():
        if k not in ("label", "status"):
            assert d2[k] == v, k                                   # every D1 parameter, value for value
    assert set(d2.params) - set(d1.params) == {"extends", "mclean", "hysteresis", "bulk_control", "alongline",
                                                    "diagnostic_iv"}
    assert "ordered_domains" not in d2.params                      # P5 stays a proposal
    assert d2["mclean"] == {"core_site_radius": 5.0, "missing_sites": "tabulated_only", "min_coverage": 0.9}
    assert d2["hysteresis"] == {"factor": 2.0, "quantities": [["Gamma", "20", "r"], ["c", "C", ""]]}
    assert d2["bulk_control"] == {"propagate_x_ref_error": True, "min_run_length_fraction": 0.95,
                                  "swap_rate_tolerance": 0.10, "smaller_cell_allowed": True,
                                  "smaller_cell_for_diagnostic_iv": False}           # rev. 3: 0.95 / 10 % for every bulk control
    assert d2["alongline"] == {"regions": ["C"], "bin_width": 5.0, "p_threshold": 0.01,
                               "n_eff_from": "block_averaging_c_C", "n_eff_without_plateau": 1}   # rev. 3 (N3)
    assert d2.status.startswith("rev. 3")
    assert d2["diagnostic_iv"] == {"not_tested_without_bulk_control": "all"}       # section 17 (Erik 25 Nov, N1)
    # D1.yaml itself is the released file (v0.1.0), byte for byte
    import hashlib
    from pathlib import Path

    import disegreg.definitions as dd
    raw = (Path(dd.__file__).parent / "D1.yaml").read_bytes()
    assert hashlib.sha256(raw).hexdigest() == d1.sha256 and yaml.safe_load(raw)["label"] == "D1"


def test_default_is_d2():
    from disegreg.definitions import DEFAULT_LABEL, label_number

    assert DEFAULT_LABEL == "D2" and load_definitions().label == "D2"
    assert label_number("D2") == 2 and label_number("D10") == 10 and label_number("P-v0.2") is None
