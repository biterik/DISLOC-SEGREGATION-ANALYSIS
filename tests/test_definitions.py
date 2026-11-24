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
    assert "definition set D1" in capsys.readouterr().out
