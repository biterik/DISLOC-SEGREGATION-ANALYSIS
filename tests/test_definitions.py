import pytest

from disegreg.cli import main
from disegreg.definitions import available, load_definitions


def test_d1_loads_with_decided_values():
    d = load_definitions("D1")
    assert d.label == "D1" and "draft" in d.status
    assert len(d.sha256) == 64
    reg = d["regions"]
    assert reg["far_field_radius"] == 60.0       # decided 22 Nov
    assert reg["surface_exclusion"] == 15.0
    assert reg["core_radius"] == 5.0 and reg["shell_radius"] == 10.0 and reg["fault_half_width"] == 4.0
    assert d["distance"]["polyline_spacing"] == 0.25
    fb = d["distance"]["fallback_axis"]
    assert fb["min_slab_fraction"] == 0.5 and fb["active_below_dxa_coverage"] == 0.5
    assert d["order"]["shell2_cutoff_a"] == 1.112
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
