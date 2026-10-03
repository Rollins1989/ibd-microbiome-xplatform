"""Regression tests for strict MAG tool-failure handling."""
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "collect_bins.py"
spec = importlib.util.spec_from_file_location("collect_bins", SCRIPT)
collect_bins = importlib.util.module_from_spec(spec)
spec.loader.exec_module(collect_bins)


def test_swallowed_binner_failure_is_reported(tmp_path):
    outdir = tmp_path / "results" / "mags" / "bins_passing"
    log = tmp_path / "results" / "logs" / "binning"
    log.mkdir(parents=True)
    (log / "S1.metabat2.log").write_text(
        "WARNING: MetaBAT2 returned non-zero (typically: no bins formed)\n"
    )

    with pytest.raises(SystemExit, match="swallowed tool failure"):
        collect_bins.validate_mag_logs(outdir, ["S1"])


def test_empty_binner_output_is_not_itself_an_error(tmp_path):
    outdir = tmp_path / "results" / "mags" / "bins_passing"
    log = tmp_path / "results" / "logs" / "binning"
    log.mkdir(parents=True)
    (log / "S1.metabat2.log").write_text("MetaBAT2 completed with no bins.\n")

    collect_bins.validate_mag_logs(outdir, ["S1"])
