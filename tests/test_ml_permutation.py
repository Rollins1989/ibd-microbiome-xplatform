"""Regression tests for the ML negative-control label permutation."""
import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "workflow" / "scripts" / "ml.py"
spec = importlib.util.spec_from_file_location("workflow_ml", SCRIPT)
ml = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ml)


def test_varying_longitudinal_labels_preserve_complete_subject_vectors():
    y = np.array([0, 1, 0, 1, 1, 1, 0, 0, 1])
    groups = np.array(["A", "A", "A", "B", "B", "B", "C", "C", "C"])
    original = sorted(tuple(y[groups == subject]) for subject in np.unique(groups))

    permuted = ml.permute_labels(y, groups, np.random.default_rng(11))
    observed = sorted(tuple(permuted[groups == subject]) for subject in np.unique(groups))

    assert observed == original
    assert not np.array_equal(permuted, y)


def test_constant_labels_are_permuted_between_subjects():
    y = np.array([0, 0, 1, 1, 1, 1])
    groups = np.array(["A", "A", "B", "B", "C", "C"])

    permuted = ml.permute_labels(y, groups, np.random.default_rng(3))

    assert sorted(permuted.tolist()) == sorted(y.tolist())
    for subject in np.unique(groups):
        assert len(set(permuted[groups == subject])) == 1
