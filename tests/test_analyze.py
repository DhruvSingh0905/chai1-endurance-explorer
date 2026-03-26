"""Validate analysis utilities (confidence parsing, RMSD computation)."""

import tempfile
from pathlib import Path

import pytest

from src.analyze import (
    compute_ligand_rmsd,
    parse_plddt_from_pdb,
    generate_summary,
)


# Minimal PDB content for testing confidence parsing
MOCK_PDB = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 85.50           N
ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 90.20           C
ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 78.30           C
HETATM    4  C1  LIG B   1       5.000   6.000   7.000  1.00 65.10           C
HETATM    5  C2  LIG B   1       6.000   7.000   8.000  1.00 70.40           C
END
"""

# Two PDBs with known ligand positions for RMSD testing
PDB_PRED = """\
HETATM    1  C1  LIG B   1       1.000   0.000   0.000  1.00 80.00           C
HETATM    2  C2  LIG B   1       0.000   1.000   0.000  1.00 80.00           C
HETATM    3  C3  LIG B   1       0.000   0.000   1.000  1.00 80.00           C
END
"""

PDB_REF = """\
HETATM    1  C1  LIG B   1       1.000   0.000   0.000  1.00 80.00           C
HETATM    2  C2  LIG B   1       0.000   1.000   0.000  1.00 80.00           C
HETATM    3  C3  LIG B   1       0.000   0.000   1.000  1.00 80.00           C
END
"""

PDB_REF_SHIFTED = """\
HETATM    1  C1  LIG B   1       2.000   0.000   0.000  1.00 80.00           C
HETATM    2  C2  LIG B   1       0.000   2.000   0.000  1.00 80.00           C
HETATM    3  C3  LIG B   1       0.000   0.000   2.000  1.00 80.00           C
END
"""


def _write_temp_pdb(content: str) -> Path:
    f = tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False)
    f.write(content)
    f.close()
    return Path(f.name)


def test_parse_plddt_from_pdb():
    pdb_path = _write_temp_pdb(MOCK_PDB)
    scores = parse_plddt_from_pdb(pdb_path)
    assert len(scores) > 0
    assert all(0 <= s <= 100 for s in scores)
    # Check specific values from our mock
    assert abs(scores[0] - 85.50) < 0.01


def test_rmsd_identical_structures():
    pred = _write_temp_pdb(PDB_PRED)
    ref = _write_temp_pdb(PDB_REF)
    rmsd = compute_ligand_rmsd(pred, ref)
    assert abs(rmsd) < 0.01, f"Identical structures should have RMSD ~0, got {rmsd}"


def test_rmsd_shifted_structures():
    pred = _write_temp_pdb(PDB_PRED)
    ref = _write_temp_pdb(PDB_REF_SHIFTED)
    rmsd = compute_ligand_rmsd(pred, ref)
    assert rmsd > 0.5, f"Shifted structures should have RMSD > 0.5, got {rmsd}"


def test_generate_summary():
    """Summary generation works with a mock results directory."""
    with tempfile.TemporaryDirectory() as tmpdir:
        results_dir = Path(tmpdir)
        # Create a mock result
        pair_dir = results_dir / "test_pair"
        pair_dir.mkdir()
        pdb_path = pair_dir / "pred_model_0.pdb"
        pdb_path.write_text(MOCK_PDB)

        summary = generate_summary(results_dir)
        assert len(summary) > 0
        assert "test_pair" in summary[0]["name"]
