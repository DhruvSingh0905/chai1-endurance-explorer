"""Validate targets.yaml schema, SMILES, and protein sequences."""

import re
from pathlib import Path

import pytest
import yaml

TARGETS_PATH = Path(__file__).parent.parent / "data" / "targets.yaml"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
EXPECTED_NAMES = {
    "ivabradine_hcn4",
    "cardarine_ppard",
    "itpp_hemoglobin",
    "roxadustat_phd2",
    "bpc157_vegfr2",
    "tb500_actin",
}


@pytest.fixture
def targets():
    with open(TARGETS_PATH) as f:
        return yaml.safe_load(f)["targets"]


def clean(seq: str) -> str:
    return "".join(seq.split())


def test_all_six_targets_present(targets):
    names = {t["name"] for t in targets}
    assert names == EXPECTED_NAMES, f"Missing or extra targets: {names ^ EXPECTED_NAMES}"


def test_required_fields(targets):
    for t in targets:
        assert "name" in t
        assert "category" in t
        assert "drug" in t
        assert "biology" in t
        drug = t["drug"]
        assert "name" in drug
        assert "kind" in drug
        assert drug["kind"] in ("smiles", "peptide_sequence")
        assert "value" in drug and len(drug["value"]) > 0


def test_smiles_strings_parse(targets):
    """Validate SMILES strings are parseable by RDKit."""
    try:
        from rdkit import Chem
    except ImportError:
        pytest.skip("rdkit not installed")

    for t in targets:
        if t["drug"]["kind"] == "smiles":
            mol = Chem.MolFromSmiles(t["drug"]["value"])
            assert mol is not None, f"Invalid SMILES for {t['drug']['name']}: {t['drug']['value']}"


def test_protein_sequences_valid(targets):
    """All protein sequences contain only valid amino acid characters."""
    for t in targets:
        tgt = t["target"]
        if "chains" in tgt:
            seqs = [(c["name"], clean(c["sequence"])) for c in tgt["chains"]]
        else:
            seqs = [(tgt["name"], clean(tgt["sequence"]))]

        for name, seq in seqs:
            invalid = set(seq) - VALID_AA
            assert not invalid, f"{name} has invalid AA chars: {invalid}"
            assert len(seq) > 10, f"{name} sequence too short: {len(seq)}"


def test_peptide_sequences_valid(targets):
    """Peptide drug sequences are valid amino acid strings."""
    for t in targets:
        if t["drug"]["kind"] == "peptide_sequence":
            seq = clean(t["drug"]["value"])
            invalid = set(seq) - VALID_AA
            assert not invalid, f"{t['drug']['name']} has invalid AA chars: {invalid}"


def test_uniprot_ids_wellformed(targets):
    """UniProt IDs match expected format (P/Q/O followed by alphanumeric)."""
    pattern = re.compile(r"^[A-Z][0-9A-Z]{5}(/[A-Z][0-9A-Z]{5})*$")
    for t in targets:
        uid = t["target"]["uniprot"]
        assert pattern.match(uid), f"Malformed UniProt ID for {t['name']}: {uid}"


def test_categories_valid(targets):
    valid = {"control", "frontier_small_molecule", "frontier_peptide"}
    for t in targets:
        assert t["category"] in valid, f"Invalid category: {t['category']}"


def test_control_has_rmsd_target(targets):
    for t in targets:
        if t["category"] == "control":
            assert t.get("rmsd_target") is not None, f"Control {t['name']} missing rmsd_target"
            assert t.get("reference_pdb") is not None, f"Control {t['name']} missing reference_pdb"
