"""Validate FASTA generation from targets.yaml."""

import tempfile
from pathlib import Path

import pytest
import yaml

from src.prepare import generate_fasta, load_targets, write_all_fastas

TARGETS_PATH = Path(__file__).parent.parent / "data" / "targets.yaml"


@pytest.fixture
def targets():
    return load_targets(TARGETS_PATH)


def test_generates_six_fasta_files():
    with tempfile.TemporaryDirectory() as tmpdir:
        paths = write_all_fastas(TARGETS_PATH, Path(tmpdir))
        assert len(paths) == 6
        for p in paths:
            assert p.exists()
            assert p.suffix == ".fasta"
            assert p.stat().st_size > 0


def test_fasta_header_format(targets):
    for t in targets:
        fasta = generate_fasta(t)
        lines = fasta.strip().split("\n")
        for line in lines:
            if line.startswith(">"):
                assert "|name=" in line, f"Missing |name= in header: {line}"
                assert line.startswith(">protein|") or line.startswith(">ligand|"), (
                    f"Invalid header type: {line}"
                )


def test_smiles_drugs_use_ligand_header(targets):
    for t in targets:
        if t["drug"]["kind"] == "smiles":
            fasta = generate_fasta(t)
            assert ">ligand|name=" in fasta, f"{t['name']} SMILES drug should use >ligand| header"


def test_peptide_drugs_use_protein_header(targets):
    for t in targets:
        if t["drug"]["kind"] == "peptide_sequence":
            fasta = generate_fasta(t)
            # Should have at least 2 >protein| headers (target + peptide)
            protein_headers = [l for l in fasta.split("\n") if l.startswith(">protein|")]
            assert len(protein_headers) >= 2, (
                f"{t['name']} peptide should have >=2 protein headers, got {len(protein_headers)}"
            )


def test_hemoglobin_has_four_chains_plus_ligand(targets):
    hb = next(t for t in targets if t["name"] == "itpp_hemoglobin")
    fasta = generate_fasta(hb)
    headers = [l for l in fasta.split("\n") if l.startswith(">")]
    protein_headers = [h for h in headers if h.startswith(">protein|")]
    ligand_headers = [h for h in headers if h.startswith(">ligand|")]
    assert len(protein_headers) == 4, f"Expected 4 protein chains, got {len(protein_headers)}"
    assert len(ligand_headers) == 1, f"Expected 1 ligand, got {len(ligand_headers)}"


def test_sequences_have_no_whitespace(targets):
    for t in targets:
        fasta = generate_fasta(t)
        lines = fasta.strip().split("\n")
        for line in lines:
            if not line.startswith(">"):
                assert " " not in line, f"Whitespace in sequence line: {line[:50]}..."
                assert "\t" not in line, f"Tab in sequence line: {line[:50]}..."
