"""Generate Chai-1 compatible FASTA files from targets.yaml."""

import argparse
from pathlib import Path

import yaml


def load_targets(targets_path: Path) -> list[dict]:
    with open(targets_path) as f:
        return yaml.safe_load(f)["targets"]


def clean_sequence(seq: str) -> str:
    """Remove whitespace/newlines from YAML multi-line sequences."""
    return "".join(seq.split())


def generate_fasta(target: dict) -> str:
    """Generate a Chai-1 compatible FASTA string for a single target pair."""
    lines = []
    drug = target["drug"]
    tgt = target["target"]

    # Handle multi-chain targets (hemoglobin tetramer)
    if "chains" in tgt:
        for chain in tgt["chains"]:
            lines.append(f">protein|name={chain['name']}")
            lines.append(clean_sequence(chain["sequence"]))
    else:
        lines.append(f">protein|name={tgt['name']}")
        lines.append(clean_sequence(tgt["sequence"]))

    # Add drug/peptide
    if drug["kind"] == "smiles":
        lines.append(f">ligand|name={drug['name']}")
        lines.append(drug["value"])
    elif drug["kind"] == "peptide_sequence":
        lines.append(f">protein|name={drug['name']}")
        lines.append(clean_sequence(drug["value"]))

    return "\n".join(lines) + "\n"


def write_all_fastas(targets_path: Path, output_dir: Path) -> list[Path]:
    """Generate FASTA files for all targets. Returns list of written paths."""
    output_dir.mkdir(parents=True, exist_ok=True)
    targets = load_targets(targets_path)
    written = []

    for target in targets:
        fasta_content = generate_fasta(target)
        out_path = output_dir / f"{target['name']}.fasta"
        out_path.write_text(fasta_content)
        written.append(out_path)
        print(f"  {out_path.name} ({len(fasta_content)} bytes)")

    return written


def main():
    parser = argparse.ArgumentParser(description="Generate FASTA files from targets.yaml")
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path(__file__).parent.parent / "data" / "targets.yaml",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).parent.parent / "data" / "fasta",
    )
    args = parser.parse_args()

    print(f"Reading targets from {args.targets}")
    paths = write_all_fastas(args.targets, args.output)
    print(f"\nGenerated {len(paths)} FASTA files in {args.output}")


if __name__ == "__main__":
    main()
