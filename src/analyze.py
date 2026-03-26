"""Analyze Chai-1 predictions: confidence scores, RMSD vs reference PDBs."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


def load_npz_scores(npz_path: Path) -> dict:
    """Load Chai-1 score file (.npz) and return key metrics."""
    data = np.load(npz_path)
    return {
        "aggregate_score": float(data["aggregate_score"][0]),
        "ptm": float(data["ptm"][0]),
        "iptm": float(data["iptm"][0]),
        "has_clashes": bool(data["has_inter_chain_clashes"][0]),
    }


def parse_plddt_from_pdb(pdb_path: Path) -> list[float]:
    """Extract per-residue pLDDT scores from B-factor column of a PDB file."""
    scores = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    bfactor = float(line[60:66].strip())
                    scores.append(bfactor)
                except (ValueError, IndexError):
                    continue
    return scores


def extract_ligand_coords(pdb_path: Path) -> np.ndarray:
    """Extract HETATM coordinates (ligand atoms) from a PDB file."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except (ValueError, IndexError):
                    continue
    return np.array(coords)


def compute_ligand_rmsd(pred_pdb: Path, ref_pdb: Path) -> float:
    """Compute RMSD between ligand atoms in predicted vs reference PDB.

    Uses simple coordinate RMSD (no alignment) — assumes same atom ordering.
    For proper analysis, use MDAnalysis or BioPython with alignment.
    """
    pred_coords = extract_ligand_coords(pred_pdb)
    ref_coords = extract_ligand_coords(ref_pdb)

    if len(pred_coords) == 0 or len(ref_coords) == 0:
        return float("nan")

    n = min(len(pred_coords), len(ref_coords))
    diff = pred_coords[:n] - ref_coords[:n]
    return float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))


def find_score_file(pair_dir: Path, model_idx: int) -> Path | None:
    """Find the NPZ score file for a given model index, handling renamed files."""
    patterns = [
        f"scores.model_idx_{model_idx}.npz",
        f"scores_model_index_{model_idx}.npz",
        f"scores_model_idx_{model_idx}.npz",
        f"score_model_index_{model_idx}.npz",
    ]
    for pattern in patterns:
        path = pair_dir / pattern
        if path.exists():
            return path
    # Fallback: search for any npz with the model index
    for f in pair_dir.glob("*scores*npz"):
        if str(model_idx) in f.name:
            return f
    return None


def analyze_pair(pair_dir: Path) -> dict:
    """Analyze all predictions for a single drug-target pair."""
    # Look for CIF files (Chai-1 output) or PDB files
    cif_files = sorted(pair_dir.glob("*.cif"))
    pdb_files = sorted(pair_dir.glob("*.pdb"))
    structure_files = cif_files or pdb_files

    if not structure_files:
        return {"name": pair_dir.name, "n_predictions": 0, "best_score": None}

    results = []
    for i, sf in enumerate(structure_files):
        # Try to load NPZ scores
        score_file = find_score_file(pair_dir, i)
        if score_file:
            scores = load_npz_scores(score_file)
        else:
            # Fallback: try to parse from B-factors (PDB only)
            if sf.suffix == ".pdb":
                plddt = parse_plddt_from_pdb(sf)
                scores = {
                    "aggregate_score": float(np.mean(plddt)) / 100 if plddt else 0.0,
                    "ptm": 0.0,
                    "iptm": 0.0,
                    "has_clashes": False,
                }
            else:
                scores = {"aggregate_score": 0.0, "ptm": 0.0, "iptm": 0.0, "has_clashes": False}

        results.append({
            "file": sf.name,
            "aggregate_score": scores["aggregate_score"],
            "ptm": scores["ptm"],
            "iptm": scores["iptm"],
            "has_clashes": scores["has_clashes"],
        })

    results.sort(key=lambda r: r["aggregate_score"], reverse=True)
    best = results[0]

    return {
        "name": pair_dir.name,
        "n_predictions": len(structure_files),
        "best_file": best["file"],
        "best_score": best["aggregate_score"],
        "best_ptm": best["ptm"],
        "best_iptm": best["iptm"],
        "has_clashes": best["has_clashes"],
        "all_results": results,
    }


def generate_summary(results_dir: Path) -> list[dict]:
    """Generate analysis summary for all pairs in results directory."""
    summaries = []
    for pair_dir in sorted(results_dir.iterdir()):
        if pair_dir.is_dir() and (
            any(pair_dir.glob("*.cif")) or any(pair_dir.glob("*.pdb"))
        ):
            summaries.append(analyze_pair(pair_dir))
    return summaries


def validate_against_reference(
    pred_pdb: Path, ref_pdb: Path, target_rmsd: float
) -> dict:
    """Validate a prediction against a known reference structure."""
    rmsd = compute_ligand_rmsd(pred_pdb, ref_pdb)
    return {
        "rmsd": rmsd,
        "target": target_rmsd,
        "passes": rmsd <= target_rmsd if not np.isnan(rmsd) else False,
    }


def main():
    parser = argparse.ArgumentParser(description="Analyze Chai-1 prediction results")
    parser.add_argument(
        "--results",
        type=Path,
        default=Path(__file__).parent.parent / "results",
    )
    parser.add_argument(
        "--targets",
        type=Path,
        default=Path(__file__).parent.parent / "data" / "targets.yaml",
    )
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    print("Analyzing Chai-1 predictions...\n")
    summaries = generate_summary(args.results)

    if not summaries:
        print("No results found. Run Chai-1 inference first (see colab/fold.ipynb).")
        return

    # Load targets for metadata
    with open(args.targets) as f:
        targets = {t["name"]: t for t in yaml.safe_load(f)["targets"]}

    # Print summary table
    print(f"{'Pair':<25} {'#':>3} {'Score':>7} {'pTM':>7} {'ipTM':>7} {'Clashes':>8}")
    print("-" * 65)

    rows = []
    for s in summaries:
        t = targets.get(s["name"], {})
        score_str = f"{s['best_score']:.4f}" if s["best_score"] is not None else "—"
        ptm_str = f"{s['best_ptm']:.4f}" if s.get("best_ptm") else "—"
        iptm_str = f"{s['best_iptm']:.4f}" if s.get("best_iptm") else "—"
        clash_str = "YES" if s.get("has_clashes") else "no"

        print(f"{s['name']:<25} {s['n_predictions']:>3} {score_str:>7} {ptm_str:>7} {iptm_str:>7} {clash_str:>8}")

        rows.append({
            "name": s["name"],
            "category": t.get("category", ""),
            "n_predictions": s["n_predictions"],
            "best_file": s.get("best_file", ""),
            "aggregate_score": s["best_score"],
            "ptm": s.get("best_ptm"),
            "iptm": s.get("best_iptm"),
            "has_clashes": s.get("has_clashes"),
        })

    # Save CSV
    out_path = args.output or (args.results / "summary.csv")
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print(f"\nSummary saved to {out_path}")

    # Interpretation
    print("\n--- Score Interpretation ---")
    print("aggregate_score = 0.8*ipTM + 0.2*pTM")
    print("pTM: protein fold confidence (0-1, >0.5 = reasonable fold)")
    print("ipTM: interface/binding confidence (0-1, >0.5 = confident binding)")
    print("For ligand docking: ipTM is the key metric.\n")

    for s in summaries:
        t = targets.get(s["name"], {})
        iptm = s.get("best_iptm", 0)
        ptm = s.get("best_ptm", 0)
        if ptm > 0.5 and iptm < 0.3:
            print(f"  {s['name']}: Protein fold OK (pTM={ptm:.3f}) but binding uncertain (ipTM={iptm:.3f})")
            print(f"    → Consider using MSAs or restraints to improve binding prediction")
        elif ptm > 0.5 and iptm > 0.3:
            print(f"  {s['name']}: Good fold (pTM={ptm:.3f}) and reasonable binding (ipTM={iptm:.3f})")
        elif ptm < 0.5:
            print(f"  {s['name']}: Low fold confidence (pTM={ptm:.3f}) — protein may be too large or disordered")


if __name__ == "__main__":
    main()
