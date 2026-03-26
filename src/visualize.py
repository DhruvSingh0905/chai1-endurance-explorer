"""3D visualization of Chai-1 predictions using py3Dmol."""

from pathlib import Path

import py3Dmol


def render_prediction(cif_path: Path, width: int = 900, height: int = 550) -> py3Dmol.view:
    """Render a structure colored by pLDDT confidence.

    Blue (>90) = high confidence, cyan (70-90) = confident,
    yellow (50-70) = uncertain, orange/red (<50) = low confidence.
    """
    cif_text = cif_path.read_text()
    view = py3Dmol.view(width=width, height=height)
    view.addModel(cif_text, "cif")

    # Protein: cartoon colored by pLDDT (B-factor column)
    view.setStyle(
        {"and": [{"not": {"resn": "LIG2"}}, {"not": {"hetflag": True}}]},
        {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 30, "max": 100}}},
    )

    # Ligand: white sticks
    view.setStyle(
        {"or": [{"resn": "LIG2"}, {"hetflag": True}]},
        {"stick": {"radius": 0.25, "colorscheme": "whiteCarbon"}},
    )

    # Transparent surface near ligand
    view.addSurface(
        py3Dmol.VDW,
        {"opacity": 0.15, "color": "white"},
        {"within": {"distance": 6, "sel": {"or": [{"resn": "LIG2"}, {"hetflag": True}]}}},
    )

    view.zoomTo()
    return view


def render_peptide(cif_path: Path, width: int = 900, height: int = 550) -> py3Dmol.view:
    """Render peptide-protein interaction. Target=pLDDT cartoon, peptide=pink."""
    cif_text = cif_path.read_text()
    view = py3Dmol.view(width=width, height=height)
    view.addModel(cif_text, "cif")

    # Chain A (target): cartoon by pLDDT
    view.setStyle(
        {"chain": "A"},
        {"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 30, "max": 100}}},
    )

    # Chain B (peptide): pink cartoon + sticks
    view.setStyle(
        {"chain": "B"},
        {"cartoon": {"color": "#FF1493", "opacity": 0.9}, "stick": {"radius": 0.2, "color": "#FF1493"}},
    )

    view.zoomTo()
    return view


def render_comparison(pred_path: Path, ref_path: Path, width: int = 900, height: int = 550) -> py3Dmol.view:
    """Overlay predicted (blue) vs reference (green) structure."""
    view = py3Dmol.view(width=width, height=height)

    pred_text = pred_path.read_text()
    fmt = "cif" if pred_path.suffix == ".cif" else "pdb"
    view.addModel(pred_text, fmt)
    view.setStyle({"model": 0}, {"cartoon": {"color": "lightblue", "opacity": 0.7}})
    view.setStyle({"model": 0, "hetflag": True}, {"stick": {"colorscheme": "blueCarbon", "radius": 0.2}})

    ref_text = ref_path.read_text()
    ref_fmt = "cif" if ref_path.suffix == ".cif" else "pdb"
    view.addModel(ref_text, ref_fmt)
    view.setStyle({"model": 1}, {"cartoon": {"color": "lightgreen", "opacity": 0.7}})
    view.setStyle({"model": 1, "hetflag": True}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}})

    view.zoomTo()
    return view
