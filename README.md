# Endurance Pharmacology Explorer

Predicting how endurance-relevant drugs bind their protein targets using [Chai-1](https://github.com/chaidiscovery/chai-lab).

## Why I Built This

I'm a competitive triathlete who studies molecular biology to understand what drives endurance performance. I wanted to see what Chai-1 could tell me about the compounds that alter aerobic capacity at the protein level — and where the model hits its limits.

This project runs 6 drug-target predictions through Chai-1: 1 control with a known crystal structure, 3 frontier small molecules, and 2 peptide-protein interactions that push the model into genuinely difficult territory.

## The Pipeline

```
targets.yaml → prepare.py → FASTA files → Chai-1 (Colab, A100) → CIF structures → analyze.py
```

Single YAML config defines all 6 pairs (SMILES, sequences, UniProt IDs, reference PDBs). One command generates Chai-1-compatible FASTA files. Inference runs on Google Colab with MSA server enabled. Analysis scores each prediction and validates against known structures where available.

## Interactive 3D Models

Every prediction is viewable as an interactive 3D structure — rotate, zoom, inspect. No install needed.

**[Browse all 6 models on GitHub Pages →](https://dhruvsingh0905.github.io/endurance-pharma-explorer/)**

For the full analysis with per-residue confidence heatmaps and binding pocket detail, run `notebooks/explorer.ipynb` locally after cloning.

## The Roster

### Control

| Drug | Target | ipTM | Purpose |
|------|--------|------|---------|
| Ivabradine | HCN4 channel | 0.128 | Pipeline validation against cryo-EM structure PDB 8OFI |

### Tier 1 — Strong Predictions

**ITPP → Hemoglobin** | ipTM: **0.849** | pTM: 0.889 | Protein pLDDT: 92.0

ITPP binds allosterically to the hemoglobin tetramer, stabilizing the T-state to right-shift the O2 dissociation curve. Chai-1 folds the tetramer near-perfectly and places the ligand with the highest binding confidence of any target.

**[View 3D Model →](https://dhruvsingh0905.github.io/endurance-pharma-explorer/models/itpp_hemoglobin.html)**

**Cardarine (GW501516) → PPARδ** | ipTM: **0.744** | Ligand pLDDT: 88.5

PPARδ agonist that upregulates fatty acid oxidation and drives mitochondrial biogenesis via PGC-1α. Ligand pLDDT of 88.5 — the model places Cardarine in the ligand-binding domain with exceptional confidence.

**[View 3D Model →](https://dhruvsingh0905.github.io/endurance-pharma-explorer/models/cardarine_ppard.html)**

### Tier 2 — Solid

**Roxadustat (FG-4592) → PHD2** | ipTM: **0.464** | Ligand pLDDT: 72.7

Competitive inhibitor of prolyl hydroxylase PHD2. Blocks the oxygen sensor that degrades HIF-1α, mimicking hypoxia at sea level — oral altitude training. Most consistent binding pose of any target (16.9A spread across 5 models).

**[View 3D Model →](https://dhruvsingh0905.github.io/endurance-pharma-explorer/models/roxadustat_phd2.html)**

### Tier 3 — The Frontier

These are where it gets interesting.

**Ivabradine → HCN4** | ipTM: 0.128 | Protein pLDDT: 86.7 | Ligand pLDDT: 18.3

Protein fold is excellent. Ligand placement is random — 42A spread across 5 models. **Why:** Ivabradine is an open-channel blocker. The binding site only exists when the HCN4 pore is in the open conformation. Static structure prediction can't capture conformational gating.

**[View 3D Model →](https://dhruvsingh0905.github.io/endurance-pharma-explorer/models/ivabradine_hcn4.html)**

**TB-500 → G-Actin** | ipTM: 0.172 | Protein pTM: **0.921**

Actin fold is the best of any target (pTM 0.92). But the 7-residue peptide lands in a different spot every time (52A spread). **Why:** A heptapeptide has too many degrees of freedom and too few intramolecular contacts for the model to anchor. The interaction surface is shallow.

**[View 3D Model →](https://dhruvsingh0905.github.io/endurance-pharma-explorer/models/tb500_actin.html)**

**BPC-157 → VEGFR2** | ipTM: 0.114 | Spread: 62.7A

15-residue peptide docking into a 745-residue extracellular domain. No known co-crystal structure — genuinely novel prediction. Largest protein + most flexible drug + no structural prior = maximum uncertainty. This is exactly the class of problem Chai Discovery is working to solve.

**[View 3D Model →](https://dhruvsingh0905.github.io/endurance-pharma-explorer/models/bpc157_vegfr2.html)**

## What This Tells Us About Chai-1

| Excels at | Struggles with |
|-----------|---------------|
| Protein fold prediction (pTM 0.63–0.92 across all 6) | Small molecule placement in dynamic pockets (ion channels) |
| Small molecule docking into defined binding pockets (PPARδ, PHD2) | Short flexible peptide docking (7–15 residues) |
| Multi-chain complex assembly (hemoglobin tetramer) | Large receptor + peptide combinations |
| Allosteric binding with structurally defined pockets (ITPP) | Conformationally dependent binding sites |

## Full Results

| Target | Category | pTM | ipTM | Prot pLDDT | Lig pLDDT | Lig Spread |
|--------|----------|-----|------|-----------|----------|-----------|
| ITPP → Hemoglobin | Allosteric | 0.889 | **0.849** | 92.0 | 44.7 | 43.4A |
| Cardarine → PPARδ | Lock-and-key | 0.734 | **0.744** | 74.8 | 88.5 | 38.3A |
| Roxadustat → PHD2 | Enzyme inhibitor | 0.628 | **0.464** | 68.3 | 72.7 | 16.9A |
| TB-500 → Actin | Peptide | 0.921 | 0.172 | 83.7 | 31.9 | 52.4A |
| Ivabradine → HCN4 | Ion channel | 0.765 | 0.128 | 86.7 | 18.3 | 42.1A |
| BPC-157 → VEGFR2 | Peptide (novel) | 0.490 | 0.114 | 79.9 | 39.9 | 62.7A |

*ipTM = interface confidence (0–1). Lig Spread = max pairwise ligand COM distance across 5 samples (lower = more consistent).*

## Reproduce

```bash
pip install -r requirements.txt
python src/prepare.py              # Generate FASTA files
pytest tests/ -v                    # Validate pipeline
# Run colab/fold.ipynb on Colab with A100
# Copy results back to results/
python src/analyze.py               # Score + validate
jupyter notebook notebooks/explorer.ipynb  # 3D visualization
```

## Structure

```
data/targets.yaml       ← All 6 pairs: SMILES, sequences, UniProt IDs, references
src/prepare.py          ← FASTA generation
src/analyze.py          ← Confidence scoring + RMSD validation
src/visualize.py        ← py3Dmol rendering
tests/                  ← Pipeline validation (17 tests)
colab/fold.ipynb        ← Chai-1 inference (Colab GPU)
notebooks/explorer.ipynb ← Interactive 3D visualization
```

## Built With

[Chai-1](https://github.com/chaidiscovery/chai-lab) | [py3Dmol](https://github.com/avirshup/py3dmol) | Sequences from [UniProt](https://www.uniprot.org/) | SMILES from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) | References from [RCSB PDB](https://www.rcsb.org/)
