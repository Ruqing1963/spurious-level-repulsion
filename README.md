# spurious-level-repulsion

**Spurious level repulsion in geophysical random-matrix analysis: declustering, seasonal gating, and the Wishart null**

Ruqing Chen — GUT Geoservice Inc., Montréal, Québec, Canada

A short methodological note. The consecutive spacing ratio ⟨r⟩ is increasingly
used to claim Wigner–Dyson (GOE) "level repulsion" in geophysical and climate
data. ⟨r⟩ is two-ended, and three mundane data structures mimic GOE without any
repulsion. This repo documents all three with controlled demonstrations and gives
the correct controls.

## The three hazards

| Hazard | Spurious result | Correct diagnostic |
|---|---|---|
| **Declustering** | pure random points → ⟨r⟩ ≈ GOE (0.39 → 0.57 as min-sep grows) | raw ⟨r⟩, no thinning; CV |
| **Seasonal gating** | bimodal intervals → ⟨r⟩ = 0.48 (permutation-invariant artifact) | shuffle test; de-seasonalized counts |
| **Wishart sampling** | iid data → bulk ⟨r⟩ ≈ GOE (0.53) for all P,T | eigenvalue density vs Marchenko–Pastur |

**Worked example (SSW).** For 41 major sudden stratospheric warmings (NOAA SSW
Compendium, JRA-55, 1958–2023): the naive interval ⟨r⟩ = 0.484 ≈ shuffle 0.475
is a seasonal artifact. The correct de-seasonalized test (per-winter event count)
shows mild **under-dispersion** (variance/mean = 0.67, 95% CI [0.49, 0.86]
excludes 1): a **weak refractory clock, not level repulsion**.

## Reproduce

```bash
pip install -r requirements.txt
cd code
python spurious_goe.py --outdir ..
```

Writes `results/summary.csv` and the two figures. All inputs are embedded
(public SSW dates) or synthetic; no data download is required.

## Layout

```
spurious-level-repulsion/
├── code/spurious_goe.py          # all three hazards + SSW de-seasonalized + figures
├── data/
│   ├── ssw_central_dates.csv     # 41 SSW central dates (NOAA Compendium; also embedded in code)
│   └── README.md                 # data provenance
├── results/summary.csv           # all reported numbers
├── figures/fig_hazards.{pdf,png} # Fig 1 (three hazards)
│         fig_ssw_counts.{pdf,png}# Fig 2 (SSW under-dispersion)
├── paper/spurious_goe_paper.{pdf,tex}
├── README.md  LICENSE (CC BY 4.0)  requirements.txt
```

## Correct practice

Before reading any ⟨r⟩ above Poisson as level repulsion: (i) no declustering;
(ii) run the permutation (shuffle) test — obs ≈ shuf means a marginal (scalar)
property, not repulsion; (iii) report CV alongside ⟨r⟩; (iv) de-seasonalize gated
events; (v) for correlation matrices, test the eigenvalue density against
Marchenko–Pastur, not the bulk spacing. Genuine Wigner–Dyson repulsion survives
all of these; the impostors do not.

## Data source

SSW central dates: NOAA SSW Compendium (Butler et al. 2017),
https://csl.noaa.gov/groups/csl8/sswcompendium/ .

## Citation

Cite the paper (`paper/spurious_goe_paper.pdf`) and this repository (Zenodo DOI on
release).

## License
CC BY 4.0 (see `LICENSE`).
