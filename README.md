# Protein Surface Roughness using SAFD

## Overview
This repository contains a Python implementation for analyzing protein surface roughness using the **Solvent Accessible Fractal Dimension (SAFD)** method. By measuring how the Solvent Accessible Surface Area (SASA) scales with varying probe radii, this project computes local fractal dimensions to identify and analyze functional clefts, binding pockets, and thermal adaptation characteristics on protein surfaces.

## Features
- **Multi-scale SASA Computation:** Calculates accessible surface areas using varying probe radii (e.g., 1.2Å to 2.1Å).
- **SAFD Extraction:** Determines the fractal dimension ($f$) for each surface atom, representing local structural roughness. 
- **Functional Site Validation:** Compares surface roughness between known binding sites and bulk non-site surfaces across Bound, Apo, and Non-cleft structures.
- **Site-size Analysis:** Investigates the correlation between the size of a site and its corresponding SAFD value.
- **Thermal Adaptation Analysis:** Analyzes differences in surface roughness among Mesophilic, Thermophilic, and Hyperthermophilic proteins.
- **Pocket Prediction:** Identifies potential novel binding pockets based on highly rough surface patches (high $f$ values).

## Project Structure
- `safd.py`: Core algorithm for extracting protein atoms, calculating multi-scale neighbor areas, and fitting SAFD values.
- `sasa.py`: Implementation for computing atomic contact areas.
- `run_full_project.py`: Main pipeline script that downloads necessary PDB files, runs full statistical validation suites, and generates comparison plots.
- `predict_pockets.py`: Script dedicated to predicting potential binding pockets from surface atom data.
- `benchmark.py`, `site_size_analysis.py`, `thermal_analysis.py`: Specialized analysis modules for distinct scientific questions.

## Setup and Installation

**Prerequisites:**
Ensure you have Python 3.7+ installed.

**1. Clone the repository**
```bash
git clone https://github.com/YourUsername/safd_project.git
cd safd_project
```

**2. Install dependencies**
This project requires several scientific Python libraries:
```bash
pip install numpy scipy matplotlib biopython requests
```

## Usage

To run the complete data analysis pipeline:

```bash
python run_full_project.py
```
This script will automatically:
1. Download all required `.pdb` files from the RCSB Protein Data Bank into the `pdbs/` directory.
2. Run validation analyses (Bound vs Apo vs Non-Cleft).
3. Generate visualization plots (`validation_results.png`, `site_size_result.png`, `thermal_adaptation_result.png`).

You can also run individual scripts (e.g., `python safd.py`) to test the SAFD calculation on a single PDB file and list its top 10 roughest surface atoms.

## Output
The pipeline produces comparative bar charts and line graphs highlighting the mean surface roughness (SAFD) of functional binding sites vs. non-site surfaces, which provides structural insights into protein interactions and stability.
