import os
import requests
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.stats import ttest_ind
from Bio.PDB import PDBParser

from safd import compute_safd

DATASET_CONFIG = {
    "bound": ["3DFR", "7HVP", "1RBP", "8ACN", "2GBP"],
    "apo": ["2ACT", "4BLM", "3EBX", "2LIV"],
    "non_cleft": ["1FXD", "1FXA", "5PAL", "1LTE"],
    "small_sites": ["3DFR", "1RBP", "8ACN", "2GBP"],
    "large_sites": ["1ACB", "7HVP"],
    "thermal": {
        "Mesophile (20-45 C)": ["1LZ1", "3DFR"],
        "Thermophile (45-80 C)": ["1V43", "1A5Z"],
        "Hyperthermophile (>80 C)": ["1B44", "1V7S"]
    }
}

os.makedirs("pdbs", exist_ok=True)
os.makedirs("results", exist_ok=True)

def ensure_downloaded(pdb_id):
    path = f"pdbs/{pdb_id}.pdb"
    if not os.path.exists(path):
        print(f"Downloading {pdb_id}...")
        res = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
        if res.status_code == 200:
            with open(path, "w") as f:
                f.write(res.text)
        else:
            print(f"Failed downloading {pdb_id}")
    return path

def get_site_atom_indices(pdb_path, atoms, contact_cutoff=4.0):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("p", pdb_path)
    het_coords = []
    for res in struct.get_residues():
        if res.id[0].startswith("H_") and res.get_resname() not in ["HOH", "WAT", "DOD"]:
            for a in res.get_atoms():
                het_coords.append(a.coord)
                
    if not het_coords:
        site_res_nums = set()
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("SITE"):
                    parts = line.split()
                    for item in parts:
                        if item.isdigit():
                            site_res_nums.add(int(item))
        if site_res_nums:
            return set(i for i, a in enumerate(atoms) if a.get_parent().id[1] in site_res_nums)
            
        chains = list(struct.get_chains())
        if len(chains) > 1:
            for c in chains[1:]:
                for a in c.get_atoms():
                    het_coords.append(a.coord)

    if het_coords:
        tree = cKDTree(het_coords)
        return set(i for i, a in enumerate(atoms) if tree.query(a.coord)[0] <= contact_cutoff)
    return set()

def analyze_dataset_group(pdb_list):
    site_safds, non_site_safds = [], []
    for pdb in pdb_list:
        path = ensure_downloaded(pdb)
        try:
            atoms, surface_indices, safd_values = compute_safd(path)
            site_indices = get_site_atom_indices(path, atoms)
            for idx in surface_indices:
                f_val = safd_values.get(idx, 2.0)
                if idx in site_indices:
                    site_safds.append(f_val)
                else:
                    non_site_safds.append(f_val)
        except Exception as e:
            print(f"Error on {pdb}: {e}")
    site_m = np.mean(site_safds) if site_safds else 0.0
    non_m = np.mean(non_site_safds) if non_site_safds else 0.0
    _, p_val = ttest_ind(site_safds, non_site_safds, equal_var=False) if site_safds else (0, 1.0)
    return site_m, non_m, p_val

def main():
    print("=" * 60)
    print("RUNNING FULL SAFD PIPELINE (OBJECTIVES 1 - 5)")
    print("=" * 60)

    # 1. Validation Suites
    print("\n[1/3] Running Validation Suites (Bound, Apo, Non-Cleft)...")
    b_site, b_non, b_p = analyze_dataset_group(DATASET_CONFIG["bound"])
    a_site, a_non, a_p = analyze_dataset_group(DATASET_CONFIG["apo"])
    nc_site, nc_non, nc_p = analyze_dataset_group(DATASET_CONFIG["non_cleft"])
    
    print(f"  Bound:      Site = {b_site:.2f} | Non-Site = {b_non:.2f} | p = {b_p:.3e}")
    print(f"  Apo:        Site = {a_site:.2f} | Non-Site = {a_non:.2f} | p = {a_p:.3e}")
    print(f"  Non-Cleft:  Site = {nc_site:.2f} | Non-Site = {nc_non:.2f} | p = {nc_p:.3e}")

    labels = ['Bound Sites', 'Apo (No Ligand)', 'Non-Cleft Sites']
    x = np.arange(len(labels))
    w = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - w/2, [b_site, a_site, nc_site], w, label='Functional Site', color='#006699')
    ax.bar(x + w/2, [b_non, a_non, nc_non], w, label='Non-Site', color='#99ccdd')
    ax.set_ylabel('Mean Roughness (SAFD)')
    ax.set_title('Validation: Site vs Non-Site Roughness')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 5)
    ax.legend()
    plt.tight_layout()
    plt.savefig("validation_results.png")

    # 2. Site-Size Analysis
    print("\n[2/3] Running Site-Size Analysis (<600 vs >600 A^2)...")
    sm_site, sm_non, sm_p = analyze_dataset_group(DATASET_CONFIG["small_sites"])
    lg_site, lg_non, lg_p = analyze_dataset_group(DATASET_CONFIG["large_sites"])
    print(f"  Small Sites (<600 A^2): Site = {sm_site:.2f} | Non-Site = {sm_non:.2f} | p = {sm_p:.3e}")
    print(f"  Large Sites (>600 A^2): Site = {lg_site:.2f} | Non-Site = {lg_non:.2f} | p = {lg_p:.3e}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar([0 - w/2, 1 - w/2], [sm_site, lg_site], w, label='Functional Site', color='#e67e22')
    ax.bar([0 + w/2, 1 + w/2], [sm_non, lg_non], w, label='Non-Site', color='#f5cba7')
    ax.set_ylabel('Mean Roughness (SAFD)')
    ax.set_title('Roughness vs. Site Size')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Small (<600 Å²)', 'Large (>600 Å²)'])
    ax.set_ylim(0, 5)
    ax.legend()
    plt.tight_layout()
    plt.savefig("site_size_result.png")

    # 3. Thermal Extension
    print("\n[3/3] Running Thermal Adaptation Extension...")
    thermal_res = {}
    for grp, pdbs in DATASET_CONFIG["thermal"].items():
        s_m, b_m, _ = analyze_dataset_group(pdbs)
        thermal_res[grp] = {"site": s_m, "bulk": b_m}
        print(f"  {grp:<26}: Site = {s_m:.2f} | Bulk = {b_m:.2f}")

    classes = list(thermal_res.keys())
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(classes, [thermal_res[c]["site"] for c in classes], marker='o', lw=2.5, color='#e67e22', label='Functional Site')
    ax.plot(classes, [thermal_res[c]["bulk"] for c in classes], marker='s', lw=2.5, color='#2980b9', label='Bulk Surface')
    ax.set_ylabel('Mean Roughness (SAFD)')
    ax.set_title('Thermal Adaptation: Site vs Bulk Surface Roughness')
    ax.set_ylim(2.0, 5.0)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend()
    plt.tight_layout()
    plt.savefig("thermal_adaptation_result.png")

    print("\n" + "=" * 60)
    print("ALL RUNS COMPLETE")
    print("Generated files: validation_results.png, site_size_result.png, thermal_adaptation_result.png")
    print("=" * 60)

if __name__ == "__main__":
    main()