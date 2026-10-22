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
    """
    Computes group roughness per Pettit & Bowie (1999) Table 1:
    Takes the median SAFD of exposed site and non-site atoms for each protein,
    then averages across all proteins in the set, quoting SEM and Student's t-test.
    """
    protein_site_meds = []
    protein_non_meds = []
    
    for pdb in pdb_list:
        path = ensure_downloaded(pdb)
        try:
            atoms, surface_indices, safd_values = compute_safd(path)
            site_indices = get_site_atom_indices(path, atoms)
            
            p_site = [safd_values.get(idx, 2.0) for idx in surface_indices if idx in site_indices]
            p_non = [safd_values.get(idx, 2.0) for idx in surface_indices if idx not in site_indices]
            
            if p_site and p_non:
                protein_site_meds.append(np.median(p_site))
                protein_non_meds.append(np.median(p_non))
        except Exception as e:
            print(f"Error on {pdb}: {e}")
            
    if protein_site_meds:
        n = len(protein_site_meds)
        site_m = float(np.mean(protein_site_meds))
        site_sem = float(np.std(protein_site_meds, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        non_m = float(np.mean(protein_non_meds))
        non_sem = float(np.std(protein_non_meds, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        _, p_val = ttest_ind(protein_site_meds, protein_non_meds, equal_var=False) if n > 1 else (0, 1.0)
    else:
        site_m, site_sem, non_m, non_sem, p_val = 0.0, 0.0, 0.0, 0.0, 1.0
        
    return site_m, site_sem, non_m, non_sem, p_val

def main():
    print("=" * 60)
    print("RUNNING FULL SAFD PIPELINE (PETTIT & BOWIE 1999 REPLICATION)")
    print("=" * 60)

    # 1. Validation Suites
    print("\n[1/3] Running Validation Suites (Bound, Apo, Non-Cleft)...")
    b_site, b_s_sem, b_non, b_ns_sem, b_p = analyze_dataset_group(DATASET_CONFIG["bound"])
    a_site, a_s_sem, a_non, a_ns_sem, a_p = analyze_dataset_group(DATASET_CONFIG["apo"])
    nc_site, nc_s_sem, nc_non, nc_ns_sem, nc_p = analyze_dataset_group(DATASET_CONFIG["non_cleft"])
    
    print(f"  Bound:      Site = {b_site:.2f} ± {b_s_sem:.2f} | Non-Site = {b_non:.2f} ± {b_ns_sem:.2f} | p = {b_p:.3e}")
    print(f"  Apo:        Site = {a_site:.2f} ± {a_s_sem:.2f} | Non-Site = {a_non:.2f} ± {a_ns_sem:.2f} | p = {a_p:.3e}")
    print(f"  Non-Cleft:  Site = {nc_site:.2f} ± {nc_s_sem:.2f} | Non-Site = {nc_non:.2f} ± {nc_ns_sem:.2f} | p = {nc_p:.3e}")

    labels = ['Bound Sites', 'Apo (No Ligand)', 'Non-Cleft Sites']
    x = np.arange(len(labels))
    w = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(x - w/2, [b_site, a_site, nc_site], w, yerr=[b_s_sem, a_s_sem, nc_s_sem], capsize=5, label='Functional Site', color='#006699', edgecolor='black', alpha=0.9)
    ax.bar(x + w/2, [b_non, a_non, nc_non], w, yerr=[b_ns_sem, a_ns_sem, nc_ns_sem], capsize=5, label='Non-Site', color='#99ccdd', edgecolor='black', alpha=0.9)
    ax.set_ylabel('Mean Roughness (SAFD)')
    ax.set_title('Validation: Site vs Non-Site Roughness (Pettit & Bowie Table 1)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 6.5)
    ax.axhline(2.0, color='gray', linestyle='--', linewidth=0.8, label='Smooth baseline (f=2)')
    ax.legend()
    plt.tight_layout()
    plt.savefig("validation_results.png", dpi=150)

    # 2. Site-Size Analysis
    print("\n[2/3] Running Site-Size Analysis (<600 vs >600 A^2)...")
    sm_site, sm_s_sem, sm_non, sm_ns_sem, sm_p = analyze_dataset_group(DATASET_CONFIG["small_sites"])
    lg_site, lg_s_sem, lg_non, lg_ns_sem, lg_p = analyze_dataset_group(DATASET_CONFIG["large_sites"])
    print(f"  Small Sites (<600 A^2): Site = {sm_site:.2f} ± {sm_s_sem:.2f} | Non-Site = {sm_non:.2f} ± {sm_ns_sem:.2f} | p = {sm_p:.3e}")
    print(f"  Large Sites (>600 A^2): Site = {lg_site:.2f} ± {lg_s_sem:.2f} | Non-Site = {lg_non:.2f} ± {lg_ns_sem:.2f} | p = {lg_p:.3e}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.bar([0 - w/2, 1 - w/2], [sm_site, lg_site], w, yerr=[sm_s_sem, lg_s_sem], capsize=5, label='Functional Site', color='#e67e22', edgecolor='black', alpha=0.9)
    ax.bar([0 + w/2, 1 + w/2], [sm_non, lg_non], w, yerr=[sm_ns_sem, lg_ns_sem], capsize=5, label='Non-Site', color='#f5cba7', edgecolor='black', alpha=0.9)
    ax.set_ylabel('Mean Roughness (SAFD)')
    ax.set_title('Roughness vs. Site Size (Pettit & Bowie Table 1)')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Small (<600 Å²)', 'Large (>600 Å²)'])
    ax.set_ylim(0, 6.5)
    ax.axhline(2.0, color='gray', linestyle='--', linewidth=0.8, label='Smooth baseline (f=2)')
    ax.legend()
    plt.tight_layout()
    plt.savefig("site_size_result.png", dpi=150)

    # 3. Thermal Extension (Novel Project Contribution)
    print("\n[3/3] Running Thermal Adaptation (Project Extension beyond Paper)...")
    thermal_res = {}
    for grp, pdbs in DATASET_CONFIG["thermal"].items():
        s_m, _, b_m, _, _ = analyze_dataset_group(pdbs)
        thermal_res[grp] = {"site": s_m, "bulk": b_m}
        print(f"  {grp:<26}: Site = {s_m:.2f} | Bulk = {b_m:.2f}")

    classes = list(thermal_res.keys())
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(classes, [thermal_res[c]["site"] for c in classes], marker='o', lw=2.5, color='#e67e22', label='Functional Site')
    ax.plot(classes, [thermal_res[c]["bulk"] for c in classes], marker='s', lw=2.5, color='#2980b9', label='Bulk Surface')
    ax.set_ylabel('Mean Roughness (SAFD)')
    ax.set_title('Thermal Adaptation: Site vs Bulk Surface Roughness (Extension)')
    ax.set_ylim(2.0, 5.5)
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend()
    plt.tight_layout()
    plt.savefig("thermal_adaptation_result.png", dpi=150)

    print("\n" + "=" * 60)
    print("ALL RUNS COMPLETE")
    print("Generated files: validation_results.png, site_size_result.png, thermal_adaptation_result.png")
    print("=" * 60)

if __name__ == "__main__":
    main()