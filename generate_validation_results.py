import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from Bio.PDB import PDBParser
from safd import compute_safd
from predict_pockets import extract_ligand_coords
from scipy.spatial import cKDTree

def get_site_atom_indices(pdb_path, atoms, contact_cutoff=4.0):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("p", pdb_path)
    
    # 1. From HETATM ligands
    het_coords = []
    for res in struct.get_residues():
        if res.id[0].startswith("H_") and res.get_resname() not in ["HOH", "WAT", "DOD"]:
            for a in res.get_atoms():
                het_coords.append(a.coord)
                
    # 2. If no ligand, check PDB SITE records (common in Apo structures)
    if len(het_coords) == 0:
        site_res_nums = set()
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("SITE"):
                    parts = line.split()
                    for item in parts:
                        if item.isdigit():
                            site_res_nums.add(int(item))
        if site_res_nums:
            site_indices = [i for i, a in enumerate(atoms) if a.get_parent().id[1] in site_res_nums]
            return set(site_indices)

    if len(het_coords) > 0:
        tree = cKDTree(het_coords)
        site_indices = []
        for i, a in enumerate(atoms):
            d, _ = tree.query(a.coord)
            if d <= contact_cutoff:
                site_indices.append(i)
        return set(site_indices)
        
    return set()

def analyze_group(pdb_list):
    """
    Computes group roughness per Pettit & Bowie (1999) Table 1:
    Takes the median SAFD of exposed site and non-site atoms for each protein,
    then averages across all proteins in the set, quoting standard error of the mean (SEM).
    """
    protein_site_medians = []
    protein_non_medians = []
    
    for pdb in pdb_list:
        path = f"pdbs/{pdb}.pdb"
        try:
            atoms, surface_indices, safd_values = compute_safd(path)
            site_indices = get_site_atom_indices(path, atoms)
            
            p_site = [safd_values.get(idx, 2.0) for idx in surface_indices if idx in site_indices]
            p_non = [safd_values.get(idx, 2.0) for idx in surface_indices if idx not in site_indices]
            
            if p_site and p_non:
                s_med = np.median(p_site)
                ns_med = np.median(p_non)
                protein_site_medians.append(s_med)
                protein_non_medians.append(ns_med)
                print(f"  {pdb}: Site Median = {s_med:.2f} | Non-Site Median = {ns_med:.2f}")
        except Exception as e:
            print(f"Skipping {pdb}: {e}")
            
    if protein_site_medians:
        n = len(protein_site_medians)
        site_mean = float(np.mean(protein_site_medians))
        site_sem = float(np.std(protein_site_medians, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        non_mean = float(np.mean(protein_non_medians))
        non_sem = float(np.std(protein_non_medians, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        _, p_val = ttest_ind(protein_site_medians, protein_non_medians, equal_var=False)
    else:
        site_mean, site_sem, non_mean, non_sem, p_val = 0.0, 0.0, 0.0, 0.0, 1.0
        
    return site_mean, site_sem, non_mean, non_sem, p_val

print("--- Calculating Table 1: Bound Structures ---")
bound_site, bound_s_sem, bound_non, bound_ns_sem, bound_p = analyze_group(["3DFR", "7HVP", "1RBP", "8ACN", "2GBP"])
print(f"Bound -> Site: {bound_site:.2f} ± {bound_s_sem:.2f}, Non-Site: {bound_non:.2f} ± {bound_ns_sem:.2f}, p-value: {bound_p:.4e}\n")

print("--- Calculating Table 1: Apo Structures ---")
apo_site, apo_s_sem, apo_non, apo_ns_sem, apo_p = analyze_group(["2ACT", "4BLM", "3EBX", "2LIV"])
print(f"Apo -> Site: {apo_site:.2f} ± {apo_s_sem:.2f}, Non-Site: {apo_non:.2f} ± {apo_ns_sem:.2f}, p-value: {apo_p:.4e}\n")

print("--- Calculating Table 1: Non-Cleft Structures ---")
nc_site, nc_s_sem, nc_non, nc_ns_sem, nc_p = analyze_group(["1FXD", "1FXA", "5PAL", "1LTE"])
print(f"Non-Cleft -> Site: {nc_site:.2f} ± {nc_s_sem:.2f}, Non-Site: {nc_non:.2f} ± {nc_ns_sem:.2f}, p-value: {nc_p:.4e}\n")

# Generate Validation Plot with SEM error bars (Pettit & Bowie format)
labels = ['Bound Sites', 'Apo (No Ligand)', 'Non-Cleft Sites']
site_means = [bound_site, apo_site, nc_site]
site_sems = [bound_s_sem, apo_s_sem, nc_s_sem]
non_site_means = [bound_non, apo_non, nc_non]
non_site_sems = [bound_ns_sem, apo_ns_sem, nc_ns_sem]

x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(x - width/2, site_means, width, yerr=site_sems, capsize=5, label='Functional Site', color='#006699', edgecolor='black', alpha=0.9)
ax.bar(x + width/2, non_site_means, width, yerr=non_site_sems, capsize=5, label='Non-Site', color='#99ccdd', edgecolor='black', alpha=0.9)

ax.set_ylabel('Mean Roughness (SAFD)')
ax.set_title('Validation: Site vs Non-Site Roughness (Pettit & Bowie 1999 Table 1)')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylim(0, 6.5)
ax.axhline(2.0, color='gray', linestyle='--', linewidth=0.8, label='Smooth baseline (f=2)')
ax.legend()
plt.tight_layout()
plt.savefig("validation_results.png", dpi=150)
print("Saved publication-grade plot with error bars to validation_results.png")