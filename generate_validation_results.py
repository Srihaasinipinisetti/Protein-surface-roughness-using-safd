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
    site_safds = []
    non_site_safds = []
    
    for pdb in pdb_list:
        path = f"pdbs/{pdb}.pdb"
        try:
            atoms, surface_indices, safd_values = compute_safd(path)
            site_indices = get_site_atom_indices(path, atoms)
            
            for idx in surface_indices:
                val = safd_values.get(idx, 2.0)
                if idx in site_indices:
                    site_safds.append(val)
                else:
                    non_site_safds.append(val)
        except Exception as e:
            print(f"Skipping {pdb}: {e}")
            
    site_mean = np.mean(site_safds) if site_safds else 0
    non_site_mean = np.mean(non_site_safds) if non_site_safds else 0
    stat, p_val = ttest_ind(site_safds, non_site_safds, equal_var=False) if site_safds else (0, 1.0)
    return site_mean, non_site_mean, p_val

print("--- Calculating Slide 14: Bound Structures ---")
bound_site, bound_non, bound_p = analyze_group(["3DFR", "7HVP", "1RBP", "8ACN", "2GBP"])
print(f"Bound - Site: {bound_site:.2f}, Non-Site: {bound_non:.2f}, p-value: {bound_p:.4e}")

print("--- Calculating Slide 15: Apo Structures ---")
apo_site, apo_non, apo_p = analyze_group(["2ACT", "4BLM", "3EBX", "2LIV"])
print(f"Apo - Site: {apo_site:.2f}, Non-Site: {apo_non:.2f}, p-value: {apo_p:.4e}")

print("--- Calculating Slide 15: Non-Cleft Structures ---")
nc_site, nc_non, nc_p = analyze_group(["1FXD", "1FXA", "5PAL", "1LTE"])
print(f"Non-Cleft - Site: {nc_site:.2f}, Non-Site: {nc_non:.2f}, p-value: {nc_p:.4e}")

# Generate Plots for Slides 14 & 15
labels = ['Bound Sites', 'Apo (No Ligand)', 'Non-Cleft Sites']
site_means = [bound_site, apo_site, nc_site]
non_site_means = [bound_non, apo_non, nc_non]

x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(8, 5))
ax.bar(x - width/2, site_means, width, label='Functional Site', color='#006699')
ax.bar(x + width/2, non_site_means, width, label='Non-Site', color='#99ccdd')

ax.set_ylabel('Mean SAFD (Roughness)')
ax.set_title('Validation: Site vs Non-Site Roughness (Slides 14 & 15)')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylim(0, 5)
ax.legend()
plt.tight_layout()
plt.savefig("validation_results.png")
print("Saved plot to validation_results.png")