import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from Bio.PDB import PDBParser
from safd import compute_safd
from scipy.spatial import cKDTree

SMALL_SITES = ["3DFR", "1RBP", "8ACN", "2GBP"]  # < 600 A^2
LARGE_SITES = ["1ACB", "7HVP"]                  # > 600 A^2 (complexes/dimers)

def evaluate_site_size_roughness(pdb_list):
    site_vals, non_site_vals = [], []
    for pdb in pdb_list:
        path = f"pdbs/{pdb}.pdb"
        try:
            atoms, surface_indices, safd_values = compute_safd(path)
            parser = PDBParser(QUIET=True)
            struct = parser.get_structure("p", path)
            partner_coords = []
            for res in struct.get_residues():
                if res.id[0].startswith("H_") and res.get_resname() not in ["HOH", "WAT", "DOD"]:
                    for a in res.get_atoms():
                        partner_coords.append(a.coord)
            if not partner_coords:
                chains = list(struct.get_chains())
                if len(chains) > 1:
                    for c in chains[1:]:
                        for a in c.get_atoms():
                            partner_coords.append(a.coord)
                            
            if partner_coords:
                tree = cKDTree(partner_coords)
                for idx in surface_indices:
                    val = safd_values.get(idx, 2.0)
                    dist, _ = tree.query(atoms[idx].coord)
                    if dist <= 4.0:
                        site_vals.append(val)
                    else:
                        non_site_vals.append(val)
        except Exception as e:
            print(f"Skipping {pdb}: {e}")
            
    mean_site = np.mean(site_vals) if site_vals else 0
    mean_non = np.mean(non_site_vals) if non_site_vals else 0
    _, p_val = ttest_ind(site_vals, non_site_vals, equal_var=False) if site_vals else (0, 1.0)
    return mean_site, mean_non, p_val

print("Evaluating Small Sites (< 600 A^2)...")
sm_site, sm_non, sm_p = evaluate_site_size_roughness(SMALL_SITES)
print(f"Small Sites: Site = {sm_site:.2f}, Non-Site = {sm_non:.2f}, p = {sm_p:.4e}")

print("\nEvaluating Large Sites (> 600 A^2)...")
lg_site, lg_non, lg_p = evaluate_site_size_roughness(LARGE_SITES)
print(f"Large Sites: Site = {lg_site:.2f}, Non-Site = {lg_non:.2f}, p = {lg_p:.4e}")

labels = ['Small Sites (<600 A^2)', 'Large Sites (>600 A^2)']
x = np.arange(len(labels))
width = 0.35

fig, ax = plt.subplots(figsize=(7, 5))
ax.bar(x - width/2, [sm_site, lg_site], width, label='Functional Site', color='#e67e22')
ax.bar(x + width/2, [sm_non, lg_non], width, label='Non-Site', color='#f5cba7')
ax.set_ylabel('Mean Roughness (SAFD)')
ax.set_title('Roughness vs. Site Size (Slide 16)')
ax.set_xticks(x)
ax.set_xticklabels(labels)
ax.set_ylim(0, 5)
ax.legend()
plt.tight_layout()
plt.savefig("site_size_result.png")
print("\nSaved plot to site_size_result.png")