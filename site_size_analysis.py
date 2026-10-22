import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from Bio.PDB import PDBParser
from safd import compute_safd
from scipy.spatial import cKDTree

SMALL_SITES = ["3DFR", "1RBP", "8ACN", "2GBP", "2LHB", "4TGL"]  # < 600 A^2 (small molecules)
LARGE_SITES = ["1ACB", "7HVP"]                                    # > 600 A^2 (complexes/dimers)

def evaluate_site_size_roughness(pdb_list):
    """
    Evaluates site vs non-site roughness and Pettit & Bowie (1999) Figure 3 Relative Roughness:
    Relative Roughness = median(f - 2)_site / median(f - 2)_all
    A ratio of 1.0 means the site is just as rough as typical bulk protein surface.
    """
    protein_site_meds = []
    protein_non_meds = []
    relative_roughness_list = []
    
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
                site_vals, non_site_vals, all_f_minus_2 = [], [], []
                
                for idx in surface_indices:
                    val = safd_values.get(idx, 2.0)
                    all_f_minus_2.append(max(val - 2.0, 0.001))
                    dist, _ = tree.query(atoms[idx].coord)
                    if dist <= 4.0:
                        site_vals.append(val)
                    else:
                        non_site_vals.append(val)
                        
                if site_vals and non_site_vals:
                    s_med = np.median(site_vals)
                    ns_med = np.median(non_site_vals)
                    protein_site_meds.append(s_med)
                    protein_non_meds.append(ns_med)
                    
                    # Pettit & Bowie Figure 3 Relative Roughness
                    site_f_minus_2 = [max(v - 2.0, 0.001) for v in site_vals]
                    rel_rough = np.median(site_f_minus_2) / np.median(all_f_minus_2)
                    relative_roughness_list.append(rel_rough)
                    print(f"  {pdb}: Site Med={s_med:.2f}, Non-Site Med={ns_med:.2f}, Rel Roughness={rel_rough:.2f}")
        except Exception as e:
            print(f"Skipping {pdb}: {e}")
            
    n = len(protein_site_meds)
    site_m = float(np.mean(protein_site_meds)) if n > 0 else 0.0
    site_sem = float(np.std(protein_site_meds, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    non_m = float(np.mean(protein_non_meds)) if n > 0 else 0.0
    non_sem = float(np.std(protein_non_meds, ddof=1) / np.sqrt(n)) if n > 1 else 0.0
    rel_m = float(np.mean(relative_roughness_list)) if relative_roughness_list else 1.0
    rel_sem = float(np.std(relative_roughness_list, ddof=1) / np.sqrt(len(relative_roughness_list))) if len(relative_roughness_list) > 1 else 0.0
    _, p_val = ttest_ind(protein_site_meds, protein_non_meds, equal_var=False) if n > 1 else (0, 1.0)
    
    return site_m, site_sem, non_m, non_sem, rel_m, rel_sem, p_val

print("=== Evaluating Small Sites (< 600 A^2) ===")
sm_site, sm_s_sem, sm_non, sm_ns_sem, sm_rel, sm_rel_sem, sm_p = evaluate_site_size_roughness(SMALL_SITES)
print(f"Small Sites: Site = {sm_site:.2f} ± {sm_s_sem:.2f}, Non-Site = {sm_non:.2f} ± {sm_ns_sem:.2f}, Rel Roughness = {sm_rel:.2f}, p = {sm_p:.4e}\n")

print("=== Evaluating Large Sites (> 600 A^2) ===")
lg_site, lg_s_sem, lg_non, lg_ns_sem, lg_rel, lg_rel_sem, lg_p = evaluate_site_size_roughness(LARGE_SITES)
print(f"Large Sites: Site = {lg_site:.2f} ± {lg_s_sem:.2f}, Non-Site = {lg_non:.2f} ± {lg_ns_sem:.2f}, Rel Roughness = {lg_rel:.2f}, p = {lg_p:.4e}\n")

# 2-Panel Plot: (A) Absolute Roughness with SEM, (B) Relative Roughness (Pettit & Bowie Fig 3)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

# Panel A
labels = ['Small (<600 Å²)', 'Large (>600 Å²)']
x = np.arange(len(labels))
width = 0.35

ax1.bar(x - width/2, [sm_site, lg_site], width, yerr=[sm_s_sem, lg_s_sem], capsize=5, label='Functional Site', color='#e67e22', edgecolor='black', alpha=0.9)
ax1.bar(x + width/2, [sm_non, lg_non], width, yerr=[sm_ns_sem, lg_ns_sem], capsize=5, label='Non-Site', color='#f5cba7', edgecolor='black', alpha=0.9)
ax1.set_ylabel('Mean Roughness (SAFD)')
ax1.set_title('(A) Roughness vs. Site Size (Table 1)')
ax1.set_xticks(x)
ax1.set_xticklabels(labels)
ax1.set_ylim(0, 7.0)
ax1.axhline(2.0, color='gray', linestyle='--', linewidth=0.8, label='Smooth (f=2)')
ax1.legend()

# Panel B: Relative Roughness
ax2.bar(labels, [sm_rel, lg_rel], width=0.5, yerr=[sm_rel_sem, lg_rel_sem], capsize=5, color=['#d35400', '#2980b9'], edgecolor='black', alpha=0.85)
ax2.axhline(1.0, color='red', linestyle='--', linewidth=1.5, label='Bulk Surface Baseline (Ratio=1.0)')
ax2.set_ylabel('Relative Roughness: median(f-2)site / median(f-2)all')
ax2.set_title('(B) Relative Roughness (Pettit & Bowie Fig 3)')
ax2.set_ylim(0, max(sm_rel + 0.8, 3.0))
ax2.legend()

plt.tight_layout()
plt.savefig("site_size_result.png", dpi=150)
print("Saved 2-panel figure to site_size_result.png")