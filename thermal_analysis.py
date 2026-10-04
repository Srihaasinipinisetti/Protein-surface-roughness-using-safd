import numpy as np
import matplotlib.pyplot as plt
from safd import compute_safd
from generate_validation_results import get_site_atom_indices

THERMAL_SETS = {
    "Mesophile (20-45 C)": ["1LZ1", "3DFR"],
    "Thermophile (45-80 C)": ["1V43", "1A5Z"],
    "Hyperthermophile (>80 C)": ["1B44", "1V7S"]
}

results = {}

for group, pdbs in THERMAL_SETS.items():
    site_vals = []
    bulk_vals = []
    for pdb in pdbs:
        path = f"pdbs/{pdb}.pdb"
        try:
            atoms, surface_indices, safd_values = compute_safd(path)
            site_indices = get_site_atom_indices(path, atoms)
            for idx in surface_indices:
                f_val = safd_values.get(idx, 2.0)
                if idx in site_indices:
                    site_vals.append(f_val)
                else:
                    bulk_vals.append(f_val)
        except Exception as e:
            print(f"Skipping {pdb}: {e}")
            
    results[group] = {
        "site": np.mean(site_vals) if site_vals else np.nan,
        "bulk": np.mean(bulk_vals) if bulk_vals else np.nan
    }

print("\nTHERMAL ADAPTATION RESULTS (Slide 17):")
print(f"{'Class':<26} {'Functional Site SAFD':<22} {'Bulk Surface SAFD'}")
print("-" * 65)
for grp, vals in results.items():
    print(f"{grp:<26} {vals['site']:<22.2f} {vals['bulk']:.2f}")

# Plotting Slide 17 trendline
classes = list(results.keys())
site_trends = [results[c]["site"] for c in classes]
bulk_trends = [results[c]["bulk"] for c in classes]

plt.figure(figsize=(8, 5))
plt.plot(classes, site_trends, marker='o', linewidth=2.5, color='#e67e22', label='Functional Site (Expected: Stays High)')
plt.plot(classes, bulk_trends, marker='s', linewidth=2.5, color='#2980b9', label='Bulk Surface (Expected: Smooths)')
plt.ylabel('Mean Roughness (SAFD)')
plt.title('Thermal Adaptation: Site vs Bulk Surface Roughness (Slide 17)')
plt.ylim(2.0, 5.0)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend()
plt.tight_layout()
plt.savefig("thermal_adaptation_result.png")
print("\nSaved plot to thermal_adaptation_result.png")