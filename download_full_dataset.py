import os
import requests

# Datasets specified across Slides 7, 8, 14, 15, and 17
DATASETS = {
    # Slide 14: Bound benchmark subset
    "bound": ["3DFR", "7HVP", "1RBP", "8ACN", "2GBP", "2LHB", "4TGL", "1MSB"],
    # Slide 15: Apo / No bound partner
    "apo": ["2ACT", "4BLM", "3EBX", "2LIV", "1LZ1", "1MVP"],
    # Slide 15: Shallow / Non-cleft structures
    "non_cleft": ["1FXD", "1FXA", "5PAL", "1LTE", "1HNE", "2PK4", "1SHA", "1CBV"],
    # Slide 17: Thermal adaptation extension
    "thermal_meso": ["1LZ1", "3DFR", "1ACB"],
    "thermal_thermo": ["1V43", "1A5Z"],
    "thermal_hyper": ["1B44", "1V7S"]
}

os.makedirs("pdbs", exist_ok=True)

for group, pdbs in DATASETS.items():
    for pdb_id in pdbs:
        path = f"pdbs/{pdb_id}.pdb"
        if not os.path.exists(path):
            res = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb")
            if res.status_code == 200:
                with open(path, "w") as f:
                    f.write(res.text)
                print(f"[{group}] Downloaded {pdb_id}")
            else:
                print(f"[{group}] Failed {pdb_id}")
        else:
            print(f"[{group}] Already present {pdb_id}")