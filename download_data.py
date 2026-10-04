import os
import requests

BENCHMARK_PDBS = ["3DFR", "7HVP", "1RBP", "1ACB"]

os.makedirs("pdbs", exist_ok=True)

for pdb_id in BENCHMARK_PDBS:
    file_path = os.path.join("pdbs", f"{pdb_id}.pdb")
    if not os.path.exists(file_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        if response.status_code == 200:
            with open(file_path, "w") as f:
                f.write(response.text)
            print(f"Downloaded: {pdb_id}")
        else:
            print(f"Failed to fetch {pdb_id}: Status {response.status_code}")
    else:
        print(f"Already exists: {pdb_id}")