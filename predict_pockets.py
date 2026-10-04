import numpy as np
from sklearn.cluster import DBSCAN
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree
from safd import compute_safd

def extract_ligand_coords(pdb_path):
    """Extracts coordinates of biologically relevant bound heteroatoms (non-water)."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    ligand_coords = []
    
    for res in structure.get_residues():
        if res.id[0].startswith("H_") and res.get_resname() not in ["HOH", "WAT", "DOD"]:
            for atom in res.get_atoms():
                ligand_coords.append(atom.coord)
                
    return np.array(ligand_coords)

def cluster_high_safd_atoms(atoms, safd_values, f_threshold=3.4, eps=4.5, min_samples=3):
    """Slide 12: DBSCAN clustering of high-roughness atoms."""
    candidate_indices = [idx for idx, f in safd_values.items() if f >= f_threshold]
    if not candidate_indices:
        return []
    
    candidate_coords = np.array([atoms[idx].coord for idx in candidate_indices])
    
    # Run DBSCAN
    db = DBSCAN(eps=eps, min_samples=min_samples).fit(candidate_coords)
    labels = db.labels_

    clusters = []
    unique_labels = set(labels) - {-1}  # -1 is noise
    for label in unique_labels:
        cluster_atom_indices = [candidate_indices[i] for i in range(len(candidate_indices)) if labels[i] == label]
        clusters.append(cluster_atom_indices)

    # Rank clusters by total SAFD score
    clusters.sort(key=lambda cl: sum(safd_values[idx] for idx in cl), reverse=True)
    return clusters

def evaluate_prediction(clusters, atoms, ligand_coords, contact_cutoff=4.0):
    """Checks if predicted clusters overlap true site atoms (within 4 A of ligand)."""
    if len(ligand_coords) == 0:
        print("No bound ligand found in PDB file to benchmark against.")
        return

    ligand_tree = cKDTree(ligand_coords)
    
    print("\n" + "=" * 55)
    print("BINDING SITE PREDICTION RESULTS")
    print("=" * 55)
    
    for rank, cluster in enumerate(clusters[:3], start=1):
        cluster_coords = np.array([atoms[idx].coord for idx in cluster])
        dists, _ = ligand_tree.query(cluster_coords)
        hits = np.sum(dists <= contact_cutoff)
        overlap_pct = (hits / len(cluster)) * 100
        
        residues = set(f"{atoms[idx].get_parent().get_resname()} {atoms[idx].get_parent().id[1]}" for idx in cluster)
        status = "SUCCESS (Hits Ligand)" if hits > 0 else "FALSE POSITIVE"
        print(f"Cluster #{rank}: {len(cluster)} atoms | {status}")
        print(f"  Overlap with site: {hits}/{len(cluster)} atoms ({overlap_pct:.1f}%)")
        print(f"  Residues: {', '.join(sorted(residues)[:8])}...")
        print("-" * 55)

if __name__ == "__main__":
    pdb_file = "pdbs/3DFR.pdb"
    atoms, surface_indices, safd_values = compute_safd(pdb_file)
    clusters = cluster_high_safd_atoms(atoms, safd_values)
    ligand_coords = extract_ligand_coords(pdb_file)
    evaluate_prediction(clusters, atoms, ligand_coords)