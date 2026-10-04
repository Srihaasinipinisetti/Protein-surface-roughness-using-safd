import os
import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser, PDBIO
from safd import compute_safd
from predict_pockets import cluster_high_safd_atoms

def extract_binding_partner_coords(pdb_path, receptor_chain_ids=None):
    """
    Extracts true site contacts. Supports both small molecule ligands (HETATM)
    and bound protein partner chains (e.g., 1ACB chain E vs chain I).
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    
    het_coords = []
    chains = list(structure.get_chains())
    
    # 1. Check for small-molecule HETATM records
    for res in structure.get_residues():
        if res.id[0].startswith("H_") and res.get_resname() not in ["HOH", "WAT", "DOD"]:
            for atom in res.get_atoms():
                het_coords.append(atom.coord)
                
    if len(het_coords) > 0:
        return np.array(het_coords)
        
    # 2. If no small molecule exists, treat the secondary chain as the ligand partner
    if len(chains) > 1:
        partner_coords = []
        # Use first chain as receptor, remaining chains as binding partner
        for chain in chains[1:]:
            for atom in chain.get_atoms():
                partner_coords.append(atom.coord)
        return np.array(partner_coords)
        
    return np.array([])

def save_safd_to_pdb(input_pdb, output_pdb, safd_values):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)
    
    atom_idx = 0
    for atom in structure.get_atoms():
        residue = atom.get_parent()
        if residue.id[0] == ' ':
            atom.set_bfactor(float(safd_values.get(atom_idx, 2.0)))
            atom_idx += 1
            
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def evaluate_paper_metrics(clusters, atoms, partner_coords, contact_cutoff=4.0):
    if len(partner_coords) == 0:
        return None

    partner_tree = cKDTree(partner_coords)
    total_surface_residues = set((a.get_parent().get_resname(), a.get_parent().id[1]) for a in atoms)
    n_total_surf = len(total_surface_residues)
    
    # Identify true functional site residues
    site_residues = set()
    for atom in atoms:
        dist, _ = partner_tree.query(atom.coord)
        if dist <= contact_cutoff:
            site_residues.add((atom.get_parent().get_resname(), atom.get_parent().id[1]))
            
    n_site_res = len(site_residues)
    if n_site_res == 0:
        return None

    # Paper rule: evaluate top two highest-scoring clusters combined
    predicted_residues = set()
    for cl in clusters[:2]:
        for idx in cl:
            predicted_residues.add((atoms[idx].get_parent().get_resname(), atoms[idx].get_parent().id[1]))

    n_pred_res = len(predicted_residues)
    overlap = predicted_residues.intersection(site_residues)
    n_overlap = len(overlap)
    
    sensitivity = (n_overlap / n_site_res) * 100 if n_site_res > 0 else 0
    specificity = (n_overlap / n_pred_res) * 100 if n_pred_res > 0 else 0
    
    # Expected random rates per Pettit & Bowie (1999) page 1380
    random_sens = (n_pred_res / n_total_surf) * 100 if n_total_surf > 0 else 1.0
    random_spec = (n_site_res / n_total_surf) * 100 if n_total_surf > 0 else 1.0
    
    # Criterion: at least 1 overlap, and sensitivity & specificity >= 2x random expectation
    success = (n_overlap > 0) and (sensitivity >= 2.0 * random_sens) and (specificity >= 2.0 * random_spec)

    return {
        "n_site_res": n_site_res,
        "n_pred_res": n_pred_res,
        "overlap": n_overlap,
        "sensitivity": sensitivity,
        "specificity": specificity,
        "random_sens": random_sens,
        "random_spec": random_spec,
        "success": success
    }

def run_benchmark(pdb_list):
    os.makedirs("results", exist_ok=True)
    
    print(f"\n{'PDB':<8} {'SiteRes':<9} {'PredRes':<9} {'Overlap':<9} {'Sens(%)':<9} {'Spec(%)':<9} {'RandomExp(Sens/Spec)':<22} {'Result'}")
    print("-" * 88)
    
    for pdb_id in pdb_list:
        pdb_path = f"pdbs/{pdb_id}.pdb"
        if not os.path.exists(pdb_path):
            continue
            
        atoms, _, safd_values = compute_safd(pdb_path)
        # Using f >= 3.6 captures the full local pocket pocket without over-clustering
        clusters = cluster_high_safd_atoms(atoms, safd_values, f_threshold=3.6, max_dist=5.0)
        partner_coords = extract_binding_partner_coords(pdb_path)
        
        save_safd_to_pdb(pdb_path, f"results/{pdb_id}_safd.pdb", safd_values)
        metrics = evaluate_paper_metrics(clusters, atoms, partner_coords)
        
        if metrics:
            status = "PASS" if metrics["success"] else "FAIL"
            rnd_str = f"{metrics['random_sens']:.1f}% / {metrics['random_spec']:.1f}%"
            print(f"{pdb_id:<8} {metrics['n_site_res']:<9} {metrics['n_pred_res']:<9} {metrics['overlap']:<9} {metrics['sensitivity']:<9.1f} {metrics['specificity']:<9.1f} {rnd_str:<22} {status}")
        else:
            print(f"{pdb_id:<8} Unable to locate bound partner/ligand.")

if __name__ == "__main__":
    run_benchmark(["3DFR", "7HVP", "1RBP", "1ACB"])