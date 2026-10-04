import numpy as np
from scipy.spatial import cKDTree
from Bio.PDB import PDBParser
from sasa import compute_atomic_contact_areas

def extract_protein_atoms(pdb_path):
    """Extracts only standard amino acid backbone/sidechain atoms."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    
    coords = []
    elements = []
    atom_records = []
    
    for atom in structure.get_atoms():
        residue = atom.get_parent()
        # Keep only standard amino acid residues (hetero flag is ' ')
        if residue.id[0] == ' ':
            coords.append(atom.coord)
            elements.append(atom.element.strip() if atom.element else atom.name[0])
            atom_records.append(atom)
            
    return np.array(coords), elements, atom_records

def compute_safd(pdb_path, probe_radii=[1.2, 1.4, 1.6, 1.8, 2.1], cutoff=5.0):
    coords, elements, atoms = extract_protein_atoms(pdb_path)
    print(f"Extracted {len(coords)} protein atoms.")
    
    # 1. Multi-scale probe calculations
    areas_by_probe = []
    for r in probe_radii:
        print(f"Calculating SASA for probe radius Rp = {r} A...")
        areas = compute_atomic_contact_areas(coords, elements, probe_radius=r)
        areas_by_probe.append(areas)
        
    areas_by_probe = np.array(areas_by_probe) # Shape: (len(radii), n_atoms)
    
    # 2. Neighborhood smoothing and fitting
    tree = cKDTree(coords)
    log_r = np.log(probe_radii)
    safd_values = {}
    
    # Filter for surface atoms (accessible at Rp = 1.4 A)
    surface_mask = areas_by_probe[1] > 0.0
    surface_indices = np.where(surface_mask)[0]
    print(f"Identified {len(surface_indices)} surface atoms.")
    
    for idx in surface_indices:
        # 5.0 Angstrom neighborhood summation
        neighbors = tree.query_ball_point(coords[idx], r=cutoff)
        summed_area = areas_by_probe[:, neighbors].sum(axis=1)
        
        valid = summed_area > 0
        if np.sum(valid) >= 3:
            slope, _ = np.polyfit(log_r[valid], np.log(summed_area[valid]), deg=1)
            safd_values[idx] = 2.0 - slope
        else:
            safd_values[idx] = 2.0
            
    return atoms, surface_indices, safd_values

if __name__ == "__main__":
    atoms, surf_idx, safd = compute_safd("pdbs/3DFR.pdb")
    
    # Preview top 10 atoms with highest SAFD (potential clefts/pockets)
    sorted_atoms = sorted(safd.items(), key=lambda x: x[1], reverse=True)
    print("\nTop 10 High-SAFD Atoms (Cleft Candidates):")
    print(f"{'Residue':<12} {'Atom':<6} {'SAFD (f)':<10}")
    print("-" * 30)
    for idx, f_val in sorted_atoms[:10]:
        atom = atoms[idx]
        res = atom.get_parent()
        res_name = f"{res.get_resname()} {res.id[1]}"
        print(f"{res_name:<12} {atom.get_name():<6} {f_val:.3f}")