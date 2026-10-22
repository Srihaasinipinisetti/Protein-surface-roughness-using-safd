import numpy as np
from scipy.spatial import cKDTree

# Standard Bondi van der Waals radii (Angstroms)
VDW_RADII = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 
    'F': 1.47, 'P': 1.80, 'S': 1.80, 'CL': 1.75
}

def generate_golden_spiral_sphere(n_points=96):
    """Generates roughly equidistant test points on a unit sphere."""
    indices = np.arange(0, n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + 5**0.5) * indices
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack((x, y, z))

def compute_atomic_contact_areas(coords, elements, probe_radius=1.4, n_test_points=96, use_contact_area=True):
    """
    Computes atomic contact area (Richards, 1977) or accessible surface area (Lee & Richards, 1971).
    Per Pettit & Bowie (1999) Eq. 2, SAFD requires contact area on the atomic van der Waals
    sphere surface (r_vdw^2) rather than accessible area ((r_vdw + Rp)^2).
    """
    n_atoms = len(coords)
    radii = np.array([VDW_RADII.get(e.upper(), 1.70) for e in elements])
    effective_radii = radii + probe_radius
    
    unit_sphere_pts = generate_golden_spiral_sphere(n_test_points)
    tree = cKDTree(coords)
    areas = np.zeros(n_atoms)
    
    for i in range(n_atoms):
        r_eff = effective_radii[i]
        r_base = radii[i] if use_contact_area else r_eff
        # Generate test points on the expanded sphere surface to test solvent probe clearance
        test_points = coords[i] + unit_sphere_pts * r_eff
        
        # Find neighboring atoms that could occlude points
        neighbor_indices = tree.query_ball_point(coords[i], r=r_eff + np.max(effective_radii))
        neighbor_indices = [idx for idx in neighbor_indices if idx != i]
        
        if not neighbor_indices:
            areas[i] = 4.0 * np.pi * (r_base**2)
            continue
            
        neigh_coords = coords[neighbor_indices]
        neigh_eff_radii = effective_radii[neighbor_indices]
        
        # Check occlusion of test points by neighboring spheres
        # Distance squared from test points to neighbor centers
        diff = test_points[:, np.newaxis, :] - neigh_coords[np.newaxis, :, :]
        dist_sq = np.sum(diff**2, axis=-1)
        
        # A point is buried if its distance to a neighbor < neighbor's effective radius
        is_occluded = np.any(dist_sq < (neigh_eff_radii**2), axis=1)
        accessible_points = np.sum(~is_occluded)
        
        areas[i] = (4.0 * np.pi * (r_base**2)) * (accessible_points / n_test_points)
        
    return areas