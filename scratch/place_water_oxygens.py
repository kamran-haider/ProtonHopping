"""
Script to place water oxygen atoms on a grid.
"""

import sys
import numpy as np
import mdtraj as md
from utils import NeighborSearch, write_watpdb_from_coords, initialize_grid



def main():
    input_file = "../data/gramicidin/prep_structures/prot_step2_out.pdb"
    lig = md.load_pdb(input_file, no_boxchk=True)
    com = np.zeros((lig.n_frames, 3))
    masses = np.ones(lig.n_atoms)
    masses /= masses.sum()
    com[0, :] = lig.xyz[0, :].astype('float64').T.dot(masses)
    grid_center = com[0, :] * 10.0
    grid_res = [0.5, 0.5, 0.5]
    grid_dims = [16.0, 16.0, 44.0]
    g = initialize_grid(grid_center, grid_res, grid_dims)
    prot_coords = lig.xyz[0, :, :] * 10.0
    valid_pts = []
    search_space = NeighborSearch(prot_coords, 2.4)
    n_wats = 0
    for pt_i in range(g.shape[0]):
        pt_coords = g[pt_i, 1:4]
        nn = search_space.query_nbrs_single_point(pt_coords)
        if len(nn) == 0:
            #print pt_i, pt_coords
            valid_pts.append(pt_coords)
            n_wats += 1
    write_watpdb_from_coords("../data/gramicidin/water_grid_placement/init_o_placement", valid_pts)
    #write_watpdb_from_coords("../data/gramicidin/water_grid_placement/o_grid", g[:, 1:4])

if __name__ == "__main__":
    sys.exit(main())
