"""
Script to place water oxygen atoms on a grid.
"""

import sys
import os
from shutil import copy2
import numpy as np
import mdtraj as md
from pymcce.utils import NeighborSearch, write_watpdb_from_coords, initialize_grid
#from pymcce.sampler import write_watpdb_from_coords, get_last_prot_at_index, generate_conformers


def main():
    curr_dir = os.getcwd()
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    water_dir = os.path.abspath("../data/gramicidin/water_placement")

    prot_structs = [os.path.join(input_dir, "%02d_step1_out.pdb" % i) for i in range(10)]
    #print(prot_structs)
    for index, st in enumerate(prot_structs):
        print("Processing %s: " % st)
        prefix = "%02d" %index
        lig = md.load_pdb(st, no_boxchk=True)
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
                valid_pts.append(pt_coords)
                n_wats += 1
        wat_o_pdb_file = os.path.join(water_dir, "%s_init_o_placement" % prefix)
        write_watpdb_from_coords(wat_o_pdb_file, valid_pts)
        wat_pdb_file = os.path.join(water_dir, "%s_init_wat_placement" % prefix)        
        """
        wat_o_pdb_file = os.path.join(water_dir, "%s_init_o_placement" % prefix)
        #write_watpdb_from_coords(wat_o_pdb_file, valid_pts)
        wat_pdb_file = os.path.join(water_dir, "%s_init_wat_placement" % prefix)
        """
        with open("tleap_type_waters.in", "w") as f2:
            f2.write("source leaprc.water.tip3p\n")
            f2.write("mol = loadpdb %s.pdb\n" % wat_o_pdb_file)
            f2.write("savepdb mol %s.pdb\n" % wat_pdb_file)
            f2.write("quit\n")
        run_command = "tleap -f tleap_type_waters.in"
        os.system(run_command)


if __name__ == "__main__":
    sys.exit(main())
