"""
Script to generate water conformers.
"""

import os
import sys
import numpy as np
import mdtraj as md
from pymcce.sampler import write_watpdb_from_coords, get_last_prot_at_index





def main():
    charge_sets = {"tip3p" : [-0.834, 0.417, 0.417], "pm6" : [-2.000, 1.000, 1.000]}
    curr_dir = os.getcwd()
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    pdb_file = os.path.join(input_dir, "step2_out.pdb")
    water_dir = os.path.abspath("../data/gramicidin/water_placement")
    charge_set = charge_sets["pm6"]
    init_waters = os.path.join(water_dir, "init_wat_placement.pdb")
    start_index = get_last_prot_at_index(pdb_file)
    #for n_conf in [1, 5, 10, 15, 20, 25]:
    for n_conf in [1]:
        n = "%02d" % n_conf
        print(n)
        coords = generate_conformers(init_waters, n_conf)
        write_filename = os.path.join(water_dir, n + "_n_conf_perwater" + "_pm6")
        write_watpdb_from_coords_ext(write_filename, coords, start_index, charge_set)
 
if __name__ == "__main__":
    sys.exit(main())
