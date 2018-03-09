"""
Script to insert gramicidin structure in a membrane (actually a low dielectric slab of neutral atoms to mimic membrane).
The input structure is the outut from generate_starting_structure.py.
"""

import os
import sys


def main():
    curr_dir = os.getcwd()
    # set directory paths
    input_dir = os.path.abspath("../data/gramicidin/prep_structures/")
    ipece_dir =  os.path.abspath("../code/ipece")
    # set file paths
    ipece_parm = os.path.join(ipece_dir, input_dir, "ipece.prm")
    input_pdb = os.path.join(ipece_dir, input_dir, "step1_out.pdb")
    output_pdb = os.path.join(ipece_dir, input_dir, "mem_step1_out.pdb")
    # assemble ipece command and run
    os.chdir(input_dir)
    run_command = ipece_dir + "/ipece " + ipece_parm + " " + input_pdb + " " + output_pdb
    os.system(run_command)
    os.chdir(curr_dir)


if __name__ == "__main__":
    sys.exit(main())
