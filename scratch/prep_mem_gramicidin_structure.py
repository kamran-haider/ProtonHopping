"""
Script to prepare gramicidin structure within membrane, by re-rerunning mcce steps 1 and 2
on the output pdb from ipece which is generated from insert_membrane.py.
"""

import os
import sys
from automated_mcce import MCCEParams


def main():
    curr_dir = os.getcwd()
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    m = MCCEParams(os.path.abspath("../code/mcce3.5"))
    m.edit_parameters(DO_PREMCCE="t", DO_ROTAMERS="t", INPDB=input_dir + "/mem_step1_out.pdb")
    m.write_runprm(input_dir + "/")
    os.chdir(input_dir)
    run_command = m.mcce_directory + "/mcce > prep_struct_run.log"
    os.system(run_command)
    os.chdir(curr_dir)


if __name__ == "__main__":
    sys.exit(main())
