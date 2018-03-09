"""
Script to pre-process input gramicidin structure to generate file ready for simulation.
This prepares only the gramicidin structure, the memebrane or water molecules are not placed.
"""

import os
import sys
from pymcce.automated_mcce import MCCEParams


def main():
    curr_dir = os.getcwd()
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    m = MCCEParams(os.path.abspath("../code/mcce3.5"))
    m.edit_parameters(DO_PREMCCE="t", DO_ROTAMERS="f", TERMINALS="t")
    os.chdir(input_dir)
    terminal_atoms = ['CY', 'OY', 'NT', 'C1', 'C2', 'OG']

    with open('noh_test_01.pdb', "r") as oldfile, open('input_struct_01.pdb', 'w') as newfile:
        for line in oldfile:
            atm_name = line[13:15]
            if not any(at == atm_name for at in terminal_atoms):
                newfile.write(line)

    m.edit_parameters(INPDB=input_dir + "/input_struct_01.pdb")
    m.write_runprm(input_dir + "/")

    run_command = m.mcce_directory + "/mcce > prep_struct_run.log"
    os.system(run_command)
    os.chdir(curr_dir)


if __name__ == "__main__":
    sys.exit(main())
