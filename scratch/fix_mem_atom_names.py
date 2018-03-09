"""
Script to fix membrane atom names in step2_out.pdb, originally generated from prep_mem_gramicidin_structure.py
"""

import os
import sys


def main():
    curr_dir = os.getcwd()
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    pdb = os.path.join(input_dir, "step2_out.pdb")
    fixed_pdb = os.path.join(input_dir, "fixed_step2_out.pdb")
    with open(pdb, "r") as f1:
        with open(fixed_pdb, "w") as f2:
            lines = f1.readlines()
            for l in lines:
                if "EM_ _" in l:
                    l = l.replace("x", "C")
                    l = l.replace("EM_ _", "MEM M")
                    f2.write(l)
                else:
                    f2.write(l)

if __name__ == "__main__":
    sys.exit(main())
