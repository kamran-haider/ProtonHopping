"""
Script to run MCCE simulation at different charges for water molecules.
"""
import numpy as np
import os
import sys
from shutil import copy2
import mdtraj as md
from pymcce.automated_mcce import MCCEParams
from pymcce.utils import write_watpdb_from_coords, get_last_prot_at_index
from pymcce.sampler import write_watpdb_from_coords_ext, generate_conformers


# define a function that assembles a step2out pdb file
def assemble_step2_out(water_pdb, other_pdb, write_pdb, last_prot_index):
    with open(water_pdb, "r") as f1:
        with open(other_pdb, "r") as f2:
            with open(write_pdb, "w") as f3:
                wat_lines = f1.readlines()
                l2 = f2.readlines()
                prot_lines = l2[:last_prot_index-1]
                mem_lines = l2[last_prot_index-1:]
                #print(prot_lines[-1])
                #print(wat_lines[0])
                #print(mem_lines[0])
                #print(mem_lines[-1])
                all_lines = prot_lines + wat_lines + mem_lines
                for l in all_lines:
                    f3.write(l)

def write_msgold(pdb_file, chain_id, resname):
    st = md.load_pdb(pdb_file, no_boxchk=True)
    with open("ms_gold", "w") as f:
        index = 1
        for res in st.topology.residues:
            if resname in res.name:
                s = "%s%s%04d\n" % (res.name, chain_id, index)
                index += 1
                f.write(s)

def main():

    # get current working directory
    curr_dir = os.getcwd()
    # set dict for using different water charges
    charge_sets = {"t3p" : [-0.834, 0.417, 0.417], "pm6" : [-2.000, 1.000, 1.000]}
    # fix number of conformers per water to 1
    n_conf = 1

    # obtain water charge set option from the command-line
    prefix = sys.argv[1]
    charge_set = charge_sets[prefix]

    # set input directories
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    run_dir = os.path.abspath("../data/gramicidin/simulations/vdw_scaling_new")
    mcce_exec_dir = os.path.abspath("../code/mcce3.5")
    water_dir = os.path.abspath("../data/gramicidin/water_placement")
    
    # set filenames
    src_run_prm = os.path.join(input_dir, "run.prm")
    src_step2out = os.path.join(input_dir, "00_step2_out.pdb")
    src_newtpl = os.path.join(input_dir, "new.tpl")
    src_ms_gold = os.path.join(input_dir, "new.tpl")

    # copy files to where MCCE will be run
    copy2(src_run_prm, run_dir)
    copy2(src_step2out, run_dir + "/temp_step2_out.pdb")
    copy2(src_newtpl, run_dir)

    os.chdir(run_dir)

    # get assembled full step2out
    copy2("../wat_charges/" + prefix + "_step2_out.pdb", "step2_out.pdb")
    copy2("../wat_charges/" + prefix + "_head3.lst", "head3.lst")
    write_msgold("step2_out.pdb", "W", "HOH")
    # set up MCCE run
    scalings = [1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.00]
    n_runs = len(scalings)

    copy2("../wat_charges/" + prefix + "_energies.opp", "energies.opp")
    for n in range(n_runs):
        with open("extra.tpl", "w") as f1:
            f1.write("SCALING  VDW0        1.00\n")
            f1.write("SCALING  VDW1        1.00\n")
            f1.write("SCALING  VDW         %s\n" % str(scalings[n]))
            f1.write("SCALING  TORS       1.0\n")

        m = MCCEParams(mcce_exec_dir)
        m.edit_parameters(INPDB="step2_out.pdb",
                          DO_PREMCCE="f", DO_ROTAMERS="f", DO_ENERGY="f",
                          DO_MONTE="t", MONTE_ADV_OPT="f", MONTE_MS="t",
                          MONTE_REDUCE="0.0", BIG_PAIRWISE="-0.5")
        m.write_runprm(run_dir + "/")
        prefix = sys.argv[1] + "_%02d" % n
        run_command = m.mcce_directory + "/mcce > " + run_dir + "/" + prefix + "_run.log"
        print("Running MCCE ...")
        print(run_command)
        os.system(run_command)
        # post-processing
        copy2("ms.dat", prefix + "_ms.dat")
        copy2("fort.38", prefix + "_fort.38")
        copy2("head3.lst", prefix + "_head3.lst")        

    os.chdir(curr_dir)


if __name__ == "__main__":
    status = sys.exit(main())
