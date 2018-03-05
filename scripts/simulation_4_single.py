"""
Script to run MCCE simulation at different charges for water molecules.
"""
import numpy as np
import os
import time
import sys
from shutil import copy2, rmtree
from subprocess import call

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
    # get current working dir
    curr_dir = os.getcwd()
    # set dict for using different water charges
    charge_sets = {"t3p" : [-0.834, 0.417, 0.417], "pm6" : [-2.000, 1.000, 1.000]}
    # fix number of conformers per water to 1
    n_conf = 5
    # set number of separate MC runs
    n_runs = 1
    # obtain water charge set option from the command-line
    prefix = sys.argv[1]
    charge_set = charge_sets[prefix]

    # set input directories
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    home_dir = os.path.abspath("../data/gramicidin/simulations/input_struct_%s" % prefix)
    mcce_exec_dir = os.path.abspath("../code/mcce3.5")
    water_dir = os.path.abspath("../data/gramicidin/water_placement")
    
    # set filenames
    src_run_prm = os.path.join(input_dir, "run.prm")
    src_newtpl = os.path.join(input_dir, "new.tpl")

    # run n different MC simulations, each on a different protein conf and coresponding water orientations
    for n in range(n_runs):
        run_dir = os.path.join(home_dir, "run_%02d" % n)
        print("Run number %d, dir %s" % (n, run_dir))

        if os.path.exists(run_dir):
            rmtree(run_dir)
        os.makedirs(run_dir)
        os.chdir(run_dir)
        # copy files to where MCCE will be run
        copy2(src_run_prm, run_dir)
        copy2(src_newtpl, run_dir)

        # obtain correct step2out
        src_step2out = os.path.join(input_dir, "%02d_final_step2_out.pdb" % n)
        write_step2out = os.path.join(run_dir, "temp_%02d_step2_out.pdb" % n)
        #print(src_step2out)
        #print(write_step2out)
        copy2(src_step2out, write_step2out)

        # obtain corresponding water positions and generate conformers
        init_waters = os.path.join(water_dir, "%02d_init_wat_placement.pdb" % n)
        start_index = get_last_prot_at_index(src_step2out)
        coords = generate_conformers(init_waters, n_conf)
        write_filename = os.path.join(run_dir, "water_confs")
        write_watpdb_from_coords_ext(write_filename, coords, start_index, charge_set)
        water_conformers = write_filename + ".pdb"

        # assemble full step2out and write ms_gold
        last_prot_at = get_last_prot_at_index(write_step2out)
        assemble_step2_out(water_conformers, write_step2out, "step2_out.pdb", last_prot_at)
        write_msgold("step2_out.pdb", "W", "HOH")
        """
        with open("head3.lst", "r") as f1, open("temp_head3.lst", "w") as f2:
            lines = f1.readlines()
            for l in lines:
                if "CTR-" in l:
                    new_line = l[:21] + "t" + l[22: ]
                    f2.write(new_line)
                else:
                    f2.write(l)
        copy2("temp_head3.lst", "head3.lst")
        """

        # set up MCCE run
        m = MCCEParams(mcce_exec_dir)
        m.edit_parameters(INPDB="step2_out.pdb",
                          DO_PREMCCE="f", DO_ROTAMERS="f", DO_ENERGY="t",
                          DO_MONTE="t", MONTE_ADV_OPT="f", MONTE_MS="t",
                          MONTE_REDUCE="0.0001", BIG_PAIRWISE="1.0", MONTE_RUNS="20")
        m.write_runprm(run_dir + "/")
        # run MCCE
        m.write_submitsh("", run_name=str(n))
        call("qsub submit.sh", shell=True)
        time.sleep(10)
        os.chdir(curr_dir)




if __name__ == "__main__":
    status = sys.exit(main())
