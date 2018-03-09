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
    curr_dir = os.getcwd()
    charge_sets = ["t3p", "pm6"]
    chg = sys.argv[1]
    if chg not in charge_sets:
        sys.exit("Charge set %s not allowed!" % chg)
    # input directories
    input_dir = os.path.abspath("../data/gramicidin/prep_structures")
    run_dir = os.path.abspath("../data/gramicidin/simulations/wat_charges")
    mcce_exec_dir = os.path.abspath("../code/mcce3.5")
    water_dir = os.path.abspath("../data/gramicidin/water_placement")
    
    # filenames
    src_run_prm = os.path.join(input_dir, "run.prm")
    src_step2out = os.path.join(input_dir, "fixed_step2_out.pdb")
    src_newtpl = os.path.join(input_dir, "new.tpl")
    src_ms_gold = os.path.join(input_dir, "new.tpl")
    # copy files to where MCCE will be run
    copy2(src_run_prm, run_dir)
    copy2(src_step2out, run_dir + "/temp_step2_out.pdb")
    copy2(src_step2out, run_dir + "/ms_gold")
    copy2(src_newtpl, run_dir)

    prefix = chg + "_"
    n_conf = "%02d" % 1
    water_conformers = os.path.join(water_dir, str(n_conf) + "_n_conf_perwater_" + chg + ".pdb")    

    # assemble full step2out
    os.chdir(run_dir)
    last_prot_at = get_last_prot_at_index("temp_step2_out.pdb")
    assemble_step2_out(water_conformers, "temp_step2_out.pdb", "step2_out.pdb", last_prot_at)
    write_msgold("step2_out.pdb", "W", "HOH")
    # set up MCCE run
    m = MCCEParams(mcce_exec_dir)
    m.edit_parameters(INPDB="step2_out.pdb",
                    DO_PREMCCE="f", DO_ROTAMERS="f", DO_ENERGY="f",
                    DO_MONTE="t", MONTE_ADV_OPT="f", MONTE_MS="t",
                    MONTE_REDUCE="0.0", BIG_PAIRWISE="-1.0")
    m.write_runprm(run_dir + "/")
    
    # run MCCE
    run_command = m.mcce_directory + "/mcce > " + run_dir + "/" + prefix + "run.log"
    print("Running MCCE ...")
    print(run_command)
    copy2(prefix + "energies.opp", "energies.opp")
    
    os.system(run_command)
    
    # post-processing
    copy2("fort.38", prefix + "fort.38")
    copy2("energies.opp", prefix + "energies.opp")
    copy2("head3.lst", prefix + "head3.lst")
    copy2("ms.dat", prefix + "ms.dat")
    copy2("step2_out.pdb", prefix + "step2_out.pdb")
    
    os.chdir(curr_dir)


if __name__ == "__main__":
    status = sys.exit(main())
