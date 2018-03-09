"""
Script to run MCCE simulation at different charges for water molecules.
"""
import os
import sys
import numpy as np
import mdtraj as md
from scipy import stats
from pymcce.automated_mcce import MCCEParams
from pymcce.mcce_simulation import Simulation
from pymcce.utils import write_watpdb_from_coords, get_last_prot_at_index


import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pylab
import seaborn as sns
sns.set(style="white")

fontsize = 10
font = {'family' : 'sans-serif',
        'weight' : 'normal',
        'size'   : fontsize}

plt.rc('font', **font)


def main():
    prefix = sys.argv[1]
    curr_dir = os.getcwd()
    n_conf = 1
    n_runs = 10
    charge_sets = ["t3p", "pm6"]
    # input files
    data_dir = os.path.abspath("../data/gramicidin/simulations/input_struct")


    mu_z = []
    for n in range(n_runs):
        prefix = sys.argv[1] + "_%02d" % n
        step2out = os.path.join(data_dir, prefix + "_step2_out.pdb")
        print(step2out)
        v = np.array([0, 0, 1])
        with open(step2out, "r") as f1:
            lines = f1.readlines()
            prot_lines = [l[31:76].split() for l in lines if "EM_" not in l and "HOH" not in l]
            n_prot_atoms = len(prot_lines)
            prot_coords = np.zeros((n_prot_atoms, 3))
            prot_charges = np.zeros((n_prot_atoms, 1))
            for index, l in enumerate(prot_lines):
                x, y , z = float(l[0]), float(l[1]), float(l[2])
                c = float(l[4])
                prot_coords[index :] += [x, y, z]
                prot_charges[index :] = c
            pos = -1.0 * prot_coords.T
            chg = prot_charges[:, 0]
            dp = pos.dot(chg)
            #proj = np.multiply(dp.dot(v) / v.dot(v), v)
            scaling_factor = dp.dot(v) / v.dot(v)
            print(scaling_factor)

    os.chdir(curr_dir)


if __name__ == "__main__":
    status = sys.exit(main())
