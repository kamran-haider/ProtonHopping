"""
Script to run MCCE simulation at different charges for water molecules.
"""
import os
import sys
import numpy as np
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
    water_charges = {"t3p" : [-0.834, 0.417, 0.417], "pm6" : [-2.000, 1.000, 1.000]}
    prefix = sys.argv[1]
    curr_dir = os.getcwd()
    n_conf = 5
    n_runs = 1
    charge_sets = ["t3p", "pm6"]
    chg = sys.argv[1]
    if chg not in charge_sets:
        sys.exit("Charge set %s not allowed!" % chg)
    charges = water_charges[chg]
    # input files
    data_dir = os.path.abspath("../data/gramicidin/simulations/input_struct_%s" % chg)
    mu = []
    amide_mu = []
    for n in range(n_runs):
        mu_z = []
        prefix = "run_%02d" % n
        print("Processing %s:" % prefix)
        msdat = os.path.join(data_dir, prefix, "ms.dat")
        head3lst = os.path.join(data_dir, prefix, "head3.lst")
        fort38 = os.path.join(data_dir, prefix, "fort.38")
        energies_opp = os.path.join(data_dir, prefix, "energies.opp")
        step2out = os.path.join(data_dir, prefix, "step2_out.pdb")

        msa = Simulation(msdat, head3lst, fort38)
        msa.parse_trajectory(sample_frequency=10)
        msa.parse_struct(step2out, n_conf)
        msa.calculate_water_dipoles(charges)

if __name__ == "__main__":
    status = sys.exit(main())
