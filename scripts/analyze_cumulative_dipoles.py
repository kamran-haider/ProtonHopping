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


def calc_amide_dipole(step2out):
    v = np.array([0, 0, 1])
    with open(step2out, "r") as f1:
        lines = f1.readlines()
        prot_lines = [l[31:76].split() for l in lines if "EM_" not in l and "HOH" not in l]
        n_prot_atoms = len(prot_lines)
        prot_coords = np.zeros((n_prot_atoms, 3))
        prot_charges = np.zeros((n_prot_atoms, 1))

        amide_lines = [l[31:76].split() for l in lines if l[13:15] == "N " or l[13:15] == "H " and "HOH" not in l]
        n_amide_atoms = len(amide_lines)
        amide_coords = np.zeros((n_amide_atoms, 3))
        amide_charges = np.zeros((n_amide_atoms, 1))

        for index, l in enumerate(amide_lines):
            x, y, z = float(l[0]), float(l[1]), float(l[2])
            c = float(l[4])
            amide_coords[index:] += [x, y, z]
            amide_charges[index:] = c
        pos = -1.0 * amide_coords.T
        chg = amide_charges[:, 0]
        dp = pos.dot(chg)
        # proj = np.multiply(dp.dot(v) / v.dot(v), v)
        scaling_factor = dp.dot(v) / v.dot(v)
        return scaling_factor


def calc_muz(n, n_conf, chg, charges, data_dir, save=False, save_freq=1000):
    mu_z = []
    prefix = "run_%02d" % n
    print("Processing %s:" % prefix)
    msdat = os.path.join(data_dir, prefix, "ms.dat")
    head3lst = os.path.join(data_dir, prefix, "head3.lst")
    fort38 = os.path.join(data_dir, prefix, "fort.38")
    energies_opp = os.path.join(data_dir, prefix, "energies.opp")
    step2out = os.path.join(data_dir, prefix, "step2_out.pdb")

    msa = Simulation(msdat, head3lst, fort38)
    msa.parse_trajectory(sample_frequency=100)
    msa.parse_struct(step2out, n_conf)
    msa.calculate_water_dipoles(charges)
    pdb_save_dir = os.path.abspath("pdbs/apolar_gram_%s/%s" % (chg, prefix))
    if not os.path.exists(pdb_save_dir):
        os.makedirs(pdb_save_dir)
    for i in range(msa.trajectory.shape[0]):
        microstate_conf_ids = msa.trajectory[i, :]
        wat_ids = []
        dps = []
        for c in microstate_conf_ids:
            dp = msa.conformer_data[msa.conf_id_name_map[c + 1]][-1]
            if dp is not None:
                dps.append(dp)
                wat_ids.extend(msa.conformer_data[msa.conf_id_name_map[c + 1]][2])
        mu_z.append(dps)
        t_dp = sum(dps)
        if save:
            if not i % save_freq:
                print("Saving config %d with dipole %g." % (i, sum(dps)))
                print("Energy = %g" % msa.energies[i])
                snapshot_st = msa.structure.atom_slice(wat_ids)
                snapshot_st.save_pdb("%s/%05d_config.pdb" % (pdb_save_dir, i))
                print("Microstate %d: Energy =  %g, Dipole = %g" % (i, msa.energies[i], t_dp))

    mu_z = np.array([np.cumsum(m) for m in mu_z])
    return mu_z

def main():
    water_charges = {"t3p" : [-0.834, 0.417, 0.417], "pm6" : [-2.000, 1.000, 1.000]}
    prefix = sys.argv[1]
    curr_dir = os.getcwd()
    n_conf = 5
    n_runs = 10
    charge_sets = ["t3p", "pm6"]
    chg = sys.argv[1]
    if chg not in charge_sets:
        sys.exit("Charge set %s not allowed!" % chg)
    charges = water_charges[chg]
    # input files
    run_type_1 = "apolar_gram"
    data_dir_1 = os.path.abspath("../data/gramicidin/simulations/%s_%s" % (run_type_1, chg))
    run_type_2 = "input_struct"
    data_dir_2 = os.path.abspath("../data/gramicidin/simulations/%s_%s" % (run_type_2, chg))

    for n in range(n_runs):
        mu_z_1 = calc_muz(n, n_conf, chg, charges, data_dir_1)
        mu_z_2 = calc_muz(n, n_conf, chg, charges, data_dir_2)

        figure = pylab.figure(figsize=(7.5, 3.5), dpi=300)
        for index, m in enumerate(mu_z_1):
            plt.plot(m, 'k-', color="red", alpha=0.01, linewidth=1.0)
        for index, m in enumerate(mu_z_2):
            plt.plot(m, 'k-', color="blue", alpha=0.01, linewidth=1.0)

        #plt.xlim(-6.0, 6.0)
        #plt.ylim(0.0, np.max(p_y_2) + 0.2)
        plt.xlabel(r'$\mu_{z,total}$', size=14)
        plt.ylabel(r'$P(\mu_{z,total})$', size=14)
        plt.tight_layout()
        plt.savefig("images/%s_cumulative_dipole_%02d" % (chg, n))
            

    os.chdir(curr_dir)

if __name__ == "__main__":
    status = sys.exit(main())
