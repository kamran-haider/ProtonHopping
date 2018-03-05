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
        test = calc_amide_dipole(step2out)
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
            #if  t_dp > -0.2 and t_dp < 0.2:
                #print("Found config %d with desired dipole %g, saving.." % (i, sum(dps)))
                #print("Energy = %g" % msa.energies[i])
                #snapshot_st = msa.structure.atom_slice(wat_ids)
                #snapshot_st.save_pdb("test_pdbs/%05d_config.pdb" % i)
                #print("Microstate %d: Energy =  %g, Dipole = %g" % (i, msa.energies[i], t_dp))

        mu_z = np.array([sum(m) for m in mu_z])
        mu.append(np.mean(mu_z))
        amide_mu.append(test)
        figure = pylab.figure(figsize=(7.5, 3.5), dpi=300)
        y, bin_edges = np.histogram(mu_z, bins=50)
        kernel = stats.gaussian_kde(mu_z)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
        plt.hist(mu_z, bins=50, normed=1, histtype='stepfilled')
        p_y = kernel.evaluate(bin_centers)
        plt.plot(bin_centers, p_y, 'k-')
        plt.xlim(-6.0, 6.0)
        plt.ylim(0.0, np.max(p_y) + 0.2)
        #plt.ylim(0.0, 1.2)
        style = dict(size=10, color='gray')
        plt.text(-5.8, 0.1, "%.2f" % test, **style)
        plt.xlabel(r'$\mu_{z,total}$', size=14)
        plt.ylabel(r'$P(\mu_{z,total})$', size=14)
        plt.tight_layout()
        plt.savefig("sim_04_%s_histogram_total_dipole.png" % prefix)
            

    os.chdir(curr_dir)
    figure = pylab.figure(figsize=(4.5, 4.5), dpi=300)
    plt.scatter(mu, amide_mu)
    plt.xlabel(r'$\mu_{z,water}$', size=14)
    plt.ylabel(r'$\mu_{z,amide}$', size=14)
    plt.tight_layout()
    plt.savefig("%s_wat_amide_dipole.png" % chg)

if __name__ == "__main__":
    status = sys.exit(main())
