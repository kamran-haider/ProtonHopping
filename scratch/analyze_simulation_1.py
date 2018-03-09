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
    charge_sets = ["t3p", "pm6"]
    chg = sys.argv[1]
    if chg not in charge_sets:
        sys.exit("Charge set %s not allowed!" % chg)
    charges = water_charges[chg]
    # input files
    data_dir = os.path.abspath("../data/gramicidin/simulations/apolar_gram_t3p/run_00")

    msdat = os.path.join(data_dir, "ms.dat")
    head3lst = os.path.join(data_dir, "head3.lst")
    fort38 = os.path.join(data_dir, "fort.38")
    energies_opp = os.path.join(data_dir, "energies.opp")
    step2out = os.path.join(data_dir, "step2_out.pdb")

    msa = Simulation(msdat, head3lst, fort38)
    msa.parse_trajectory(sample_frequency=100)

    print(msa.trajectory.shape)
    #print(msa.conformer_data)
    msa.parse_struct(step2out, n_conf)
    msa.calculate_water_dipoles(charges)

    mu_z = []
    save_freq = 1000
    for i in range(msa.trajectory.shape[0]):
        microstate_conf_ids = msa.trajectory[i, :]
        #print(microstate_conf_ids)
        #sys.exit()
        wat_ids = []
        dps = []
        for c in microstate_conf_ids:
            dp = msa.conformer_data[msa.conf_id_name_map[c + 1]][-1]
            if dp is not None:
                dps.append(dp)
                wat_ids.extend(msa.conformer_data[msa.conf_id_name_map[c + 1]][2])
        mu_z.append(dps)
        t_dp = sum(dps)
        if not i % save_freq:
            snapshot_st = msa.structure.atom_slice(wat_ids)
            snapshot_st.save_pdb("test_pdbs/%05d_config.pdb" % i)
            print("Microstate %d: Energy =  %g, Dipole = %g" % (i, msa.energies[i], t_dp))
            #break

    mu_z = np.array([sum(m) for m in mu_z])
    figure = pylab.figure(figsize=(7.5, 3.5), dpi=300)    
    y, bin_edges = np.histogram(mu_z, bins=50)
    kernel = stats.gaussian_kde(mu_z)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    plt.hist(mu_z, bins=50, normed=1, histtype='stepfilled')
    p_y = kernel.evaluate(bin_centers)
    plt.plot(bin_centers, p_y, 'k-')
    #plt.xlim(-6.0, 10.0)
    #plt.ylim(0.0, np.max(p_x) + 0.1)
    #plt.ylim(0.0, np.max(p_x) + 0.1)

    plt.xlabel(r'$\mu_{z,total}$', size=14)
    plt.ylabel(r'$P(\mu_{z,total})$', size=14)
    plt.tight_layout()
    plt.savefig("sim_01_%s_histogram_total_dipole.png" % str(chg))


    os.chdir(curr_dir)


if __name__ == "__main__":
    status = sys.exit(main())