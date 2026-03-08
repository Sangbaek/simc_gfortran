import numpy as np 
import pandas as pd 
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib
from glob import glob
import uproot


rcParams = {
    "font.size": 24,
    "axes.titlesize": 26,
    "axes.labelsize": 24,
    "xtick.labelsize": 22,
    "ytick.labelsize": 22,
    "legend.fontsize": 22,
    "figure.figsize": (10, 6),
    }
plt.rcParams.update(rcParams)

kine_settings = ["1", "2", "3", "4", "5", "6", "7"]
prefix        = ["a", "b", "b", "c", "d"]
kine_with_prefix = [k + p for p in prefix for k in kine_settings]
def beam_energy_by_prefix(prefix):
  if prefix in ["a", "b"]:
    return 850
  else:
    return 1645

def Q2_by_prefix(prefix):
  if prefix == "a":
    return 0.026
  elif prefix == "b":
    return 0.030
  elif prefix == "c":
    return 0.040
  elif prefix == "d":
    return 0.050

beam_energies = [beam_energy_by_prefix(p) for p in prefix for k in kine_settings]
Q2s           = [Q2_by_prefix(p) for p in prefix for k in kine_settings]
theta_pqs_nominal_cm     = [7,15,15,35,55,75,85,7,15,21,35,55,75,85,7,15,28,35,55,75,85,7,15,33,33,55,75,85]
phi_pqs_nominal_cm       = [0,0,180,180,180,180,180,0,180,0,180,180,180,180,0,180,0,180,180,180,180,0,180,0,180,180,180,180]

mpi0 = 0.135  # GeV/c^2
Ngen = 10000
selected_keys = ["Weight", "sigcm", "thcm", "phicm", "W", "Q2", "mmnuc"]
def main():
  for i in range(len(kine_with_prefix)):
    kine_setting = kine_with_prefix[i]
    beam_energy  = beam_energies[i]

    normfacs = []
    successful_random_seeds = []
    for random_seed in range(1000, 1001):
      filename = "nd26_{}_e{}_rs{}".format(kine_setting, beam_energy, random_seed)
      metafile = "outfiles/{}.hist".format(filename)
      if len(glob(metafile)) == 0:
        continue
      successful_random_seeds.append(random_seed)
      with open(metafile, "r") as f:
        lines = f.readlines()
        for line in lines:
          if "normfac" in line:
            normfacs.append(float(line.split()[-1]))
      
    normfac = np.mean(normfacs)
    dfs_kinematics = pd.DataFrame()

    #print("Kine setting: {}, Beam energy: {}, Successful random seeds: {}".format(kine_setting, beam_energy, successful_random_seeds))

    for random_seed in successful_random_seeds:
      filename = "nd26_{}_e{}_rs{}".format(kine_setting, beam_energy, random_seed)
      rootfile = uproot.open("worksim/{}.root".format(filename))
      tree     = rootfile["h666"]

      df_kinematics = pd.DataFrame()
      df_Weight  = ak.to_dataframe(tree["Weight"].arrays())
      df_W       = ak.to_dataframe(tree["W"].arrays())
      df_Q2    = ak.to_dataframe(tree["Q2"].arrays())
      df_thcm = ak.to_dataframe(tree["thcm"].arrays())
      df_mmnuc  = ak.to_dataframe(tree["mmnuc"].arrays())

      assert len(df_Weight) == len(df_W) == len(df_Q2) == len(df_thcm) == Ngen
      for key in tree.keys():
        df_key = ak.to_dataframe(tree[key].arrays())
        df_kinematics[key] = df_key
      df_kinematics.reset_index(drop=True, inplace=True)
      dfs_kinematics = pd.concat([dfs_kinematics, df_kinematics], ignore_index=True)
    
    dfs_kinematics.loc[:, "Weight"] = normfac * dfs_kinematics.Weight / len(dfs_kinematics)

    Q2_cond       = (dfs_kinematics.Q2 -Q2s[i]).abs()  < 0.005
    thcm_cond  = (dfs_kinematics.thcm - theta_pqs_nominal_cm[i]).abs() < 2.5
    W_cond        = (dfs_kinematics.W - 1.230).abs() < 0.01

    cond = Q2_cond & thcm_cond & W_cond

    print(r"{:.1f}E{:.0f}".format(dfs_kinematics.loc[cond].Weight.sum()/10**int(np.log10(dfs_kinematics.loc[cond].Weight.sum())), int(np.log10(dfs_kinematics.loc[cond].Weight.sum())))) 
    if i%7 ==6:
      print("")
    fig, axs = plt.subplots(2, 5, figsize = (25,10))

    Q2_bins      = np.linspace(Q2s[i]-0.01, Q2s[i]+0.01, 101)
    thcm_bins = np.linspace(theta_pqs_nominal_cm[i]-10.0, theta_pqs_nominal_cm[i]+10.0, 101)
    W_bins       = np.linspace(1.2, 1.3, 101)
    mmnuc_bins   = np.linspace(0.1, 0.3, 101)


    axs[0,0].hist(dfs_kinematics.Q2, bins=Q2_bins, weights=dfs_kinematics.Weight, histtype='step', color='k')
    axs[0,0].axvline(Q2s[i], color='r', linestyle='--')
    axs[0,0].fill_betweenx([0, axs[0,0].get_ylim()[1]], Q2s[i]-0.005, Q2s[i]+0.005, color='r', alpha=0.5)
    axs[0,0].set_xlabel(r"$Q^2$ (GeV/c)$^2$")
    axs[0,0].set_ylabel(r"Counts before cuts" + "\n" + r"(1/0.002 $(GeV/c)^2$)")
    axs[0,1].hist(dfs_kinematics.thcm, bins=thcm_bins, weights=dfs_kinematics.Weight, histtype='step', color='k')
    axs[0,1].axvline(theta_pqs_nominal_cm[i], color='r', linestyle='--') 
    axs[0,1].fill_betweenx([0, axs[0,1].get_ylim()[1]], theta_pqs_nominal_cm[i]-2.5, theta_pqs_nominal_cm[i]+2.5, color='k', alpha=0.5)
    axs[0,1].set_xlabel(r"$\theta_{pq, CM}$ (${}^{\circ}$)")
    axs[0,1].set_ylabel(r"Counts before cuts" + "\n" + r"(1/0.1 ${}^{\circ}$)")
    axs[0,2].hist(dfs_kinematics.W, bins=W_bins, weights=dfs_kinematics.Weight, histtype='step', color='k')
    axs[0,2].axvline(1.230, color='r', linestyle='--')
    axs[0,2].fill_betweenx([0, axs[0,2].get_ylim()[1]], 1.220, 1.240, color='r', alpha=0.5)
    axs[0,2].set_xlabel(r"$W$ (GeV)")
    axs[0,2].set_ylabel(r"Counts before cuts" + "\n" + r"(1/0.001 GeV)")
    axs[0,3].hist(dfs_kinematics.mmnuc, bins=mmnuc_bins, weights=dfs_kinematics.Weight, histtype='step', color='k')
    axs[0,3].axvline(mpi0, color='r', linestyle='--')
    axs[0,3].set_xlabel(r"$M_{miss}$ (GeV/c$^2$)")
    axs[0,3].set_ylabel(r"Counts before cuts" + "\n" + r"(1/0.002 GeV/c$^2$)")
    axs[0,4].hist(dfs_kinematics.phicm, bins=np.linspace(0, 360, 361), weights=dfs_kinematics.Weight, histtype='step', color='k')
    axs[0,4].axvline(phi_pqs_nominal_cm[i], color='r', linestyle='--') 
    axs[0,4].axvline(360 - phi_pqs_nominal_cm[i], color='r', linestyle='--') 
    axs[0,4].set_xlabel(r"$\phi_{CM}$ (${}^{\circ}$)")
    axs[0,4].set_ylabel(r"Counts before cuts" + "\n" + r"(1/1 ${}^{\circ}$)")
    axs[0,4].set_xlim(0, 360)
    axs[0,4].set_xticks([0, 90, 180, 270, 360])

    axs[1,0].hist(dfs_kinematics.loc[cond].Q2, bins=Q2_bins, weights=dfs_kinematics.loc[cond].Weight, histtype='step', color='k')
    axs[1,0].axvline(Q2s[i], color='r', linestyle='--')
    axs[1,0].set_xlabel(r"$Q^2$ (GeV/c)$^2$")
    axs[1,0].set_ylabel(r"Counts after cuts" + "\n" + r"(1/0.002 $(GeV/c)^2$)")
    axs[1,1].hist(dfs_kinematics.loc[cond].thcm, bins=thcm_bins, weights=dfs_kinematics.loc[cond].Weight, histtype='step', color='k')
    axs[1,1].axvline(theta_pqs_nominal_cm[i], color='r', linestyle='--') 
    axs[1,1].set_xlabel(r"$\theta_{pq, CM}$ (${}^{\circ}$)")
    axs[1,1].set_ylabel(r"Counts after cuts" + "\n" + r"(1/0.1 ${}^{\circ}$)")
    axs[1,2].hist(dfs_kinematics.loc[cond].W, bins=W_bins, weights=dfs_kinematics.loc[cond].Weight, histtype='step', color='k')
    axs[1,2].axvline(1.230, color='r', linestyle='--')
    axs[1,2].set_xlabel(r"$W$ (GeV)")
    axs[1,2].set_ylabel(r"Counts after cuts" + "\n" + r"(1/0.001 GeV)")
    axs[1,3].hist(dfs_kinematics.loc[cond].mmnuc, bins=mmnuc_bins, weights=dfs_kinematics.loc[cond].Weight, histtype='step', color='k')
    axs[1,3].axvline(mpi0, color='r', linestyle='--')
    axs[1,3].set_xlabel(r"$M_{miss}$ (GeV/c$^2$)")
    axs[1,3].set_ylabel(r"Counts after cuts" + "\n" + r"(1/0.002 GeV/c$^2$)")
    axs[1,4].hist(dfs_kinematics.loc[cond].phicm, bins=np.linspace(0, 360, 361), weights=dfs_kinematics.loc[cond].Weight, histtype='step', color='k')
    axs[1,4].axvline(phi_pqs_nominal_cm[i], color='r', linestyle='--') 
    axs[1,4].axvline(360 - phi_pqs_nominal_cm[i], color='r', linestyle='--') 
    axs[1,4].set_xlabel(r"$\phi_{CM}$ (${}^{\circ}$)")
    axs[1,4].set_ylabel(r"Counts after  cuts" + "\n" + r"(1/1 ${}^{\circ}$)")
    axs[1,4].set_xlim(0, 360)
    axs[1,4].set_xticks([0, 90, 180, 270, 360])

    fig.text(0.5, 0.95, r"Kin {}, Beam E {} MeV, $Q^2$ {} GeV$^2/c^2$, Estimated events after cut: ${{{:.1f}}}\times 10^{{{:.0f}}}$".format(kine_setting, beam_energy, Q2s[i], dfs_kinematics.loc[cond].Weight.sum()/10**int(np.log10(dfs_kinematics.loc[cond].Weight.sum())), int(np.log10(dfs_kinematics.loc[cond].Weight.sum()))), fontsize=30, ha='center')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.savefig("plots/nd26_{}.pdf".format(kine_setting), bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
  main()
