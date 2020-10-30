from mofapy2.run.entry_point import entry_point
import pandas as pd
from mofapy2.simulate import simulate_mofa as simmofa
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as ss
import sys

if sys.argv[1] == "time_aware":
    time_aware = True
else:
    time_aware = False

# simulate data
# set simulation parameters
G = 2
K = 3
N = 50
s_sim = [1, 1, 0]
ls_sim = [0.2, 0.2, 0.0]
sharedness = [True, False, False] # one shared, one non-shared, one non-smooth

# simulate data
sim = simmofa.simulate_data(N=N, seed=2020, views=["0", "1", "2"], D=[50, 50, 50],
                            K=K, G=G, lscales=ls_sim, noise_level=0.01, scales=s_sim,
                            shared=sharedness, plot=False)
data = sim['data']
sample_cov = sim['sample_cov']

# import matplotlib.pyplot as plt
# for k in range(K):
#     plt.figure(k)
#     for g in range(G):
#         plt.scatter(sample_cov[g], sim['Z'][g][:,k])
#     plt.savefig("factor%s.png" % k)


# mask data
data = simmofa.mask_samples(sim, 0.3, perc_all_views=0.3)

# prepare MOFA model
ent = entry_point()
ent.set_data_options(scale_views=False)

views_names = ["view_1", "view_2", "view_3"]
grp_names = ["group_1", "group_2"]
ent.set_data_matrix(data,
                    views_names=views_names,
                    groups_names=grp_names
                    )
ent.set_model_options(factors=K)
ent.set_train_options()
# for time-aware multi-modal FA with GP model add covariates and smooth_options
if time_aware:
    ent.set_covariates(sim['sample_cov'], covariates_names = ["time"])
    ent.set_smooth_options(model_groups = True)


ent.build()
ent.run()

if time_aware:
    x = np.linspace(0, 1, 50)
    x = np.array(x)
    x = x[:, None]
    ent.predict_factor(new_covariates=x)

# Save the model
if time_aware:
    outfile = "out/toy_model.hdf5"
else:
    outfile = "out/toy_model_noGP.hdf5"
ent.save(outfile, save_data=True)