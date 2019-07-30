import SPTCata as spt
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, sys
from collections import defaultdict
import itertools

# for aesthetics
import matplotlib
matplotlib.rcParams['text.usetex'] = True

#####################################
# PARAMETERS OF ANALYSIS
#####################################

# the quality threshold
quality = 50.0

# temporal and spatial scale
dt = 0.015             # in seconds (15 ms)
scale_l = 1.0          # in microns
cutoff = 12            # number of values to consider for the MSD fit
N = 200                # number of bins in diffusion model
m = 0.0                # minimum displacement value for diffusion model fit
M = 0.05               # maximum displacement value for diffusion model fit

#####################################
# LOAD DATA AND PERFORM ANALYSIS
#####################################

catadata_rootdir = '%s/work/CRG/projects/catadata'%(os.getenv('HOME'))

# check for proper invocation
if len(sys.argv) < 3 :
    print("Usage: %s <experiment_label> <msd_outdir> <datadir(s)>"%(sys.argv[0]))
    sys.exit(1)

# get command line arguments
experiment_label = sys.argv[1]
msd_outdir = sys.argv[2]
datadirs = sys.argv[3:]

# init data structures that will hold the data
experiments = []
track_length = []
tracks = []
angles = []
r2 = []

# cycle through the subdirectories of the treatments
for datadir in datadirs :
    for subdir in os.listdir(datadir) :

        # load the experiment data and append it to the `experiments` dictionary
        full_dir_name = '%s/%s'%(datadir, subdir)
        experiment = spt.SPT(full_dir_name, scale_l, links=False, quality=quality)
        experiments.append(experiment)

        # extract all the tracks from this experiment and put them into a
        # separate list
        tracks.extend([t for t in experiment.trajectory_spots])

        # calculate the square displacement of each track in this experiment and
        # put it into a dictionary
        r2.extend(np.array(experiment.displacement_spots)**2)

        # perform the angle analysis for the given treatment
        experiment.angle_analysis()
        angles.extend(experiment.angles)

# perform the MSD analysis associated to this treatment
msd = spt.msd_t(tracks)
msd_fit = spt.fit_msd(msd, cutoff, dt)

# save output of MSD... it will be plotted later
msd_fname = '%s/%s-msd.dat'%(msd_outdir, experiment_label)
np.savetxt(msd_fname, msd)
msd_fit_fname = '%s/%s-msd_fit.dat'%(msd_outdir, experiment_label)
np.savetxt(msd_fit_fname, msd_fit)

# do the angle analysis for the entire treatment
angles = np.array(angles)
angles = angles[~np.isnan(angles)]

#####################################
# FITTING AND PLOTTING
#####################################

# plot of fit of diffusion coefficient
x, y, coeff_1s, var_matrix_1s, coeff_2s, var_matrix_2s =\
spt.fit_diffusion_models(r2, m, M, N, dt)

fig = plt.figure(figsize=(6,4))
plt.plot(x, y, 'o', markersize=5, label='Data')
plt.plot(x, spt.diffusion_model_1s((x, dt), coeff_1s[0]), label='One species fit')
plt.plot(x, spt.diffusion_model_2s((x, dt), coeff_2s[0],
                                  coeff_2s[1],
                                  coeff_2s[2]), label='Two species fit')
plt.title(experiment_label, fontsize=24)
plt.axhline(y = 1.0, linestyle = '--', color = 'k', linewidth = 0.75)
plt.xlabel("Displacement [$\mu$m]")
plt.ylabel("Cumulative distribution")
plt.xlim(0., 0.05)
plt.xticks(np.arange(0., 0.0205, 0.005))
plt.legend(loc='lower right', fontsize = 16)

plt.text(0.007, 0.8, "One species fit: D = %.3f $\mu m^2/s$"%(coeff_1s[0]), fontsize = 14)
plt.text(0.007, 0.5, "Two species fit: D1 = %.3f $\mu m^2/s$\n"
                      "                D2 = %.3f $\mu m^2/s$\n"
                      "                f1 = %.3f"%(coeff_2s[0],
                                                            coeff_2s[1],
                                                            coeff_2s[2]), fontsize = 14)   

# save figure
fig.tight_layout()
fname = '%s/figures/2019-05-15-Cata_LM_MB/r2-distribution-%s.pdf'%(catadata_rootdir,
                                                                   experiment_label)
fig.savefig(fname)

# angles plot
fname = '%s/figures/2019-05-15-Cata_LM_MB/angles-%s.pdf'%(catadata_rootdir,
                                                                   experiment_label)
fig, ax = spt.plot_angles(angles, 24, experiment_label)
fig.tight_layout()
fig.savefig(fname)
