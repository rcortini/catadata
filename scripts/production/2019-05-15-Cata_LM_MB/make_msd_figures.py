import os, sys
import numpy as np
import matplotlib.pyplot as plt
import itertools

# parameters
dt = 0.015             # in seconds (15 ms)
cutoff = 12

# directories
catadata_root = '%s/work/CRG/projects/catadata'%(os.getenv('HOME'))
tannic_acid_datadir = '%s/data/msd/TannicAcid'%(catadata_root)
olaparib_datadir = '%s/data/msd/Olaparib'%(catadata_root)

# for aesthetics
import matplotlib
matplotlib.rcParams['text.usetex'] = True

def load_msd_data(data, msd_data_dir) :
    msd = {}
    msd_fits = {}
    for name, datafile in data.items() :
        msd_fname = '%s/%s-msd.dat'%(msd_data_dir, datafile)
        msd_fit_fname = '%s/%s-msd_fit.dat'%(msd_data_dir, datafile)
        msd[name] = np.loadtxt(msd_fname)
        msd_fits[name] = np.loadtxt(msd_fit_fname)
    return msd, msd_fits


def plot_msd(title, msd, msd_fits, out_fname) :
    # prepare figure
    fig = plt.figure(figsize=(8,6))
    ax = plt.subplot(111)

    # markers and colors
    marker = itertools.cycle(('+', '.', '^'))
    colors = itertools.cycle(('r', 'g', 'b'))

    # plot for all the treatment conditions
    for treatment, m in msd.items() :
        
        # t axis
        t = dt*np.arange(cutoff)
        
        # get the fit values
        a, da, D, dD = msd_fits[treatment]
        
        # the diffusion model: 2D diffusion
        yfit = 4 * D * t ** a
        
        # plot the data and the fit to the model
        color = next(colors)
        ax.loglog(dt*np.arange(m.shape[0]), m, next(marker),
              color = color, markersize = 8, label = r'%s: $\alpha$ = %.3f +/-'
                  '%.3f\par D = %.3f +/- %.3f $\mu m^2/s$'%(treatment, a, da, D, dD))
        ax.loglog(t, yfit, '--', color = color)

    # finish up with plot style
    ax.legend(loc='lower left', fontsize = 14)
    ax.set_xlabel(r'$log[\Delta t]$ [s]', fontsize = 18)
    ax.set_ylabel(r'$\log [MSD (\Delta t)]$', fontsize = 18)
    ax.set_title(title, fontsize=24)

    # and save
    fig.tight_layout()
    fig.savefig(out_fname)


########################
# Tannic Acid
########################

data = {
    'Tannic Acid' : 'Tannic-Acid+R5020',
    'Tannic Acid Control' : 'Tannic-Acid-Control',
}

msd, msd_fits = load_msd_data(data, tannic_acid_datadir)
out_fname = '%s/figures/2019-05-15-Cata_LM_MB/Tannic-Acid-msd.pdf'%(catadata_root)
plot_msd('Tannic Acid', msd, msd_fits, out_fname)

########################
# Tannic Acid Subset
########################

data = {
    'Tannic Acid' : 'Tannic-Acid-Subset+R5020',
    'Tannic Acid Control' : 'Tannic-Acid-Subset-R5020-Control',
}

msd, msd_fits = load_msd_data(data, tannic_acid_datadir)
out_fname = '%s/figures/2019-05-15-Cata_LM_MB/TannicAcid-Subset-msd.pdf'%(catadata_root)
plot_msd('Tannic Acid Subset', msd, msd_fits, out_fname)

########################
# Olaparib
########################

data = {
    'Olaparib' : 'Olaparib+R5020',
    'Olaparib Control' : 'Olaparib-R5020-Control',
    'Olaparib No Hormone' : 'EtOH-NoHormone-Control',
}

msd, msd_fits = load_msd_data(data, olaparib_datadir)
out_fname = '%s/figures/2019-05-15-Cata_LM_MB/Olaparib-msd.pdf'%(catadata_root)
plot_msd('Olaparib', msd, msd_fits, out_fname)

########################
# Olaparib Subset
########################

data = {
    'Olaparib' : 'Olaparib-Subset+R5020',
    'Olaparib Control' : 'Olaparib-Subset-R5020-Control',
    'Olaparib No Hormone' : 'Olaparib-Subset-EtOH-Control',
}

msd, msd_fits = load_msd_data(data, olaparib_datadir)
out_fname = '%s/figures/2019-05-15-Cata_LM_MB/Olaparib-Subset-msd.pdf'%(catadata_root)
plot_msd('Olaparib Subset', msd, msd_fits, out_fname)

################################
# Olaparib without Olaparib
################################

data = {
    'Olaparib Control' : 'Olaparib-R5020-Control',
    'Olaparib No Hormone' : 'EtOH-NoHormone-Control',
}

msd, msd_fits = load_msd_data(data, olaparib_datadir)
out_fname = '%s/figures/2019-05-15-Cata_LM_MB/Olaparib-no-Olaparib-msd.pdf'%(catadata_root)
plot_msd('R5020 and Control', msd, msd_fits, out_fname)
