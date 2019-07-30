import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys, os
import SPTCata as spt
from collections import defaultdict

def idx_grid(val, minval, delta) :
    return int((val-minval)//delta)

def make_grid(experiment, Nx = 200, Ny = 200, qmin = 30.0) :
    # the file containing the map
    map_fname = '%s/Spots in tracks statistics_MAP.txt'%(experiment.datadir)
    
    # don't attempt to do anything if the file doesn't exist
    if not os.path.exists(map_fname) :
        return False

    # parse the file
    density_map = pd.read_csv(map_fname, sep='\t')
    
    # load all the positions in the tracks
    X = np.array(density_map['POSITION_X'])
    Y = np.array(density_map['POSITION_Y'])
    Q = np.array(density_map['QUALITY'])
    
    # make the grid
    dx = (X.max()-X.min())/Nx
    dy = (Y.max()-Y.min())/Ny
    grid = np.zeros((Nx+1, Ny+1), dtype=np.int32)

    # fill the grid
    experiment.xm = X.min()
    experiment.ym = Y.min()
    for x, y, q in zip(X, Y, Q) :
        if q >= qmin :
            i = idx_grid(x, experiment.xm, dx)
            j = idx_grid(y, experiment.ym, dy)
            grid[i, j] += 1
            
    # put the grid into the experiment
    experiment.grid = grid
    experiment.dx = dx
    experiment.dy = dy
    
    return True

def grid_analysis(experiment) :
            
    # do the calculations
    density = [[] for i in range(experiment.grid.max())]
    counter = 0
    for n, trajectory in enumerate(experiment.trajectory_spots) :
        for x, y in trajectory[1:] :

            # get the current value of the displacement
            displacement = experiment.displacement_spots[counter]

            # get the grid indices
            i = idx_grid(x, experiment.xm, experiment.dx)
            j = idx_grid(y, experiment.ym, experiment.dy)

            # add the current value to the value of the density map
            try :
                density[experiment.grid[i,j]-1].append(displacement)
            except IndexError :
                print(experiment.datadir, n)
                break
            counter += 1
    
    # perform the averages and return
    dvd = np.zeros(experiment.grid.max())
    for i, d in enumerate(density) :
        dvd[i] = np.mean(density[i])
        
    # put the results in
    experiment.dvd = dvd

######################################

catadata_root = '%s/work/CRG/projects/catadata'%(os.getenv('HOME'))
datadirs = {
    'Olaparib_R5020' : '%s/data/1_Olaparib_R5020'%(catadata_root),
    'R5020' : '%s/data/2_DMSO_R5020_Control'%(catadata_root),
    'EtOH' : '%s/data/3_EtOH_Nohormone_Control'%(catadata_root)
}

# the quality
quality = 50.0

# temporal and spatial scale
dt = 0.015             # in seconds (15 ms)
scale_l = 0.160        # in microns (160 um)

# cycle through all the directories and do the analysis
experiments = defaultdict(list)
for treatment, datadir in datadirs.items() :
    for subdir in os.listdir(datadir) :
        full_dir_name = '%s/%s'%(datadir, subdir)
        
        # in this part we load the map
        map_fname = '%s/Spots in tracks statistics_MAP.txt'%(full_dir_name)
        if os.path.exists(map_fname) :
            experiment = spt.SPT(full_dir_name, 1.0, links=False, quality=quality)
            make_grid(experiment)
            grid_analysis(experiment)
            experiments[treatment].append(experiment)

figures_dir = '../../figures/grid_trajectories/'
for treatment, experiment_batch in experiments.items() :
    for i, experiment in enumerate(experiment_batch) :
        fig = plt.figure(figsize=(10,10))
        plt.imshow(np.log(1+experiment.grid), cmap=plt.cm.Greys, aspect='equal', origin='lower')
        for trajectory in experiment.trajectory_spots :
            plt.plot(trajectory[:, 1]/experiment.dy, trajectory[:, 0]/experiment.dx, alpha = 0.5)
            
        plt.title(experiment.datadir.split('/')[-2] + " " + experiment.datadir.split('/')[-1])
        fig.tight_layout()
        fig.savefig('%s/%s-%d.png'%(figures_dir, treatment, i))

for treatment, experiment_batch in experiments.items() :
    fig = plt.figure(figsize=(6,4))
    for i, experiment in enumerate(experiment_batch) :
        plt.plot(experiment.dvd)
    plt.title(treatment, fontsize=24)
    plt.xlabel("Number of passages in map", fontsize=18)
    plt.ylabel("Average displacement", fontsize=18)
    plt.xlim(0, 40)
    fig.tight_layout()
    fig.savefig("%s/%s-dvd.pdf"%(figures_dir, treatment))
