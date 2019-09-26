import numpy as np
import pandas as pd
import sys, os, re
from collections import defaultdict

def load_spt(datadir, quality = None) :
    # load the SPT data file
    spt_fname = '%s/Spots in tracks statistics.txt'%(datadir)
    spt = pd.read_csv(spt_fname, sep='\t')
    
    # group by TRACK_ID
    spt_by_track = spt.groupby('TRACK_ID')

    # get quality of tracks and list of excluded tracks
    if quality is not None :
        track_statistics = pd.read_csv('%s/Track statistics.csv'%(datadir))
        idx_tracks_to_exclude = track_statistics.TRACK_MEAN_QUALITY < quality
        tracks_to_exclude = track_statistics[idx_tracks_to_exclude].TRACK_ID
    else :
        tracks_to_exclude = []
    
    # extract trajectories
    trajectories = []
    for track_id, track in spt_by_track :

        # skip tracks that did not have sufficiently high average quality
        if track_id in tracks_to_exclude.values :
            continue

        # extract x and y from the trajectory
        x = track['POSITION_X']
        y = track['POSITION_Y']

        # finally, append the current trajectory to the list of trajectories
        trajectories.append(np.array([x, y]).T)
    
    return trajectories

def load_experiments(root_datadir, exp_names, bad_cells=[], quality=None) :
    """
    """

    # init the output data structure
    experiments = defaultdict(list)

    # iterate over the subdirectories. Notice that the variable `subdir` will contain only
    # the name of the subdirectory, not the full path of it.
    for exp_name in exp_names :
        datadir = '{}/{}'.format(root_datadir, exp_name)

        for d, sds, fs in os.walk(datadir) :

            # here we test the name of the directory. We use a regular expression to check
            # whether the name of the subdirectory contains a format "StackN_CellN".
            for sd in sds :
                if re.match('Stack[0-9]+_Cell[0-9]+', sd) is None :
                    continue

                # exclude the cells that are in the list of bad cells
                exp_id = '{}/{}'.format(exp_name, sd)
                if exp_id not in bad_cells :
                    full_dir = '%s/%s'%(d, sd)
                    experiments[datadir].append(load_spt(full_dir, quality=quality))
    return experiments
