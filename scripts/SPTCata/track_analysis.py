import pandas as pd
import numpy as np

def trajectories(data, scale_l) :
    
    # extract the tracks
    tracks = data.groupby('TRACK_ID')
    trajs = []
    for name, track in tracks :
        try :
            x = track['EDGE_X_LOCATION']
            y = track['EDGE_Y_LOCATION']
        except KeyError :
            x = track['POSITION_X']
            y = track['POSITION_Y']
        trajs.append(scale_l * np.array([x, y]).T)
    return trajs

def angle(v1, v2) :
    d = np.dot(v1, v2)
    v = v1[0]*v2[1] - v1[1]*v2[0]
    nv1 = np.linalg.norm(v1)
    nv2 = np.linalg.norm(v2)
    N = d/(nv1*nv2)
    return np.sign(v)*np.arccos(N)


def angle_analysis(trajectories) :

    # iterate through the trajectories
    angles = []
    for traj in trajectories :
        
        # calculate the difference in x values and y values
        delta_X = np.diff(traj[:, 0], axis=0)[1:]
        delta_Y = np.diff(traj[:, 1], axis=0)[1:]

        # calculate the normalized tangent vectors
        T = np.array([delta_X, delta_Y]).T
        
        # angle analysis
        for i in range(1, len(T)) :
            angles.append(angle(T[i-1], T[i]))
    
    return angles

class SPT :
    
    def __init__(self, datadir, scale_l, links = False, quality = None) :
        
        # data directory
        self.datadir = datadir
        self.spots_fname = '%s/Spots in tracks statistics.txt'%(self.datadir)
        
        # load data
        self.spots = pd.read_csv(self.spots_fname, sep = '\t')

        # process tracks
        tracks_spots = self.spots.groupby('TRACK_ID')
        
        # get quality of tracks and list of excluded tracks
        if quality is not None :
            track_statistics = pd.read_csv('%s/Track statistics.csv'%(datadir))
            idx_tracks_to_keep = track_statistics.TRACK_MEAN_QUALITY > quality
            tracks_to_keep = track_statistics[idx_tracks_to_keep].TRACK_ID
            self.spots = tracks_spots.filter(lambda x : x.name in
                                             tracks_to_keep.values)

        # process trajectories
        self.trajectory_spots = trajectories(self.spots, scale_l)
        
        # process displacement
        self.displacement_spots = []
        for traj in self.trajectory_spots :
            self.displacement_spots.extend(np.linalg.norm(np.diff(traj, axis = 0), axis = 1))

        # if user requests to process links, do it
        if links :
            self.links_fname = '%s/Links in tracks statistics.csv'%(self.datadir)
            self.links = pd.read_csv(self.links_fname)
            tracks_links = self.links.groupby('TRACK_ID')
            self.trajectory_links = trajectories(self.links, scale_l)
            self.displacement_links = scale_l * self.links['DISPLACEMENT']

    def angle_analysis(self) :
        # process angles
        self.angles = angle_analysis(self.trajectory_spots)

