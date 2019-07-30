import numpy as np
import matplotlib.pyplot as plt

def plot_angles(angles, bins, title) :
    counts, angle_edges = np.histogram(angles, bins = bins)
    angle_centers = angle_edges[1:] - np.ediff1d(angle_edges)
    fig = plt.figure(figsize = (5,5))
    ax = plt.subplot(111, projection = 'polar')
    bars = ax.bar(angle_centers, counts, width = 2*np.pi/bins, edgecolor = 'k')
    ax.set_title(title, y = 1.1, fontsize = 18)
    ax.set_yticklabels([])
    return fig, ax
