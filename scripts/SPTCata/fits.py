import numpy as np
from scipy.optimize import curve_fit
from scipy import stats

def msd_track(track) :
    """
    This function calculates the MSD of a single track. `track` needs to be
    an np.array with dimensions (T, 2)
    """
    # get the length of the track
    T = track.shape[0]
    
    # initialize a matrix that will allow us to calculate the MSD: the matrix
    # has T-1 rows, corresponding to each time point except the last one (because
    # for the last one we cannot calculate any difference in time); it has then
    # T-1 columns, so that it stores all the possible time intervals from 1 to T
    delta = np.zeros((T-1, T-1))
    delta.fill(np.nan)
    
    # calculate the MSD for the single track
    for t0 in range(T-1) :
        for delay in range(1,T-t0) :
            t1 = t0 + delay
            pos1 = track[t1,:]
            pos0 = track[t0,:]
            delta[t0, delay-1] = np.sum((pos1-pos0)**2)
    
    # now we can return the mean of this track
    return np.nanmean(delta, axis=0)

# this is the function that we will use to calculate the MSD
def msd_t(tracks) :
    
    # init the variables
    ntracks = len(tracks)
    lengths = [t.shape[0] for t in tracks]
    max_t = int(max(lengths)) - 1
    
    # We init a matrix that will contain the MSD calculated for each of the
    # tracks. It will have `ntracks` rows and `max_t - 1` columns.
    MSD = np.zeros((ntracks, max_t))
    MSD.fill(np.nan)

    # iterate through the all the tracks
    for i, track in enumerate(tracks) :
        m = msd_track(track)
        MSD[i, 0:len(m)] = m
    
    # return the MSD
    return np.nanmean(MSD, axis=0)

def diffusion_model_1s(X, *p):
    """
    One species diffusion model.
    """
    x = X[0]
    dt = X[1]
    D = p[0]
    return 1-np.exp(-x/(4*D*dt))


def diffusion_model_2s(X, *p) :
    """
    Two species diffusion model.
    """
    x = X[0]
    dt = X[1]
    D1, D2, f1 = p
    return 1-f1*np.exp(-x/(4*D1*dt))-(1-f1)*np.exp(-x/(4*D2*dt))

def fit_diffusion_models(r2_tracks, m, M, N, dt,
                         p0_1s = 10.0, p0_2s = [10.0, 1.0, 0.1]) :
    
    # calculate the histogram of square displacements of the tracks
    bins = np.linspace(m, M, N+1)
    h, bins = np.histogram(r2_tracks, bins=bins)
    
    # prepare x and y for fitting
    delta = bins[1]-bins[0]
    x = np.zeros(N+1)
    y = np.zeros(N+1)
    x = bins + delta/2.
    y[1:] = np.cumsum(h)/sum(h)
    X = (x, dt)
    
    # one species diffusion model fit
    coeff_1s, var_matrix_1s = curve_fit(diffusion_model_1s, X, y, p0=p0_1s)
    
    # two species diffusion model fit
    coeff_2s, var_matrix_2s = curve_fit(diffusion_model_2s, X, y, p0=p0_2s,
                                       bounds = (0,
                                                 [np.inf, np.inf, 1.0]))
    
    return x, y, coeff_1s, var_matrix_1s, coeff_2s, var_matrix_2s

def linear_regression (x,y,prob) :
    """
    Fit (x,y) to a linear function, using maximum likelihood estimation of the
    confidence intervals on the coefficients, given the user-supplied
    probability *prob*
    """
    n = len(x)
    xy = x*y
    xx = x*x
    # estimates
    xmean = x.mean()
    ymean = y.mean()
    xxmean = xx.mean()
    xymean = xy.mean()
    b1 = (xymean-xmean*ymean) / (xxmean-xmean**2)
    b0 = ymean-b1*xmean
    s2 = 1./n * sum([(y[i] - b0 - b1 * x[i])**2 for i in range(n)])
    #confidence intervals
    alpha = 1 - prob
    c1 = stats.chi2.ppf(alpha/2.,n-2)
    c2 = stats.chi2.ppf(1-alpha/2.,n-2)
    # get results and return
    c = -1 * stats.t.ppf(alpha/2.,n-2)
    bb1 = c * (s2 / ((n-2) * (xxmean - (xmean)**2)))**.5
    bb0 = c * ((s2 / (n-2)) * (1 + (xmean)**2 / (xxmean - xmean**2)))**.5
    return b0,b1,bb0,bb1

def fit_msd (msd,cutoff,delta_t) :
    # prepare the values to fit: exclude the first value because it is zero
    t = np.arange(1, msd.size+1)*delta_t
    x = np.log(t[:cutoff])
    y = np.log(msd[:cutoff])
    # perform fit to y = ax + b with their errors
    b,a,db,da = linear_regression (x,y,0.99)
    # now convert the value of b into a diffusion coefficient
    D = np.exp(b)/4
    dD = np.exp(db)/4
    return a,da,D,dD
