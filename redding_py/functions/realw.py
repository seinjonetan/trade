import numpy as np
from scipy.special import gamma

def realw(param, fund, L, w, tradesh, dist, nobs, sigma):
    # parameters
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # fund[:,0]=a; fund[:,1]=b; fund[:,2]=H; fund[:,3]=Iwest; fund[:,4]=Ieast; 
    a = fund[:,0]
    b = fund[:,1]
    H = fund[:,2]
    Iwest = fund[:,3]
    Ieast = fund[:,4]

    # delta function
    deltaf = gamma((epsilon - 1) / epsilon)
    # gamma function
    gammaf = gamma((theta + 1 - sigma) / theta)

    # domestic trade share
    dtradesh = np.diag(tradesh)

    # real wage
    realwage = ((a / dtradesh) ** (alpha / theta)) * ((L / H) ** (-(1 - alpha)))
    realwage = realwage / (alpha * (gammaf ** alpha) * (((1 - alpha) / alpha) ** (1 - alpha)))

    return realwage