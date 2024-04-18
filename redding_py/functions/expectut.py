import numpy as np
from scipy.special import gamma

def expectut(param, fund, L, w, tradesh, dist, nobs, sigma):
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

    # expected utility
    EU = b * (gammaf ** (-alpha * epsilon)) * (alpha ** -epsilon) * (((1 - alpha) / alpha) ** (-epsilon * (1 - alpha)))
    EU = ((a / dtradesh) ** (alpha * epsilon / theta)) * ((L / H) ** (-epsilon * (1 - alpha))) * EU
    EU[Iwest == 1] = deltaf * (np.sum(EU[Iwest == 1]) ** (1 / epsilon))
    EU[Ieast == 1] = deltaf * (np.sum(EU[Ieast == 1]) ** (1 / epsilon))

    return EU