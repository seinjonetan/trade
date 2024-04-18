import numpy as np
from scipy.special import gamma

def pindex(param, fund, w, dtradesh, nobs, sigma):
    # parameters
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # fund[:,0]=a; fund[:,1]=b; fund[:,2]=H; 
    a = fund[:,0]
    b = fund[:,1]
    H = fund[:,2]

    # gamma function
    gammaf = gamma((theta + 1 - sigma) / theta)

    # price index
    P = ((gammaf ** -theta) * a * (w ** -theta) / dtradesh) ** (-1 / theta)

    return P