import numpy as np

def landprice(param, fund, L, w, dist, nobs):
    # parameters
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # fund[:,0]=a; fund[:,1]=b; fund[:,2]=H; 
    a = fund[:,0]
    b = fund[:,1]
    H = fund[:,2]

    r = ((1 - alpha) / alpha) * ((w * L) / H)

    return r