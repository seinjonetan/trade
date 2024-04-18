import numpy as np
from scipy.special import gamma

def welfaregains(param, Ctradesh, tradesh, CL, L, nobs):
    # parameters;
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # delta function;
    deltaf = gamma((epsilon - 1) / epsilon)
    # gamma function;
    gammaf = gamma((theta + 1 - sigma) / theta)

    # domestic trade share;
    dtradesh = np.diag(tradesh)
    Cdtradesh = np.diag(Ctradesh)

    # welfare gains;
    welfgain = ((dtradesh / Cdtradesh) ** (alpha / theta)) * ((L / CL) ** ((1 / epsilon) + (1 - alpha)))

    return welfgain