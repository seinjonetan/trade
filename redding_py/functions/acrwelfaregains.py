import numpy as np

def acrwelfaregains(param, Ctradesh, tradesh, CL, L, nobs):
    # parameters
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # domestic trade share
    dtradesh = np.diag(tradesh)
    Cdtradesh = np.diag(Ctradesh)

    # welfare gains
    acrwelfgain = (dtradesh / Cdtradesh) ** (alpha / theta)

    return acrwelfgain