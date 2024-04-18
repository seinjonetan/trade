import numpy as np

def mobwelfaregains(param, Ctradesh, tradesh, CL, L, nobs):
    # parameters
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # domestic trade share
    dtradesh = np.diag(tradesh)
    Cdtradesh = np.diag(Ctradesh)

    # welfare gains
    welfgain = (dtradesh / Cdtradesh) ** (alpha / theta) * (L / CL) ** (1 - alpha)

    return welfgain