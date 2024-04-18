import numpy as np
from scipy.special import gamma

def welfare(param, fund, L, w, tradesh, dist, nobs, sigma, LLwest, LLeast):
    # parameters;
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # fund[:,0]=a; fund[:,1]=b; fund[:,2]=H; fund[:,3]=Iwest; fund[:,4]=Ieast;
    a = fund[:,0]
    b = fund[:,1]
    H = fund[:,2]
    Iwest = fund[:,3]
    Ieast = fund[:,4]

    # delta function;
    deltaf = gamma((epsilon - 1) / epsilon)
    # gamma function;
    gammaf = gamma((theta + 1 - sigma) / theta)

    # domestic trade share;
    dtradesh = np.diag(tradesh)

    # welfare;
    welf = deltaf * (b ** (1 / epsilon)) * ((a / dtradesh) ** (alpha / theta)) * (H ** (1 - alpha)) * (L ** (-((1 / epsilon) + (1 - alpha))))
    welf[Iwest == 1] = welf[Iwest == 1] / (alpha * (((1 - alpha) / alpha) ** (1 - alpha)) * (gammaf ** alpha) * (LLwest ** (-1 / epsilon)))
    welf[Ieast == 1] = welf[Ieast == 1] / (alpha * (((1 - alpha) / alpha) ** (1 - alpha)) * (gammaf ** alpha) * (LLeast ** (-1 / epsilon)))

    return welf
