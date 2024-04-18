import numpy as np
from scipy.stats import gmean

def solveabCtyOpen(param, observe, dwght, dist, nobs):
    # param=[alpha theta epsilon]
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # observe[:,0]=L; observe[:,1]=w; observe[:,2]=H; observe[:,3]=Iwest; observe[:,4]=Ieast
    L = observe[:,0]
    w = observe[:,1]
    H = observe[:,2]
    Iwest = observe[:,3]
    Ieast = observe[:,4]

    # convergence indicator
    aconverge = 0
    bconverge = 0

    # Initializations
    a_i = np.ones(nobs)
    b_i = np.ones(nobs)

    # trade costs
    dd = np.power(dist, -theta)
    dd = dwght * dd

    print('>>>> Start productivity and amenities Convergence <<<<')

    # Start outer loop for amenities
    xx = 1
    while xx < 2000:
        # Start inner loop for productivity
        x = 1
        while x < 2000:
            # Trade share
            pwmat = (a_i * (w ** -theta))[:, None]
            nummat = dd * pwmat
            denom = nummat.sum(axis=0)
            tradesh = nummat / denom

            # Income equals expenditure
            income = w * L
            expend = tradesh @ income

            # Convergence criterion
            income_r = np.round(income * 10**6)
            expend_r = np.round(expend * 10**6)

            # Update loop
            if np.all(income_r == expend_r):
                x = 10000
                aconverge = 1
            else:  
                a_e = a_i * (income / expend)
                a_i = 0.25 * a_e + 0.75 * a_i
                # Normalization
                a_i[Iwest == 1] /= gmean(a_i[Iwest == 1])
                aconverge = 0
                x += 1

        # Domestic trade share
        dtradesh = np.diag(tradesh)

        # Population
        num = b_i * ((a_i / dtradesh) ** (alpha * epsilon / theta)) * ((L / H) ** (-epsilon * (1 - alpha)))
        L_e = np.zeros(nobs)
        L_e[Iwest == 1] = num[Iwest == 1] / num[Iwest == 1].sum()
        L_e[Ieast == 1] = num[Ieast == 1] / num[Ieast == 1].sum()
        L_e[Iwest == 1] *= LLwest
        L_e[Ieast == 1] *= LLeast

        # Convergence criterion
        L_r = np.round(L * 10**6)
        L_e_r = np.round(L_e * 10**6)
        gap = np.abs(L_e_r - L_r).max()

        # Update loop
        if gap == 0:
            xx = 10000
            bconverge = 1
        else:  
            b_e = b_i * (L / L_e)
            b_i = 0.25 * b_e + 0.75 * b_i
            # Normalization
            b_i[Iwest == 1] /= gmean(b_i[Iwest == 1])
            b_i[Ieast == 1] /= gmean(b_i[Ieast == 1])
            bconverge = 0
            xx += 1

    return a_i, b_i