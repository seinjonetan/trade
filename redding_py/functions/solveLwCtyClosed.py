import numpy as np
import time

def solveLwCtyClosed(param, fund, dwght, dist, nobs, LL, LLwest, LLeast):
    # Start timer
    xtic = time.time()

    # param=[alpha theta epsilon]
    alpha, theta, epsilon = param

    # fund[:,0]=a; fund[:,1]=b; fund[:,2]=H; fund[:,3]=Iwest; fund[:,4]=Ieast
    a, b, H, Iwest, Ieast = fund.T

    # convergence indicator
    wconverge = 0
    Lconverge = 0

    # Initializations
    L_i = np.ones(nobs) * (LL / nobs)
    w_i = np.ones(nobs)

    # trade costs
    dd = dist ** -theta
    dd = dwght * dd

    print('>>>> Start Wage and Population Convergence <<<<')

    xx = 1
    while xx < 2000:

        # Start inner loop for wages
        x = 1
        while x < 2000:

            # Trade share
            pwmat = a * (w_i ** -theta) * np.ones(nobs)
            nummat = dd * pwmat
            denom = np.sum(nummat)
            denommat = np.ones((nobs, 1)) * denom
            tradesh = nummat / denommat

            # Income equals expenditure
            income = w_i * L_i
            expend = np.dot(tradesh, income)

            # Convergence criterion
            income_r = np.round(income * 10**6)
            expend_r = np.round(expend * 10**6)

            # Update loop
            if np.all(income_r == expend_r):
                x = 10000
                wconverge = 1
            else:
                w_e = w_i * (expend / income) ** (1 / theta)
                w_i = 0.25 * w_e + 0.75 * w_i
                # Normalization
                w_i[Iwest == 1] /= np.mean(w_i[Iwest == 1])
                w_i[Ieast == 1] /= np.mean(w_i[Ieast == 1])
                wconverge = 0
                x += 1

        # End inner loop for wages

        # Domestic trade share
        dtradesh = np.diag(tradesh)

        # Population
        num = b * (a / dtradesh) ** (alpha * epsilon / theta) * (L_i / H) ** (-epsilon * (1 - alpha))
        L_e = np.zeros(nobs)
        L_e[Iwest == 1] = num[Iwest == 1] / np.sum(num[Iwest == 1])
        L_e[Ieast == 1] = num[Ieast == 1] / np.sum(num[Ieast == 1])
        L_e[Iwest == 1] *= LLwest
        L_e[Ieast == 1] *= LLeast

        # Convergence criterion
        L_i_r = np.round(L_i * 10**6)
        L_e_r = np.round(L_e * 10**6)

        # Update loop
        if np.all(L_i_r == L_e_r):
            xx = 10000
            Lconverge = 1
        else:
            L_e = L_i * (L_e / L_i) ** (1 / (epsilon * (1 - alpha)))
            L_i = 0.25 * L_e + 0.75 * L_i
            Lconverge = 0
            xx += 1
    
    # End outer loop for population
            
    return w_i, L_i, tradesh, dtradesh, Lconverge, wconverge, xtic