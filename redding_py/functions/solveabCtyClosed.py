import numpy as np
import time
from scipy.stats.mstats import gmean

def solveabCtyClosed(param, observe, dwght, dist, nobs, LLwest, LLeast):
    # start timer
    xtic = time.time()

    # parameters
    alpha = param[0]
    theta = param[1]
    epsilon = param[2]

    # observations
    L = observe[:, 0]
    w = observe[:, 1]
    H = observe[:, 2]
    Iwest = observe[:, 3]
    Ieast = observe[:, 4]

    # convergence indicators
    aconverge = 0
    bconverge = 0

    # initializations
    a_i = np.ones(nobs)
    b_i = np.ones(nobs)

    # trade costs
    dd = dist ** -theta
    dd = dwght * dd

    print('>>>> Start productivity and amenities Convergence <<<<')

    xx = 1
    while xx < 2000:
        # Start inner loop for productivity
        x = 1
        while x < 2000:
            # Trade share
            pwmat = (a_i * (w ** -theta))[:, None]
            nummat = dd * pwmat
            denom = np.sum(nummat)
            denommat = np.ones((nobs, 1)) * denom
            tradesh = nummat / denommat
            # test
            test = np.sum(tradesh)
            mntest = np.mean(test)
            # Income equals expenditure
            income = w * L
            expend = tradesh @ income
            # Convergence criterion
            income_r = np.round(income * (10 ** 6))
            expend_r = np.round(expend * (10 ** 6))
            x = np.max(np.abs(income_r - expend_r))

            # Update loop
            if np.array_equal(income_r, expend_r):
                # print('>>>> Productivity Convergence Achieved <<<<')
                x = 10000
                aconverge = 1
            else:
                a_e = a_i * (income / expend)
                a_i = 0.25 * a_e + 0.75 * a_i
                # Normalization
                # Separate countries
                a_i[Iwest == 1] = a_i[Iwest == 1] / gmean(a_i[Iwest == 1])
                a_i[Ieast == 1] = a_i[Ieast == 1] / gmean(a_i[Ieast == 1])
                aconverge = 0
                x += 1

        # domestic trade share
        dtradesh = np.diag(tradesh)

        # population
        num = b_i * ((a_i / dtradesh) ** (alpha * epsilon / theta)) * ((L / H) ** (-epsilon * (1 - alpha)))
        L_e = np.zeros(nobs)
        L_e[Iwest == 1] = num[Iwest == 1] / np.sum(num[Iwest == 1])
        L_e[Ieast == 1] = num[Ieast == 1] / np.sum(num[Ieast == 1])
        L_e[Iwest == 1] = L_e[Iwest == 1] * LLwest
        L_e[Ieast == 1] = L_e[Ieast == 1] * LLeast

        # Convergence criterion
        L_r = np.round(L * (10 ** 6))
        L_e_r = np.round(L_e * (10 ** 6))
        gap = np.max(np.abs(L_e_r - L_r))

        # Update loop
        if gap == 0:
            # print('>>>> Population Convergence Achieved <<<<')
            xx = 10000
            bconverge = 1
        else:
            b_e = b_i * (L / L_e)
            b_i = 0.25 * b_e + 0.75 * b_i
            # Normalization
            # Separate countries
            b_i[Iwest == 1] = b_i[Iwest == 1] / gmean(b_i[Iwest == 1])
            b_i[Ieast == 1] = b_i[Ieast == 1] / gmean(b_i[Ieast == 1])
            bconverge = 0
            xx += 1

    xtic = time.time() - xtic

    return a_i, b_i, dtradesh, aconverge, bconverge, xtic