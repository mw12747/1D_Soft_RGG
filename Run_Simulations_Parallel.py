from Simulations import run_sims
from multiprocessing import Pool
import numpy as np
from scipy.special import gamma

if __name__ == "__main__":

    conn_fun = "Waxman"
    eta = 1.0

    new = True
    Torus = False

    L = 10000

    rt = 1

    mean_degree_temp = np.arange(15.6, 20, 0.3)
    mu_list_overall = [((2 / x)*gamma((eta + 1)/eta))**eta for x in mean_degree_temp]

    # mu_list = mu_list_overall[0:5]
    # mu_list = mu_list_overall[5:9]
    # mu_list = mu_list_overall[9:12]
    mu_list = mu_list_overall[12:15]


    # with Pool(processes=4) as pool:
    #     ntlist = [10 + i for i in range(25)]
    #
    #     data = [(ntrials, conn_fun, eta, L, Torus, mu, rt) for ntrials in ntlist
    #             for mu in mu_list]
    #
    #     results_1 = [pool.apply_async(run_sims, args=(data[i])) for i in range(len(data))]
    #
    #     pool.close()
    #     pool.join()

    for mu in mu_list:
        run_sims(500, conn_fun, eta, L, Torus, mu, new, rt)

