from Simulations_Bulk import run_sims
from multiprocessing import Pool
import numpy as np

if __name__ == "__main__":

    conn_fun = "Waxman"
    eta = 1.0

    Torus = False

    L = 1000

    rt = 1

    mu_list_overall = np.arange(0.15, 0.4, 0.01)

    # mu_list = mu_list_overall[0:4]
    # mu_list = mu_list_overall[4:9]
    # mu_list = mu_list_overall[9:16]
    mu_list = mu_list_overall[16:25]


    for mu in mu_list:
        run_sims(1000, conn_fun, eta, L, Torus, mu, rt)

