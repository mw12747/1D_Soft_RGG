from Simulations_Beta import run_sims
import numpy as np

conn_fun = "Rayleigh"
eta = 1001.0

Torus = False

L = 1001
ntrials = 500

mu_list_temp = np.arange(0.2, 2.0, 0.1)

# mu_list = mu_list_temp[0:6]
# mu_list = mu_list_temp[6:12]
mu_list = mu_list_temp[12:18]

for mu in mu_list:
    print(mu)
    run_sims(ntrials, conn_fun, eta, L, Torus, mu)
