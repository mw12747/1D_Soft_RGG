from Simulated_Data import import_data_bulk, import_data
import matplotlib.pyplot as plt
import numpy as np
from math import log

def poisson_approx(*kargs):
    if len(kargs) == 1:
        return [1 - np.exp(-x) for x in kargs[0]]
    else:
        return [1 - np.exp(-x - y) for (x, y) in zip(kargs[0], kargs[1])]


conn_fun = "Waxman"
eta = 1.0
L = 1000
Torus = False
data_type = "Overall"

x_list = np.arange(0.0001, 0.89, 0.0001)
poi = [1 - np.exp(-L*np.exp(-2 / x)) for x in x_list]

y_list = np.arange(0, 1, 0.0001)
mu_crit = [2 / log(L) for i in y_list]


mu_bulk, P_fc_bulk, P_iso_bulk, P_gap_bulk, P_split_bulk, \
mean_degree_bulk, n_iso_bulk, n_gaps_bulk = import_data_bulk(conn_fun, eta, L, Torus, data_type)

mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps = import_data(conn_fun, eta, L, Torus, data_type)

plt.figure()

plt.title("Simulations for L = %s, %s conn. fun, $\eta$ =  %s" % (L, conn_fun, eta))

plt.xlabel("$\mu$")
plt.ylabel("Prob. of event")

# plt.plot(mu, P_fc, 'o', label="P_fc")
plt.plot(mu_bulk, P_iso_bulk, '.', label="P_iso_bulk")
# plt.plot(mu_bulk, P_gap_bulk, '.', label="P_gap_bulk")
plt.plot(mu, P_iso, 'x', label="P_iso")
# plt.plot(mu, P_gap, 'x', label="P_gap")

# plt.plot(mu, poisson_approx(n_iso, n_gaps), 'x', label="Poisson Approx. P_fc")
# plt.plot(mu, poisson_approx(n_iso), 'x', label="Poisson Approx. P_iso")
# plt.plot(mu, poisson_approx(n_gaps), 'x', label="Poisson Approx. P_gap")

plt.plot(x_list, poi)
# plt.plot(mu_crit, y_list)

plt.legend()

plt.xlim([0.1, 0.4])

# plt.yscale('log')

plt.show()
