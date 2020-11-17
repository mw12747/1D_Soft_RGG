import matplotlib.pyplot as plt
import numpy as np
import json

def import_data(conn_fun, eta, L):
    fname = f"Simulated_Data/{conn_fun}/Eta_{eta}/L_{L}/Line/Overall.txt"

    data = json.load(open(fname))
    muvalues = list(data.keys())
    mu = [float(i) for i in muvalues]
    P_fc = np.sort([1 - data[mu][1] / data[mu][0] for mu in muvalues])
    P_iso = [data[mu][2] / data[mu][0] for mu in muvalues]
    P_gap = np.sort([data[mu][3] / data[mu][0] for mu in muvalues])
    P_split = np.sort([data[mu][4] / data[mu][0] for mu in muvalues])
    mean_degree_temp = np.sort([data[mu][5] / data[mu][0] for mu in muvalues])
    mean_degree = mean_degree_temp[::-1]
    n_iso = [data[mu][6] / data[mu][0] for mu in muvalues]
    n_gaps = np.sort([data[mu][7] / data[mu][0] for mu in muvalues])

    return mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps


conn_fun = "Rayleigh"

mu_1000, P_fc_1000, P_iso_1000, P_gap_1000, P_split_1000, mean_degree_1000, n_iso_1000, n_gaps_1000 \
    = import_data(conn_fun, 1000.0, 1000)

mu_1001, P_fc_1001, P_iso_1001, P_gap_1001, P_split_1001, mean_degree_1001, n_iso_1001, n_gaps_1001 \
    = import_data(conn_fun, 1001.0, 1001)

mu_10000, P_fc_10000, P_iso_10000, P_gap_10000, P_split_10000, mean_degree_10000, n_iso_10000, n_gaps_10000 \
    = import_data(conn_fun, 10000.0, 10000)

x_list = np.arange(0.1, 6, 0.001)
y_list = [1/x for x in x_list]

y_list_poi = [1 - np.exp(-x) for x in y_list]

plt.figure()

plt.rcParams.update({'font.size': 12})

plt.grid(color='whitesmoke', which='major', linestyle='-', linewidth='1')
plt.minorticks_on()
plt.grid(color='whitesmoke', which='minor', linestyle='--', linewidth='1')

plt.plot(x_list, y_list_poi, color='k', label=r"$f(\beta) = 1 - e^{-\beta}$", linewidth=1)

plt.plot(mu_1000, P_iso_1000, '.', ms=7, color='silver', label=r"$\bar{L}$ = 1000")
# plt.plot(mu_1001, P_iso_1001, 'o', label=r"$\bar{L}$ = 1000, Rayleigh")
plt.plot(mu_10000, P_iso_10000, '.', ms=7, color='gray', label=r"$\bar{L}$ = 10000")



# plt.title(f"Waxman conn. fun, $\eta$ =  1")

plt.xlabel(r"$\beta$")
plt.ylabel(r"$\mathbb{P}(N_{iso} \geq 1)$")

# plt.xlim([3.4, 20.1])
# plt.ylim([-0.03, 1.03])

plt.legend()
plt.show()

