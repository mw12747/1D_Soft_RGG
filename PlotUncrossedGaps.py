import numpy as np
import matplotlib.pyplot as plt
import json
from Simulated_Data import import_data

########################################################################################################################

# def f(mu, w):
#     return np.exp(-mu * w) / (1 - np.exp(-1/mu))**2 * \
#            ( 1/(mu**2) * (1 - np.exp(-1/mu) / 2)**2 +
#                                   2 / mu * (1 - np.exp(-1/mu) / 2) * (1 - np.exp(-1/mu)) +
#                                   (1 - np.exp(-1/mu))**2)
#
# w_list = np.arange(0, 30, 0.001)
# E_edges = [f(mu, w) for w in w_list]
# E_edges_og = [np.exp(-mu * w)*(1/mu**2 + 2/mu + 1) for w in w_list]


########################################################################################################################

def EdgesAcross_NonIso_Sims(rt, mu, eta):
    fname = "Simulated_Data/EdgesAcrossGaps/Non_Iso/Edges_Across_NonIso_rt_%s_mu_%s_eta_%s.txt" % (rt, mu, eta)

    with open(fname, "r") as my_file:
        data = json.load(my_file)

    w_values_temp = list(data.keys())
    w_sims_list = np.sort([float(i) for i in w_values_temp])
    w_values = ['%s' %i for i in w_sims_list]

    n_edges_across = [data[w][1] / data[w][0] for w in w_values]
    n_edges_across_2 = [data[w][2]/ data[w][0] for w in w_values]
    n_edges_across_3 = [data[w][3] / data[w][0] for w in w_values]
    n_edges_across_4 = [data[w][4] / data[w][0] for w in w_values]
    n_edges_across_5 = [data[w][5] / data[w][0] for w in w_values]
    n_edges_across_6 = [data[w][6] / data[w][0] for w in w_values]
    n_edges_across_7 = [data[w][7] / data[w][0] for w in w_values]
    n_edges_across_8 = [data[w][8] / data[w][0] for w in w_values]
    n_edges_across_9 = [data[w][9] / data[w][0] for w in w_values]
    n_edges_across_10 = [data[w][10] / data[w][0] for w in w_values]
    n_edges_across_11 = [data[w][11] / data[w][0] for w in w_values]
    n_edges_across_12 = [data[w][12] / data[w][0] for w in w_values]

    p_NC = [data[w][13] / data[w][0] for w in w_values]

    return w_sims_list, n_edges_across, n_edges_across_2, n_edges_across_3, n_edges_across_4, n_edges_across_5, n_edges_across_6, n_edges_across_7, n_edges_across_8, n_edges_across_9, n_edges_across_10, n_edges_across_11, n_edges_across_12, p_NC

def N_ucg(rt, mu, eta, L):
    summand = [EdgesAcross_NonIso_Sims(rt, mu, eta)[13][i] * np.exp(-EdgesAcross_NonIso_Sims(rt, mu, eta)[0][i]) * 0.1 for i in range(len(EdgesAcross_NonIso_Sims(rt, mu, eta)[13]))]
    return L * sum(summand)


########################################################################################################################



def EdgesAcross_No_Condition(rt, mu):
    fname = "Simulated_Data/EdgesAcrossGaps/No_Condition/EdgesAcrossGaps_rt_%s.txt" % rt

    with open(fname, "r") as my_file:
        data = json.load(my_file)

    mu = str(mu)

    mu_list_no_condition = [float(i) for i in data.keys()]

    w_sims_list = data[mu][1]
    n_edges_across = [data[mu][2][i]/data[mu][0] for i in range(len(data[mu][2]))]

    return(w_sims_list, n_edges_across, mu_list_no_condition)

def N_ucg_no_condition(rt, mu, L):
    summand = [EdgesAcross_No_Condition(rt, mu)[1][i] * np.exp(-EdgesAcross_No_Condition(rt, mu)[0][i]) * 0.1 for i in range(len(EdgesAcross_No_Condition(rt, mu)[0]))]
    return L * sum(summand)

########################################################################################################################

conn_fun = "Waxman"
rt = 1
eta = 1.0
L = 1000
Torus = False
data_type = "Overall"

########################################################################################################################

mu_list_2, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps = import_data(conn_fun, eta, L, Torus, data_type)

########################################################################################################################

# mu_list_no_condition = EdgesAcross_No_Condition(rt, 0.2)[2]
# mu_list = [0.1, 0.15, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.35, 0.4, 0.45]
mu_list = [round(x, 4) for x in np.arange(0.15, 0.44, 0.01)] + [round(x, 4) for x in np.arange(0.45, 0.95, 0.05)]

print(mu_list)

########################################################################################################################

eta = 1

P_ucg_sim = [1 - np.exp(-N_ucg(rt, mu, eta, L)) for mu in mu_list]
N_ucg_sim = [N_ucg(rt, mu, eta, L) for mu in mu_list]
# P_ucg_sim_no_condition = [1 - np.exp(-N_ucg_no_condition(rt, x, L)) for x in mu_list_no_condition]

# print(P_ucg_sim_no_condition)

plt.figure()

plt.title("Simulations for L = %s, %s conn. fun, $\eta$ =  %s" % (L, conn_fun, eta))

plt.xlabel("$\mu$")
plt.ylabel("Prob. of event")

plt.plot(mu_list, P_ucg_sim, 'x', label="Sims of sims P_ucg")
plt.plot(mu_list_2, P_gap, label="P_gap")
# plt.plot(mu_list_no_condition, P_ucg_sim_no_condition, label='P_gap_no_condition')

plt.xlim([0.09, 1])

plt.legend()


plt.figure()

plt.title("Simulations for L = %s, %s conn. fun, $\eta$ =  %s" % (L, conn_fun, eta))

plt.xlabel("$\mu$")
plt.ylabel("Expected Number")

plt.plot(mu_list, N_ucg_sim, 'x', label="Sims of sims N_ucg")
plt.plot(mu_list_2, n_gaps, label="N_gap")

plt.show()
