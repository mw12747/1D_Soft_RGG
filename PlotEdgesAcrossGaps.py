import numpy as np
import matplotlib.pyplot as plt
import json
from math import factorial

def E_cond(mu, w):
    return np.exp(-mu * w) / (1 - np.exp(-1/mu))**2 * \
           ( 1/(mu**2) * (1 - np.exp(-1/mu) / 2)**2 +
                                  2 / mu * (1 - np.exp(-1/mu) / 2) * (1 - np.exp(-1/mu)) +
                                  (1 - np.exp(-1/mu))**2)

def e_h(mu):
    return np.exp(-1/mu)

def p_z(mu):
    return (1 - e_h(mu))**2



# def Var_cond(mu, w):
#     return E_cond(mu, w) + 1/(1 - np.exp(-1/mu))**4 * \
#            (((1 - np.exp(-1/mu)) / (2*mu) + np.exp(-1/mu) / (3*mu))
#              * ((1 - np.exp(-1/mu))*np.exp(-mu*w)/mu + np.exp(-1/mu)*np.exp(-mu*w)/(2*mu))**2 +
#              ((1 - np.exp(-1 / mu)) / mu + np.exp(-1 / mu) / (2 * mu))**2
#              * ((1-np.exp(-1/mu))*np.exp(-mu*w)/(2*mu)+np.exp(-1/mu)*np.exp(-mu*w)/(3*mu)) +
#             2 * (((1-np.exp(-1/mu))/mu+np.exp(-1/mu)/(2 * mu)) *
#                   ((1-np.exp(-1/mu))*np.exp(-2*mu*w)/(2*mu)+np.exp(-1/mu)*np.exp(-2*mu*w)/(3*mu))
#                   + ((1-np.exp(-1/mu))/(2*mu)+np.exp(-1/mu)/(3 * mu))*
#                   ((1-np.exp(-1/mu))*np.exp(-2*mu*w)/mu+np.exp(-1/mu)*np.exp(-2*mu*w)/(2*mu)))
#             ) - np.exp(-4 * mu * w)

def Var_cond(mu, w):
    return E_cond(mu, w) + np.exp(-2 * mu * w) * (1 / (p_z(mu) ** (3/2) * mu ** 3) * (1 - e_h(mu) / 2) ** 2 * (1 - e_h(mu) / 2) +
                                                  1 / (p_z(mu) ** 2 * mu ** 2) * (1 - e_h(mu) / 2) * (1 - e_h(mu) / 3)
                                                                - 1 / p_z(mu)**4)
            # + np.exp(-2 * mu * w) * (p_z(mu) - 1) / p_z(mu) ** 2 * (2 * (1 - e_h(mu) / 2) ** 2 / mu ** 2 +
            #                                                         2 * p_z(mu) ** (1/2) * (1 - e_h(mu) / 2) ** 3 / mu ** 3 +
            #                                                         (1 - e_h(mu) / 2) ** 4 / mu ** 4)

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

def EdgesAcross_No_Condition(rt, mu):
    fname = "Simulated_Data/EdgesAcrossGaps/No_Condition/EdgesAcrossGaps_rt_%s.txt" % rt

    with open(fname, "r") as my_file:
        data = json.load(my_file)

    mu = str(mu)

    w_sims_list = data[mu][1]
    n_edges_across = [data[mu][2][i]/data[mu][0] for i in range(len(data[mu][2]))]

    return(w_sims_list, n_edges_across)

rt = 1
mu = 0.7
eta = 1
w_list = np.arange(0, 30, 0.001)


E_edges = [E_cond(mu, w) for w in w_list]
E_edges_og = [np.exp(-mu * w)*(1/mu**2 + 2/mu + 1) for w in w_list]


Var_edges = [Var_cond(mu, w) for w in w_list]

Var_edges_og = [1 / mu**2 * (np.exp(-mu * w) + np.exp(-2 * mu * w)) +
                1 / mu**3 * np.exp(-2 * mu * w) + 2 / mu * np.exp(-mu * w) +
                np.exp(-mu * w) - np.exp(-2 * mu * w) for w in w_list]

########################################################################################################################

w_sims_list, n_edges_across, n_edges_across_2, n_edges_across_3, n_edges_across_4, n_edges_across_5, n_edges_across_6, n_edges_across_7, n_edges_across_8, n_edges_across_9, n_edges_across_10, n_edges_across_11, n_edges_across_12, p_NC= EdgesAcross_NonIso_Sims(rt, mu, eta)

# w_sims_list_no_con, n_edges_across_no_con = EdgesAcross_No_Condition(rt, mu)

poisson_approx_sims = [np.exp(-n_edges_across[i]) for i in range(len(n_edges_across))]
poisson_approx = [np.exp(-E_edges[i]) for i in range(len(E_edges))]

var_sim = [n_edges_across_2[w] - n_edges_across[w] ** 2 for w in range(len(n_edges_across))]

########################################################################################################################

approx = [1 - (n_edges_across[i]
               - 1/2 * (n_edges_across_2[i] - n_edges_across[i])
               + 1/factorial(3) * (n_edges_across_3[i] - 3*n_edges_across_2[i] + 2*n_edges_across[i])
               - 1/factorial(4) * (n_edges_across_4[i] - 6*n_edges_across_3[i] + 11*n_edges_across_2[i] - 6*n_edges_across[i])
               + 1/factorial(5) * (n_edges_across_5[i] - 10*n_edges_across_4[i] + 35*n_edges_across_3[i] - 50*n_edges_across_2[i] +
                                   24*n_edges_across[i])
               - 1/factorial(6) * (n_edges_across_6[i] - 15*n_edges_across_5[i] + 85*n_edges_across_4[i] - 225*n_edges_across_3[i] +
                                   274*n_edges_across_2[i] - 120*n_edges_across[i])
               + 1/factorial(7) * (n_edges_across_7[i] - 21*n_edges_across_6[i] + 175*n_edges_across_5[i] - 735*n_edges_across_4[i] +
                                   1624*n_edges_across_3[i] - 1764*n_edges_across_2[i] + 720*n_edges_across[i])
               - 1/factorial(8) * (n_edges_across_8[i] - 28*n_edges_across_7[i] + 322*n_edges_across_6[i] - 1960*n_edges_across_5[i] +
                                   6769*n_edges_across_4[i] - 13132*n_edges_across_3[i] + 13068*n_edges_across_2[i] - 5040*n_edges_across[i])
               + 1/factorial(9) * (n_edges_across_9[i] - 36*n_edges_across_8[i] + 546*n_edges_across_7[i] - 4536*n_edges_across_6[i] +
                                   22449*n_edges_across_5[i] - 67284*n_edges_across_4[i] + 118124*n_edges_across_3[i] -
                                   109584*n_edges_across_2[i] + 40320*n_edges_across[i])
               - 1/factorial(10) * (n_edges_across_10[i] - 45*n_edges_across_9[i] + 870*n_edges_across_8[i] - 9450*n_edges_across_7[i] +
                                    63273*n_edges_across_6[i] - 269325*n_edges_across_5[i] + 723680*n_edges_across_4[i] -
                                    1172700*n_edges_across_3[i] + 1026576*n_edges_across_2[i] - 362880*n_edges_across[i])
               + 1/factorial(11) * (n_edges_across_11[i] - 55*n_edges_across_10[i] + 1320*n_edges_across_9[i] - 18150*n_edges_across_8[i] +
                                    157773*n_edges_across_7[i] - 902055*n_edges_across_6[i] + 3416930*n_edges_across_5[i] -
                                    8409500*n_edges_across_4[i] + 12753576*n_edges_across_3[i] - 10628640*n_edges_across_2[i] +
                                    3628800*n_edges_across[i])
               - 1/factorial(12) * (n_edges_across_12[i] - 66*n_edges_across_11[i] + 1925*n_edges_across_10[i] - 32670*n_edges_across_9[i] +
                                    357423*n_edges_across_8[i] - 2637558*n_edges_across_7[i] + 13339535*n_edges_across_6[i] -
                                    45995730*n_edges_across_5[i] + 105258076*n_edges_across_4[i] - 150917976*n_edges_across_3[i] +
                                    120543840*n_edges_across_2[i] - 39916800*n_edges_across[i])
               )
          for i in range(len(w_sims_list))]


########################################################################################################################

Markov = [max(0, 1 - i) for i in n_edges_across]
Chebyshev = [min(1, var_sim[i] / n_edges_across[i]**2) for i in range(len(var_sim)) if n_edges_across[i] != 0]

Markov_Analytic = [max(0, 1 - i) for i in E_edges]

Chebyshev_Analytic = [min(1, i / j**2) for i, j in zip(E_edges, Var_edges)]

########################################################################################################################

plt.figure()

plt.plot(w_sims_list, n_edges_across, 'o', label='Sims - Conditional')
plt.plot(w_list, E_edges, '-', label='Analytic - Conditional')

# plt.plot(w_sims_list_no_con, n_edges_across_no_con, 'x', label='Sims - No Condition')
plt.plot(w_list, E_edges_og, '-', label='Analytic - No Condition')

plt.title('$\mu$ = %s' %mu)

plt.xlim([-0.1, 15])

plt.xlabel('Gap size')
plt.ylabel('No. of conns across')

plt.legend()

# plt.figure()
#
# plt.plot(w_sims_list_no_con, n_edges_across_no_con, 'x', label='Sims E[X] No Condition')
# plt.plot(w_list, E_edges_og, '-', label='Analytic Original')
#
# plt.title('No condition, $\mu$ = %s' %mu)
#
# plt.xlabel('Gap size')
# plt.ylabel('No. of conns across')
#
# plt.legend()
#

########################################################################################################################

plt.figure()

plt.plot(w_sims_list, p_NC, 'x', label='Sims $P_{nc}$')
# plt.plot(w_sims_list, poisson_approx_sims, 'o', label='Poisson Approx Sims')
plt.plot(w_list, poisson_approx, '-', label='Poisson Approx Analytic')

plt.plot(w_sims_list, Markov, '+', label='Sims Markov bound')
plt.plot(w_sims_list[0:len(Chebyshev)], Chebyshev, '+', label='Sims Chebyshev bound')

plt.plot(w_list, Markov_Analytic, '-', label='Analytic Markov bound')
plt.plot(w_list, Chebyshev_Analytic, '-', label='Analytic Chebyshev bound')

plt.plot(w_sims_list, approx, '.', label='Sims approx')

plt.title('12 Moments, $\mu$ = %s' %mu)

plt.xlabel('Gap size')
plt.ylabel('No. of conns across')

plt.ylim([-0.1, 1.1])

plt.legend()

########################################################################################################################

plt.figure()

plt.plot(w_list, Var_edges_og, label="Variance og Analytic")
plt.plot(w_list, Var_edges, label="Variance Analytic")
plt.plot(w_sims_list, var_sim, 'x', label="Variance Sims")

plt.legend()

plt.show()