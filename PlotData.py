from Simulated_Data import import_data
import matplotlib.pyplot as plt
import numpy as np
from math import log, sqrt, pi, gamma
from scipy.special import gamma, erf, expi
from scipy import integrate
import json

def poisson_approx(*kargs):
    if len(kargs) == 1:
        return [1 - np.exp(-x) for x in kargs[0]]
    else:
        return [1 - np.exp(-x - y) for (x, y) in zip(kargs[0], kargs[1])]


def iso_crit(eta, L):
    return (2 * gamma((eta+1) / eta) / log(L))**eta


def gap_crit(eta, L, c):
    if eta == 1:
        ret = 2*log(log(L))**(1/2)/c * log(L)/log(log(L))
    else:
        ret = 2*log(log(L))**(1 /4)/c * log(L) / log(log(L)**(1/eta))
    return ret


data_type = "Overall"

conn_fun = "Rayleigh"
eta = 2.0
Torus = False
new = True
pois_approx = False
md = True
crit = True

x_list = np.arange(0.0001, 0.8, 0.0001)
y_list = np.arange(-0.3, 1.3, 0.0001)

for L in [10000]:
    plt.figure()

    plt.rcParams.update({'font.size': 12})

    plt.grid(color='whitesmoke', which='major', linestyle='-', linewidth='1')
    plt.minorticks_on()
    plt.grid(color='whitesmoke', which='minor', linestyle='--', linewidth='1')

    ### Import Data ####################################################################################################

    if new == True:
        mu, P_fc, P_iso, P_gap, P_isoUucg, P_split, mean_degree, n_iso, n_gaps = import_data(conn_fun, eta, L, Torus,
                                                                                         data_type, new=True)
    else:
        mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps = import_data(conn_fun, eta, L, Torus, data_type)

    ####################################################################################################################
    y_test = [1 - np.exp(-L * np.exp(- 1 / x**2 + 1 - np.exp(-1 / x))) for x in x_list]
    ### Plots ##########################################################################################

    if conn_fun == "Waxman":
        mu = [mu[i] for i in range(len(mu)) if i%2 == 0]
        mean_degree = [mean_degree[i] for i in range(len(mean_degree)) if i%2==0]
        P_fc = [P_fc[i] for i in range(len(P_fc)) if i%2 == 0]
        P_iso = [P_iso[i] for i in range(len(P_iso)) if i%2 == 0]
        P_gap = [P_gap[i] for i in range(len(P_gap)) if i%2 == 0]
        P_isoUucg = [P_isoUucg[i] for i in range(len(P_isoUucg)) if i%2 == 0]

    if md:
        x_axis = mean_degree
        x_label = "Mean Degree"
        x_lim = [3.4, 20.1]
    else:
        x_axis = mu
        x_label = r"$\mu$"
        x_lim = [0.1, 0.6]

    if pois_approx == True:
        if md == True:
            crit_iso = [log(L) for x in y_list]
            plt.plot(crit_iso, y_list, '--', color='silver', label="$\ln L$")

            x_list_poi = np.arange(3.4, 20.1, 0.001)
            y_list_poi = [1 - np.exp(-L * np.exp(-x)) for x in x_list_poi]

            plt.plot(x_list_poi, y_list_poi, '-', color='k', label="Poisson Approximation")
            plt.plot(mean_degree, P_iso, '.', ms=8,  color='royalblue', label="$P_{iso}$ (Waxman)")

            # plt.title(f"Mean no. of nodes = {L}, {conn_fun} conn. fun, $\eta$ =  {eta}")

            plt.xlabel("Mean Degree")
            plt.ylabel("Prob. of event")

            plt.xlim([3.4, 20.1])
            plt.ylim([-0.03, 1.03])

            plt.legend()

        else:
            if conn_fun == "Waxman":
                y_poi_temp = [
                    integrate.quad(lambda y: np.exp(-1 / x * (2 - np.exp(-x * y) - np.exp(-x * (L - y)))), 0, L)[0]
                    for x in x_list]
            else:
                y_poi_temp = [integrate.quad(lambda y: np.exp(-1 / 2 * sqrt(pi / x) * (erf(L - y) + erf(y))), 0, L)[0]
                              for x in x_list]
            y_poi = [1 - np.exp(-y) for y in y_poi_temp]

            plt.plot(x_list, y_poi, '-', color='k', label='Poisson Approximation')
            plt.plot(mu, P_iso, '.', ms=7, color='gray', label="$P_{iso}$")

            plt.title(f"Mean no. of nodes = {L}, {conn_fun} conn. fun, $\eta$ =  {eta}")

            plt.xlabel("$\mu$")
            plt.ylabel("Prob. of event")

            plt.xlim([0, 0.8])
            plt.ylim([-0.03, 1.03])

            plt.legend()

    elif pois_approx == False and new == True:

        plt.plot(x_axis, P_fc, 'o', ms=6, color='k', label=r"$P_{dis}$")
        plt.plot(x_axis, P_gap, 's', ms=4, color='forestgreen', label="$P_{ucg}$")
        plt.plot(x_axis, P_iso, '.', ms=8,  color='royalblue', label="$P_{iso}$")
        plt.plot(x_axis, P_isoUucg, 'p', ms=4, color='gainsboro', label="$P_{iso \cup ucg}$")

        # L_temp = '{:,}'.format(L).replace(',', ' ')

        plt.title(f"{conn_fun} conn. fun., $\eta$ =  {eta} \n Mean no. of nodes = {L}")

        plt.xlabel(x_label)
        plt.ylabel("Prob. of event")

        plt.xlim(x_lim)
        plt.ylim([-0.03, 1.03])

        plt.legend()

    else:

        plt.plot(x_axis, P_fc, 'o', ms=6, color='k', label=r"$\bar{P}_{fc}$")
        plt.plot(x_axis, P_iso, '.', ms=7, color='dimgrey', label="$P_{iso}$")
        plt.plot(x_axis, P_gap, 's', ms=5, color='darkgrey', label="$P_{ucg}$")

        # plt.title(f"Mean no. of nodes = {L}, {conn_fun} conn. fun, $\eta$ =  {eta}")

        plt.xlabel(x_label)
        plt.ylabel("Prob. of event")

        plt.xlim(x_lim)
        plt.ylim([-0.03, 1.03])

        plt.legend()

    # if crit == True:
    #     crit_iso = [log(L) for x in y_list]
    #     # crit_gap = [gap_crit(eta, L) for x in y_list]
    #
    #     plt.plot(crit_iso, y_list, color='k', label="$\ln L$")
    #     # plt.plot(crit_gap, y_list, color='darkgrey', label="Gap-Crit")
    #
    #     plt.title(f"Mean no. of nodes = {L}, {conn_fun} conn. fun, $\eta$ =  {eta}")
    #
    #     plt.xlabel("Mean Degree")
    #     plt.ylabel("Prob. of event")
    #
    #     plt.xlim([3.4, 20.1])
    #     plt.ylim([-0.03, 1.03])
    #
    #     plt.legend()

    ####################################################################################################################

conn_fun = "Rayleigh"
eta = 2.0
for L in [1000]:
    new = True
    ### Import Data ####################################################################################################

    if new == True:
        mu, P_fc, P_iso, P_gap, P_isoUucg, P_split, mean_degree, n_iso, n_gaps = import_data(conn_fun, eta, L, Torus,
                                                                                         data_type, new=True)
    else:
        mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps = import_data(conn_fun, eta, L, Torus, data_type)

    ####################################################################################################################
    y_test = [1 - np.exp(-L * np.exp(- 1 / x**2 + 1 - np.exp(-1 / x))) for x in x_list]
    ### Plots ##########################################################################################

    if pois_approx == True:
        if md == True:
            plt.plot(mean_degree, P_iso, '.', ms=8,  color='firebrick', label="$P_{iso}$ (Rayleigh)")

            # plt.title(f"Mean no. of nodes = {L}")

            plt.xlabel("Mean Degree")
            plt.ylabel("Prob. of event")

            plt.xlim([3.4, 20.1])
            plt.ylim([-0.03, 1.03])

            plt.legend()

        else:
            break

    else:
        break

    ####################################################################################################################
# em = 0.57721566490153286060651209008240243104215933593992
# a = 0.5
#  L = 1000
# x_list_md = np.arange(0.0001, 20, 0.0001)
# y_list_test = [1 - np.exp(-L * np.exp(-x/2 * (log(x/2) - expi(-x/2) + em) - 1 + np.exp(-x/2))) for x in x_list_md]
# y_list_test_new = [1 - np.exp(-L*np.exp(-(1 + log(L) * (2/x)**(-1/eta) * (log(log(L)**a))**(1/eta) + (log(L)**a)))) for x in x_list_md]
#
# # plt.plot(x_list_md, y_list_test_new, label='UCG upper bound')
#
# test_data = {"0.15": [1100, 1096, 4, 3, 4, 1, 14556.106294611398, 4, 3], "0.2": [1100, 1018, 66, 28, 80, 56, 10936.618227901497, 67, 31], "0.25": [1100, 614, 332, 235, 476, 299, 8757.554572945917, 430, 295], "0.3": [1101, 91, 806, 717, 1002, 755, 7312.237977933036, 1523, 1409], "0.35": [1100, 2, 1049, 1031, 1097, 1032, 6254.014732023808, 3711, 4281], "0.5": [2, 1, 1, 1, 1, 1, 8.250592885375495, 5, 4], "0.6": [5, 0, 5, 5, 5, 5, 18.98993667714778, 14, 21]}
# x_test_temp = list(test_data.keys())
# x_test = [float(i) for i in x_test_temp]
# y_test = [test_data[x][3] / test_data[x][0] for x in x_test_temp]
# plt.plot(x_test, y_test, 'x', label='new data')
# plt.figure()
# plt.plot(mu, P_gap, 'x')
# plt.plot(x_list, y_test)
# crit_iso = [2*log(10000) for x in y_list]
# plt.plot(crit_iso, y_list, '--', color='silver', label="$\ln L$")

gap_test = [gap_crit(2, 10000, 1.043) for y in y_list]
plt.plot(gap_test, y_list, label='GapCritTest')

plt.show()
