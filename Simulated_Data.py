import numpy as np
import json
import platform

"""
SimulatedData is a file used to access all of the simulated data that has been accumulated.
"""

def import_data(conn_fun, eta, L, Torus, data_type, new=False):

    """
    This function imports the simulated data as a dictionary.

    :param conn_fun: Input the connection function as a string, choices are:  Exponential, Polynomial, Rayleigh, Waxman
    :param eta: Choose the value of eta for the connection function, choices are: 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0
    :param L: Choose the system size
    :param Torus: Boolean input to choose whether you are looking at a line or a torus
    :param data_type: Define the type of data you wish to import, choices are: Overall, GapSize, GapError, IsoLoc
    :return: Simulated data as a dictionary: mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps
    """

    if Torus == True:
        if data_type == "Overall":
            fname = "Simulated_Data/%s/Eta_%s/L_%s/Torus/Overall.txt" % (conn_fun, eta, L)

            data = json.load(open(fname))
            muvalues = list(data.keys())
            mu = np.sort([float(i) for i in muvalues])
            P_fc = np.sort([1 - data[mu][1] / data[mu][0] for mu in muvalues])
            P_iso = np.sort([data[mu][2] / data[mu][0] for mu in muvalues])
            P_gap = np.sort([data[mu][3] / data[mu][0] for mu in muvalues])
            P_split = np.sort([data[mu][4] / data[mu][0] for mu in muvalues])
            mean_degree_temp = np.sort([data[mu][5] / data[mu][0] for mu in muvalues])
            mean_degree = mean_degree_temp[::-1]
            n_iso = np.sort([data[mu][6] / data[mu][0] for mu in muvalues])
            n_gaps = np.sort([data[mu][7] / data[mu][0] for mu in muvalues])

            return mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps

    elif new == True:
        if data_type == "Overall":
            fname = f"Simulated_Data/{conn_fun}/Eta_{eta}/L_{L}/UnionTest/Overall.txt"

            data = json.load(open(fname))
            muvalues = list(data.keys())
            mu = np.sort([float(i) for i in muvalues])
            P_fc = np.sort([1 - data[mu][1] / data[mu][0] for mu in muvalues])
            P_iso = np.sort([data[mu][2] / data[mu][0] for mu in muvalues])
            P_gap = np.sort([data[mu][3] / data[mu][0] for mu in muvalues])
            P_isoUucg = np.sort([data[mu][4] / data[mu][0] for mu in muvalues])
            P_split = np.sort([data[mu][5] / data[mu][0] for mu in muvalues])
            mean_degree_temp = np.sort([data[mu][6] / data[mu][0] for mu in muvalues])
            mean_degree = mean_degree_temp[::-1]
            n_iso = np.sort([data[mu][7] / data[mu][0] for mu in muvalues])
            n_gaps = np.sort([data[mu][8] / data[mu][0] for mu in muvalues])

            return mu, P_fc, P_iso, P_gap, P_isoUucg, P_split, mean_degree, n_iso, n_gaps

    else:
        if data_type == "Overall":
            fname = "Simulated_Data/%s/Eta_%s/L_%s/Line/Overall.txt" % (conn_fun, eta, L)

            data = json.load(open(fname))
            muvalues = list(data.keys())
            mu = np.sort([float(i) for i in muvalues])
            P_fc = np.sort([1 - data[mu][1] / data[mu][0] for mu in muvalues])
            P_iso = np.sort([data[mu][2] / data[mu][0] for mu in muvalues])
            P_gap = np.sort([data[mu][3] / data[mu][0] for mu in muvalues])
            P_split = np.sort([data[mu][4] / data[mu][0] for mu in muvalues])
            mean_degree_temp = np.sort([data[mu][5] / data[mu][0] for mu in muvalues])
            mean_degree = mean_degree_temp[::-1]
            n_iso = np.sort([data[mu][6] / data[mu][0] for mu in muvalues])
            n_gaps = np.sort([data[mu][7] / data[mu][0] for mu in muvalues])

            return mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps

            # data = json.load(open(fname))
            # muvalues = list(data.keys())
            # mu = [float(i) for i in muvalues]
            # P_fc = [1 - data[mu][1] / data[mu][0] for mu in muvalues]
            # P_iso = [data[mu][2] / data[mu][0] for mu in muvalues]
            # P_gap = [data[mu][3] / data[mu][0] for mu in muvalues]
            # P_split = [data[mu][4] / data[mu][0] for mu in muvalues]
            # mean_degree_temp = [data[mu][5] / data[mu][0] for mu in muvalues]
            # mean_degree = mean_degree_temp[::-1]
            # n_iso = [data[mu][6] / data[mu][0] for mu in muvalues]
            # n_gaps = [data[mu][7] / data[mu][0] for mu in muvalues]
            #
            # return mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps

def import_data_bulk(conn_fun, eta, L, Torus, data_type):

    """
    This function imports the simulated data as a dictionary.

    :param conn_fun: Input the connection function as a string, choices are:  Exponential, Polynomial, Rayleigh, Waxman
    :param eta: Choose the value of eta for the connection function, choices are: 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0
    :param L: Choose the system size
    :param Torus: Boolean input to choose whether you are looking at a line or a torus
    :param data_type: Define the type of data you wish to import, choices are: Overall, GapSize, GapError, IsoLoc
    :return: Simulated data as a dictionary: mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps
    """

    if Torus == True:
        if data_type == "Overall":
            fname = "Simulated_Data/%s/Eta_%s/L_%s/Torus/Overall.txt" % (conn_fun, eta, L)

            data = json.load(open(fname))
            muvalues = list(data.keys())
            mu = np.sort([float(i) for i in muvalues])
            P_fc = np.sort([1 - data[mu][1] / data[mu][0] for mu in muvalues])
            P_iso = np.sort([data[mu][2] / data[mu][0] for mu in muvalues])
            P_gap = np.sort([data[mu][3] / data[mu][0] for mu in muvalues])
            P_split = np.sort([data[mu][4] / data[mu][0] for mu in muvalues])
            mean_degree_temp = np.sort([data[mu][5] / data[mu][0] for mu in muvalues])
            mean_degree = mean_degree_temp[::-1]
            n_iso = np.sort([data[mu][6] / data[mu][0] for mu in muvalues])
            n_gaps = np.sort([data[mu][7] / data[mu][0] for mu in muvalues])

            return mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps

    else:
        if data_type == "Overall":
            fname = "Simulated_Data/Remove_Boundary_Effects/%s/Eta_%s/L_%s/Line/Overall.txt" % (conn_fun, eta, L)

            # data = json.load(open(fname))
            # muvalues = list(data.keys())
            # mu = np.sort([float(i) for i in muvalues])
            # P_fc = np.sort([1 - data[mu][1] / data[mu][0] for mu in muvalues])
            # P_iso = np.sort([data[mu][2] / data[mu][0] for mu in muvalues])
            # P_gap = np.sort([data[mu][3] / data[mu][0] for mu in muvalues])
            # P_split = np.sort([data[mu][4] / data[mu][0] for mu in muvalues])
            # mean_degree_temp = np.sort([data[mu][5] / data[mu][0] for mu in muvalues])
            # mean_degree = mean_degree_temp[::-1]
            # n_iso = np.sort([data[mu][6] / data[mu][0] for mu in muvalues])
            # n_gaps = np.sort([data[mu][7] / data[mu][0] for mu in muvalues])
            #
            # return mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps

            data = json.load(open(fname))
            muvalues = list(data.keys())
            mu = [float(i) for i in muvalues]
            P_fc = [1 - data[mu][1] / data[mu][0] for mu in muvalues]
            P_iso = [data[mu][2] / data[mu][0] for mu in muvalues]
            P_gap = [data[mu][3] / data[mu][0] for mu in muvalues]
            P_split = [data[mu][4] / data[mu][0] for mu in muvalues]
            mean_degree_temp = [data[mu][5] / data[mu][0] for mu in muvalues]
            mean_degree = mean_degree_temp[::-1]
            n_iso = [data[mu][6] / data[mu][0] for mu in muvalues]
            n_gaps = [data[mu][7] / data[mu][0] for mu in muvalues]

            return mu, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps
