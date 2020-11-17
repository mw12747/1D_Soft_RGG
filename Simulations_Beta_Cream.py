# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 18:00 2019

@author: mw12747
"""
# Code to simulate a 1D Soft RGG on a PPP

import networkx as nx
import json
from multiprocessing import Pool
from numpy import exp, random, arange, var
import os
from math import log
import sys
import time
import platform


## Define initial functions ##

def conn_fun(x, y, beta, L):
    """
    Defines the connection function that we wish to look at
    :param x: location of first node
    :param y: location of second node
    :param mu: inverse of the "connection range"
    :return: Returns the probability of nodes x and y being connected
    """

    rayleigh = exp(- 2 * (abs(x - y)) / log(beta * L))
    return rayleigh


def PPP(rt, L):
    """
    1-dimensional Poisson point process of rate rt
    :param rt: Gives the rate of the PPP
    :param L: length of line segment
    :return:
    [0] - Returns a list of the locations of nodes for a 1D PPP
    [1] - no. of nodes
    """
    L = float(L)
    n_nodes = random.poisson(rt * L)  # Poisson number of points
    node_locations = sorted(random.uniform(0, L, n_nodes))  # Uniformly distributed along the line

    return node_locations, n_nodes


###############################

## Create node and edge lists ##

def edges(rt, L, mu, eta):
    """
    Gives an edge and node list for our 1D soft RGG
    :param rt: Input the rate of the PPP
    :param L: Give the length of the segment we want to look at
    :param mu: Input the mu for the connection function
    :return: Outputs a list of edges and a list of the nodes in our graph along with node positions
    """

    pos_list = PPP(rt, L)[0]
    n_nodes = len(pos_list)
    node_list = [i for i in range(n_nodes)]
    edge_list = []

    for i in range(n_nodes):
        for j in range(i + 1, n_nodes, 1):
            if conn_fun(pos_list[i], pos_list[j], mu, eta) > 10 ** (-5):  # Checks that nodes aren't too far apart, if they are then break
                if random.uniform(0, 1) < conn_fun(pos_list[i], pos_list[j], mu, eta):
                    edge_list.append((i, j))
            else:
                break

    return edge_list, node_list, pos_list


################################

## Generate graph ##

def graph(edge_list, node_list, pos_list):
    """
    Creates a 1-dimensional soft RGG
    :param edge_list: list of edges
    :param node_list: list of node locations
    :param pos_list: list of node positions
    :return: nx.Graph object
    """
    G = nx.Graph()
    G.add_edges_from(edge_list)
    G.add_nodes_from(node_list)
    for n in range(len(node_list)):
        G.nodes[n]['pos'] = pos_list[n]

    return G


####################

## Define metrics ##

def metrics(G):
    """
    Gives the relevant quantities we want to test on G
    :param G: Input is a networkx graph object
    :return:
    [0] - is the graph fully connected
    [1] - are there isolated nodes
    [2] - is there a gap
    [3] - is there a split
    [4] - mean degree of the network
    [5] - number of isolated nodes
    [6] - number of gaps
    [7] - locations of isolated nodes
    [8] - locations of gaps
    [9] - size of gaps
    [10] - gaps error
    [11] - no. of clusters (w/o isolated nodes)
    [12] - cluster size (w/o isolated nodes)
    """
    P_fc = 0  # Give starting values for all the parameters we wish to return
    P_iso = 0
    P_gap = 0
    P_split = 0
    mean_degree = sum(list(dict(G.degree()).values())) / len(G)
                                                # Overall degree of the network divided by the number of nodes
    n_iso = 0
    n_gaps = 0
    loc_iso = []
    loc_gaps =[]
    size_gaps = []
    gaps_error = []
    n_clusters = 0
    # cluster_size = 0

    if nx.is_connected(G):
        P_fc = 1  # Test full connectivity, if true can stop
        n_clusters += 1
    else:
        if list(nx.isolates(G)):
            P_iso = 1  # Test for isolated nodes

        isolated = [i for i in nx.isolates(G)]
        n_iso += len(isolated)  # Record number of isolated nodes
        if P_iso == 1:
            loc_iso = [G.nodes[i]['pos'] for i in isolated]  # Record location of isolated nodes

        G.remove_nodes_from(isolated)  # need to remove isolated nodes for later tests

        clusters = list(nx.connected_components(G))  # Gives a list of all of the connected clusters in our graph G
        n_clusters += len(clusters)

        cc = sorted([[min(clusters[i]), max(clusters[i])] for i in range(len(clusters))])

        k = len(cc)  # gives the number of clusters we are looking at

        if k > 1:

            i = 0  # gives a counter value
            while i < k - 1:
                if cc[i + 1][0] - cc[i][1] == 1:
                    P_gap = 1
                    n_gaps += 1  # Record number of gaps
                    loc_gaps += [G.nodes[cc[i][1]]['pos']]  # Record locations of gaps
                    size_gaps += [G.nodes[cc[i + 1][0]]['pos'] - G.nodes[cc[i][1]]['pos']]  # Record size of gaps
                else:
                    P_split = 1  # If no iso and no gap then must be split

                # elif cc[i + 1][0] < cc[i][1]: # Test for splits
                #     if cc[i + 1][1] > cc[i][1]:
                #         P_split = 1
                #     else:  # If split is within another cluster need to remove
                #         P_split = 1
                #         cc.remove(cc[i+1])
                #         k -= 1
                #         i -= 1

                i += 1

        gaps_error += [n_gaps]

    return P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps, loc_iso, loc_gaps, size_gaps, gaps_error, n_clusters


####################

def metrics_Torus(G):
    """
    Gives the relevant quantities we want to test on G
    :param G: Input is a networkx graph object
    :return:
    [0] - is the graph fully connected
    [1] - are there isolated nodes
    [2] - is there a gap
    [3] - is there a split
    [4] - mean degree of the network
    [5] - number of isolated nodes
    [6] - number of gaps
    [7] - locations of isolated nodes
    [8] - locations of gaps
    [9] - size of gaps
    [10] - gaps error
    [11] - no. of clusters (w/o isolated nodes)
    [12] - cluster size (w/o isolated nodes)
    """
    P_fc = 0  # Give starting values for all the parameters we wish to return
    P_iso = 0
    P_gap = 0
    P_split = 0
    mean_degree = sum(list(dict(G.degree()).values())) / len(G)  # Overall degree of the network divided by the number of nodes
    n_iso = 0
    n_gaps = 0
    loc_iso = []
    loc_gaps =[]
    size_gaps = []
    gaps_error = []
    n_clusters = 0
    # cluster_size = 0

    if nx.is_connected(G):
        P_fc = 1  # Test full connectivity, if true can stop
        n_clusters += 1
    else:
        if list(nx.isolates(G)):
            P_iso = 1  # Test for isolated nodes

        isolated = [i for i in nx.isolates(G)]
        n_iso += len(isolated)  # Record number of isolated nodes
        if P_iso == 1:
            loc_iso = [G.nodes[i]['pos'] for i in isolated]  # Record location of isolated nodes

        G.remove_nodes_from(isolated)  # need to remove isolated nodes for later tests

        clusters = list(nx.connected_components(G))  # Gives a list of all of the connected clusters in our graph G
        n_clusters += len(clusters)

        cc = sorted([[min(clusters[i]), max(clusters[i])] for i in range(len(clusters))])

        k = len(cc)  # gives the number of clusters we are looking at

        if k > 1:
            i = 0  # gives a counter value
            while i < k - 1:
                if cc[i + 1][0] - cc[i][1] == 1:
                    P_gap = 1
                    n_gaps += 1  # Record number of gaps
                    loc_gaps += [G.nodes[cc[i][1]]['pos']]  # Record locations of gaps
                    size_gaps += [G.nodes[cc[i + 1][0]]['pos'] - G.nodes[cc[i][1]]['pos']]  # Record size of gaps
                else:
                    P_split = 1  # If no iso and no gap then must be split

                # elif cc[i + 1][0] < cc[i][1]: # Test for splits
                #     if cc[i + 1][1] > cc[i][1]:
                #         P_split = 1
                #     else:  # If split is within another cluster need to remove
                #         P_split = 1
                #         cc.remove(cc[i+1])
                #         k -= 1
                #         i -= 1

                i += 1

        gaps_error += [n_gaps]

    return P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps, loc_iso, loc_gaps, size_gaps, gaps_error, n_clusters


####################

def metrics_new(G):
    """
    Gives the relevant quantities we want to test on G
    :param G: Input is a networkx graph object
    :return:
    [0] - is the graph fully connected
    [1] - are there isolated nodes
    [2] - is there a gap
    [3] - is there a split
    [4] - mean degree of the network
    [5] - number of isolated nodes
    [6] - number of gaps
    [7] - locations of isolated nodes
    [8] - locations of gaps
    [9] - size of gaps
    [10] - gaps error
    [11] - no. of clusters (w/o isolated nodes)
    [12] - cluster size (w/o isolated nodes)
    """
    P_fc = 0  # Give starting values for all the parameters we wish to return
    P_iso = 0
    P_gap = 0
    P_split = 0
    P_isoUucg = 0
    mean_degree = sum(list(dict(G.degree()).values())) / len(G)
                                                # Overall degree of the network divided by the number of nodes
    n_iso = 0
    n_gaps = 0
    loc_iso = []
    loc_gaps =[]
    size_gaps = []
    gaps_error = []
    n_clusters = 0
    # cluster_size = 0

    if nx.is_connected(G):
        P_fc = 1  # Test full connectivity, if true can stop
        n_clusters += 1
    else:
        if list(nx.isolates(G)):
            P_iso = 1  # Test for isolated nodes

        isolated = [i for i in nx.isolates(G)]
        n_iso += len(isolated)  # Record number of isolated nodes
        if P_iso == 1:
            loc_iso = [G.nodes[i]['pos'] for i in isolated]  # Record location of isolated nodes

        G.remove_nodes_from(isolated)  # need to remove isolated nodes for later tests

        clusters = list(nx.connected_components(G))  # Gives a list of all of the connected clusters in our graph G
        n_clusters += len(clusters)

        cc = sorted([[min(clusters[i]), max(clusters[i])] for i in range(len(clusters))])

        k = len(cc)  # gives the number of clusters we are looking at

        if k > 1:

            i = 0  # gives a counter value
            while i < k - 1:
                if cc[i + 1][0] - cc[i][1] == 1:
                    P_gap = 1
                    n_gaps += 1  # Record number of gaps
                    loc_gaps += [G.nodes[cc[i][1]]['pos']]  # Record locations of gaps
                    size_gaps += [G.nodes[cc[i + 1][0]]['pos'] - G.nodes[cc[i][1]]['pos']]  # Record size of gaps
                else:
                    P_split = 1  # If no iso and no gap then must be split

                # elif cc[i + 1][0] < cc[i][1]: # Test for splits
                #     if cc[i + 1][1] > cc[i][1]:
                #         P_split = 1
                #     else:  # If split is within another cluster need to remove
                #         P_split = 1
                #         cc.remove(cc[i+1])
                #         k -= 1
                #         i -= 1

                i += 1
        if P_iso == 1 or P_gap == 1:
            P_isoUucg = 1

        gaps_error += [n_gaps]

    return P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps, loc_iso, loc_gaps, size_gaps, gaps_error, n_clusters, P_isoUucg

####################

def run_sims(ntrials, conn_fun, eta, L, Torus, mu, new=False, rt=1):
    """
    This function will run the simulation
    :param ntrials: Requires as an input the number of trials you wish to run
    :param L: Requires as an input the length of the line segment
    :param mu: Requires as an input the value of mu for the connection function
    :param rt: Requires as an input the rate of the PPP
    :param eta: Requires as an input the value of eta for the connection function
    :return: Will then write the data to files
    """

    count = 0
    P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps = 0, 0, 0, 0, 0, 0, 0
    loc_iso, loc_gaps, size_gaps, gaps_error = [], [], [], []

    if new == True:
        P_isoUucg = 0

    mu = round(mu, 4)
    print("mu = %s" %mu)

    while count < ntrials:
        if (count / ntrials * 100) % 10 == 0:
            print(count / ntrials * 100, "% complete")
        if conn_fun == "Waxman":
            eta = 1.0
        else:
            eta = eta

        edge_list, node_list, pos_list = edges(rt, L, mu, eta)
        G = graph(edge_list, node_list, pos_list)

        if Torus == True:
            vals = metrics_Torus(G)

        elif new == True:
            vals = metrics_new(G)
            P_isoUucg += vals[12]

        else:
            vals = metrics(G)
        P_fc += vals[0]
        P_iso += vals[1]
        P_gap += vals[2]
        P_split += vals[3]
        mean_degree += vals[4]
        n_iso += vals[5]
        n_gaps += vals[6]

        loc_iso += vals[7]
        loc_gaps += vals[8]
        size_gaps += vals[9]
        gaps_error += vals[10]

        count += 1

    if Torus == True:
        fname = "Simulated_Data/%s/Eta_%s/L_%s/Torus/Overall.txt" % (conn_fun, eta, L)
        path = fname[:-12]

    elif new == True:
        fname = "Simulated_Data/%s/Eta_%s/L_%s/UnionTest/Overall.txt" % (conn_fun, eta, L)
        path = fname[:-12]

    else:
        fname = "Simulated_Data/%s/Eta_%s/L_%s/Line/Overall.txt" % (conn_fun, eta, L)
        path = fname[:-12]


    mu = round(mu, 4)  # makes sure that we have a readable mu

    if new == True:
        if not os.path.exists(path):
            os.makedirs(path)

        if not os.path.isfile(fname):

            with open(fname, "w") as my_file:
                data = {mu: [ntrials, P_fc, P_iso, P_gap, P_isoUucg, P_split, mean_degree, n_iso, n_gaps]}
                json.dump(data, my_file)

        else:
            with open(fname, "r") as my_file:
                data = json.load(my_file)
            with open(fname, "w") as my_file:
                mulabel = str(mu)
                if mulabel in list(data.keys()):
                    lastcount, last_pfc, last_iso, last_gap, last_isoUucg, last_split, last_degree, last_n_iso, last_n_gaps = data[
                        mulabel]
                    data[mulabel] = [lastcount + ntrials,
                                     last_pfc + P_fc,
                                     last_iso + P_iso,
                                     last_gap + P_gap,
                                     last_isoUucg + P_isoUucg,
                                     last_split + P_split,
                                     last_degree + mean_degree,
                                     last_n_iso + n_iso,
                                     last_n_gaps + n_gaps]
                    json.dump(data, my_file)
                else:
                    data[mu] = [ntrials, P_fc, P_iso, P_gap, P_isoUucg, P_split, mean_degree, n_iso, n_gaps]
                    json.dump(data, my_file)

    else:
        if not os.path.exists(path):
            os.makedirs(path)

        if not os.path.isfile(fname):

            with open(fname, "w") as my_file:
                data = {mu: [ntrials, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps]}
                json.dump(data, my_file)

        else:
            with open(fname, "r") as my_file:
                data = json.load(my_file)
            with open(fname, "w") as my_file:
                mulabel = str(mu)
                if mulabel in list(data.keys()):
                    lastcount, last_pfc, last_iso, last_gap, last_split, last_degree, last_n_iso, last_n_gaps = data[mulabel]
                    data[mulabel] = [lastcount + ntrials,
                                     last_pfc + P_fc,
                                     last_iso + P_iso,
                                     last_gap + P_gap,
                                     last_split + P_split,
                                     last_degree + mean_degree,
                                     last_n_iso + n_iso,
                                     last_n_gaps + n_gaps]
                    json.dump(data, my_file)
                else:
                    data[mu] = [ntrials, P_fc, P_iso, P_gap, P_split, mean_degree, n_iso, n_gaps]
                    json.dump(data, my_file)

    # if Torus == True:
    #     iso_loc_fname = "Simulated_Data/%s/Eta_%s/L_%s/Torus/IsoLoc.txt" % (conn_fun, eta, L)
    #
    # else:
    #     iso_loc_fname = "Simulated_Data/%s/Eta_%s/L_%s/Line/IsoLoc.txt" % (conn_fun, eta, L)
    #
    # if not os.path.isfile(iso_loc_fname):
    #
    #     with open(iso_loc_fname, "w") as my_file:
    #         data = {mu: loc_iso}
    #         json.dump(data, my_file)
    #
    # else:
    #     with open(iso_loc_fname, "r") as my_file:
    #         data = json.load(my_file)
    #     with open(iso_loc_fname, "w") as my_file:
    #         mulabel = str(mu)
    #         if mulabel in list(data.keys()):
    #             last_loc_iso = data[mulabel]
    #             data[mulabel] = last_loc_iso + loc_iso
    #             json.dump(data, my_file)
    #         else:
    #             data[mu] = loc_iso
    #             json.dump(data, my_file)
    #
    # if Torus == True:
    #     gap_loc_fname = "Simulated_Data/%s/Eta_%s/L_%s/Torus/GapLoc.txt" % (conn_fun, eta, L)
    #
    # else:
    #     gap_loc_fname = "Simulated_Data/%s/Eta_%s/L_%s/Line/GapLoc.txt" % (conn_fun, eta, L)
    #
    # if not os.path.isfile(gap_loc_fname):
    #
    #     with open(gap_loc_fname, "w") as my_file:
    #         data = {mu: loc_gaps}
    #         json.dump(data, my_file)
    #
    # else:
    #     with open(gap_loc_fname, "r") as my_file:
    #         data = json.load(my_file)
    #     with open(gap_loc_fname, "w") as my_file:
    #         mulabel = str(mu)
    #         if mulabel in list(data.keys()):
    #             last_loc_gaps = data[mulabel]
    #             data[mulabel] = last_loc_gaps + loc_gaps
    #             json.dump(data, my_file)
    #         else:
    #             data[mu] = loc_gaps
    #             json.dump(data, my_file)
    #
    # if Torus == True:
    #     gap_size_fname = "Simulated_Data/%s/Eta_%s/L_%s/Torus/GapSize.txt" % (conn_fun, eta, L)
    #
    # else:
    #     gap_size_fname = "Simulated_Data/%s/Eta_%s/L_%s/Line/GapSize.txt" % (conn_fun, eta, L)
    #
    # if not os.path.isfile(gap_size_fname):
    #
    #     with open(gap_size_fname, "w") as my_file:
    #         data = {mu: size_gaps}
    #         json.dump(data, my_file)
    #
    # else:
    #     with open(gap_size_fname, "r") as my_file:
    #         data = json.load(my_file)
    #     with open(gap_size_fname, "w") as my_file:
    #         mulabel = str(mu)
    #         if mulabel in list(data.keys()):
    #             last_gap_size = data[mulabel]
    #             data[mulabel] = last_gap_size + size_gaps
    #             json.dump(data, my_file)
    #         else:
    #             data[mu] = size_gaps
    #             json.dump(data, my_file)
    #
    # if Torus == True:
    #     gap_error_fname = "Simulated_Data/%s/Eta_%s/L_%s/Torus/GapError.txt" % (conn_fun, eta, L)
    #
    # else:
    #     gap_error_fname = "Simulated_Data/%s/Eta_%s/L_%s/Line/GapError.txt" % (conn_fun, eta, L)
    #
    # if not os.path.isfile(gap_error_fname):
    #
    #     with open(gap_error_fname, "w") as my_file:
    #         data = {mu: gaps_error}
    #         json.dump(data, my_file)
    #
    # else:
    #     with open(gap_error_fname, "r") as my_file:
    #         data = json.load(my_file)
    #     with open(gap_error_fname, "w") as my_file:
    #         mulabel = str(mu)
    #         if mulabel in list(data.keys()):
    #             last_gap_error = data[mulabel]
    #             data[mulabel] = last_gap_error + gaps_error
    #             json.dump(data, my_file)
    #         else:
    #             data[mu] = gaps_error
    #             json.dump(data, my_file)
