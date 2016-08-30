#! /usr/bin/env/python

"""
contains code that is used to cluster the fingerprint functions and then plot the clusters or refactor them in the periodic table.
"""
import numpy as np
from sklearn.cluster import KMeans

def clust_finger(finger_array, num_clusts=15, get_object=True):
	with KMeans(n_clusters=num_clusts,random_state=42) as Km:
		clusts=Km.fit_predict(finger_array)
		if get_object:
			return (clusts, Km)
		else:
			return (clusts)


def get_scatterplot_summed(ref_array, Z_tag):
    """
    Return the arrays required to create scatter plots for the summed array

    :param ref_array: The array of clustered elements
    :param Z_tag: the atomic number of element you want to get arrays for
    :return: (x,y,z,dz) array for plotting 3d bar plot
    """
    count = 0
    xarr = np.zeros(8464)
    yarr = np.zeros(8464)
    zarr = np.zeros(8464)
    dzarr = np.zeros(8464)
    ref_summed=np.sum(ref_array,axis=3)
    for (a, b), d in np.ndenumerate(ref_summed[Z_tag - 1]):
        if d > 0:
            xarr[count] = a + 1
            yarr[count] = b + 1
            dzarr[count] = d
            count = count + 1
    xarr = xarr[0:count]
    yarr = yarr[0:count]
    zarr = zarr[0:count]
    dzarr = dzarr[0:count]
    return (xarr, yarr, zarr, dzarr)

def get_scatterplot_summed_refactored(ref_array, Z_tag, Zlist):
    """
    get the refactored arrays for creating scatter plots for elements relabeled by Zlist.
    :param ref_array: The array of refactored clustered elements
    :param Z_tag: the atomic number of element you want to get arrays for
    :param Zlist: the list containing numeric indices for the different elements after regrouping
    :return: (x,y,z,dz) array for plotting 3d bar plot
    """
    count = 0
    xarr = np.zeros(8464)
    yarr = np.zeros(8464)
    zarr = np.zeros(8464)
    ref_summed = np.sum(ref_array, axis=3)
    for (a, b), d in np.ndenumerate(ref_summed[Z_tag - 1]):
        if d > 0:
            xarr[count] = Zlist.index(a + 1)
            yarr[count] = Zlist.index(b + 1)
            zarr[count] = d
            count = count + 1
    xarr = xarr[0:count]
    yarr = yarr[0:count]
    zarr = zarr[0:count]
    return (xarr, yarr, zarr)


def get_restructured_array(Formulas_s,clusters_s):
    """
    Get refactored array with clustering information reorganized by the order of elements defined by their atomic numbers
    :param Formulas_s: The Formulas of the compounds
    :param clusters_s: The clusters of the compounds
    :return: The refactored array containing the clusters as indexed by Zlist
    """
    from pymatgen import Composition
    from itertools import permutations
    Comps_s = [Composition(i) for i in Formulas_s]
    Z_max = max([E.Z for comp in Comps_s for E in comp.elements])
    Z_array = np.array([[E.Z for E in comp.elements] for comp in Comps_s])
    refactor_array = np.zeros((Z_max, Z_max, Z_max, np.amax(clusters_s) + 1))
    for i in range(len(clusters_s)):
        for index in set(permutations(Z_array[i])):
            refactor_array[index[0] - 1, index[1] - 1, index[2] - 1][clusters_s[i]] += 1
    return refactor_array


def get_refactored_array(Formulas_s, clusters_s,Z_list):
    """
    Get refactored array with clustering information reorganized by the order of elements defined in Z_list
    :param Formulas_s: The Formulas of the compounds
    :param clusters_s: The clusters of the compounds
    :param Z_list: the list containing the reordered atomic numbers
    :return: The refactored array containing the clusters as indexed by Zlist
    """
    from pymatgen import Composition
    from itertools import permutations
    Comps_s = [Composition(i) for i in Formulas_s]
    Z_max = max([E.Z for comp in Comps_s for E in comp.elements])
    Z_array = np.array([[E.Z for E in comp.elements] for comp in Comps_s])
    refactor_array = np.zeros((Z_max, Z_max, Z_max, np.amax(clusters_s) + 1))
    for i in range(len(clusters_s)):
        for index in set(permutations(Z_array[i])):
            refactor_array[Z_list.index([index[0] - 1]), Z_list.index([index[1] - 1]),Z_list.index([index[2] - 1])][clusters_s[i]] += 1
    return refactor_array




	
	
