import numpy as np
import pickle
import matplotlib.pyplot as plt
from collections import Counter
import scipy.integrate


def read_pickle(filename):
    """ Reads in a pickle file and returns the object

    :param filename:  string containing name of filename
    :return: The object stored in the pickle file
    """
    return pickle.load(open(filename, 'rb'))

def createKernel(size,sigma,delta):
    """
    Creates a Kernel for approximating a delta function by a gaussian

    :param size: the number of array elements that the Gaussian spans on each side of the maximum
    :param sigma: The sigma parameter of the gaussian
    :param delta: The size of each array element in angstroms (mesh size so to speak)
    :return: an array of length 2*size+1 which contains the gaussian
    """
    prefac = 1.0 / (np.sqrt(2 * np.pi * sigma * sigma))
    kern = np.array([prefac * (np.exp(-(i * delta) * (i * delta) / (2 * sigma * sigma))) for i in range(-1 * size, size + 1)])
    return kern


def getz(specie):
    """
    Return the atomic number of the specie

    :param specie: the pymatgen elemental Specie
    :return: The atomic number of the specie
    """
    return specie.Z


def getchi(specie):
    """
    Return the atomic number of the specie

    :param specie: the pymatgen elemental Specie
    :return: The electronegativity number of the specie
    """
    return specie.X

def get_ones(specie):
    return 1.0

def create_obs_dict(obser,specie_list):
    func_map= {"Z":getz,"Chi":getchi,"ones":get_ones}
    obs_dict=dict()
    for spec in specie_list:
        obs_dict[spec]=func_map[obser](spec)
    return obs_dict


def get_phi(struct,obser='ones',rmax=15,delta=0.05,sigma=0.075,kernsize=10,tol=0.1):

    """
    Get the fingerprint phi for a given structure
    :param struct: The pymatgen structure you want to use
    :param rmax: The maximum r in Angstroms you want to create the fingerprint for
    :param obser: The observable you want to calculate the phi for
    :param delta: The size of each array element in angstroms (mesh size so to speak)
    :param sigma: The sigma parameter of the gaussian used to approximate delta function
    :param kernsize: the number of array elements that the Gaussian spans on each side of the maximum
    :param tol: The tolerance upto which two sites are assumed to be the same in space
    """
    spec = struct.species
    kern=createKernel(kernsize,sigma,delta)
    n_bins = rmax / delta
    count = Counter(spec)
    els = [a for a in count.keys()]
    obs=create_obs_dict(obser,els)
    finger = np.zeros((int(n_bins + kernsize)), dtype=float)
    nb_unit=sum([count[entry1]*obs[entry1]*count[entry2]*obs[entry2] for entry1 in els for entry2 in els])
    prefac_outer=struct.volume/(nb_unit)
    for i, s in enumerate(spec):
        b_i=obs[s]
        list_at = struct.get_sites_in_sphere(struct.sites[i].coords, rmax)
        dists = np.array([entry[-1] for entry in list_at if entry[-1] > tol])
        names = [entry[0].specie for entry in list_at if entry[-1] > tol]
        for j, dist in enumerate(dists):
            b_j=obs[names[j]]
            prefac_inner=b_i*b_j/(4 * np.pi * dist * dist)
            poin = np.floor(dist / delta)
            finger[int(poin - kernsize):int(poin + kernsize + 1)] += ((prefac_inner) * kern)

    finger=(finger*prefac_outer)-1
    print "Average of last twenty bins", np.mean(finger[-1 * (kernsize + 20):-1 * kernsize])
    return finger[0:-kernsize]

def get_phi_scaled(struct,obser='ones',n_bins=100,delta=0.05,sigma=0.075,kernsize=12,tol=1e-3,debug=False):

    """
    Get the fingerprint scaled phi for a given structure
    :param struct: The pymatgen structure you want to use
    :param rmax: The maximum r in Angstroms you want to create the fingerprint for
    :param obser: The observable you want to calculate the phi for
    :param delta: The size of each array element in units scale factor (mesh size so to speak)
    :param sigma: The sigma parameter of the gaussian used to approximate delta function in units of scale factor
    :param kernsize: the number of array elements that the Gaussian spans on each side of the maximum
    :param tol: The tolerance upto which two sites are assumed to be the same in space
    """
    spec = struct.species
    # sf is the calculated scaling factor
    sf=(struct.volume/len(struct.species))**(0.33)
    kern=createKernel(kernsize,sigma*sf,delta*sf)
    rmax=n_bins*sf*delta
    count = Counter(spec)
    els = [a for a in count.keys()]
    obs=create_obs_dict(obser,els)
    finger = np.zeros((int(n_bins + 2*kernsize)), dtype=float)
    nb_unit=sum([count[entry1]*obs[entry1]*count[entry2]*obs[entry2] for entry1 in els for entry2 in els])
    prefac_outer=(struct.volume)/(4*np.pi*nb_unit*(sf*sf))
    for i, s in enumerate(spec):
        b_i=obs[s]
        list_at = struct.get_sites_in_sphere(struct.sites[i].coords, rmax)
        dists = np.array([entry[-1]/sf for entry in list_at if entry[-1] > tol])
        names = [entry[0].specie for entry in list_at if entry[-1] > tol]
        for j, dist in enumerate(dists):
            b_j=obs[names[j]]
            prefac_inner=b_i*b_j/(dist * dist)
            poin = np.floor(dist / delta)+kernsize
            finger[int(poin - kernsize):int(poin + kernsize + 1)] += ((prefac_inner) * kern)
		    
    finger=(finger*prefac_outer)-1
    if debug:
	    print "Scale factor=",sf
	    print "Average of last twenty bins", np.mean(finger[-1 * (kernsize + 20):-1 * kernsize])
    return finger[kernsize:-kernsize]


def check_sumrule(phi,struct,delta=0.05):
    R = np.array([(i) * delta for i in range(len(phi))])
    y=np.array([(R[i]**2)*phi[i] for i in range(len(phi))])
    plt.plot(R,y)
    N=len(struct.species)
    plt.show()
    return scipy.integrate.simps(4*np.pi*y*N/struct.volume,dx=delta)

def main():
    structs = read_pickle('struct_all.pickle')
    phi_1=get_phi(structs[0],'ones')
    plt.figure(figsize=(10, 5))
    plt.plot(phi_1[0])
    plt.show()


if __name__=="__main__":
    main()



