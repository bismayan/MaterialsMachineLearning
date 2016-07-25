#! /usr/bin/env/python

import numpy as np
import pandas as pd
from pymatgen.core import Element, Composition
import csv

def read_data(filename):
    """
    :argument
    filename - The filename of the csv file to read from
    :returns
    DataFrame - Pandas Dataframe containing the formatted parsed data
    """
    data=pd.read_csv(filename)
    uniq_data=get_unique_comps(data)
    space_groups=get_spacegroups(uniq_data)
    (stoich_coeffs,at_nos,eneg)=get_comp_data(unique_data)
    DataFrame=pd.DataFrame({"Z":at_nos,"Stoichiometry Coefficients":stoich_coeffs,"Electronegativities":eneg,"Space Groups":space_groups})
    return DataFrame

def get_unique_comps(csv_data,max_num=100,constraint="integer"):
    unique_data = []
    found_comps = []
    for line in data:
        form = Composition(line[2]).formula
        if form not in found_comps:
            unique_data.append(line)
            found_comps.append(form)
    return unique_data

def get_spacegroups(un_data):
    for row1 in un_data:
        row1[1] = row1[1].replace(' ', '')
    list_space = [row1[1].rstrip('Z').rstrip('S').rstrip("H").rstrip('R') for row1 in un_data]
    with open("ICSD/spacegroups.dat", 'r') as f:
        dat = csv.reader(f, dialect='excel-tab', quoting=csv.QUOTE_NONE)
        list_dat = [element.strip() for row in dat for element in row]
        list1 = [[int(list_dat[i * 2]), list_dat[i * 2 + 1]] for i in range(int(len(list_dat) / 2))]
    dict_space = {}
    for i in range(len(list1)):
        dict_space[list1[i][1]] = list1[i][0]
    with open('ICSD/spacegroups_2.dat', 'r') as f1:
        f = f1.readlines()
        for line in f:
            data2 = [element.strip() for element in line.split()]
            if data2[1] not in dict_space.keys():
                dict_space[data2[1]] = int(data2[0])
    with open('ICSD/spacegroups_3.dat', 'r') as f1:
        f = f1.readlines()
        for line in f:
            data3 = [element.strip() for element in line.split()]
            if data3[0] not in dict_space.keys():
                dict_space[data3[0]] = int(data3[1])
    space_num_array = np.zeros(len(list_space), dtype=float)
    for i, s in enumerate(list_space):
        space_num_array[i] = dict_space[s]
    return space_num_array

    
def get_comp_data(un_data):
    element_universe = [str(e) for e in Element]
    dict_element = {}
    for i, j in enumerate(element_universe):
        dict_element[str(j)] = i
    stoich_array = np.zeros((len(un_data), 3), dtype=float)
    at_num_array = np.zeros((len(un_data), 3), dtype=float)
    electroneg_array = np.zeros((len(un_data), 3), dtype=float)
    for index, entry in enumerate(un_data):
        comp = Composition(entry[2])
        temp_dict = dict(comp.get_el_amt_dict())
        for count, key in enumerate(temp_dict.keys()):
            stoich_array[index][count] = temp_dict[key]
            if key not in ['D', 'T']:
                at_num_array[index][count] = Element(key).Z
                electroneg_array[index][count] = Element(key).X
            else:
                at_num_array[index][count] = Element('H').Z
                electroneg_array[index][count] = Element('H').X
    return (stoich_array,at_num_array,electroneg_array)












