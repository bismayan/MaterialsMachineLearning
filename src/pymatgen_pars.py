import pymatgen as mg
from pymatgen.matproj.rest import MPRester
import pandas as pd
from pymatgen import Element,Composition
import multiprocessing as mp
import pickle
import json
from monty.json import MontyEncoder,MontyDecoder
import numpy as np


def ret_struct_obj(i):
    return mg.Structure.from_str(i,fmt="cif")

def return_struct_list():
        with open("../data/ternaries_from_mg.pickle",'r') as f:
            temp_list=pickle.load(f)
        p=mp.Pool(4)
        struct_lis=p.map(ret_struct_obj,temp_list)
        return struct_lis

def read_ternaries():
    with MPRester() as m:
        ternaries1 = m.query(criteria={"nelements": 3}, properties=['icsd_ids', 'pretty_formula', 'cif'])
        list_cif = [i['cif'] for i in ternaries1]
        outfile=open("../data/ternaries_from_mg.pickle",'w')
        pickle.dump(list_cif,outfile)
        del(list_cif)
        outfile.close()

def read_unique_data(filename):
    with open(filename,'r') as f:
        structs=json.load(f,cls=MontyDecoder)
    return structs

def get_space_groups(strts):
    sgroups=np.array([a.get_spacegroup_info()[1] for a in strts])
    return sgroups

def read_data(filename):
    """
    :argument
    filename - The filename of the json: file to read from
    :returns
    DataFrame - Pandas Dataframe containing the formatted parsed data
    """
    uniq_data=read_unique_data(filename)
    space_groups=get_space_groups(uniq_data)
    (comps,stoich_coeffs,at_nos,eneg)=get_comp_data(uniq_data)
    DataFrame = pd.DataFrame({"Z1": at_nos[:, 0]})
    DataFrame["Z2"] = at_nos[:, 1]
    DataFrame["Z3"] = at_nos[:, 2]
    DataFrame["St_coeff1"] = stoich_coeffs[:, 0]
    DataFrame["St_coeff2"] = stoich_coeffs[:, 1]
    DataFrame["St_coeff3"] = stoich_coeffs[:, 2]
    DataFrame["Eneg1"] = eneg[:, 0]
    DataFrame["Eneg2"] = eneg[:, 1]
    DataFrame["Eneg3"] = eneg[:, 2]
    DataFrame["Space Group"] = space_groups
    DataFrame["Formulae"] = [c.formula for c in comps]
    return DataFrame

def get_comp_data(un_data):
    element_universe = [str(e) for e in Element]
    dict_element = {}
    for i, j in enumerate(element_universe):
        dict_element[str(j)] = i
    stoich_array = np.zeros((len(un_data), 3), dtype=float)
    at_num_array = np.zeros((len(un_data), 3), dtype=int)
    electroneg_array = np.zeros((len(un_data), 3), dtype=float)
    comp_array=[a.composition for a in un_data]
    temp_dict_list = [dict(comp.get_el_amt_dict()) for comp in comp_array]
    for index,temp_dict in enumerate(temp_dict_list):
        for count, key in enumerate(temp_dict.keys()):
            stoich_array[index][count] = temp_dict[key]
            if key not in ['D', 'T']:
                at_num_array[index][count] = Element(key).Z
                electroneg_array[index][count] = Element(key).X
            else:
                at_num_array[index][count] = Element('H').Z
                electroneg_array[index][count] = Element('H').X
    del(dict_element)
    del(temp_dict_list)
    return (comp_array,stoich_array,at_num_array,electroneg_array)

#if __name__=="__main__":


