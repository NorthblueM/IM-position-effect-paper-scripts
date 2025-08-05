#-*- coding: UTF-8 -*-

# @Date     : 14th Sep, 2021
# @Author   : Northblue

"""get amino acid and element mass
"""

import os

AA_FILE_DIR = os.path.dirname(os.path.abspath(__file__))
ELE_FILE_DIR = os.path.dirname(os.path.abspath(__file__))

def read_aa2ele(fpath=os.path.join(AA_FILE_DIR, 'aa.ini')):
    """read aa.ini
    """

    aa2comp = {}
    with open(fpath, 'r') as fin:
        for line in fin:
            if line[0] == 'R':
                aa, comp = line.split('=')[1].strip().split('|')[0:2]
                ele2n = {}
                for ele in comp.strip()[:-1].split(')'):
                    e, n = ele.split('(')
                    if e not in ele2n: # if the same elements written separately
                        ele2n[e] = int(n)
                    else:
                        ele2n[e] += int(n)
                aa2comp[aa] = ele2n
    return aa2comp

def read_ele2mass(fpath=os.path.join(ELE_FILE_DIR, 'element.ini')):
    """read element.ini
    """

    ele2mass = {}
    with open(fpath, 'r') as fin:
        for line in fin:
            if line[0] == 'E':
                ele, mass, prob = line.split('=')[1].strip().split('|')[0:3]
                ms=[float(i) for i in mass.split(',')[:-1]]
                ps=[float(i) for i in prob.split(',')[:-1]]
                comb=list(zip(ms, ps))
                mass_prob_high = sorted(comb, key=lambda x:x[1])[-1][0]
                ele2mass[ele] = mass_prob_high
    return ele2mass

def label_aa2ele(aa2comp, mp_label):
    """update labeled aa components
    """

    if '*' in mp_label:
        for aa in aa2comp:
            mp_label[aa] = mp_label['*']

    aa2comp_label = {}
    for aa, comp in aa2comp.items():
        if aa in mp_label:
            aa2comp_label[aa] = {}
            for e, n in comp.items():
                if e in mp_label[aa]:
                    aa2comp_label[aa][mp_label[aa][e]] = n
                else:
                    aa2comp_label[aa][e] = n
        else:
            aa2comp_label[aa] = comp

    return aa2comp_label


def read_aa2mass(fpath_aa=os.path.join(AA_FILE_DIR, 'aa.ini'), fpath_ele=os.path.join(ELE_FILE_DIR, 'element.ini'), mp_label={}):
    """get aa monoisotopic mass
    Args:
        mp_label: {aa:{ele:ele_labeled,},}, aa=* means all aa 
    """

    aa2comp = read_aa2ele(fpath_aa)
    ele2mass = read_ele2mass(fpath_ele)

    aa2comp = label_aa2ele(aa2comp, mp_label)

    aa2mass = {}
    for aa, comp in aa2comp.items():
        mass = 0.0
        for e, n in comp.items():
            mass += (ele2mass[e] * n)
        aa2mass[aa] = mass
    
    return aa2mass


def get_seq_mass(seq):
    mass = 0.0
    aa2mass = read_aa2mass()
    for aa in seq:
        mass += aa2mass[aa]
    return mass


def get_molecule_comp(molecule_str):
    """
    Args:
        molecule_str: C(16)H(17)N(5)O(3)
    Returns:
        {'C':16, 'H':17, 'N':5, 'O':3}
    """

    comp = molecule_str
    ele2n = {}
    for ele in comp.strip()[:-1].split(')'):
        e, n = ele.split('(')
        if e not in ele2n: # if the same elements written separately
            ele2n[e] = int(n)
        else:
            ele2n[e] += int(n)
    return ele2n


def get_comp_mass(comp, ele2mass):
    """
    Args:
        comp: {'C':16, 'H':17, 'N':5, 'O':3}
    Returns:
        mass
    """

    mass = 0.0
    for e, n in comp.items():
        mass += (ele2mass[e] * n)
    return mass


def get_molecule_mass(molecule_str, ele2mass):

    comp = get_molecule_comp(molecule_str)
    mass = get_comp_mass(comp, ele2mass)
    return mass        
          
                    
def _main():
    aa2mass = read_aa2mass()
    print(aa2mass)


if __name__ == '__main__':
    _main()
