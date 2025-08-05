# -*- coding: UTF-8 -*-
# @Date     : 15th Jul, 2022
# @Author   : Northblue

"""Singleton pattern of ion matching parameters
"""
import os
from functools import wraps

import aa_ele
import linker_file
import modification_file
from pathlib import Path

PARAM_FILE_DIR = os.path.dirname(os.path.abspath(__file__))

# # multiprocessing
# # Can't pickle <class 'ion_match_param.Param'>: it's not the same object as ion_match_param.Param
# https://stackoverflow.com/questions/52185507/pickle-and-decorated-classes-picklingerror-not-the-same-object
# def singleton(cls):
#     _instances = {}
#     @wraps(cls)
#     def getinstance(*args, **kwargs):
#         if cls not in _instances:
#             _instances[cls] = cls(*args, **kwargs)
#         return _instances[cls]
#     return getinstance


# @singleton
class Param(object):

    _init_flag = False

    def __init__(self):
        if self._init_flag:
            return

        # print('__init__')
        self._init_flag = True

        # # ======param
        self.fpath_res = r'F:\DataSet\pLink2_test\Leiker_15N_labeling\mpz\pLink_task_2025.01.09.12.12.44\reports\uniprot-ecoli-20171023_2025.01.09.filtered_cross-linked_spectra.csv'
        self.dpath_mgf = r'F:\DataSet\pLink2_test\Leiker_15N_labeling\data'
        self.is_get_ms1_inten = True
        self.n_process = 2 # if multiple processes, need a lot of memory
        self.progress_visual = True

        self.ms2_match_tolerance = 20 # ppm

        # limit the charge of the fragment ion
        self.matching_fragment_charge_lower_precursor = True

        # charge_mode
        self.charge_mode = '+' # '-'

        # charges considered for the MS-Cleavable crosslinked characteristic ion
        self.cleav_characteristic_charge = [1, 2, 3]
        self.is_dsso_long_add_h2o = True

        # param path
        # self.dpath_plink_bin = r'C:\ProgramPfind\pFindStudio\pLink\2.3.11\bin'
        self.dpath_plink_bin = PARAM_FILE_DIR

        self.fpath_aa = os.path.join(self.dpath_plink_bin, 'aa.ini')
        self.fpath_ele = os.path.join(self.dpath_plink_bin, 'element.ini')
        self.fpath_mod = str(Path(self.dpath_plink_bin)/'modification.ini')
        self.fpath_linker = os.path.join(self.dpath_plink_bin, 'xlink.ini')

        # mass list
        self.aa2mass = aa_ele.read_aa2mass(self.fpath_aa, self.fpath_ele)
        # self.mod2mass = {'Carbamidomethyl[C]':57.021464, 'Oxidation[M]':15.994915}
        self.mod2mass = modification_file.read_mod2mass(self.fpath_mod)
        self.mod2mass['TDS_OH'] = 522.170
        self.df_linker = self.get_df_linker(self.fpath_linker)

        # isotopic labeling
        self.is_labeled = False
        self.mp_label_15n = {'R':{'*':{'N':'15N'}}}
        self.plink_labelid = {1:{'R':{}}, 2:self.mp_label_15n}

        # labeled mass list
        self.aa2mass_labeled = {} # mass
        for label_id, mp_label in self.plink_labelid.items():
            self.aa2mass_labeled[label_id] = aa_ele.read_aa2mass(self.fpath_aa, self.fpath_ele, mp_label['R'])
        
        self.aa2ele_labeled = {} # element
        for label_id, mp_label in self.plink_labelid.items():
            self.aa2ele_labeled[label_id] = aa_ele.label_aa2ele(aa_ele.read_aa2ele(self.fpath_aa), mp_label['R'])

        # # ======default param
        self.ms2_match_wind = 0.5 # Da # 0.5Da must more than 20ppm

        # # ======constant
        # self.precursor_charge_max = 100
        
        # mass
        self.mass_H = 1.00782503214
        self.mass_C = 12.00000000000
        self.mass_N = 14.00307400524
        self.mass_O = 15.9949146221
        self.mass_P = 30.973761998
        self.mass_S = 31.972070690
        self.mass_electron = 0.0005485799
        self.mass_proton = 1.007276466621
        self.mass_H2O = self.mass_O + 2*self.mass_H
        self.mass_NH3 = self.mass_N + 3*self.mass_H


    def get_df_linker(self, fpath):
        linker = linker_file.Linker()
        linker.read_plink_xlink(fpath)
        return linker.df_linker


    # # singleton
    # # https://www.cnblogs.com/huchong/p/8244279.html
    # def __new__(cls, *args, **kwargs):
        # if not hasattr(cls, '_instance'):
            # print('------not hasattr')
            # cls._instance = super().__new__(cls, *args, **kwargs)
        # return cls._instance

    # # singleton, init once, set _init_flag
    # # https://cloud.tencent.com/developer/article/1564189
    # # https://blog.csdn.net/qq_26442553/article/details/94393191
    _instance = None
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super().__new__(cls,*args,**kwargs)
        return cls._instance