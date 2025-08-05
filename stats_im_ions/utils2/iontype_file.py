# -*- coding: UTF-8 -*-
# @Date     : 15th Jul, 2022
# @Author   : Northblue

"""iontype
"""

import os
import pandas as pd

from ion_match_param import Param

IONTYPE_FILE_DIR = os.path.dirname(os.path.abspath(__file__))
FPATH_IONTYPE_TEMPLATE = os.path.join(IONTYPE_FILE_DIR, 'iontype_ini_template.csv')
FPATH_IONTYPE_USED = os.path.join(IONTYPE_FILE_DIR, 'iontype_ini_used.csv')

class IonType(object):

    def __init__(self, param=Param()):


        self.param = param

        # ['sym', 'term', 'mass_shift', 'charge', 'link_mode', 'link_mass', 'cleav_sym', 'cleav_mass', 'characteristic_ion']
        # add ['iontype_str']
        self.df_iontype = pd.DataFrame()
        
        self.fpath = ''

        """
        columns:
            term:
                0, n-term ion
                1, c-term ion
            mass_shift: 
                increase in mass summation of amino acid residues
            link_mode:
                0, without linker;
                1, with non-cleav linker;
                3, with cleav linker
            characteristic_ion: 
                0, fragment ion; 
                1, cleavable characteristic ion
        """

    def set_id(self):
        """id start at 0
        """
        self.df_iontype['id'] = list(range(self.df_iontype.shape[0]))


    def reset_index(self):
        self.df_iontype = self.df_iontype.reset_index(drop=True)


    def read_iontype(self, fpath=FPATH_IONTYPE_TEMPLATE):
        if os.path.exists(fpath):
            self.df_iontype = pd.read_csv(fpath)
        self.reset_index()


    def write_iontype(self, fpath=FPATH_IONTYPE_USED, linker_name=''):
        if linker_name:
            fpath_out = fpath[:-4] + '_' + linker_name + '.csv'
        else:
            fpath_out = fpath
        self.df_iontype.to_csv(fpath_out, index=False)
        
        self.fpath = fpath_out


    def set_linker_mass(self, ser):
        """ser: linker info
        """
        param = self.param

        self.df_iontype['link_mass'] = ser['linker_mass']
        df_tmp = self.df_iontype.copy(deep=True)

        # corsslink ion with non-cleavable linker
        # mass1 = sum(aa)1 + mass_shift + link_mass + mass_pep2
        df_mode1 = df_tmp.copy(deep=True)
        df_mode1['link_mode'] = 1
        self.df_iontype = pd.concat([self.df_iontype, df_mode1])
        
        # corsslink ion with cleavable linker
        # mass1 = sum(aa)1 + mass_shift + cleav_mass
        if ser['ms_cleavable'] == 1:
            df_long = df_tmp.copy(deep=True)
            df_long['link_mode'] = 3
            df_long['cleav_sym'] = 'L'
            df_long['cleav_mass'] = ser['long_mass']

            df_short = df_tmp.copy(deep=True)
            df_short['link_mode'] = 3
            df_short['cleav_sym'] = 'S'
            df_short['cleav_mass'] = ser['short_mass']

            self.df_iontype = pd.concat([self.df_iontype, df_long, df_short])

            if param.is_dsso_long_add_h2o and ('DSSO' in ser['name']):
                df_long_h2o = df_tmp.copy(deep=True)
                df_long_h2o['link_mode'] = 3
                df_long_h2o['cleav_sym'] = 'L+H2O'
                df_long_h2o['cleav_mass'] = ser['long_mass']+param.mass_H2O
                self.df_iontype = pd.concat([self.df_iontype, df_long_h2o])

            # cleavable characteristic ion
            # based on summation of amino acid residues, not include H2O
            # mass1 = sum(aa)1 + mass_shift + cleav_mass
            if param.cleav_characteristic_charge:
                s = self.df_iontype.iloc[0].copy(deep=True)
                s['characteristic_ion'] = 1
                s['sym'] = '-'
                s['mass_shift'] = param.mass_H2O
                s['link_mode'] = 3
                s['link_mass'] = ser['linker_mass']

                for charge in param.cleav_characteristic_charge:
                    s['charge'] = charge
                    s['cleav_sym'] = 'L'
                    s['cleav_mass'] = ser['long_mass']
                    self.df_iontype = pd.concat([self.df_iontype, s.to_frame().T])

                    s['cleav_sym'] = 'S'
                    s['cleav_mass'] = ser['short_mass']
                    self.df_iontype = pd.concat([self.df_iontype, s.to_frame().T])

                    if param.is_dsso_long_add_h2o and ('DSSO' in ser['name']):
                        s['cleav_sym'] = 'L+H2O'
                        s['cleav_mass'] = ser['long_mass']+param.mass_H2O
                        self.df_iontype = pd.concat([self.df_iontype, s.to_frame().T])

        self.reset_index()


    def set_xl_ion(self, ser, fpath=FPATH_IONTYPE_TEMPLATE):
        """ser: linker info
        """
        self.read_iontype(fpath)
        self.set_linker_mass(ser)
        self.reset_index()
        self.write_iontype(linker_name=ser['name'])


    def set_df_iontype_str(self):
        if 'iontype_str' not in  self.df_iontype.columns:
            mp_iontype_id2str = get_iontype_id2str(self.df_iontype)
            self.df_iontype['iontype_str'] = list(mp_iontype_id2str.values())



def get_iontype_str(ser, param=Param()):
    """ser: one iontype
    return example:
        b|-|+
        y-H2O|L+H2O|++
    """

    charge_mode = param.charge_mode

    s = []
    if ser['characteristic_ion'] == 0: # fragment ion
        s.append(ser['sym'])
        if ser['link_mode'] == 0:
            s.append('-')
        if ser['link_mode'] == 1: # contain linker and another complete peptide
            s.append('X')
        elif ser['link_mode'] == 3: # contain breaking linker
            s.append(ser['cleav_sym'])

    elif ser['characteristic_ion'] == 1: # ms-cleavable characteristic ion
        s.append(ser['sym'])
        s.append(ser['cleav_sym'])
    
    elif ser['characteristic_ion'] == 2: # small molecule characteristic ion
        s.append(ser['sym'])
        s.append(':')

    s.append(str(ser['charge'])+charge_mode)
    return '|'.join(s)


def get_iontype_id2str(df):
    """
    Returns:
        dict = {1:'b|-|+'}
    """

    mp_id2str = {}
    for i in range(df.shape[0]):
        mp_id2str[i] = get_iontype_str(df.iloc[i])

    return mp_id2str

