# -*- coding: UTF-8 -*-
# @Date     : 15th Jul, 2022
# @Author   : Northblue

"""Chemical crosslinking agent file
"""

import pandas as pd

class Linker(object):

    def __init__(self):

        # columns see read_plink_xlink()
        self.df_linker = pd.DataFrame()

    def set_id(self):
        """id start at 0
        """
        self.df_linker['id'] = list(range(self.df_linker.shape[0]))

    def reset_index(self):
        self.df_linker = self.df_linker.reset_index(drop=True)

    def read_plink_xlink(self, fpath):
        """read ./plink/bin/xlink.ini
        """

        cols = ['name', 'alpha_sites', 'beta_sites', 'linker_mass', 'mono_mass', 'linker_composition', 'mono_composition', 'ms_cleavable', 'long_mass', 'short_mass']
        ls = []
        name = ''
        with open(fpath, 'r') as f:
            for line in f:
                if line[:4] == 'name':
                    name = line.split('=')[1].strip()
                if name and line[:len(name)] == name:
                    segs = line.split('=')[1].strip().split()
                    alpha_sites = segs[0]
                    beta_sites = segs[1]
                    linker_mass = float(segs[2])
                    mono_mass = float(segs[4])
                    linker_composition = segs[6]
                    mono_composition = segs[7]
                    ms_cleavable = int(segs[8])
                    long_mass = float(segs[9])
                    short_mass = float(segs[10])

                    linker_one = [name, alpha_sites, beta_sites, linker_mass, mono_mass, linker_composition, mono_composition, ms_cleavable, long_mass, short_mass]

                    ls.append(linker_one)

        self.df_linker = pd.DataFrame(ls, columns=cols)
        self.reset_index()