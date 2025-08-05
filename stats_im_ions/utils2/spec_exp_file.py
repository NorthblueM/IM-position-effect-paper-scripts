# -*- coding: UTF-8 -*-
# @Date     : 15th Jul, 2022
# @Author   : Northblue

"""experimental spectra
"""

import sys
import pandas as pd
import numpy as np

from ion_match_param import Param

import linecache

class SpecExp(object):
    def __init__(self):
        self.title = ''
        self.charge = 0
        self.rt = 0.0 # units are seconds
        self.pepmass = 0.0 # m/z
        self.inten_max = sys.float_info.max
        self.inten_sum = sys.float_info.max
        
        # columns ['mz', 'inten', 'inten_relat']
        self.df_peak = pd.DataFrame()


    def set_spec_mgf(self, fpath, title, idx_title2line={}):
        """load one spec according to title from mgf file
        fpath: fpath_mgf
        idx_title = {title:[start_i, end_i]}
        """

        ls_peak = []
        flag_begin = False # Whether read 'BEGIN IONS'
        flag_end = False # Whether read 'END IONS'

        # # # read too slow # read num=10 cProfile=32s
        # with open(fpath, 'r') as f:
        #     for i, line in enumerate(f):
        #         if i < idx_title2line[title][0]:
        #             continue
        #         if i > idx_title2line[title][1]:
        #             break

        # # # fast read # read num=10 cProfile<1 s
        if not linecache.getline(fpath, 0): # first line is blank
            add = 1
        else:
            add = 0
        start_i = idx_title2line[title][0] + add
        end_i = idx_title2line[title][1] + 1 + add
        for i in list(range(start_i, end_i)):
            line = linecache.getline(fpath, i)
            if True:

                # print(line.strip())
                if 'BEGIN IONS' in line:
                    flag_begin = True
                    ls_peak = [] # initialize
                    continue
                if 'END IONS' in line:
                    flag_end = True
                    self.df_peak = pd.DataFrame(ls_peak, columns=['mz', 'inten'])
                    continue
                if 'TITLE=' in line:
                    self.title = line.split('=')[1].strip()
                    continue
                if 'CHARGE=' in line:
                    self.charge = int(line.split('=')[1].strip()[:-1])
                    continue
                if 'RTINSECONDS=' in line:
                    self.rt = float(line.split('=')[1].strip())
                    continue
                if 'PEPMASS=' in line: # m/z
                    self.pepmass = float(line.split('=')[1].strip())
                    continue
                if not line.isspace() and line.split()[0][0].isdigit():
                    segs = line.strip().split()
                    mz = float(segs[0])
                    inten = float(segs[1])
                    ls_peak.append([mz, inten])
                    continue

        assert flag_begin, 'not read BEGIN IONS'
        assert flag_end, 'not read END IONS'

        self.inten_max = self.df_peak['inten'].max()
        self.inten_sum = self.df_peak['inten'].sum()
    

    def set_spec_mgf_offset(self, fpath, title, idx_title2offset={}):
        """load one spec according to title from mgf file
        fpath: fpath_mgf
        idx_title = {title:[offset_i, line_total_i]}

        use offset but linecache.getline(), because linecache.getline() load full in memory
        """

        ls_peak = []
        flag_begin = False # Whether read 'BEGIN IONS'
        flag_end = False # Whether read 'END IONS'

        offset_i = idx_title2offset[title][0]
        line_total_i = idx_title2offset[title][1]
        with open(fpath, 'r') as f:
            f.seek(offset_i)
            for _ in range(line_total_i):
                line = f.readline()
                # print(line.strip())
                if 'BEGIN IONS' in line:
                    flag_begin = True
                    ls_peak = [] # initialize
                    continue
                if 'END IONS' in line:
                    flag_end = True
                    self.df_peak = pd.DataFrame(ls_peak, columns=['mz', 'inten'])
                    continue
                if 'TITLE=' in line:
                    self.title = line.split('=')[1].strip()
                    continue
                if 'CHARGE=' in line:
                    s = line.split('=')[1].strip()
                    if s.endswith('+') or s.endswith('-'):
                        self.charge = int(s[:-1])
                    continue
                if 'RTINSECONDS=' in line:
                    self.rt = float(line.split('=')[1].strip())
                    continue
                if 'PEPMASS=' in line: # m/z
                    self.pepmass = float(line.split('=')[1].strip())
                    continue
                if not line.isspace() and line.split()[0][0].isdigit():
                    segs = line.strip().split()
                    mz = float(segs[0])
                    inten = float(segs[1])
                    ls_peak.append([mz, inten])
                    continue

        assert flag_begin, 'not read BEGIN IONS'
        assert flag_end, 'not read END IONS'

        self.inten_max = self.df_peak['inten'].max()
        self.inten_sum = self.df_peak['inten'].sum()


    def set_inten_relat(self, mode='max'):
        """add column 'intensity of relative' of 'df_peak'
        """

        if mode == 'max':
            self.inten_max = self.df_peak['inten'].max()
            self.df_peak['inten_relat'] = self.df_peak['inten']/self.inten_max
        if mode == 'sum':
            self.inten_sum = self.df_peak['inten'].sum()
            self.df_peak['inten_relat'] = self.df_peak['inten']/self.inten_sum
        if mode == 'max_sqrt':
            self.inten_max_sqrt = np.sqrt(self.df_peak['inten'].max())
            self.df_peak['inten_relat'] = np.sqrt(self.df_peak['inten'])/self.inten_max_sqrt


def get_mgf_title2line_idx(fpath):
    """
    Return: {title:[start_i, end_i]}
        start_i: i at 'BEGIN IONS'
        end_i: i at 'END IONS'
    """

    idx_title2line={}
    with open(fpath, 'r') as f:
        for i, line in enumerate(f):
            if 'BEGIN IONS' in line:
                start_i = i
                continue
            if 'END IONS' in line:
                end_i = i
                idx_title2line[title] = [start_i, end_i]
                continue
            if 'TITLE=' in line:
                title = line.split('=')[1].strip()
                continue

    return idx_title2line


def get_mgf_title2offset_idx(fpath):
    """
    Return: {title:[offset_i, line_total_i]}
        offset_i: i at 'BEGIN IONS'
        line_total_i: read line total, end at 'END IONS'
    """

    idx_title2offset = {}
    line_total = 0 # the total of line including 'BEGIN IONS' and 'END IONS'
    offset = 0
    with open(fpath, 'r') as f:
        for i, line in enumerate(f):
            offset_line = len(line)
            if line.endswith('\n'): # +1 is '\n','\r\n'
                offset_line += 1
            offset += offset_line
            line_total += 1

            if 'BEGIN IONS' in line:
                offset_i = offset - offset_line
                line_total = 1
                continue
            if 'END IONS' in line:
                idx_title2offset[title] = [offset_i, line_total]
                continue
            if 'TITLE=' in line:
                title = line.split('=')[1].strip()
                continue
    
    return idx_title2offset


def get_mgf_order2offset_idx(fpath):
    """
    Return: {title:[offset_i, line_total_i]}
        offset_i: i at 'BEGIN IONS'
        line_total_i: read line total, end at 'END IONS'
    """

    idx_title2offset = {}
    line_total = 0 # the total of line including 'BEGIN IONS' and 'END IONS'
    offset = 0
    order = 0
    with open(fpath, 'r') as f:
        for i, line in enumerate(f):
            offset_line = len(line) 
            # if line.endswith('\n'): # +1 is '\n','\r\n'
                # print(line)
                # offset_line += 1
            offset += offset_line
            line_total += 1

            if 'BEGIN IONS' in line:
                offset_i = offset - offset_line
                line_total = 1
                continue
            if 'END IONS' in line:
                idx_title2offset[order] = [offset_i, line_total]
                order += 1
                continue
    
    return idx_title2offset
