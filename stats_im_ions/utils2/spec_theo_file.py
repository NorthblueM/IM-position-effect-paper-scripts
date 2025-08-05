# -*- coding: UTF-8 -*-
# @Date     : 15th Jul, 2022
# @Author   : Northblue

"""calculate the theoretical spectrum
"""

import pandas as pd
from ion_match_param import Param

import iontype_file

# https://stackoverflow.com/questions/78813431/numpy-numbers-output-with-type-cast-around-them-which-is-a-bug-in-my-program
# 写出tuple(ser[['iontype_str', 'pepid', 'site1', 'site2']])不带np.int64()
import numpy as np
np.set_printoptions(legacy="1.25")


IL_LOSS = {'I': 29.039123, 'L': 43.0547722}

class SpecTheo(object):

    def __init__(self, param=Param()):

        # ['iontype_id', 'pepid', 'site1', 'site2', 'mz', 'iontype_str']
        # see get_df_ion_col()
        self.df_ion = pd.DataFrame()

        self.df_iontype = pd.DataFrame() # depend on the linker
        self.df_linker = pd.DataFrame()
        self.aa2mass = {}
        self.mod2mass = {}

        ## Initialization parameters, include mass
        self.config_param = param
        self.set_from_param()


    def set_from_param(self):
        self.aa2mass = self.config_param.aa2mass
        self.mod2mass = self.config_param.mod2mass
        self.df_linker = self.config_param.df_linker


    def update_mass_labeled(self, label_id):
        # to do: 如果label_id一样，不更新
        self.aa2mass = self.config_param.aa2mass_labeled[label_id]


    def get_aa_mass_list_nterm(self, seq, mods):
        """ modifictaion site start at 0
        ls: from nterm to cterm, 0->len(seq)
        """
        ls = []
        mass = 0.0
        for i, aa in enumerate(seq):
            if i in mods:
                mass += self.aa2mass[aa]+self.mod2mass[mods[i]]
            else:
                mass += self.aa2mass[aa]
            ls.append(mass)
        return ls


    def get_aa_mass_list_cterm(self, seq, mods):
        """cterm ion, for x/y/z ions
        ls: from nterm to cterm, have reversed
        """
        ls = []
        mass = 0.0
        for i in range(len(seq)-1, -1, -1):
            if i in mods:
                mass += self.aa2mass[seq[i]]+self.mod2mass[mods[i]]
            else:
                mass += self.aa2mass[seq[i]]
            ls.append(mass)
        # reverse, from cterm to nterm -> from nterm to cterm
        ls.reverse()
        return ls
    

    def get_aa_mass_ls_nterm_modnc(self, seq, mods):
        """ nterm ion, for a/b/c ions
        
        modifictaion site start at 1
        n-term modifictaion site at 0
        c-term modifictaion site at len(seq)+1

        ls: from nterm to cterm, 0->len(seq)
        """
        ls = []
        mass = 0.0

        if 0 in mods:
            mass += self.mod2mass[mods[0]]

        for i, aa in enumerate(seq):
            mass += self.aa2mass[aa]
            s = i + 1
            if s in mods: mass += self.mod2mass[mods[s]]
            ls.append(mass)
        
        if len(seq)+1 in mods:
            ls[-1] += self.mod2mass[mods[len(seq)+1]]

        return ls


    def get_aa_mass_ls_cterm_modnc(self, seq, mods):
        """ cterm ion, for x/y/z ions

        ls: from nterm to cterm, have reversed
        """
        ls = []
        mass = 0.0

        if len(seq)+1 in mods:
            mass += self.mod2mass[mods[len(seq)+1]]

        for i in range(len(seq)-1, -1, -1):
            mass += self.aa2mass[seq[i]]
            s = i + 1
            if s in mods: mass += self.mod2mass[mods[s]]
            ls.append(mass)

        if 0 in mods:
            ls[-1] += self.mod2mass[mods[0]]

        # reverse, from cterm to nterm -> from nterm to cterm
        ls.reverse()
        return ls


    def get_residues_mass(self, seq, mods):
        """not include H2O
        """
        mass = 0.0
        for aa in seq:
            mass += self.aa2mass[aa]
        for mod in mods.values():
            mass += self.mod2mass[mod]
        return mass


    def get_pep_mass(self, seq, mods):
        """include H2O
        """
        mass = self.get_residues_mass(seq, mods)
        mass += self.config_param.mass_H2O
        return mass


    # def cal_mz(self, mass, charge, charge_mode='+'):
    def cal_mz(self, mass, charge):
        """
        Args:
            mass: uncharged mass
        """

        # Positive or negative mode
        charge_mode = self.config_param.charge_mode

        mass /= charge
        if charge_mode == '+':
            mass += self.config_param.mass_proton
        if charge_mode == '-':
            mass -= self.config_param.mass_proton
        return mass


    def get_df_ion_col(self):
        """self.df_ion.columns
        """
        # dissociation pep_id: 
            # 1, alpha;
            # 2, beta;
            # 3, dissociative twice in two pep, site1 in alpha, site2 in beta
        # site:
            # start at 0
            # default, -1
            # if site1!=-1 & site2!=-1, dissociative twice
        
            # c-term frag ion site alse from n to c
            # n-term site align with c-term
        
        # ---linear peptide---
            # pepid:
                # 0, small molecule
                # 1, peptide fragment
            # frag site:
                # start at 0 from n to c;
                # b/y frag site align with c-term
                # site1=-1 & site2=-1, precursor
        
        return ['iontype_id', 'pepid', 'site1', 'site2', 'mz']


    def set_df_iontype_str(self):
        """add column 'iontype_str'
        """

        # # reload iontype_str, slow, num=100, cProfile=17s/40s
        # mp_iontype_id2str = iontype_file.get_iontype_id2str(self.df_iontype)
        # self.df_ion['iontype_str'] = self.df_ion['iontype_id'].apply(lambda x:mp_iontype_id2str[x])

        # # fast opt
        if 'iontype_str' not in self.df_iontype.columns:
            mp_iontype_id2str = iontype_file.get_iontype_id2str(self.df_iontype)
            self.df_iontype['iontype_str'] = list(mp_iontype_id2str.values())
        self.df_ion['iontype_str'] = self.df_ion['iontype_id'].apply(lambda x:self.df_iontype.iloc[x]['iontype_str'])


    def set_df_iontype_charge(self):
        """add column 'charge'
        """

        self.df_ion['charge'] = self.df_ion['iontype_id'].apply(lambda x:self.df_iontype.iloc[x]['charge'])



class SpecTheoXl(SpecTheo):

    def __init__(self):
        super().__init__()
        self.seq1 = ''
        self.seq2 = ''
        self.link1 = 0 # start at 0 
        self.link2 = 0
        self.mod1 = {} # {mod_site: mod_name} # site start at 0
        self.mod2 = {}
        self.linker = ''
        self.ser_linker = pd.Series()

        self.mass1_nterm = []
        self.mass2_nterm = []
        self.mass1_cterm = []
        self.mass2_cterm = []

        self.mass_pep1 = 0.0
        self.mass_pep2 = 0.0
        self.mass_pep = 0.0

    def link_mod_site_subtract1(self):
        """plink format start at 1, site-1, link_site/mod_site start at 0
        """
        if self.link1 > 0:
            self.link1 -= 1
        if self.link2 > 0:
            self.link2 -= 1
        
        mod_new = {}
        for site, name in self.mod1.items():
            if site > 0:
                mod_new[site-1] = name
        self.mod1 = mod_new

        mod_new = {}
        for site, name in self.mod2.items():
            if site > 0:
                mod_new[site-1] = name
        self.mod2 = mod_new


    def read_df_iontype(self, fpath):
        """custom ion types
        """

        iontype = iontype_file.IonType()
        iontype.read_iontype(fpath)
        self.df_iontype = iontype.df_iontype

    def update_linker(self, linker):
        """更新交联剂, 主要是交联剂质量"""
        self.linker = linker
        ser_linker = self.df_linker[self.df_linker['name']==linker]
        ser_linker = pd.Series(list(ser_linker.iloc[0]), index=ser_linker.columns)

        self.ser_linker = ser_linker


    def update_df_iontype(self, linker):
        """plink format linker, based on ion type templates
        """
        ser_linker = self.df_linker[self.df_linker['name']==linker]
        ser_linker = pd.Series(list(ser_linker.iloc[0]), index=ser_linker.columns)
        self.ser_linker = ser_linker
        self.linker = linker

        iontype = iontype_file.IonType()
        iontype.set_xl_ion(ser_linker)
        self.df_iontype = iontype.df_iontype


    def set_pep_linker_iontype_plink(self, pep, mods, linker):
        """plink output format pepitde & linker, set ion types"""

        self.set_pep_plink(pep, mods)

        if linker != self.linker:
            self.linker = linker
            self.update_df_iontype(self.linker)


    def set_pep_plink(self, pep, mods):
        """plink output format, initialize the peptide"""

        self.seq1 = pep.split('-')[0].split('(')[0]
        self.seq2 = pep.split('-')[1].split('(')[0]
        
        self.link1 = int(pep.split('-')[0].split('(')[1][:-1])
        self.link2 = int(pep.split('-')[1].split('(')[1][:-1])

        self.mod1 = {} # must be initialized, otherwise the previous modification will not be cleared
        self.mod2 = {}
        if not pd.isnull(mods):
            segs = mods.split(';')
            for seg in segs:
                mod_name = seg.split('(')[0]
                mod_site = int(seg.split('(')[1][:-1])
                if mod_site > len(self.seq1):
                    self.mod2[mod_site-3-len(self.seq1)] = mod_name
                else:
                    self.mod1[mod_site] = mod_name

        self.link_mod_site_subtract1() 


    def cal_theo_pep_mass_xl(self):
        """计算交联肽质量"""

        self.mass_pep1 = self.get_pep_mass(self.seq1, self.mod1)
        self.mass_pep2 = self.get_pep_mass(self.seq2, self.mod2)

        self.mass_pep = self.mass_pep1 + self.mass_pep2 + self.ser_linker['linker_mass']



    def cal_theo_mz(self):

        # initialize
        self.mass1_nterm = self.get_aa_mass_list_nterm(self.seq1, self.mod1)
        self.mass1_cterm = self.get_aa_mass_list_cterm(self.seq1, self.mod1)
        self.mass2_nterm = self.get_aa_mass_list_nterm(self.seq2, self.mod2)
        self.mass2_cterm = self.get_aa_mass_list_cterm(self.seq2, self.mod2)
        self.mass_pep1 = self.get_pep_mass(self.seq1, self.mod1)
        self.mass_pep2 = self.get_pep_mass(self.seq2, self.mod2)

        ls_ion = []
        ls_ion.extend(self.cal_theo_mz_frag1())
        ls_ion.extend(self.cal_thero_mz_characteristic1())

        self.df_ion = pd.DataFrame(ls_ion, columns=self.get_df_ion_col())


    def cal_theo_mz_frag1(self):
        """ions that dissociate once

        not included ions of complete peptides
        """
        ls = []

        for i in range(self.df_iontype.shape[0]):
            ser_iontype = self.df_iontype.iloc[i]
            if ser_iontype['characteristic_ion'] != 0:
                continue
            link_mode = ser_iontype['link_mode']
            ion_term = ser_iontype['term']
            mass_shift = ser_iontype['mass_shift']
            charge = ser_iontype['charge']
            link_mass = ser_iontype['link_mass']
            cleav_mass = ser_iontype['cleav_mass']
            
            if link_mode == 0: # without linker
                if ion_term == 0: # n-term
                    # alpha
                    for j in range(0, self.link1):
                        # seq: ABC, mass_nterm idx=0 between AB
                        mass = self.mass1_nterm[j] + mass_shift
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 1, j, -1, mz]
                        ls.append(ion)
                    # beta
                    for j in range(0, self.link2):
                        mass = self.mass2_nterm[j] + mass_shift
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 2, j, -1, mz]
                        ls.append(ion)
                elif ion_term == 1: # c-term
                    # alpha
                    for j in range(self.link1+1, len(self.seq1)):
                        # seq: ABC, mass_cterm idx=0 before A, idx=1 between AB
                        mass = self.mass1_cterm[j] + mass_shift
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 1, j-1, -1, mz]
                        ls.append(ion)
                    # beta
                    for j in range(self.link2+1, len(self.seq2)):
                        mass = self.mass2_cterm[j] + mass_shift
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 2, j-1, -1, mz]
                        ls.append(ion)

            elif link_mode == 1: # with non-cleav linker;
                if ion_term == 0: # n-term
                    # alpha
                    for j in range(self.link1, len(self.seq1)-1):
                        # seq: ABC, mass_nterm idx=len(seq)-1 after C
                        mass = self.mass1_nterm[j] + mass_shift + link_mass + self.mass_pep2
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 1, j, -1, mz]
                        ls.append(ion)
                    # beta
                    for j in range(self.link2, len(self.seq2)-1):
                        mass = self.mass2_nterm[j] + mass_shift + link_mass + self.mass_pep1
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 2, j, -1, mz]
                        ls.append(ion)
                elif ion_term == 1: # c-term
                    # alpha
                    for j in range(1, self.link1+1):
                        # seq: ABC, mass_cterm idx=0 before A, idx=1 between AB
                        mass = self.mass1_cterm[j] + mass_shift + link_mass + self.mass_pep2
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 1, j-1, -1, mz]
                        ls.append(ion)
                    # beta
                    for j in range(1, self.link2+1):
                        mass = self.mass2_cterm[j] + mass_shift + link_mass + self.mass_pep1
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 2, j-1, -1, mz]
                        ls.append(ion)

            elif link_mode == 3: # with cleav linker;
                if ion_term == 0: # n-term
                    # alpha
                    for j in range(self.link1, len(self.seq1)-1):
                        # seq: ABC, mass_nterm idx=len(seq)-1 after C
                        mass = self.mass1_nterm[j] + mass_shift + cleav_mass
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 1, j, -1, mz]
                        ls.append(ion)
                    # beta
                    for j in range(self.link2, len(self.seq2)-1):
                        mass = self.mass2_nterm[j] + mass_shift + cleav_mass
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 2, j, -1, mz]
                        ls.append(ion)
                elif ion_term == 1: # c-term
                    # alpha
                    for j in range(1, self.link1+1):
                        # seq: ABC, mass_cterm idx=0 before A, idx=1 between AB
                        mass = self.mass1_cterm[j] + mass_shift + cleav_mass
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 1, j-1, -1, mz]
                        ls.append(ion)
                    # beta
                    for j in range(1, self.link2+1):
                        mass = self.mass2_cterm[j] + mass_shift + cleav_mass
                        mz = self.cal_mz(mass, charge)
                        ion = [i, 2, j-1, -1, mz]
                        ls.append(ion)
        
        return ls
    
    def cal_thero_mz_characteristic1(self):
        """cleavable characteristic ion
        """

        ls = []

        for i in range(self.df_iontype.shape[0]):
            ser_iontype = self.df_iontype.iloc[i]
            if ser_iontype['characteristic_ion'] != 1:
                continue
            mass_shift = ser_iontype['mass_shift']
            charge = ser_iontype['charge']
            cleav_mass = ser_iontype['cleav_mass']
            
            # alpha
            mass = self.mass_pep1 + cleav_mass
            mz = self.cal_mz(mass, charge)
            ion = [i, 1, -1, -1, mz]
            ls.append(ion)

            # beta
            mass = self.mass_pep2 + cleav_mass
            mz = self.cal_mz(mass, charge)
            ion = [i, 2, -1, -1, mz]
            ls.append(ion)

        return ls



class SpecTheoLinearIL(SpecTheo):

    def __init__(self):
        super().__init__()
        self.seq = ''
        self.mod = {} # {mod_site: mod_name} # site start at 0, pFind修饰N端为0, C端修饰为长度+1

        self.mass_nterm = []
        self.mass_cterm = []

        self.mass_pep = 0.0


    def read_df_iontype(self, fpath):
        """custom ion types
        """

        iontype = iontype_file.IonType()
        iontype.read_iontype(fpath)
        self.df_iontype = iontype.df_iontype


    def set_pep_pfind(self, pep, mods):
        """pfind output format, initialize the peptide"""

        self.seq = pep

        self.mod = {} # must be initialized, otherwise the previous modification will not be cleared
        if (pd.isnull(mods)) or (mods=='null') or (not mods):
            pass
        else:
            segs = mods.split(';')[:-1]
            for seg in segs:
                ss = seg.strip().split(',')
                mod_name = ss[1]
                mod_site = int(ss[0])
                self.mod[mod_site] = mod_name


    def get_df_ion_col_il(self):
        """区分i和L的ion计算column name
        """

        return  ['iontype_id', 'pepid', 'site1', 'site2', 'mz', 'loss_str']


    def cal_theo_mz(self):

        # initialize
        self.mass_nterm = self.get_aa_mass_ls_nterm_modnc(self.seq, self.mod)
        self.mass_cterm = self.get_aa_mass_ls_cterm_modnc(self.seq, self.mod)
        self.mass_pep = self.get_pep_mass(self.seq, self.mod)

        ls_ion = []
        ls_ion.extend(self.cal_theo_mz_il_loss())

        self.df_ion = pd.DataFrame(ls_ion, columns=self.get_df_ion_col_il())

    def comb_il(self, num):
        """num个I和L的组合
        """

        ls = [{'I':0, 'L':0}] # 默认添加中性丢失I=0,L=0
        for i in range(1, num+1): # I+L总数
            for j in range(i+1): 
                ls.append({'I':j, 'L':i-j})
        return ls


    def cal_theo_mz_il_loss(self):
        """仅包含I、L的离子

        还有他们的中性丢失离子组合
        """
        ls = []

        for i in range(self.df_iontype.shape[0]):
            ser_iontype = self.df_iontype.iloc[i]

            ion_term = ser_iontype['term']
            mass_shift = ser_iontype['mass_shift']
            charge = ser_iontype['charge']

            if ion_term == 0: # n-term
                assert False, 'not support n-term ion'

            elif ion_term == 1: # c-term, 碎裂位点0~(seq-2)


                # 一般匹配，从1开始, 是和b离子的位点对齐，碎裂位点从0计数，最大为len(seq)-2

                # 碎裂位置0表示母离子的z离子，根据机理不会形成
                # 但看到文章有含两个IL的离子，分别有I和L的丢失谱峰，所以还是计算

                # 0表示母离子，是y型，不是z型
                # for j in range(1, len(self.seq)):
                for j in range(0, len(self.seq)):
                    
                    
                    # C端质量表反转了
                    # print('\n')
                    # print(self.seq, len(self.seq), j)

                    seq_sub = self.seq[j:]
                    num_il = seq_sub.count('I') + seq_sub.count('L')

                    # print(seq_sub, len(seq_sub), num_il)
                    if num_il == 0:
                        continue

                    # # 根据机理只有z离子上会有IL的丢失
                    # if seq_sub[-1] not in ['I', 'L']:
                        # continue
                    # num_il = 1 # 只考虑断裂C端，即右位置的I/L

                    mass = self.mass_cterm[j] + mass_shift
                    for il2num in self.comb_il(num_il):
                        
                        # 断裂位点
                        # s = j - 1
                        s = j # 0表示y型母离子

                        m = mass - IL_LOSS['I']*il2num['I'] - IL_LOSS['L']*il2num['L']

                        # if s == 0: # 母离子不会形成z离子，当y离子（母离子）处理
                            # m += 18.010565 - 1.9918409 + 1.00782503214 # 得电子
                        # 0号位点，即肽N端，即母离子，也假设会形成z离子，统一处理

                        mz = self.cal_mz(m, charge)
                        ion = [i, 1, s, -1, mz, f"I-{il2num['I']}_L-{il2num['L']}"]
                        ls.append(ion)
        return ls


class SpecTheoLinear(SpecTheo):

    def __init__(self):
        super().__init__()
        self.seq = ''
        self.mod = {} # {mod_site: mod_name} # site start at 0, pFind修饰N端为0, C端修饰为长度+1

        self.mass_nterm = []
        self.mass_cterm = []

        self.mass_pep = 0.0 # 残基和，带修饰，带H2O


    def read_df_iontype(self, fpath):
        """custom ion types
        """

        iontype = iontype_file.IonType()
        iontype.read_iontype(fpath)
        self.df_iontype = iontype.df_iontype


    def set_pep_pfind(self, pep, mods):
        """pfind output format, initialize the peptide"""

        self.seq = pep

        self.mod = {} # must be initialized, otherwise the previous modification will not be cleared
        if (pd.isnull(mods)) or (mods=='null') or (not mods):
            pass
        else:
            segs = mods.split(';')[:-1]
            for seg in segs:
                ss = seg.strip().split(',')
                mod_name = ss[1]
                mod_site = int(ss[0])
                self.mod[mod_site] = mod_name

        
    def cal_theo_mz(self):

        # initialize
        self.mass_nterm = self.get_aa_mass_ls_nterm_modnc(self.seq, self.mod)
        self.mass_cterm = self.get_aa_mass_ls_cterm_modnc(self.seq, self.mod)
        self.mass_pep = self.get_pep_mass(self.seq, self.mod)

        ls_ion = []
        ls_ion.extend(self.cal_theo_mz_frag1())
        ls_ion.extend(self.cal_theo_mz_characteristic_small())

        self.df_ion = pd.DataFrame(ls_ion, columns=self.get_df_ion_col())


    def cal_theo_mz_frag1(self):
        """ions that dissociate once

        not included ions of complete peptides
        """
        ls = []

        for i in range(self.df_iontype.shape[0]):
            ser_iontype = self.df_iontype.iloc[i]
            if ser_iontype['characteristic_ion'] != 0:
                continue
            ion_term = ser_iontype['term']
            mass_shift = ser_iontype['mass_shift']
            charge = ser_iontype['charge']
            
            if ion_term == 0: # n-term, b离子
                for j in range(0, len(self.seq)-1):
                    mass = self.mass_nterm[j] + mass_shift
                    mz = self.cal_mz(mass, charge)
                    ion = [i, 1, j, -1, mz]
                    ls.append(ion)
            
            # 0表示母离子, 断裂位点记录是-1
            elif ion_term == 1: # c-term, y离子
                for j in range(0, len(self.seq)):
                    mass = self.mass_cterm[j] + mass_shift
                    mz = self.cal_mz(mass, charge)
                    ion = [i, 1, j-1, -1, mz]
                    ls.append(ion)

        return ls
    
    def cal_theo_mz_characteristic_small(self):
        """small molecule characteristic ion, not include peptide
        """

        ls = []

        for i in range(self.df_iontype.shape[0]):
            ser_iontype = self.df_iontype.iloc[i]
            if ser_iontype['characteristic_ion'] != 2: # small molecule
                continue
            mass_shift = ser_iontype['mass_shift']
            charge = ser_iontype['charge']

            mass = mass_shift
            mz = self.cal_mz(mass, charge)
            ion = [i, 0, -1, -1, mz]
            ls.append(ion)

        return ls




def get_iontheo_str(ser, iontype_str=False, match_type=''):
    """ser: one theo ion
    return example:
        ('iontype_id', 'pepid', 'site1', 'site2')

    Args:
        match_type: 匹配的类型
            ETD-IL: ETD谱图用于区分IL, 只匹配含有IL的子串、以及IL的丢失组合
    """

    if iontype_str and match_type == 'ETD-IL':
        return tuple(ser[['iontype_str', 'pepid', 'site1', 'site2', 'loss_str']].tolist())

    if iontype_str:
        # return list(ser[['iontype_id', 'iontype_str', 'pepid', 'site1', 'site2']].tolist())
        return tuple(ser[['iontype_str', 'pepid', 'site1', 'site2']].tolist())
    else:
        return tuple(ser[['pepid', 'site1', 'site2']].tolist())