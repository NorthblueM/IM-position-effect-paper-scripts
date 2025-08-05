# -*- coding: UTF-8 -*-
# Date      : 22th Mar, 2021
# Author    : Northblue

""""操作plink结果CSV文件的一些工具函数"""

import os
import sys
sys.path.append(os.path.dirname(os.getcwd()))

import pandas as pd
import numpy as np
from utils import fasta_file

import ahocorasick

from tools import df_apply_mutli as dfm


def get_title_rawscan(x, mode='pParse'):
    """获取title中raw和scan的字符串"""

    return '.'.join(x.split('.')[:-3])

def get_title_pparseid(x):
    """获取title中pparse id, int"""

    return int(x.split('.')[-2])


def get_seq_i2l(seq, is_i2l=True):
    """单肽序列I换为L"""

    if not is_i2l:
        return seq
    else:
        return seq.replace('I', 'L')
        # seq_ls = []
        # for aa in seq:
        #     if aa == 'I':
        #         seq_ls.append('L')
        #     else:
        #         seq_ls.append(aa)
        # return ''.join(seq_ls)



def get_xl_seq1(x, is_i2l=False):
    """获得交联的seq1
    Args:
        x: ['Peptide']
    """
    seq = x.split('-')[0].split('(')[0]
    return get_seq_i2l(seq, is_i2l=is_i2l)



def get_xl_seq2(x, is_i2l=False):
    """获得交联的seq2
    Args:
        x: ['Peptide']
    """
    seq = x.split('-')[1].split('(')[0]
    return get_seq_i2l(seq, is_i2l=is_i2l)



def get_xl_link1(x):
    """获得交联的link1，从1开始
    Args:
        x: ['Peptide']
    """
    return int(x.split('-')[0].split('(')[1][:-1])



def get_xl_link2(x):
    """获得交联的link2
    Args:
        x: ['Peptide']
    """
    return int(x.split('-')[1].split('(')[1][:-1])


def get_xl_link_char(x, seg_i, is_i2l=False):
    """获得交联位点的字母
    没有区分蛋白与肽段N端、C端

    Args:
        x: ['Peptide']
    """
    link = int(x.split('-')[seg_i].split('(')[1][:-1])
    seq = x.split('-')[seg_i].split('(')[0]
    seq = get_seq_i2l(seq, is_i2l=is_i2l)

    return seq[link-1]


def get_xl_link1_char(x, is_i2l=False):
    """获得交联位点的字母
    没有区分蛋白与肽段N端、C端

    Args:
        x: ['Peptide']
    """
    return get_xl_link_char(x, 0, is_i2l=is_i2l)


def get_xl_link2_char(x, is_i2l=False):
    """获得交联位点的字母
    没有区分蛋白与肽段N端、C端

    Args:
        x: ['Peptide']
    """
    return get_xl_link_char(x, 1, is_i2l=is_i2l)


def get_xl_seq_seq_str(x, is_i2l=True):
    """获得序列对字符串表示，序列、没有位点
    序列：按照字母表排序
    Returns:
        seq1-seq2
    """

    pep = x['Peptide']

    seq1 = get_xl_seq1(pep, is_i2l=is_i2l)
    seq2 = get_xl_seq2(pep, is_i2l=is_i2l)

    return '-'.join(sorted([seq1, seq2]))



def get_xl_seq_str(x, is_i2l=True):
    """获得序列对字符串表示，序列、位点
    序列：按照字母表排序
    Returns:
        seq1(link1)-seq2(link2)
    """

    pep = x['Peptide']

    seq1 = get_xl_seq1(pep, is_i2l=is_i2l)
    seq2 = get_xl_seq2(pep, is_i2l=is_i2l)
    link1 = get_xl_link1(pep)
    link2 = get_xl_link2(pep)

    # 序列按照字母表排序
    # fix bug 如果序列相同，位点不同将被覆盖
    if seq1 != seq2:
        seq2link = {seq1:link1, seq2:link2}
        seq1, seq2 = sorted([seq1, seq2])
        link1 = seq2link[seq1]
        link2 = seq2link[seq2]
    elif link1 != link2:
        link1, link2 = sorted([link1, link2])

    seq_str = seq1+'('+str(link1)+')-'+seq2+'('+str(link2)+')'

    return seq_str


def get_xl_mod_dict(seq1, mods_str):
    """获取交联的seq对应mod

    Args:
        seq1: 为了获取alpha肽长度
        mods_str: csv中存储的modifications字符串
    Returns:
        mod:修饰的字典
    """

    mod1 = {} # 位点：修饰名称，位点从1开始，N端修饰为0
    mod2 = {} # 
    if not pd.isnull(mods_str) and (not (isinstance(mods_str, (str)) and mods_str=='null')):
        segs = mods_str.split(';')
        for seg in segs:

            # 取最后一个括号切分，避免情况Label_13C(6)[R](10)
            mod_name = '('.join(seg.split('(')[:-1])
            mod_site = int(seg.split('(')[-1][:-1])

            if mod_site > len(seq1):
                mod2[mod_site-3-len(seq1)] = mod_name
            else:
                mod1[mod_site] = mod_name

    return mod1, mod2


def get_xl_mod(x, is_i2l=False):
    """获取交联的seq对应mod

    Returns:
        mod:修饰的字典
    """

    pep = x['Peptide']
    mods = x['Modifications']

    seq1 = get_xl_seq1(pep, is_i2l=is_i2l)
    seq2 = get_xl_seq2(pep, is_i2l=is_i2l)

    mod1, mod2 = get_xl_mod_dict(seq1, mods)

    return mod1, mod2


def set_xl_mod(df):
    """加入两列mod1和mod2, 分别表示seq1和seq2的上的修饰, 字典形式
    """

    df['mod1'] = df.apply(lambda x:get_xl_mod(x)[0], axis=1, result_type='reduce')
    df['mod2'] = df.apply(lambda x:get_xl_mod(x)[1], axis=1, result_type='reduce')

    return df


def get_modcnt_str(site2mod):
    if not site2mod:
        mod_str = '0*null'
    else:
        mod2cnt = {}
        for mod_site, mod_name in site2mod.items():
            if mod_name in mod2cnt:
                mod2cnt[mod_name] += 1
            else:
                mod2cnt[mod_name] = 1
        # 按照修饰名称排序
        mod2cnt = dict(sorted(mod2cnt.items(), key=lambda x:x[0],reverse=False))
        mod_str = ';'.join([str(mod_cnt)+'*'+mod_name for mod_name, mod_cnt in mod2cnt.items()])
    return mod_str


def get_xl_seq_seqmod_str(x, is_i2l=True):
    """获得序列+修饰对字符串表示，修饰、序列
    修饰不包含位点, 只包括种类和个数, alpha和beta分开统计
    序列不包含交联位点
    序列：按照字母表排序
    修饰：先按肽排序，再按照修饰名称排序
    Returns:
        seq1-seq2_mod1(cnt1)-mod2(cnt2)
    """

    pep = x['Peptide']
    mods = x['Modifications']

    seq1 = get_xl_seq1(pep, is_i2l=is_i2l)
    seq2 = get_xl_seq2(pep, is_i2l=is_i2l)

    mod1, mod2 = get_xl_mod_dict(seq1, mods)

    # 序列按照字母表排序
    if seq1 != seq2: # 序列不同，按照序列字母排序
        seq2mod = {seq1:mod1, seq2:mod2}
        seq1, seq2 = sorted([seq1, seq2])
        mod1 = seq2mod[seq1]
        mod2 = seq2mod[seq2]
    else: # 序列和交联位点均相同，按照修饰种类和个数字符串排序
        mod1_str = get_modcnt_str(mod1)
        mod2_str = get_modcnt_str(mod2)
        if mod1_str != mod2_str:
            mod_str2mod = {mod1_str:mod1, mod2_str:mod2}
            mod1_str, mod2_str = sorted([mod1_str, mod2_str])
            mod1 = mod_str2mod[mod1_str]
            mod2 = mod_str2mod[mod2_str]

    mod1_str = get_modcnt_str(mod1)
    mod2_str = get_modcnt_str(mod2)

    seq_str = seq1+'-'+seq2
    mod_str = '-'.join([mod1_str, mod2_str])
    pep_str = '_'.join([seq_str, mod_str])

    return pep_str


def get_xl_seqmod_str(x, is_i2l=True):
    """获得序列+修饰对字符串表示，修饰、序列、位点
    修饰不包含位点, 只包括种类和个数, alpha和beta分开统计
    序列：按照字母表排序
    修饰：先按肽排序，再按照修饰名称排序
    Returns:
        seq1(link1)-seq2(link2)_mod1(cnt1)-mod2(cnt2)
    """

    pep = x['Peptide']
    mods = x['Modifications']

    seq1 = get_xl_seq1(pep, is_i2l=is_i2l)
    seq2 = get_xl_seq2(pep, is_i2l=is_i2l)
    link1 = get_xl_link1(pep)
    link2 = get_xl_link2(pep)

    mod1, mod2 = get_xl_mod_dict(seq1, mods)

    # 序列按照字母表排序
    if seq1 != seq2: # 序列不同，按照序列字母排序
        seq2link = {seq1:link1, seq2:link2}
        seq2mod = {seq1:mod1, seq2:mod2}
        seq1, seq2 = sorted([seq1, seq2])
        link1 = seq2link[seq1]
        link2 = seq2link[seq2]
        mod1 = seq2mod[seq1]
        mod2 = seq2mod[seq2]
    elif link1 != link2: # 序列相同，按照交联位点大小排序
        link2mod = {link1:mod1, link2:mod2}
        link1, link2 = sorted([link1, link2])
        mod1 = link2mod[link1]
        mod2 = link2mod[link2]
    else: # 序列和交联位点均相同，按照修饰种类和个数字符串排序
        mod1_str = get_modcnt_str(mod1)
        mod2_str = get_modcnt_str(mod2)
        if mod1_str != mod2_str:
            mod_str2mod = {mod1_str:mod1, mod2_str:mod2}
            mod1_str, mod2_str = sorted([mod1_str, mod2_str])
            mod1 = mod_str2mod[mod1_str]
            mod2 = mod_str2mod[mod2_str]

    mod1_str = get_modcnt_str(mod1)
    mod2_str = get_modcnt_str(mod2)

    seq_str = seq1+'('+str(link1)+')-'+seq2+'('+str(link2)+')'
    mod_str = '-'.join([mod1_str, mod2_str])
    pep_str = '_'.join([seq_str, mod_str])

    return pep_str


def get_mod_str(site2mod):
    site2mod = dict(sorted(site2mod.items(), key=lambda x:x[0],reverse=False))
    if not site2mod:
        mod_str = 'null'
    else:
        mod_ls = []
        for mod_site, mod_name in site2mod.items():
            mod_ls.append(mod_name + '(' + str(mod_site) + ')')
        mod_str = ';'.join(mod_ls)
    return mod_str

def get_xl_pep_str(x, is_i2l=True):
    """获得肽段对字符串表示，修饰、序列、位点
    序列：按照字母表排序
    修饰：先按肽排序，再按照位点排序
    Returns:
        seq1(link1)-seq2(link2)_mod1(site1);mod2(site2)
    """

    pep = x['Peptide']
    mods = x['Modifications']

    seq1 = get_xl_seq1(pep, is_i2l=is_i2l)
    seq2 = get_xl_seq2(pep, is_i2l=is_i2l)
    link1 = get_xl_link1(pep)
    link2 = get_xl_link2(pep)

    mod1, mod2 = get_xl_mod_dict(seq1, mods)

    # 序列按照字母表排序
    # fix Bug: 同序列时，字典的键重复，可能会互相覆盖
    if seq1 != seq2: # 序列不同，按照序列字母排序
        seq2link = {seq1:link1, seq2:link2}
        seq2mod = {seq1:mod1, seq2:mod2}
        seq1, seq2 = sorted([seq1, seq2])
        link1 = seq2link[seq1]
        link2 = seq2link[seq2]
        mod1 = seq2mod[seq1]
        mod2 = seq2mod[seq2]
    elif link1 != link2: # 序列相同，按照交联位点大小排序
        link2mod = {link1:mod1, link2:mod2}
        link1, link2 = sorted([link1, link2])
        mod1 = link2mod[link1]
        mod2 = link2mod[link2]
    else:
        mod1 = dict(sorted(mod1.items(), key=lambda x:x[0],reverse=False))
        mod2 = dict(sorted(mod2.items(), key=lambda x:x[0],reverse=False))
        if mod1 != mod2: # 序列位点均相同，按照修饰排序
            mod1_str = get_mod_str(mod1)
            mod2_str = get_mod_str(mod2)
            mod_str2mod = {mod1_str:mod1, mod2_str:mod2}
            mod1_str, mod2_str = sorted([mod1_str, mod2_str])
            mod1 = mod_str2mod[mod1_str]
            mod2 = mod_str2mod[mod2_str]

    site2mod = {}
    for mod_site, mod_name in mod1.items():
        site2mod[mod_site] = mod_name
    for mod_site, mod_name in mod2.items():
        site2mod[mod_site+3+len(seq1)] = mod_name

    mod_str = get_mod_str(site2mod)
    seq_str = seq1+'('+str(link1)+')-'+seq2+'('+str(link2)+')'
    pep_str = '_'.join([seq_str, mod_str])

    return pep_str



def get_xl_prec_str(x, is_i2l=True):
    """获得母离子字符串表示，电荷、修饰、序列、位点
    序列：按照字母表排序
    修饰：先按肽排序，再按照位点排序
    Returns:
        seq1(link1)-seq2(link2)_mod1(site1);mod2(site2)_charge
    """

    charge = int(x['Charge'])

    pep_str = get_xl_pep_str(x, is_i2l=is_i2l)
    prec_str = '_'.join([pep_str, str(charge)])

    return prec_str


def get_xl_prec_seq_seqmod_str(x, is_i2l=True):
    """
    不带修饰位点，只修饰数目
    不带交联位点，只交联序列
    """

    charge = int(x['Charge'])
    seq_seqmod_str = get_xl_seq_seqmod_str(x, is_i2l=is_i2l)
    prec_str = '_'.join([seq_seqmod_str, str(charge)])

    return prec_str


def get_xl_prec_run_str(x, is_i2l=True):
    """不同raw的母离子区分对待
    """

    title = x['Title']
    run_str = '.'.join(title.split('.')[:-5])

    prec_str = get_xl_prec_str(x, is_i2l=is_i2l)
    csm_scan_str = '_'.join([prec_str, run_str])

    return csm_scan_str


def get_xl_csm_str(x, is_i2l=True):
    """获得CSM字符串表示，Title、电荷、修饰、序列、位点
    序列：按照字母表排序
    修饰：先按肽排序，再按照位点排序
    Returns:
        seq1(link1)-seq2(link2)_mod1(site1);mod2(site2)_charge_title
    """

    title = x['Title']

    prec_str = get_xl_prec_str(x, is_i2l=is_i2l)
    csm_str = '_'.join([prec_str, title])

    return csm_str


def get_xl_csm_scan_str(x, is_i2l=True):
    """获得CSM字符串表示，Title中去掉pParse编号
    Title中去掉电荷，因为prec_str中包含电荷
    """

    title = x['Title']
    scan_str = '.'.join(title.split('.')[:-3])

    prec_str = get_xl_prec_str(x, is_i2l=is_i2l)
    csm_scan_str = '_'.join([prec_str, scan_str])

    return csm_scan_str


def split_xl_pep_str(df, is_i2l=False, is_linkchar=True):
    """拆分['Peptide']到seq1/seq2/link1/link2/link1_char/link2_char"""

    df['seq1'] = df['Peptide'].apply(lambda x:get_xl_seq1(x, is_i2l=is_i2l))
    df['seq2'] = df['Peptide'].apply(lambda x:get_xl_seq2(x, is_i2l=is_i2l))
    df['link1'] = df['Peptide'].apply(get_xl_link1)
    df['link2'] = df['Peptide'].apply(get_xl_link2)
    if is_linkchar:
        df['link1_char'] = df['Peptide'].apply(lambda x:get_xl_link1_char(x, is_i2l=is_i2l))
        df['link2_char'] = df['Peptide'].apply(lambda x:get_xl_link2_char(x, is_i2l=is_i2l))
    return df

def set_xl_is_same_seq(df):
    """标记是否同序列交联"""

    df['is_sameseq'] = (df['seq1'] == df['seq2']).astype(int)

    return df



def set_xl_level_str(df, is_i2l=True, n_process=0, is_csmscan=False, is_modcnt=False):
    """获得序列对、肽段对、母离子、CSM"""

    if n_process:
        df['seq_seq_str'] = dfm.apply_raw(df, get_xl_seq_seq_str, kwargs={'is_i2l':is_i2l}, n_process=n_process)
        df['seq_str'] = dfm.apply_raw(df, get_xl_seq_str, kwargs={'is_i2l':is_i2l}, n_process=n_process)
        df['pep_str'] = dfm.apply_raw(df, get_xl_pep_str, kwargs={'is_i2l':is_i2l}, n_process=n_process)
        df['prec_str'] = dfm.apply_raw(df, get_xl_prec_str, kwargs={'is_i2l':is_i2l}, n_process=n_process)
        df['csm_str'] = dfm.apply_raw(df, get_xl_csm_str, kwargs={'is_i2l':is_i2l}, n_process=n_process)
        if is_csmscan:
            df['csm_scan_str'] = dfm.apply_raw(df, get_xl_csm_scan_str, kwargs={'is_i2l':is_i2l}, n_process=n_process)
    else:
        df['seq_seq_str'] = df.apply(lambda x:get_xl_seq_seq_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
        df['seq_str'] = df.apply(lambda x:get_xl_seq_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
        df['pep_str'] = df.apply(lambda x:get_xl_pep_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
        df['prec_str'] = df.apply(lambda x:get_xl_prec_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
        df['csm_str'] = df.apply(lambda x:get_xl_csm_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
        if is_csmscan:
            df['csm_scan_str'] = df.apply(lambda x:get_xl_csm_scan_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
        if is_modcnt:
            df['seq_seqmod_str'] = df.apply(lambda x:get_xl_seq_seqmod_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
            df['prec_seq_seqmod_str'] = df.apply(lambda x:get_xl_prec_seq_seqmod_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
    
    return df

def set_xl_level_count(df, levels_at=['seq_seq_str', 'pep_str', 'prec_str'], levels_in=['prec_str', 'csm_str']):
    """设置, 某层次, 在某层次的count
    
    Args:
        levels_at=['prts_str', 'prtsites_str', 'seq_seq_str', 'pep_str', 'prec_str']
        levels_in=['prtsites_str', 'seq_seq_str', 'pep_str', 'prec_str', 'csm_str']
    Returns:
        seq_seq-csm_cnt: seq_seq在csm层次的支撑数目
    
    """

    for lvl_at in levels_at:
        for lvl_in in levels_in:
            col = '@'.join([lvl_in[:-4]+'_cnt', lvl_at[:-4]])
            # col = '-'.join([lvl_at[:-4], lvl_in[:-4]+'_cnt'])

            _df = df[[lvl_at, lvl_in]].drop_duplicates(subset=[lvl_in], keep='first') # lvl_in只保留一行
            _ser = _df[[lvl_at]].apply(pd.Series.value_counts) # df格式

            # df[col] = df[[lvl_at]].replace(_ser) # 访问必须[lvl]，否则'Series' object has no attribute '_replace_columnwise' # 速度太慢了
            _ser = _ser[lvl_at].to_dict() # 转为字典, 速度快
            df[col] = df[lvl_at].transform(lambda x:_ser[x])

    return df


def set_xl_count_ppi_seq_seq2(df):
    """
    
    支持PPI有≥2完全不一样的序列对, 即序列a-b/c-d, 而非a-b/a-c
    简易版, return只有1和2

    准确时需要涉及图论算法, 建图, 去掉最少的点数, 使得图中各点度数为0
    一种算法, 贪心, 从度最大开始, 度相同, 选邻居度和少的(未验证), 例如M形式为2

    除seq为图上的点, 还可以有带电荷和修饰的版本prec_seq_seqmod2, seq_seqmod2
    """

    def get_cnt(x):
        if x.shape[0] == 1:
            return 1
        x = x.drop_duplicates(subset=['seq_seq_str'], keep='first')[['seq1', 'seq2']] #加快速度
        if x.shape[0] == 1:
            return 1
        for i, row1 in enumerate(x.itertuples()):
            for j, row2 in enumerate(x.itertuples()):
                if j <= i:
                    continue
                if row1.seq1 != row2.seq1 and row1.seq1 != row2.seq2 and row1.seq2 != row2.seq1 and row1.seq2 != row2.seq2:
                    return 2
        return 1

    # 给每一行
    dct_cnt = df.groupby(['prts_str']).apply(get_cnt).to_dict()
    df['seq_seq2_cnt@prts'] = df['prts_str'].transform(lambda x:dct_cnt[x])

    return df





def set_rp_is_unique(df):
    """判断RP层, 蛋白质推断是否Unique
    
        PPI是unique, RP不一定, 例如同一个蛋白两个位点
    
    """

    df['is_unique_rp'] = df.apply(lambda x:((len(x['prtsites1'])==1) & (len(x['prtsites2'])==1)), axis=1, result_type='reduce')

    return df

def set_ppi_is_unique(df):
    """判断PPI层, 蛋白质推断是否Unique"""

    df['is_unique_ppi'] = df.apply(lambda x:((len(x['prts1'])==1) & (len(x['prts2'])==1)), axis=1, result_type='reduce')

    return df


def set_xl_is_unique(df):

    df = set_rp_is_unique(df)
    df = set_ppi_is_unique(df)

    return df


def stat_all_level_num_fdr(df, col_true='is_true', method=1, protein_type_ls=[-1], level_ls=[], is_unique_ls=[], prts_level=True):
    """统计结果的各层的数目和计算FDR
    
    Args:
        col_true: 默认正确与否的标注列为is_true
        method: 1,高层结果有一个对即算对
        protein_type_ls = [-1, 1, 2], -1=全集、1=intra、2=inter
    """

    def get_1_level_num_fdr(df, level, is_unique_col='', protein_type=-1):

        is_unique = 0
        num = 0
        num_true = 0
        num_false = 0
        fdr_calc = 0.0

        _df = df.copy(deep=True)

        if protein_type >= 0:
            if protein_type == 1:
                _df = _df[_df['Protein_Type']=='Intra-Protein']
            elif protein_type == 2:
                _df = _df[_df['Protein_Type']=='Inter-Protein']

        if is_unique_col:
            _df = _df[_df[is_unique_col] == 1]
            is_unique = 1
    
        gp = _df.groupby(level)
        num = gp.ngroups
        if col_true in _df.columns:
            if method == 1:
                if num != 0: # 否则
                    num_true = sum(gp.apply(lambda x:(sum(x[col_true])>=1)))
            num_false = num - num_true
            if num != 0: fdr_calc = num_false/num

            one = [level, is_unique, is_unique_col, protein_type, num, num_true, fdr_calc, num_false]
        else:
            one = [level, is_unique, is_unique_col, protein_type, num]

        return one
    
    if prts_level:
        level_ls = ['csm', 'prec', 'pep', 'seq', 'seq_seq', 'prtsites', 'prtsites', 'prts', 'prts']
        is_unique_ls = ['', '', '', '', '', '', 'is_unique_rp', '', 'is_unique_ppi']
    else:
        level_ls = ['csm', 'prec', 'pep', 'seq', 'seq_seq']
        is_unique_ls = ['', '', '', '', '']
    
    level_ls = [x+'_str' for x in level_ls]
    param_ls = list(zip(level_ls, is_unique_ls)) #必须list, zip只能遍历一遍

    # protein_type_ls = [-1, 1, 2]

    ls_ls = []
    for protein_type in protein_type_ls:
        for level, is_unique_col in param_ls:
            one = get_1_level_num_fdr(df, level, is_unique_col, protein_type)
            ls_ls.append(one)
    
    if col_true in df.columns:
        cols = ['level', 'is_unique', 'unique_level', 'protein_type', 'num', 'num_true', 'fdr_calc', 'num_false']
    else:
        cols = ['level', 'is_unique', 'unique_level', 'protein_type', 'num']

    df_stat = pd.DataFrame(ls_ls, columns=cols)

    return df_stat


def stat_all_level_num(df, protein_type_ls=[-1], level_ls=[], is_unique_ls=[], prts_level=True):
    """统计结果的各层的数目
    
    Args:
        col_true: 默认正确与否的标注列为is_true
        protein_type_ls = [-1, 1, 2], -1=全集、1=intra、2=inter
    """

    def get_1_level_num(df, level, is_unique_col='', protein_type=-1):

        is_unique = 0
        num = 0

        _df = df.copy(deep=True)

        if protein_type >= 0:
            if protein_type == 1:
                _df = _df[_df['Protein_Type']=='Intra-Protein']
            elif protein_type == 2:
                _df = _df[_df['Protein_Type']=='Inter-Protein']

        if is_unique_col:
            _df = _df[_df[is_unique_col] == 1]
            is_unique = 1
    
        gp = _df.groupby(level)
        num = gp.ngroups

        one = [level, is_unique, is_unique_col, protein_type, num]

        return one
    
    if not level_ls:
        if prts_level:
            level_ls = ['csm', 'prec', 'pep', 'seq', 'seq_seq', 'prtsites', 'prtsites', 'prts', 'prts']
            is_unique_ls = ['', '', '', '', '', '', 'is_unique_rp', '', 'is_unique_ppi']
        else:
            level_ls = ['csm', 'prec', 'pep', 'seq', 'seq_seq']
            is_unique_ls = ['', '', '', '', '']
    level_ls = [x+'_str' for x in level_ls]
    param_ls = list(zip(level_ls, is_unique_ls))

    # protein_type_ls = [-1, 1, 2]

    ls_ls = []
    for protein_type in protein_type_ls:
        for level, is_unique_col in param_ls:
            one = get_1_level_num(df, level, is_unique_col, protein_type)
            ls_ls.append(one)
    
    df_stat = pd.DataFrame(ls_ls, columns=['level', 'is_unique', 'unique_level', 'protein_type', 'num'])

    return df_stat


def get_linear_seq(x, is_i2l=False):
    """获得单肽的seq
    Args:
        x: ['Peptide']
    """
    seq = x.strip()
    return get_seq_i2l(seq, is_i2l=is_i2l)


def get_mono_seq(x, is_i2l=False):
    """获得mono的seq
    Args:
        x: ['Peptide']
    """
    seq = x.split('(')[0]
    return get_seq_i2l(seq, is_i2l=is_i2l)


def get_loop_seq(x, is_i2l=False):
    """获得mono的seq
    Args:
        x: ['Peptide']
    """
    seq = x.split('(')[0]
    return get_seq_i2l(seq, is_i2l=is_i2l)


def set_line_seq(df, is_i2l=False):
    """获得单肽/mono/loop的序列"""

    def get_line_seq(x, is_i2l=is_i2l):
        pep_type = int(x['Peptide_Type'])
        pep = x['Peptide']
        if pep_type == 0:
            return get_linear_seq(pep, is_i2l=is_i2l)
        elif pep_type == 1:
            return get_mono_seq(pep, is_i2l=is_i2l)
        elif pep_type == 2:
            return get_loop_seq(pep, is_i2l=is_i2l)
    
    df['seq'] = df.apply(lambda x:get_line_seq(x, is_i2l=is_i2l), axis=1, result_type='reduce')
    return df


def legal_digest_n(seq, prtac2seq, prtac, start_idx, is_digest=True, digest_n={}, digest_c={}):
    """肽N端是特异酶切
    C端特异酶, 上一个肽末尾符合C端酶切
    N端特异酶, 该肽第一个氨基酸符合N端酶切

    Args:
        is_digest: 是否限制是特异酶切, 不限制直接合法
        start_idx: 肽段在蛋白上的开始未知, 从0开始
        diest_c: {'K','R'}, trypsin
    """

    if not is_digest:
        return 1
    
    # 肽在蛋白N端
    if start_idx == 0 or (start_idx == 1 and prtac2seq[prtac][0]=='M'):
        return 1
    
    # N特异酶切
    if seq[0] in digest_n:
        return 1
    
    # C特异酶切
    if prtac2seq[prtac][start_idx - 1] in digest_c:
        return 1
    
    return 0


def legal_digest_c(seq, prtac2seq, prtac, start_idx, is_digest=True, digest_n={}, digest_c={}):
    """肽C端是特异酶切
    C端特异酶, 该肽最后一个氨基酸符合C端酶切
    N端特异酶, 下一个肽开始符合N端酶切

    Args:
        is_digest: 是否限制是特异酶切, 不限制直接合法
        start_idx: 肽段在蛋白上的开始未知, 从0开始
        diest_c: {'K','R'}, trypsin
    """

    if not is_digest:
        return 1


    # 肽在蛋白C端
    if start_idx + len(seq) == len(prtac2seq[prtac]):
        return 1

    # 符合N特异酶切
    if prtac2seq[prtac][start_idx + len(seq)] in digest_n:
        return 1

    # 符合C特异酶切
    if seq[-1] in digest_c:
        return 1
    
    return 0


def legal_digest_2item(seq, prtac2seq, prtac, start_idx, is_digest_n=True, is_digest_c=True, digest_n={}, digest_c={}):
    """肽N端和C端都符合特异酶切
    
    Args:
        is_digest_n: 肽n端是否限制是特异酶切, 不限制直接合法
    """

    # 不酶切
    if not digest_n and not digest_c:
        return 1

    # 两端都不限制，非特异, 即肽开始结束氨基酸哪个都可以
    if not is_digest_n and not is_digest_c:
        return 1

    # 两端特异
    if is_digest_n:
        if not legal_digest_n(seq, prtac2seq, prtac, start_idx, is_digest=True, digest_n=digest_n, digest_c=digest_c):
            return 0

    if is_digest_c:
        if not legal_digest_c(seq, prtac2seq, prtac, start_idx, is_digest=True, digest_n=digest_n, digest_c=digest_c):
            return 0
    
    return 1


def legal_digest_1item(seq, prtac2seq, prtac, start_idx, is_digest_n=True, is_digest_c=True, digest_n={}, digest_c={}):
    """肽N端和C端, 一端符合酶切，一端一定不符合, 半特异
    
    Args:
        is_digest_n: 肽n端是否限制是特异酶切, 半特异为True
        is_digest_c: 肽c端是否限制是特异酶切, 半特异为True
    """

    flag = 0
    
    if is_digest_n:
        if not legal_digest_n(seq, prtac2seq, prtac, start_idx, is_digest=True, digest_n=digest_n, digest_c=digest_c):
            flag += 1

    if is_digest_c:
        if not legal_digest_c(seq, prtac2seq, prtac, start_idx, is_digest=True, digest_n=digest_n, digest_c=digest_c):
            flag += 1
    
    if flag == 1:
        return 1
    
    return 0


def legal_digest_type(seq, prtseq, start_idx, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """肽酶切是否合法, 根据不同酶切类型

    Args:
        seq: 肽序列
        prtseq: 蛋白序列
        start_idx: 肽在蛋白上的起始位点, 从0开始
        digest: 酶切规则, [[N端酶切],[[C端酶切AA],[C端酶切blocked AA]]]
        digest_type: 0, 两端特异; 1, 一端特异; 2, 非特异; 3, 蛋白不酶切
        digest_part: [距离N端长度, 距离C端长度], 即切掉N端或C端部分, 当做蛋白新的N端或C端, 主要考虑信号肽和合成肽
    """

    if digest_type == 3:
        if start_idx == 0 and len(seq) == len(prtseq):
            return 1
        return 0
    
    if digest_type == 2:
        return 1
    
    end_idx = start_idx + len(seq) - 1 # 肽在蛋白上的结束位点索引
    
    if digest_type == 0:

        # 肽N端
        if start_idx > 0 and (seq[0] not in digest[0][0] or prtseq[start_idx-1] in digest[0][1]) and (prtseq[start_idx-1] not in digest[1][0] or seq[0] in digest[1][1]):

            if (not (start_idx == 1 and prtseq[0]=='M')) and (not (start_idx <= digest_part[0])): # start_idx从0开始，start_idx前(不包括start_idx位点本身)允许截断n个氨基酸
                return 0
        

        # 肽C端
        if end_idx < len(prtseq)-1 and (prtseq[end_idx+1] not in digest[0][0] or seq[-1] in digest[0][1]) and (seq[-1] not in digest[1][0] or prtseq[end_idx+1] in digest[1][1]):

            if not (len(prtseq) - (end_idx+1) <= digest_part[1]): ## end_idx+1表示截断点前有多少个氨基酸, 允许截断C端n个氨基酸
                return 0
        
        return 1
    

    if digest_type == 1:
        # 肽N端
        if start_idx > 0 and (seq[0] not in digest[0][0] or prtseq[start_idx-1] in digest[0][1]) and (prtseq[start_idx-1] not in digest[1][0] or seq[0] in digest[1][1]):

            # 肽C端
            if end_idx < len(prtseq)-1 and (prtseq[end_idx+1] not in digest[0][0] or seq[-1] in digest[0][1]) and (seq[-1] not in digest[1][0] or prtseq[end_idx+1] in digest[1][1]):

                if (not (start_idx == 1 and prtseq[0]=='M')) and (not (start_idx <= digest_part[0])):

                    if not (len(prtseq) - end_idx+1 <= digest_part[1]):
                        return 0

        return 1



def infer_protein(seq_ls, fpath_prtac2seq, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """序列推断蛋白

    Args:
        fpath_prtac2seq: dict类型为prtac2seq, 否则为fpath_fasta
    Returns:
        {肽序列: [(蛋白ac, 开始位点(从1开始) ), ...], }
    """

    if not seq_ls: # 没有要推断的序列
        return {}

    if not isinstance(fpath_prtac2seq, dict):
        prtac2seq = fasta_file.read_fasta_i2l(fpath_prtac2seq, is_i2l=is_i2l, is_shortac=is_shortac)
    else:
        prtac2seq = fpath_prtac2seq
    
    seq_ls = set(list(seq_ls))

    auto = ahocorasick.Automaton()
    for idx, key in enumerate(seq_ls):
        auto.add_word(key, (idx, key))
    auto.make_automaton()

    infer_all = [] # (蛋白ac，开始位点（从0开始），结束位点，肽序列编号（从0开始），肽序列)
    for prtac, prtseq in prtac2seq.items():
        for end_idx, (insert_order, original_value) in auto.iter(prtseq):
            start_idx = end_idx - len(original_value) + 1

            # 肽限制两端氨基酸符合特异酶切
            if legal_digest_type(original_value, prtac2seq[prtac], start_idx, digest=digest, digest_type=digest_type, digest_part=digest_part):
                infer_all.append((prtac, start_idx, end_idx, insert_order, original_value))

    seq2prtac_idx = {} # 肽序列: [(蛋白ac，开始位点（从1开始）)，...]
    for infer_one in infer_all:
        prtac = infer_one[0]
        prt_idx = infer_one[1]+1
        seq = infer_one[4]
        if seq in seq2prtac_idx:
            seq2prtac_idx[seq].append((prtac, prt_idx))
        else:
            seq2prtac_idx[seq] = [(prtac, prt_idx)]

    return seq2prtac_idx



def infer_protein_x_subseq(seq_ls, fpath_prtac2seq, min_len=5, set_x={'X'}, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """包含未知氨基酸序列, 即fasta是X可以是任何氨基酸

    蛋白序列推断, 采用推断子串的方法, 可以支持一个肽中包含多个X, 只要肽足够长
    若枚举替换氨基酸的方法, 一个蛋白上有多个X时, 20^n空间太大
    
    Args:
        fpath_prtac2seq: dict类型为prtac2seq, 否则为fpath_fasta
        set_x: 未知氨基酸
    Returns:
        {肽序列: [(蛋白ac, 开始位点(从1开始) ), ...], }
    """
    if not seq_ls: # 没有要推断的序列
        return {}

    # 包含X的蛋白质序列, X表示不确定氨基酸, MS Annika软件支持鉴定, Ecoli
    # DLAETALYQLQVTTKSPMGALVYGSGGLLIDNGXLR: Ecoli蛋白序列
    # DLAETALYQLQVTTKSPMGALVYGSGGLLLDNGYLR: MS Annika鉴定序列

    # 包含X的蛋白质序列的蛋白质推断
    # 将肽段系列枚举切割，最短查询为5，所以要求序列长度大于等于10（不小于4*2+1）
    # 例如AAAAXAAAA，为9


    def align_seq_and_protein_has_x(prtac2seq, prtac, start_idx, seq, set_x):
        """
        蛋白序列包含X, 肽段为被推断

        重新推断后, 核对蛋白序列上是否包含X

        且未对齐氨基酸只能是X
        防止特殊情况, 鉴定序列ABCDEFG, 蛋白序列ABCXEFG、ABCDEFH, 后者推断是错的
        
        举例: Ecoli蛋白库序列: 
        MR,LGNGMEVEPWLKAAVRKEFVDDNR,VK
        MR,LGNGMEXEPWLKAAVRKEFVDDNR,VK

        Args:
            start_idx: 从0开始
            seq: 鉴定序列
        """

        prt_seq_str = prtac2seq[prtac][start_idx:start_idx+len(seq)]
        
        if not set(list(prt_seq_str)).intersection(set(set_x)):
            return False
        
        for i in range(len(seq)):
            if (seq[i] != prt_seq_str[i]) and (prt_seq_str[i] not in set_x):
                return False

        return True


    if not isinstance(fpath_prtac2seq, dict):
        prtac2seq = fasta_file.read_fasta_i2l(fpath_prtac2seq, is_i2l=is_i2l, is_shortac=is_shortac)
    else:
        prtac2seq = fpath_prtac2seq


    seq2prtac_idx_x = {}

    for seq in seq_ls:
        assert len(seq) > 2*(min_len-1)+1, '没有推断到蛋白质的序列\'%s\'小于可推断子串长度%d'%(seq, min_len)

        if seq in seq2prtac_idx_x:
            continue
        
        # 切割序列，{子串:肽段上起始位点从0开始，}
        seq_x_tuple = []

        # 有重复子串, 虽会覆盖dict, 但也能work, 最后取的是子串推断的最大交集, 只是两个重复子串只记录了一次
        # 当然最好是改成元组, 例如ABCDEGABCDE, 对应蛋白(1)ABCDEXABCDE、 (2)ABCDXGABCDE, 会漏掉(1)

        # 滑动窗口
        for i in range(0, len(seq)-min_len+1): # range(x,y)取不到y
            seq_x_tuple.append((seq[i:i+min_len], i))


        # # ------前后缀方法, 多个X肽两端情况可能查询不到
        # # 前缀,包括全长
        # for i in range(min_len, len(seq)+1): # range(x,y)取不到y
        #     seq_x_dict[seq[:i]] = 0 # 取前i个
        # # 后缀,不包括全长
        # for i in range(1, len(seq)-min_len+1): # range(x,y)取不到y
        #     seq_x_dict[seq[i:]] = i # 抛弃前i个

        # # 子串推断
        # seq_x_ls_nterm = [] #N端要符合酶切
        # seq_x_ls_cterm = [] #C端要符合酶切
        # for k, v in seq_x_dict.items():
        #     if v > 0:
        #         seq_x_ls_cterm.append(k)
        #     else:
        #         seq_x_ls_nterm.append(k)

        
        # # 因为后面要限制全肽的两端酶切，此处两端可都False
        # tmp_is_digest = False

        # # 仅N端要限制酶切特异性, C端不要
        # seq2prtac_idx_i_nterm = infer_protein(seq_x_ls_nterm, prtac2seq, is_i2l=is_i2l, is_shortac=is_shortac, is_digest_n=tmp_is_digest, is_digest_c=False, digest_n=digest_n, digest_c=digest_c)

        # # 仅C端要限制酶切特异性，N端不要
        # seq2prtac_idx_i_cterm = infer_protein(seq_x_ls_cterm, prtac2seq, is_i2l=is_i2l, is_shortac=is_shortac, is_digest_n=False, is_digest_c=tmp_is_digest, digest_n=digest_n, digest_c=digest_c)

        # seq2prtac_idx_i = {**seq2prtac_idx_i_nterm, **seq2prtac_idx_i_cterm}
        # # ------

        seq_x_ls = list(set([x[0] for x in seq_x_tuple])) # 去重
        
        # 因为后面要限制全肽的两端酶切，此处两端可都False
        seq2prtac_idx_i = infer_protein(seq_x_ls, prtac2seq, is_i2l=is_i2l, is_shortac=is_shortac, digest_type=2) # 非特异酶切


        # 子串推断蛋白合并，蛋白位点的最大交集
        prtsites = []
        for subseq, sub_i in seq_x_tuple:
            if subseq not in seq2prtac_idx_i:
                continue
            for prtac, idx in seq2prtac_idx_i[subseq]:
                start_idx = idx - sub_i # 从1开始！！！！！！

                # 根据蛋白上的定位, 全肽, 限制两端氨基酸符合特异酶切
                # 鉴定的seq是不包含X的, 可以判断酶切特异性

                # TODO: 可能seq在蛋白上的前后有X, 会影响判断酶切特异性
                # 是X则认为酶切合法

                if start_idx-1 < 0 or start_idx-1+len(seq) > len(prtac2seq[prtac]):
                    continue
                
                if legal_digest_type(seq, prtac2seq[prtac], start_idx-1, digest=digest, digest_type=digest_type, digest_part=digest_part):
                    if align_seq_and_protein_has_x(prtac2seq, prtac, start_idx-1, seq, set_x): # 包含X, 且除X外其他氨基酸一致
                        prtsites.append(prtac+'('+str(start_idx)+')')
        ser = pd.Series(prtsites)
        s = ser.value_counts()
        prtsites = list(s[s == s.max()].index)

        seq2prtac_idx_x[seq] = []
        for prtsite in prtsites:
            prtac = prtsite.split('(')[0]
            prt_idx = int(prtsite.split('(')[1][:-1])
            seq2prtac_idx_x[seq].append((prtac, prt_idx))
        
        assert seq2prtac_idx_x[seq], '序列\'%s\'没有推断到蛋白'%(seq)
    return seq2prtac_idx_x


def infer_protein_include_subset(seq_ls, fpath_fasta, min_len=5, set_x={'X'}, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """序列推断蛋白, 包含X的序列用子序列推断
    """

    seq2prtac_idx = infer_protein(seq_ls, fpath_fasta, is_i2l=is_i2l, is_shortac=is_shortac, digest=digest, digest_type=digest_type, digest_part=digest_part)

    # ================================================
    # 以NOTINPROTEIN_i蛋白替代
    # i = 0
    # for seq in seq_ls:
    #     if seq not in seq2prtac_idx:
    #         i+=1
    #         seq2prtac_idx[seq] = [('NOTINPROTEIN_'+str(i), 1)]
    

    seqs_x = [x for x in seq_ls if x not in seq2prtac_idx]
    seq2prtac_idx_x = infer_protein_x_subseq(seqs_x, fpath_fasta, min_len=min_len, set_x=set_x, is_i2l=is_i2l, is_shortac=is_shortac, digest=digest, digest_type=digest_type, digest_part=digest_part)

    # 在和含X的序列推断，合并前统计
    seq_ls_miss = [x for x in seq_ls if x not in seq2prtac_idx]
    
    for seq, v in seq2prtac_idx_x.items():
        seq2prtac_idx[seq] = v
    # ================================================
        
    return seq2prtac_idx, seq_ls_miss



def set_xl_infer_protein_notin_fasta(df, seq_ls_miss):
    """设置蛋白是否在fasta中找不见"""

    seq_ls_miss = set(seq_ls_miss)

    def is_notin_fasta(x):
        """肽段序列是否可以被直接推断"""

        flag = 0
        if x['seq1'] in seq_ls_miss:
            flag = flag | 1
        if x['seq2']: # 有可能是单肽, seq2为空
            if x['seq2'] in seq_ls_miss:
                flag = flag | 2
        return flag


    df['notin_fasta'] = df.apply(lambda x:is_notin_fasta(x), axis=1, result_type='reduce')
    num_miss = sum([0 if x==0 else 1 for x in list(df['notin_fasta'])])

    print(f'存在fasta中找不见蛋白的肽: {num_miss}/{df.shape[0]}')
    # assert not num_miss,'存在fasta中找不见蛋白的肽'

    return df



def set_xl_prtsites(df, fpath_fasta, min_len=5, set_x={'X'}, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """交联两条肽序列单独推断
    prtsites1格式: [蛋白AC(交联位点), ], 位点从1开始
    prtsites1/2和seq1/2一一对应

    Args:
        x_min_len: 如果fasta存在X可以代表任意氨基酸, 查询最短长度
    """

    seq_ls = []
    df = split_xl_pep_str(df, is_i2l=is_i2l) # Peptide切分
    seq_ls.extend(list(df['seq1']))
    seq_ls.extend(list(df['seq2']))
    seq_ls = list(set(seq_ls))


    seq2prtac_idx, seq_ls_miss = infer_protein_include_subset(seq_ls, fpath_fasta, min_len=min_len, set_x=set_x, is_i2l=is_i2l, is_shortac=is_shortac, digest=digest, digest_type=digest_type, digest_part=digest_part)

    df = set_xl_infer_protein_notin_fasta(df, seq_ls_miss)

    def get_seq_prtsite(x):
        prtsite1 = []
        prtsite2 = []
        # 从tuple转为括号字符串格式 # 肽段和蛋白位点均从1开始计数
        for prtac, idx in seq2prtac_idx[x['seq1']]:
            prtsite1.append(prtac+'('+str(idx+x['link1']-1)+')')
        for prtac, idx in seq2prtac_idx[x['seq2']]:
            prtsite2.append(prtac+'('+str(idx+x['link2']-1)+')')
        
        prtsite1 = sorted(list(set(prtsite1)))
        prtsite2 = sorted(list(set(prtsite2)))
        return prtsite1, prtsite2

    df['prtsites1'] = df.apply(lambda x:get_seq_prtsite(x)[0], axis=1, result_type='reduce')
    df['prtsites2'] = df.apply(lambda x:get_seq_prtsite(x)[1], axis=1, result_type='reduce')
    
    return df



def set_xl_prts(df):
    """
    肽段对应蛋白质
    存在多个蛋白位点对应一个蛋白的情况, list表示
    
    需要先获得蛋白交联位点列column[prtsites1、prtsites2]

    Returns:
        [蛋白名称, 蛋白质名称]
    """

    def prtsite2prt(x):
        prt = []
        for prtsite in x:
            prt.append(prtsite.split('(')[0])
        return sorted(list(set(prt)))

    df['prts1'] = df['prtsites1'].apply(lambda x:prtsite2prt(x))
    df['prts2'] = df['prtsites2'].apply(lambda x:prtsite2prt(x))
    return df


def get_xl_prtsites_str(x):
    """获得蛋白位点对字符串表示
    可能对应多个蛋白位点，用$$分割
    
    Returns:
        蛋白质名称(位点)$$蛋白质名称(位点)-蛋白质名称(位点)
    """

    prtsites1 = '$$'.join(x['prtsites1'])
    prtsites2 = '$$'.join(x['prtsites2'])

    # 不可以用肽序列顺序排列
    # 可能肽序列顺序不同，但对应一个PPI
    
    prtsites_str = '-'.join(sorted([prtsites1, prtsites2]))
    return prtsites_str


def get_xl_prts_str(x):
    """获得蛋白对字符串表示
    可能对应多个蛋白，用$$分割
    
    Returns:
        蛋白质名称$$蛋白质名称-蛋白质名称
    """

    prts1 = '$$'.join(x['prts1'])
    prts2 = '$$'.join(x['prts2'])

    # 不可以用肽序列顺序排列
    # 可能肽序列顺序不同，但对应一个PPI

    prts_str = '-'.join(sorted([prts1, prts2]))
    return prts_str


def set_xl_rp_ppi_str(df):
    """
    设置rp与ppi层的唯一字符串表示
    """

    df['prtsites_str'] = df.apply(lambda x:get_xl_prtsites_str(x), axis=1, result_type='reduce')
    df['prts_str'] = df.apply(lambda x:get_xl_prts_str(x), axis=1, result_type='reduce')
    return df


def set_reinfer_xl_rp_ppi_str(df, fpath_fasta, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """
    蛋白重新推断
    设置rp与ppi层的唯一字符串表示
    """

    df = set_xl_prtsites(df, fpath_fasta, is_i2l=is_i2l, is_shortac=is_shortac, digest=digest, digest_type=digest_type, digest_part=digest_part)
    df = set_xl_prts(df)
    df = set_xl_rp_ppi_str(df)

    return df


def set_xl_prts_link(df):
    """按照prtsites的顺序, 单独的蛋白质层面位点list, 不包括蛋白名称, 位点从1开始
    Returns:
        [蛋白位点1, ...]
    """

    def get_prt_link(x):
        ls = []
        for xx in x:
            ls.append(int(xx.split('(')[-1][:-1]))
        # return list(set(ls))
        return ls

    df['prts_link1'] = df['prtsites1'].apply(get_prt_link)
    df['prts_link2'] = df['prtsites2'].apply(get_prt_link)

    return df


def set_xl_prt_link_char(df, is_i2l=False):
    """获得交联位点的字母
    考虑蛋白质N端
    """

    df = set_xl_prts_link(df)

    def get_xl_prt_link_char(seq, link, prts_link, is_i2l):
        if 1 in prts_link:
            return '['
        elif 2 in prts_link and link == 1:
            return '[' # '交联在蛋白质序列第二位，第一位是M（G），也算蛋白N端'
        else:
            seq = get_seq_i2l(seq, is_i2l=is_i2l)
            return seq[link-1]

    df['prt_link1_char'] = df.apply(lambda x:get_xl_prt_link_char(x['seq1'], x['link1'], x['prts_link1'], is_i2l=is_i2l), axis=1, result_type='reduce')
    df['prt_link2_char'] = df.apply(lambda x:get_xl_prt_link_char(x['seq2'], x['link2'], x['prts_link2'], is_i2l=is_i2l), axis=1, result_type='reduce')
    
    return df

def set_xl_maybe_monollink(df):
    """交联可能是monolink"""

    def get_ls_xl_distance_prt(x):
        """蛋白的交联位点距离"""
        ls_xl_prt = [] # 是同蛋白交联, 位点距离（正值，大减小）
        for prtsite1 in x['prtsites1']:
            for prtsite2 in x['prtsites2']:
                prts1 = prtsite1.split('(')[0]
                linksite1 = int(prtsite1.split('(')[1][:-1])
                prts2 = prtsite2.split('(')[0]
                linksite2 = int(prtsite2.split('(')[1][:-1])
                if prts1 == prts2:
                    ls_xl_prt.append(abs(linksite1-linksite2))
        return ls_xl_prt
    
    def get_ls_xl_distance_pep(x):
        """肽段的交联位点距离"""
        len1 = len(x['seq1'])
        len2 = len(x['seq2'])
        link1 = x['link1']
        link2 = x['link2']
        dis1 = len1-link1+link2
        dis2 = len2-link2+link1
        ls_xl_pep = [dis1, dis2] # 两条肽交联位点，在蛋白层面的差值
        return ls_xl_pep
    
    def maybe_monolink(x):
        # if not set(x['prts1']).intersection(set(x['prts2'])): # 不是Intra
            # return 0
        if set(x['ls_xl_dis_prt']).intersection(set(x['ls_xl_dis_pep'])):
            return 1
        return 0
    
    df['ls_xl_dis_prt'] = df.apply(lambda x:get_ls_xl_distance_prt(x), axis=1, result_type='reduce')
    df['ls_xl_dis_pep'] = df.apply(lambda x:get_ls_xl_distance_pep(x), axis=1, result_type='reduce')
    df['maybe_monolink'] = df.apply(lambda x:maybe_monolink(x), axis=1, result_type='reduce')
    return df




def set_line_prts(df, fpath_fasta, is_i2l=True, is_shortac=True):
    """获得line型肽段（单肽/mono/loop）推断蛋白
    """

    df = set_line_seq(df, is_i2l=is_i2l)
    seq_ls = list(df['seq'])
    seq_ls = list(set(seq_ls))

    seq2prtac_idx = infer_protein(seq_ls, fpath_fasta, is_i2l=is_i2l, is_shortac=is_shortac)

    def get_seq_prts(x):
        prts = []
        for prtac, idx in seq2prtac_idx[x['seq']]:
            prts.append(prtac)
        prts = sorted(prts)
        return prts
    
    df['prts'] = df.apply(lambda x:get_seq_prts(x), axis=1, result_type='reduce')

    return df


def get_pfind_seq_str(x, is_i2l=True):
    """获得pFind的序列字符串表示

    Returns:
        seq
    """
    seq = x['Sequence']

    seq = seq.strip()
    return get_seq_i2l(seq, is_i2l=is_i2l)


def get_pfind_mod_dict(mods_str):
    """获取pFind的seq对应mod

    6,Carbamidomethyl[C];13,Xlink_Leiker_cAL_monoH2O[K];
    11,Xlink_Leiker_cAL_monoH2O[K];

    Args:
        mods_str: csv中存储的Modification字符串
    Returns:
        mod:修饰的字典
    """

    mod = {} # 位点：修饰名称，修饰位点从1计数, N端修饰为0, C端修饰为序列长度+1

    if (pd.isnull(mods_str)) or (mods_str=='null') or (not mods_str):
        return mod
    else:
        segs = mods_str.split(';')[:-1]
        for seg in segs:
            ss = seg.strip().split(',')
            mod_name = ss[1]
            mod_site = int(ss[0])
            mod[mod_site] = mod_name

    return mod

def get_pfind_pep_str(x, is_i2l=True):
    """获得pFind的肽段字符串表示

    Returns:
        seq_mod1(site1);
    """

    seq = get_pfind_seq_str(x, is_i2l=is_i2l)
    stie2mod = get_pfind_mod_dict(x['Modification'])
    mod_str = get_mod_str(stie2mod)

    return '_'.join([seq, mod_str])


def get_pfind_seqmod_str(x, is_i2l=True):
    """获得pFind的序列修饰字符串表示

    只包含修饰的种类和个数不包括修饰的位点
    """

    seq = get_pfind_seq_str(x, is_i2l=is_i2l)
    stie2mod = get_pfind_mod_dict(x['Modification'])
    
    mod2cnt = {}
    for mod_site, mod_name in stie2mod.items():
        if mod_name in mod2cnt:
            mod2cnt[mod_name] += 1
        else:
            mod2cnt[mod_name] = 1
    
    if mod2cnt:
        # 按照修饰名称排序
        mod2cnt = dict(sorted(mod2cnt.items(), key=lambda x:x[0],reverse=False))
        mod_str = ';'.join([str(mod_cnt)+'*'+mod_name for mod_name, mod_cnt in mod2cnt.items()])
    else:
        mod_str = '0*null'

    return '_'.join([seq, mod_str])


def get_pfind_prec_str(x, is_i2l=True):
    """获得pFind的肽段的母离子字符串表示

    Returns:
        pep-str_charge
    """

    pep_str = get_pfind_pep_str(x, is_i2l=is_i2l)
    charge = x['Charge']
    return '_'.join([pep_str, str(charge)])


def get_pfind_psm_scan_str(x, is_i2l=True):
    """获得pFind的PSM的Scan字符串表示

    Returns:
        prec-str_rawname.scan.scan
    """

    prec_str = get_pfind_prec_str(x, is_i2l=is_i2l)
    rawscan_str = '.'.join([x['rawname'], str(x['scannum']), str(x['scannum'])])


    return '_'.join([prec_str, rawscan_str])


def get_pfind_psm_str(x, is_i2l=True):
    """获得pFind的PSM的字符串表示

    Returns:
        prec-str_title
    """

    title = x['Title']

    prec_str = get_pfind_prec_str(x, is_i2l=is_i2l)
    psm_str = '_'.join([prec_str, title])

    return psm_str


def set_pfind_level_str(df, is_i2l=True, is_title=True):
    """设置pFind的肽段的字符串表示
    """

    df['seq_str'] = df.apply(lambda x:get_pfind_seq_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
    df['seq_mod_str'] = df.apply(lambda x:get_pfind_seqmod_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
    df['pep_str'] = df.apply(lambda x:get_pfind_pep_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
    df['prec_str'] = df.apply(lambda x:get_pfind_prec_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
    if is_title:
        df['psm_str'] = df.apply(lambda x:get_pfind_psm_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')

    return df


def get_plink_line_seq(x, is_i2l=False):
    """获得regular、monolink、looplink的肽段序列
    Args:
        x: ['Peptide']
    """
    seq = x.split('(')[0]
    return get_seq_i2l(seq, is_i2l=is_i2l)


def get_plink_monolink_link(x):
    """获得monolink的交联位点
    Args:
        x: ['Peptide']
    """
    return int(x.split('(')[1][:-1])


def get_plink_looplink_link2(x):
    """获得交联的link2
    Args:
        x: ['Peptide']
    """
    return int(x.split('(')[2][:-1])


def get_plink_line_seq_seq_str(x, peptype=0, is_i2l=False):
    """获得regular、monolink、looplink的肽段序列字符串(不包括交联位点)表示
    """

    pep = x['Peptide']
    seq = get_plink_line_seq(pep, is_i2l=is_i2l)
    return seq


def get_plink_line_seq_seqmod_str(x, peptype=0, is_i2l=False):
    """获得regular、monolink、looplink的肽段序列字符串(不包括交联位点)表示
    修饰不包含位点, 只包括种类和个数
    """

    pep = x['Peptide']
    mods = x['Modifications']

    seq_str = get_plink_line_seq_seq_str(x, peptype=peptype, is_i2l=is_i2l)

    seq = get_plink_line_seq(pep, is_i2l=is_i2l)
    mod1, mod2 = get_xl_mod_dict(seq, mods)
    if mod2: assert False, 'get_plink_pep_line_str has mod2'

    mod_str = get_modcnt_str(mod1) # 内部有修饰名称排序

    return '_'.join([seq_str, mod_str])


def get_plink_line_seq_str(x, peptype=0, is_i2l=False):
    """获得regular、monolink、looplink的肽段序列字符串(包括交联位点)表示

    Args:
        peptype: 0-regular, 1-monolink, 2-looplink
    
    Returns:
        regular: seq
        monolink: seq(link)
        looplink: seq1(link1)(link2)
    """

    pep = x['Peptide']
    seq = get_plink_line_seq(pep, is_i2l=is_i2l)
    if peptype == 0:
        return '_'.join([seq])
    if peptype == 1:
        link = get_plink_monolink_link(pep)
        return '_'.join([seq+'('+str(link)+')'])
    if peptype == 2:
        link1 = get_plink_monolink_link(pep)
        link2 = get_plink_looplink_link2(pep)
        if link2 < link1:
            link1, link2 = link2, link1
        return '_'.join([seq+'('+str(link1)+')'+'('+str(link2)+')'])


def get_plink_line_pep_str(x, peptype=0, is_i2l=True):
    """获得regular、monolink、looplink的肽段字符串表示

    Args:
        peptype: 0-regular, 1-monolink, 2-looplink
    
    Returns:
        regular: seq_mod1(site1);seq_mod2(site2)
        monolink: seq(link)_mod1(site1);mod2(site2)
        looplink: seq1(link1)(link2)_mod1(site1);mod2(site2), link2>link1
    """

    pep = x['Peptide']
    mods = x['Modifications']

    seq = get_plink_line_seq(pep, is_i2l=is_i2l)
    mod1, mod2 = get_xl_mod_dict(seq, mods)
    
    if mod2: assert False, 'get_plink_pep_line_str has mod2'
    
    mod_str = get_mod_str(mod1) # 内部带了修饰根据位点排序

    seq_str = get_plink_line_seq_str(x, peptype=peptype, is_i2l=is_i2l)

    return '_'.join([seq_str, mod_str])



def get_plink_line_prec_str(x, peptype=0, is_i2l=True):
    """获得regular、monolink、looplink的母离子字符串表示

    Args:
        peptype: 0-regular, 1-monolink, 2-looplink
    
    Returns:
        pep-str_charge
    """

    pep_str = get_plink_line_pep_str(x, peptype=peptype, is_i2l=is_i2l)
    charge = x['Charge']
    return '_'.join([pep_str, str(charge)])



def get_plink_line_csm_str(x, peptype=0, is_i2l=True):
    """获得regular、monolink、looplink的PSM的字符串表示

    Args:
        peptype: 0-regular, 1-monolink, 2-looplink
    
    Returns:
        prec-str_title
    """

    title = x['Title']

    prec_str = get_plink_line_prec_str(x, peptype=peptype, is_i2l=is_i2l)
    psm_str = '_'.join([prec_str, title])

    return psm_str


def set_plink_line_level_str(df, peptype=0, is_i2l=True, is_modcnt=True):

    df['seq_seq_str'] = df.apply(lambda x:get_plink_line_seq_seq_str(x, peptype=peptype, is_i2l=is_i2l), axis=1, result_type='reduce')
    df['seq_str'] = df.apply(lambda x:get_plink_line_seq_str(x, is_i2l=is_i2l), axis=1, result_type='reduce')
    df['pep_str'] = df.apply(lambda x:get_plink_line_pep_str(x, peptype=peptype, is_i2l=is_i2l), axis=1, result_type='reduce')
    df['prec_str'] = df.apply(lambda x:get_plink_line_prec_str(x, peptype=peptype, is_i2l=is_i2l), axis=1, result_type='reduce')
    df['csm_str'] = df.apply(lambda x:get_plink_line_csm_str(x, peptype=peptype, is_i2l=is_i2l), axis=1, result_type='reduce')
    if is_modcnt:
        df['seq_seqmod_str'] = df.apply(lambda x:get_plink_line_seq_seqmod_str(x, peptype=peptype, is_i2l=is_i2l), axis=1, result_type='reduce')

    return df


def split_plink_looplink_str(df, is_i2l=False, is_linkchar=True):
    """拆分['Peptide']到seq1/seq2/link1/link2/link1_char/link2_char"""

    def get_looplink_link(x):
        """link1 < link2"""
        x = x['Peptide']

        link1 = int(x.split('(')[1].split(')')[0])
        link2 = int(x.split('(')[2][:-1])
        if link1 < link2:
            return [link1, link2]
        else:
            return [link2, link1]
    
    def get_looplink_link_char(x):
        """link1 < link2"""
        x = x['Peptide']
        link1 = int(x.split('(')[1].split(')')[0])
        link2 = int(x.split('(')[2][:-1])
        if not (link1 < link2):
            link1, link2 = link2, link1
        
        seq = get_seq_i2l(x.split('(')[0], is_i2l=is_i2l)
        return seq[link1-1], seq[link2-1]

    df['seq1'] = df['Peptide'].apply(lambda x:get_seq_i2l(x.split('(')[0], is_i2l=is_i2l))
    df['seq2'] = '' # 单肽seq2为空
    df[['link1', 'link2']] = (df.apply(lambda x: get_looplink_link(x), axis=1, result_type='expand') if df.shape[0] > 0 else pd.DataFrame(columns=['link1', 'link2']))

    if is_linkchar:
        df[['link1_char', 'link2_char']] = (df.apply(lambda x: get_looplink_link_char(x), axis=1, result_type='expand') if df.shape[0] > 0 else pd.DataFrame(columns=['link1_char', 'link2_char']))
    return df



def set_plink_looplink_prtsites(df, fpath_fasta, min_len=5, set_x={'X'}, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """looplink肽序列推断蛋白
    prtsites1格式: [蛋白AC(交联位点1)(交联位点2), ], 位点从1开始
    prtsites2格式: []空列表

    Args:
        x_min_len: 如果fasta存在X可以代表任意氨基酸, 查询最短长度
    """

    seq_ls = []
    df = split_plink_looplink_str(df, is_i2l=is_i2l) # Peptide切分
    seq_ls.extend(list(df['seq1']))
    seq_ls = list(set(seq_ls))

    seq2prtac_idx, seq_ls_miss = infer_protein_include_subset(seq_ls, fpath_fasta, min_len=min_len, set_x=set_x, is_i2l=is_i2l, is_shortac=is_shortac, digest=digest, digest_type=digest_type, digest_part=digest_part)

    df = set_xl_infer_protein_notin_fasta(df, seq_ls_miss)


    def get_seq_prtsite(x):
        prtsite1 = []
        prtsite2 = []
        # 从tuple转为括号字符串格式 # 肽段和蛋白位点均从1开始计数
        for prtac, idx in seq2prtac_idx[x['seq1']]:
            prtsite1.append(prtac+'('+str(idx+x['link1']-1)+')'+'('+str(idx+x['link2']-1)+')')
        
        prtsite1 = sorted(list(set(prtsite1)))
        prtsite2 = sorted(list(set(prtsite2)))
        return prtsite1, prtsite2

    df['prtsites1'] = df.apply(lambda x:get_seq_prtsite(x)[0], axis=1, result_type='reduce')
    df['prtsites2'] = df.apply(lambda x:get_seq_prtsite(x)[1], axis=1, result_type='reduce')
    
    return df



def split_plink_monolink_str(df, is_i2l=False, is_linkchar=True):
    """拆分['Peptide']到seq1/seq2/link1/link2/link1_char/link2_char"""

    df['seq1'] = df['Peptide'].apply(lambda x:get_seq_i2l(x.split('(')[0], is_i2l=is_i2l))
    df['seq2'] = '' # 单肽seq2为空

    df['link1'] = df['Peptide'].apply(lambda x:int(x.split('(')[1][:-1]))
    df['link2'] = -1 # 单肽seq2为空
    
    if is_linkchar:
        df['link2_char'] = 'B'

        df['link1_char'] = df['Peptide'].apply(lambda x:get_seq_i2l(x.split('(')[0], is_i2l=is_i2l)[int(x.split('(')[1][:-1])-1])
    return df



def set_plink_monolink_prtsites(df, fpath_fasta, min_len=5, set_x={'X'}, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """monolink肽序列推断蛋白
    prtsites1格式: [蛋白AC(交联位点), ], 位点从1开始
    prtsites2格式: []空列表

    Args:
        x_min_len: 如果fasta存在X可以代表任意氨基酸, 查询最短长度
    """

    seq_ls = []
    df = split_plink_monolink_str(df, is_i2l=is_i2l) # Peptide切分
    seq_ls.extend(list(df['seq1']))
    seq_ls = list(set(seq_ls))

    seq2prtac_idx, seq_ls_miss = infer_protein_include_subset(seq_ls, fpath_fasta, min_len=min_len, set_x=set_x, is_i2l=is_i2l, is_shortac=is_shortac, digest=digest, digest_type=digest_type, digest_part=digest_part)

    df = set_xl_infer_protein_notin_fasta(df, seq_ls_miss)


    def get_seq_prtsite(x):
        prtsite1 = []
        prtsite2 = []
        # 从tuple转为括号字符串格式 # 肽段和蛋白位点均从1开始计数
        for prtac, idx in seq2prtac_idx[x['seq1']]:
            prtsite1.append(prtac+'('+str(idx+x['link1']-1)+')')
        
        prtsite1 = sorted(list(set(prtsite1)))
        prtsite2 = sorted(list(set(prtsite2)))
        return prtsite1, prtsite2

    df['prtsites1'] = df.apply(lambda x:get_seq_prtsite(x)[0], axis=1, result_type='reduce')
    df['prtsites2'] = df.apply(lambda x:get_seq_prtsite(x)[1], axis=1, result_type='reduce')
    
    return df


def split_plink_regular_str(df, is_i2l=False, is_linkchar=True):
    """拆分['Peptide']到seq1/seq2/link1/link2/link1_char/link2_char"""

    df['seq1'] = df['Peptide'].apply(lambda x:get_seq_i2l(x.split('(')[0], is_i2l=is_i2l))
    df['seq2'] = '' # 单肽seq2为空

    df['link1'] = -1 # 单肽seq2为空
    df['link2'] = -1 # 单肽seq2为空

    if is_linkchar:
        df['link1_char'] = 'B'
        df['link2_char'] = 'B'
    return df



def set_plink_regular_prtsites(df, fpath_fasta, min_len=5, set_x={'X'}, is_i2l=True, is_shortac=True, digest=[[[],[]],[['K','R'],[]]], digest_type=0, digest_part=[0, 0]):
    """looplink肽序列推断蛋白
    prtsites1格式: [蛋白AC(交联位点1)(交联位点2), ], 位点从1开始
    prtsites2格式: []空列表

    Args:
        x_min_len: 如果fasta存在X可以代表任意氨基酸, 查询最短长度
    """

    seq_ls = []
    df = split_plink_regular_str(df, is_i2l=is_i2l) # Peptide切分
    seq_ls.extend(list(df['seq1']))
    seq_ls = list(set(seq_ls))

    seq2prtac_idx, seq_ls_miss = infer_protein_include_subset(seq_ls, fpath_fasta, min_len=min_len, set_x=set_x, is_i2l=is_i2l, is_shortac=is_shortac, digest=digest, digest_type=digest_type, digest_part=digest_part)

    df = set_xl_infer_protein_notin_fasta(df, seq_ls_miss)


    def get_seq_prtsite(x):
        prtsite1 = []
        prtsite2 = []
        # 从tuple转为括号字符串格式 # 肽段和蛋白位点均从1开始计数
        for prtac, idx in seq2prtac_idx[x['seq1']]:
            prtsite1.append(prtac)
        
        prtsite1 = sorted(list(set(prtsite1)))
        prtsite2 = sorted(list(set(prtsite2)))
        return prtsite1, prtsite2

    df['prtsites1'] = df.apply(lambda x:get_seq_prtsite(x)[0], axis=1, result_type='reduce')
    df['prtsites2'] = df.apply(lambda x:get_seq_prtsite(x)[1], axis=1, result_type='reduce')
    
    return df


def get_plink_line_prtsites_str(x):
    """获得蛋白位点对字符串表示
    可能对应多个蛋白位点，用$$分割
    
    Returns:
        蛋白质名称(位点)$$蛋白质名称(位点)
    """

    prtsites1 = '$$'.join(x['prtsites1'])
    
    prtsites_str = prtsites1
    return prtsites_str


def get_plink_line_prts_str(x):
    """获得蛋白位点对字符串表示
    可能对应多个蛋白位点，用$$分割
    
    Returns:
        蛋白质名称(位点)$$蛋白质名称(位点)
    """

    prts1 = '$$'.join(x['prts1'])
    
    prts_str = prts1
    return prts_str


def set_plink_line_rp_ppi_str(df):

    df['prtsites_str'] = df.apply(lambda x:get_plink_line_prtsites_str(x), axis=1, result_type='reduce')
    df['prts_str'] = df.apply(lambda x:get_plink_line_prts_str(x), axis=1, result_type='reduce')
    return df


def set_reinfer_plink_res_rp_ppi_str(df, fpath_fasta, peptype=3, **kwargs):
    """
    蛋白重新推断
    设置rp与ppi层的唯一字符串表示


    Args:
        peptype: 0:regular, 1:monolink, 2:looplink, 3:xl
        **kwargs: 蛋白质推断相关字典
    """

    if peptype == 0:
        df = set_plink_regular_prtsites(df, fpath_fasta, **kwargs)
        df = set_xl_prts(df)
        df = set_plink_line_rp_ppi_str(df)
    elif peptype == 1:
        df = set_plink_monolink_prtsites(df, fpath_fasta, **kwargs)
        df = set_xl_prts(df)
        df = set_plink_line_rp_ppi_str(df)
    elif peptype == 2:
        df = set_plink_looplink_prtsites(df, fpath_fasta, **kwargs)
        df = set_xl_prts(df)
        df = set_plink_line_rp_ppi_str(df)
    elif peptype == 3:
        df = set_reinfer_xl_rp_ppi_str(df, fpath_fasta, **kwargs)

    return df
