# -*- coding: UTF-8 -*-
# Date      : 21th March, 2025
# Author    : Peng-Zhi Mao


import pandas as pd
import numpy as np
from pathlib import Path
import os
import sys
# sys.path.append(os.path.dirname(os.getcwd()))
# sys.path.append(os.path.dirname(os.path.dirname(os.getcwd())))
# sys.path.append(str(Path.cwd().parent.parent/'utils2'))

sys.path.append(str(Path.cwd()/'utils2'))

from general_link import util_plink_result as util_res
import util_pfind

import ion_match_bat_pfind_linear as ion_match_bat
from ion_match_param import Param

import importlib
importlib.reload(ion_match_bat)

from tqdm import tqdm



def generate_iontype(fpath_iontype_temp, fpath_iontype, dct_chrc):
    ## 生成离子类型表

    # 读取碎片离子类型模板
    df_temp = pd.read_csv(fpath_iontype_temp)

    ls_ls = []
    for ion_name, ion_mass in dct_chrc.items():
        ls_ls.append([ion_name, 0, ion_mass, 1, 0, 0, '-', 0, 2])
    
    df = pd.concat([df_temp, pd.DataFrame(ls_ls, columns=df_temp.columns)], axis=0)

    df.to_csv(fpath_iontype, index=False)

    return df




def filter_csv(fpath_res, fpath_out, modi_one, modi_ls, aa_one, enz_type, pep_level=False):
    ## 筛选csv行

    df = pd.read_csv(fpath_res, sep='\t')

    if pep_level:
        print('------CSM num:', df.shape[0])
        df['pep_str'] = df.apply(lambda x:util_res.get_pfind_pep_str(x, is_i2l=True), axis=1, result_type='reduce')
        # 一个肽只保留一个
        df = df.drop_duplicates('pep_str', keep='first')
        print('------Peptide num:', df.shape[0])

    print('---Enzyme Specificity: \n', df['Specificity'].value_counts())
    df = df[df['Specificity'].isin(enz_type)]

    df['modi_dct'] = df['Modification'].apply(util_res.get_pfind_mod_dict)

    def is_modi_one(modi_dct):

        # 修饰存在，且仅存在一个
        ls = [x for x in modi_dct.values() if x in modi_one]

        if int(len(ls) == 1):
            if modi_ls is None:
                return 1
            # 除目标修饰外，仅考虑的修饰
            if len(set(modi_dct.values()).difference(set(modi_one)).difference(set(modi_ls))) == 0:
                return 1
        return 0
    
    def is_aa_one(x):

        modi_dct = x['modi_dct']
        seq = x['Sequence']

        # 氨基酸存在，且仅存在一个
        ls = [x for x in seq if x in aa_one]

        if int(len(ls) == 1):
            if modi_ls is None:
                return 1
            # 除目标修饰外，仅考虑的修饰
            if len(set(modi_dct.values()).difference(set(modi_one)).difference(set(modi_ls))) == 0:
                return 1
        return 0
    
    if modi_one:
        df['is_modi'] = df['modi_dct'].apply(is_modi_one)
    elif aa_one:
        df['is_modi'] = df.apply(is_aa_one, axis=1, result_type='reduce')
    else:
        assert False, '请设置目标修饰或氨基酸'

    num_before = df.shape[0]

    df = df[df['is_modi'] == 1]

    num_after = df.shape[0]
    print(f'------Before filter: {num_before}, after filter: {num_after}')
    
    if fpath_out:
        df.to_csv(fpath_out, sep='\t', index=False)

    return df



def match_frag(dpath_mgf, fpath_res, fpath_iontype, fpath_out='', progress_visual=True, n_process=2, match_ppm=20, head=None):
    ## 肽谱匹配

    print('-'*12+'Matching fragment ions'+'-'*12)

    param = Param()
    param.ms2_match_tolerance = match_ppm

    # 开头需要sys.path.append(str(Path.cwd()/'utils2'))，不能从util2.ion_match_param import Param

    # 需要重载，MS2误差写在ion_match的头中
    importlib.reload(ion_match_bat.ion_match)

    print('---match_ppm:', param.ms2_match_tolerance, ion_match_bat.ion_match.Param().ms2_match_tolerance, ion_match_bat.ion_match.MS2_MATCH_TOL)

    # df_match的res_id是df_res的行号
    # df_res的行号和索引要一致，多进程用的是索引
    df_match = ion_match_bat.match_pfind_linear_iontype(dpath_mgf, fpath_res, fpath_iontype, fpath_out=fpath_out, progress_visual=progress_visual, n_process=n_process) #, head=10)

    # df_match和df_res需要res_id对齐
    df_res = pd.read_csv(fpath_res, sep='\t')
    df_res['res_id'] = list(range(df_res.shape[0]))
    df_res.to_csv(fpath_res, sep='\t', index=False)

    print('Matched ions:', df_match.shape[0], f'in {df_match['res_id'].nunique()} spectra')

    return df_match


def calc_coverage(fpath_res, fpath_match, ls_filter_cov=[], pep_level=False):
    """计算覆盖度"""

    print('-'*12+'calc_coverage'+'-'*12)

    df = pd.read_csv(fpath_res, sep='\t')
    df_match = pd.read_csv(fpath_match)

    df_match['match_info'] = df_match['match_info'].apply(eval)

    df['pep_len'] = df['Sequence'].apply(len)

    def calc_cov(x):

        res_id = x['res_id']
        pep_len = x['pep_len']

        ls_b = [0]*(pep_len-1)
        ls_y = [0]*(pep_len-1)

        _df = df_match[df_match['res_id'] == res_id]
        if _df.shape[0] == 0:
            return [0, 0, 0]
        
        for i in range(_df.shape[0]):
            ser = _df.iloc[i]
            for match_one in ser['match_info']:
                sym = match_one[0][0].split('|')[0]
                frag_site = int(match_one[0][2])
                if frag_site < 0: continue # 母离子
                if sym == 'b':
                    ls_b[frag_site] = 1
                elif sym == 'y':
                    ls_y[frag_site] = 1
        cov_b = sum(ls_b)/len(ls_b)
        cov_y = sum(ls_y)/len(ls_y)
        ls_by = [x or y for x, y in zip(ls_b, ls_y)]
        cov_by = sum(ls_by)/len(ls_by)
        return [cov_b, cov_y, cov_by]

    ls_ls = []
    for i in tqdm(range(df.shape[0])):
        ls_ls.append(calc_cov(df.iloc[i]))
    
    df[['cov_b', 'cov_y', 'cov_by']] = pd.DataFrame(ls_ls, columns=['cov_b', 'cov_y', 'cov_by'])

    if ls_filter_cov:
        b = (df['cov_b'] > ls_filter_cov[0])
        b = b & (df['cov_y'] > ls_filter_cov[1])
        b = b & (df['cov_by'] > ls_filter_cov[2])
        n1 = df.shape[0]
        df = df[b]
        n2 = df.shape[0]
        print(f'---Coverage {ls_filter_cov} filter, before: {n1}, after: {n2}')

    if pep_level:
        print('------CSM num:', df.shape[0])
        df['pep_str'] = df.apply(lambda x:util_res.get_pfind_pep_str(x, is_i2l=True), axis=1, result_type='reduce')
        # 一个肽只保留一个
        df = df.drop_duplicates('pep_str', keep='first')
        print('------Peptide num:', df.shape[0])

    df.to_csv(fpath_res, sep='\t', index=False)

    return df


def label_res(fpath_res, fpath_match, modi_one, aa_one, dct_chrc):
    ## 结果中增加信息

    print('-'*12+'Labeling characteristic ions'+'-'*12)

    df = pd.read_csv(fpath_res, sep='\t')
    df_match = pd.read_csv(fpath_match)

    df['pep_len'] = df['Sequence'].apply(len)


    df['modi_dct'] = df['Modification'].apply(util_res.get_pfind_mod_dict)

    def get_modi_site(modi_dct):
        ls = [x for x in modi_dct.values() if x in modi_one]
        if len(ls) == 1:
            return list(modi_dct.keys())[list(modi_dct.values()).index(ls[0])]
        assert False, '修饰多于一个'
    
    def get_aa_site(seq):
        """获取氨基酸位置"""
        ls = [x for x in seq if x in aa_one]
        if len(ls) == 1:
            aa = ls[0]
            return seq.index(aa) + 1 # 计数从1开始
        assert False, '氨基酸多于一个'

    if modi_one:
        df['modi_site'] = df['modi_dct'].apply(get_modi_site)  # 计数从1开始
    elif aa_one:
        df['modi_site'] = df['Sequence'].apply(get_aa_site)
    else:
        assert False, '请设置目标修饰或氨基酸'

    def get_modi_site_3quantile(x):
        # 修饰位置的3分位数
        len31 = np.ceil(x['pep_len']/3) # 向上取整
        len32 = np.ceil(x['pep_len']/3*2) # 向上取整

        if x['modi_site'] <= len31:
            return 1
        elif x['modi_site'] > len32:
            return 3
        else:
            return 2
    
    df['modi_site_3quantile'] = df.apply(get_modi_site_3quantile, axis=1)


    # 离子类型符号
    df_match['match_info'] = df_match['match_info'].apply(eval)
    def get_match_info_sym(x):
        ls = []
        for match_one in x:
            ls.append(match_one[0][0].split('|')[0])
        return ls
    df_match['sym'] = df_match['match_info'].apply(get_match_info_sym)


    cols = []
    for chrc_name in dct_chrc.keys():
        cols.append(f'is_{chrc_name}')
        cols.append(f'inten_{chrc_name}')
        cols.append(f'inten_{chrc_name}_match_median')
    cols.append('match_median_inten')

    # 判断特征离子是否存在，强度是多少
    def get_chrc_ion_intensity(x):
        """
        Args:
            x: res_id
        """

        ls = []
        for chrc_name in dct_chrc.keys():
            ls.append(0)
            ls.append(0) # 没匹配强度为0，之前是-1
            ls.append(0) # 没匹配强度为0，之前是-1

        ls_by_inten = []
        _df = df_match[df_match['res_id'] == x]
        for i in range(_df.shape[0]):
            ser = _df.iloc[i]
            if not set(ser['sym']).intersection(set(dct_chrc.keys())):
                ls_by_inten.append(ser['inten_relat'])
            else:
                for j, chrc_name in enumerate(dct_chrc.keys()):
                    if chrc_name in ser['sym']:
                        ls[3*j] = 1
                        if ser['inten_relat'] > ls[3*j+1]: # 多根峰匹配,取强度最大
                            ls[3*j+1] = ser['inten_relat']
        match_median = np.median(ls_by_inten)
        for j, chrc_name in enumerate(dct_chrc.keys()):
            if ls[3*j+1] > 0:
                ls[3*j+2] = ls[3*j+1]/match_median

        ls.append(match_median)
        return ls

    ls_ls = []
    for i in tqdm(list(df['res_id'])):
        ls_ls.append(get_chrc_ion_intensity(i))
    df[cols] = pd.DataFrame(ls_ls, columns=cols)

    df.to_csv(fpath_res, sep='\t', index=False)

    return df


def stats_res(fpath_res, dct_chrc):
    """统计特征离子出现的比例"""

    df = pd.read_csv(fpath_res, sep='\t')

    def stat_num(df):
        num = df.shape[0]
        print('---'*3, f'All num {num}')

        for c in dct_chrc.keys():
            print('---'*3, c)

            num = df[f'is_{c}'].sum()
            ratio = num/df.shape[0]
            print(f'---num {num}, ratio {ratio:.4f}')

            print('------stie n r-total, r-occur, \t inten: q1 q2 q3 \t q1m q2m q3m \t\t q1m*r q2m*r q3m*r')
            for i in range(1, 4):

                # 只考虑特征离子存在的, 之前不存在强度为-1
                # _df1 = df[df[f'is_{c}']==1]

                # 包括特征离子不存在的，强度为0
                _df1 = df.copy()
                
                _df2 = _df1[_df1['modi_site_3quantile'] == i]
                r = _df2.shape[0]/_df1.shape[0] # 1/3该位置占的总数
                r_occur = _df2[f'is_{c}'].sum()/_df2.shape[0] # 该位置特征离子出现比例

                q1 = _df2[f'inten_{c}'].quantile(0.25)
                q2 = _df2[f'inten_{c}'].quantile(0.5)
                q3 = _df2[f'inten_{c}'].quantile(0.75)
            
                q1m = _df2[f'inten_{c}_match_median'].quantile(0.25)
                q2m = _df2[f'inten_{c}_match_median'].quantile(0.5)
                q3m = _df2[f'inten_{c}_match_median'].quantile(0.75)

                # print(f'---{i}/3 {_df2.shape[0]} {r:.4f} \t {q1:.4f} {q2:.4f} {q3:.4f} \t {q1m:.4f} {q2m:.4f} {q3m:.4f}')
                print(f'---{i}/3 {_df2.shape[0]} {r:.4f} {r_occur:.4f} \t {q1:.4f} {q2:.4f} {q3:.4f} \t {q1m:.4f} {q2m:.4f} {q3m:.4f} \t {q1m*r:.4f} {q2m*r:.4f} {q3m*r:.4f}')

    print('==='*6, 'Stats CSM level', '==='*6)
    stat_num(df)

    df['pep_str'] = df.apply(lambda x:util_res.get_pfind_pep_str(x, is_i2l=True), axis=1, result_type='reduce')
    # 一个肽只保留一个
    df = df.drop_duplicates('pep_str', keep='first')

    print('==='*6, 'Stats Peptide level', '==='*6)
    stat_num(df)


def flow_1(fpath_iontype_temp, fpath_iontype, dct_chrc, fpath_res, fpath_out, modi_one, modi_ls, aa_one, enz_type, pep_level, dpath_mgf, fpath_out_match, n_process, frag_ppm, ls_filter_cov):

    df_iontype = generate_iontype(fpath_iontype_temp, fpath_iontype, dct_chrc)

    df_res = filter_csv(fpath_res, fpath_out, modi_one, modi_ls, aa_one, enz_type)

    df_match = match_frag(dpath_mgf, fpath_out, fpath_iontype, fpath_out=fpath_out_match, progress_visual=True, n_process=n_process, match_ppm=frag_ppm)

    df_label = calc_coverage(fpath_out, fpath_out_match, ls_filter_cov, pep_level)

    df_label = label_res(fpath_out, fpath_out_match, modi_one, aa_one, dct_chrc)

    stats_res(fpath_out, dct_chrc)


def _main():
    ## 路径

    # 是否肽段水平
    pep_level = True

    # mgf文件夹
    dpath_mgf = r'F:\pFindCooperation\202503_imide_position\20250303_Propionyl\mgf'

    # pfind结果路径
    fpath_res = r'F:\pFindCooperation\202503_imide_position\20250303_Propionyl\pFind-Filtered.spectra'

    # 输出路径（非必须修改）
    fpath_out = fpath_res.replace('.spectra', f'_pep{int(pep_level)}_label.spectra')

    # 输出碎片离子匹配路径（非必须修改）
    fpath_out_match = str(Path(fpath_out).parent/(Path(fpath_out).stem+f'_fragMatch.csv'))

    # 目标修饰，一条肽仅一个
    # modi_one = set(['Propionyl[K](Delta_H(4)C(3)O(1)[K])'])

    # 目标氨基酸，一条肽仅一个
    # 注意：如果设置了目标修饰，则不需要设置目标氨基酸
    modi_one = set([])
    aa_one_ls = ['H', 'F', 'Y', 'W']
    aa2chrc = {'H': 'H-Im', 'F': 'F-Im', 'Y': 'Y-Im', 'W': 'W-Im'} # 与亚胺离子对照

    # 除目标修饰外，仅考虑的修饰
    # 若考虑全部开放式修饰, 则modi_ls = None
    modi_ls = set(['Carbamidomethyl[C]', 'Oxidation[M]'])
    # modi_ls = None

    # 酶切类型
    enz_type = [3] # 仅特异
    # enz_type = [0, 1, 2, 3] # 0非特异 1仅C端特异，2仅N端特异，3特异

    # 序列覆盖度过滤
    ls_filter_cov = [] # 不过滤
    ls_filter_cov = [0, 0, 0.7] # b离子单独、y离子单独、by离子合并

    # 碎片离子匹配误差
    frag_ppm = 20 # ppm

    # 匹配是使用进程数, 速度慢可以调大一点, 但不要超过cpu总线程数
    n_process = 10

    # 特征离子信息Characteristic ion
    # 注意：质量不带电荷，即不加质子质量
    # dct_chrc = {'LinIm':156.1262541, 'CycIm':139.0997071}
    # 修饰：Propionyl[K](Delta_H(4)C(3)O(1)[K]), LinIm ion：C8H17N2O, MW =157.13409, CycIm ion：C8H14NO, MW =140.10754  # 一价+氢原子

    dct_chrc = {'H-Im':109.0639918, 'F-Im':119.0734946, 'Y-Im':135.0684087, 'W-Im':158.0843924}
    # H:C5H8N3 m/z=110.0718188, F:C8H10N m/z=120.08132, Y:C8H10NO m/z=136.076235, W:C10H11N2 m/z=159.0922186  # 一价+氢原子

    # 碎片离子类型模板路径（非必须修改）
    fpath_iontype_temp = str(Path.cwd()/'utils2'/'iontype_ini_used_linear_template.csv')

    # 碎片离子类型生成路径（非必须修改）
    fpath_iontype = fpath_iontype_temp.replace('_template.csv', '_use.csv')

    # ----------------------------------------------------
    if modi_one:
        flow_1(
            fpath_iontype_temp=fpath_iontype_temp,
            fpath_iontype=fpath_iontype,
            dct_chrc=dct_chrc,
            fpath_res=fpath_res,
            fpath_out=fpath_out,
            modi_one=modi_one,
            modi_ls=modi_ls,
            aa_one=set([]),
            enz_type=enz_type,
            pep_level=pep_level,
            dpath_mgf=dpath_mgf,
            fpath_out_match=fpath_out_match,
            n_process=n_process,
            frag_ppm=frag_ppm,
            ls_filter_cov=ls_filter_cov
        )
    elif aa_one_ls:
        for aa_name in aa_one_ls:
            print('='*24+f'Processing {aa_name}'+'='*24)
            aa_one = set([aa_name])
            
            flow_1(
                fpath_iontype_temp=fpath_iontype_temp,
                fpath_iontype=fpath_iontype.replace('_use.csv', f'_use_{aa_name}.csv'),
                dct_chrc={aa2chrc[aa_name]: dct_chrc[aa2chrc[aa_name]]},
                fpath_res=fpath_res,
                fpath_out=fpath_out.replace('_label.spectra', f'_label_{aa_name}.spectra'),
                modi_one=set([]),
                modi_ls=modi_ls,
                aa_one=aa_one,
                enz_type=enz_type,
                pep_level=pep_level,
                dpath_mgf=dpath_mgf,
                fpath_out_match=fpath_out_match.replace('_fragMatch.csv', f'_fragMatch_{aa_name}.csv'),
                n_process=n_process,
                frag_ppm=frag_ppm,
                ls_filter_cov=ls_filter_cov
            )
    else:
        assert False, '请设置目标修饰或氨基酸'


if __name__ == '__main__':
    _main()
