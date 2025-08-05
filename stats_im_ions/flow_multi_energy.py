# -*- coding: UTF-8 -*-
# Date      : 7th April, 2025
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


def get_nce_str(x):
    """NCE表示为字符串
    阶梯能量用,分割"""

    if isinstance(x, float):
        # 可能是nan
        if np.isnan(x):
            return '0'
        x = str(int(x))
    
    nce_str = ','.join([str(int(float(x))) for x in x.split(',')])
    return nce_str


def generate_nce_align(dpath_mgf, fpath_res, fpath_out, nce_trigger=35, nce_ls =[]):
    """生成NCE对齐的结果文件
    
    
    确认NCE枚举的能量都采集到, 例如一共采集了7个
    触发的能量和枚举的能量, 鉴定肽段若不一致, 进行投票, 相同投票选打分高的
    同一个肽段鉴定, 广播给所有能量
    """

    # 读取pfc, 记录了NCE、Master Scan Number
    df_pfc_ls = []
    for fp in [x for x in Path(dpath_mgf).iterdir() if x.name.endswith('.pfc')]:
        df_pfc = pd.read_csv(fp, sep='\t')
        df_pfc['rawname'] = fp.stem
        df_pfc_ls.append(df_pfc)
    df_pfc = pd.concat(df_pfc_ls)
    df_pfc = df_pfc.reset_index(drop=True)
    
    # 只保留MS2
    df_pfc = df_pfc[df_pfc['SpectrumType'] == 'MS2']
    df_pfc['NCE'] = df_pfc['NCE'].apply(get_nce_str)
    df_pfc['PrecursorScan'] = df_pfc['PrecursorScan'].astype(int) # 一定要int, 否则字符串转换时会带小数点
    df_pfc['ScanNo'] = df_pfc['ScanNo'].astype(int)

    # rawscan
    df_pfc['scannum'] = df_pfc['ScanNo']
    df_pfc = util_pfind.set_rawscan(df_pfc)

    # ---读取_HCDFT_precursor_score.csv, 确认存在的title
    df_prec_ls = []
    for f in [x for x in Path(dpath_mgf).iterdir() if x.name.endswith('_HCDFT_precursor_score.csv')]:
        df_prec = pd.read_csv(f)
        df_prec_ls.append(df_prec)
    df_prec = pd.concat(df_prec_ls)
    df_prec = df_prec.reset_index(drop=True)

    df_prec = util_pfind.set_split_tilte(df_prec, col_title='title', mode='pParse')
    df_prec = util_pfind.set_rawscan(df_prec)

    # 必须pParse有导出谱图
    _n_rawscan_before = df_pfc['rawscan'].nunique()
    _set_rawscan = set(df_prec['rawscan'].unique())
    df_pfc = df_pfc[df_pfc['rawscan'].isin(_set_rawscan)]
    _n_rawscan_after = df_pfc['rawscan'].nunique()
    print(f"---filter pParse extract #rawscan: {_n_rawscan_before} -> {_n_rawscan_after}")


    # ---过滤, 触发的能量必须全都有
    def get_master_scan(x):
        """不同raw, 返回rawname+scan"""

        # 触发的MS2为原MS2 scan
        if x['NCE'] == nce_trigger:
            return '.'.join([x['rawname'], str(x['ScanNo'])])
        return '.'.join([x['rawname'], str(x['PrecursorScan'])])
    # df_pfc = df_pfc[df_pfc['NCE'] != nce_trigger]
    df_pfc['rawscan_master'] = df_pfc.apply(get_master_scan, axis=1)

    print(f"{'-'*6}acquisition-trigger spec number (include base-trigger)")
    print(df_pfc['rawscan_master'].value_counts().value_counts())

    # df_prec增加NCE信息
    _dct = dict(zip(df_pfc['rawscan'], df_pfc['NCE']))
    df_prec['NCE'] = df_prec['rawscan'].map(_dct)
    _dct = dict(zip(df_pfc['rawscan'], df_pfc['rawscan_master']))
    df_prec['rawscan_master'] = df_prec['rawscan'].map(_dct)
    _dct = dict(zip(df_pfc['rawscan'], df_pfc['PrecursorScan']))
    df_prec['PrecursorScan'] = df_prec['rawscan'].map(_dct)

    # 过滤，采集的触发的能量必须全都有
    col_grp = 'rawscan_master'
    _dct = df_pfc.groupby([col_grp])['NCE'].agg(lambda x: len(set(x))).to_dict()
    df_pfc['trigger_num'] = df_pfc[col_grp].map(_dct)
    df_pfc = df_pfc[df_pfc['trigger_num'] == len(nce_ls)+1]
    print(df_pfc['trigger_num'].value_counts())



    # ---------读取pfind结果
    df_res = pd.read_csv(fpath_res, sep='\t')

    df_res = util_pfind.set_split_tilte(df_res, col_title='File_Name', mode='pParse')
    df_res = util_pfind.set_rawscan(df_res)
    print(f"---pfind result flt: {df_res.shape[0]}")

    df_res['Title'] = df_res['File_Name']
    df_res = util_res.set_pfind_level_str(df_res, is_i2l=True, is_title=True)

    # 过滤，trigger采集时是全的
    _set_rawscan = set(df_pfc['rawscan'].unique())
    df_res = df_res[df_res['rawscan'].isin(_set_rawscan)]
    print(f"---pfind result flt acquisition-trigger-all-nce: spec {df_res.shape[0]}, peptide {df_res['pep_str'].nunique()}")

    scan2master = dict(zip(df_pfc['rawscan'], df_pfc['rawscan_master']))
    df_res['rawscan_master'] = df_res['rawscan'].map(scan2master)

    scan2nce = dict(zip(df_pfc['rawscan'], df_pfc['NCE']))
    df_res['NCE'] = df_res['rawscan'].map(scan2nce)
    print(df_res['NCE'].value_counts())

    # 过滤，触发母NCE一定被鉴定
    _grp = df_res.groupby(['rawscan_master'])['NCE']
    print(f"---1 master scan, ID NCE num")
    print(_grp.apply(lambda x: x.nunique()).value_counts())

    # _dct = _grp.apply(lambda x: 1 if nce_trigger in set(x) else 0).to_dict()
    # df_res['is_trigger'] = df_res['rawscan_master'].map(_dct)
    # df_res = df_res[df_res['is_trigger'] == 1]
    # print(f"---pfind result flt NCE-{nce_trigger} ID, spec: {df_res.shape[0]}, peptide: {df_res['pep_str'].nunique()}")


    # 每个master scan, 投票选出一个peptide，包括全扫描的NCE
    def get_master_pep(x):
        
        pep_vcnt = x['pep_str'].value_counts()
        # index排序
        pep_vcnt = pep_vcnt.sort_index(ascending=False)
        # 选出投票最多的肽段
        pep_max = pep_vcnt.index[0]
        pep_cnt = pep_vcnt.values[0]

        # return [pep_str, pep_cnt]
        return {'pep_max': pep_max, 'pep_cnt': pep_cnt}
    df_pep = df_res.groupby(['rawscan_master']).apply(get_master_pep, include_groups=False).apply(pd.Series).reset_index(drop=False)
    print(f"{'-'*6} select 1trigger-1pep pep-max, spec: {df_pep.shape[0]}, pep: {df_pep['pep_max'].nunique()}")

    df_pep = df_pep.set_index('rawscan_master')
    df_res[['pep_max', 'pep_cnt']] = df_res['rawscan_master'].apply(lambda x:df_pep.loc[x][['pep_max', 'pep_cnt']])

    # df_res过滤，一个trigger选出投票最多的肽段
    df_res = df_res[df_res['pep_str'] == df_res['pep_max']]
    print(f"{'-'*6} 1pep-1master before, spec: {df_res.shape[0]}, master: {df_res['rawscan_master'].nunique()}, pep: {df_res['pep_str'].nunique()}")

    # 一个肽段被触发多次, 选不同NCE鉴定更多的; 鉴定一样多，选分数高的
    _df_drop = df_res.sort_values(['pep_cnt', 'Final_Score', 'Raw_Score'], ascending=[False, True, False])
    _df_drop = _df_drop.drop_duplicates(subset=['pep_str'], keep='first')
    _set_master = set(_df_drop['rawscan_master'].unique())
    df_res = df_res[df_res['rawscan_master'].isin(_set_master)]
    print(f"{'-'*6} 1pep-1master filter, spec: {df_res.shape[0]}, master: {df_res['rawscan_master'].nunique()}, pep: {df_res['pep_str'].nunique()}")

    # ------trigger挑选肽段，广播到所有NCE
    # 肽谱匹配读取的母离子电荷是ser['Charge']

    # 只包括触发的NCE
    _set_master = set(df_res['rawscan_master'].unique())
    _df_prec = df_prec[df_prec['rawscan_master'].isin(_set_master)]

    # {master scan: {nce: {charge: [title]}}}
    master2nce2title = {}
    for i in range(_df_prec.shape[0]):
        row = _df_prec.iloc[i]
        rawscan_master = row['rawscan_master']
        nce = row['NCE']
        title = row['title']
        charge = row['charge']

        if rawscan_master not in master2nce2title:
            master2nce2title[rawscan_master] = {}
        if nce not in master2nce2title[rawscan_master]:
            master2nce2title[rawscan_master][nce] = {}
        if charge not in master2nce2title[rawscan_master][nce]:
            master2nce2title[rawscan_master][nce][charge] = []
        master2nce2title[rawscan_master][nce][charge].append(title)
    
    nce_all = [nce_trigger]
    nce_all.extend(nce_ls)
    def broadcast_nce(x):
        """广播到所有NCE"""

        nce_id = list(x['NCE'])
        nce_broadcast = list(set(nce_all).difference(set(nce_id)))

        ser_temp = x.iloc[0]
        rawscan_master = ser_temp['rawscan_master']
        charge = ser_temp['Charge']

        ls_ls = []
        # 第一行复制n份
        for nce in nce_broadcast:
            _ser = ser_temp.copy()
            _ser['NCE'] = nce
            _ser['Charge'] = charge

            # 尽量找charge一致的title
            if charge in master2nce2title[rawscan_master][nce]:
                _title = master2nce2title[rawscan_master][nce][charge][0]
            else:
                # charge不一致，随机找一个
                c = list(master2nce2title[rawscan_master][nce].keys())[0]
                _title = master2nce2title[rawscan_master][nce][c][0]

                # raw中没有记录电荷, pfc中默认upperCharge=3,lowerCharge=2, pfb中只有2+，pParse2Plus合并，只合并了2+，主scan可能pParse导出了3+被鉴定
                # print(f"---Warning: {rawscan_master} {nce} charge {charge} not found, use {c} instead")

            _ser['File_Name'] = _title
            
            ls_ls.append(_ser.to_list())
        
        if ls_ls:
            cols = list(x.columns)
            df_new = pd.DataFrame(ls_ls, columns=cols)
            x = pd.concat([x, df_new], ignore_index=True)

        return x
    
    cols = list(df_res.columns)
    df_nce = df_res.groupby(['rawscan_master'])[cols].apply(lambda x: broadcast_nce(x)).reset_index(drop=True)
    
    # 重新设置rawscan
    df_nce['Title'] = df_nce['File_Name']
    df_nce = util_pfind.set_split_tilte(df_nce, is_print=False, col_title='File_Name', mode='pParse')
    df_nce = util_res.set_pfind_level_str(df_nce, is_i2l=True, is_title=True)

    print(f"{'-'*6} broadcast all nce, spec: {df_nce.shape[0]}, master: {df_nce['rawscan_master'].nunique()}, pep: {df_nce['pep_str'].nunique()}")

    # Charge和charge_title电荷是否相等
    df_nce['charge_equal'] = (df_nce['Charge'] == df_nce['charge_title']).astype(int)
    print(df_nce['charge_equal'].value_counts())

    if fpath_out:
        df_nce.to_csv(fpath_out, sep='\t', index=False)

    return df_nce


def calc_coverage(fpath_res, fpath_match, dpath_mgf, ls_filter_cov=[]):
    """计算覆盖度"""


    # ---读取_HCDFT_head.csv, 谱图强度加和
    df_head_ls = []
    for f in [x for x in Path(dpath_mgf).iterdir() if x.name.endswith('_HCDFT_head.csv')]:
        df_head = pd.read_csv(f)
        df_head_ls.append(df_head)
    df_head = pd.concat(df_head_ls)
    df_head = df_head.reset_index(drop=True)
    title2intensum = dict(zip(df_head['title'], df_head['inten_sum']))


    print('-'*12+'calc_coverage'+'-'*12)

    df = pd.read_csv(fpath_res, sep='\t')
    df_match = pd.read_csv(fpath_match)

    df_match['match_info'] = df_match['match_info'].apply(eval)

    df['pep_len'] = df['Sequence'].apply(len)

    def calc_cov(x):

        res_id = x['res_id']
        pep_len = x['pep_len']

        # 平均每位点匹配谱峰数
        n_b = 0
        n_y = 0

        # 匹配谱峰强度占比
        # TODO: 匹配谱峰数占比
        inten_b = 0
        inten_y = 0

        # 序列覆盖率
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
                    n_b += 1
                    # n_b += ser['inten_relat']
                    inten_b += ser['inten']
                elif sym == 'y':
                    ls_y[frag_site] = 1
                    n_y += 1
                    # n_y += ser['inten_relat']
                    inten_y += ser['inten']
        cov_b = sum(ls_b)/len(ls_b)
        cov_y = sum(ls_y)/len(ls_y)
        ls_by = [x or y for x, y in zip(ls_b, ls_y)]
        cov_by = sum(ls_by)/len(ls_by)
        
        n_by = n_b + n_y
        n_b /= (pep_len-1)
        n_y /= (pep_len-1)
        n_by /= (pep_len-1)

        title = x['File_Name']
        inten_by = inten_b + inten_y
        inten_b = inten_b/title2intensum[title]
        inten_y = inten_y/title2intensum[title]
        inten_by = inten_by/title2intensum[title]
        
        cov_inten_b = cov_b * inten_b
        cov_inten_y = cov_y * inten_y
        cov_inten_by = cov_by * inten_by
        
        return [cov_b, cov_y, cov_by, n_b, n_y, n_by
                , inten_b, inten_y, inten_by
                , cov_inten_b, cov_inten_y, cov_inten_by]

    ls_ls = []
    for i in tqdm(range(df.shape[0])):
        ls_ls.append(calc_cov(df.iloc[i]))

    cols = ['cov_b', 'cov_y', 'cov_by', 'n_b', 'n_y', 'n_by'
            , 'inten_b', 'inten_y', 'inten_by'
            , 'cov_inten_b', 'cov_inten_y', 'cov_inten_by']
    df[cols] = pd.DataFrame(ls_ls, columns=cols)

    if ls_filter_cov:
        b = (df['cov_b'] > ls_filter_cov[0])
        b = b & (df['cov_y'] > ls_filter_cov[1])
        b = b & (df['cov_by'] > ls_filter_cov[2])
        n1 = df.shape[0]
        df = df[b]
        n2 = df.shape[0]
        print(f'---Coverage {ls_filter_cov} filter, before: {n1}, after: {n2}')

    df.to_csv(fpath_res, sep='\t', index=False)

    return df


def filter_csv(fpath_res, fpath_out, modi_one, modi_ls, enz_type):
    ## 筛选csv行

    df = pd.read_csv(fpath_res, sep='\t')

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
    df['is_modi'] = df['modi_dct'].apply(is_modi_one)

    num_before = df.shape[0]

    df = df[df['is_modi'] == 1]

    num_after = df.shape[0]
    print(f'------Before filter: {num_before}, after filter: {num_after}')
    
    if fpath_out:
        df.to_csv(fpath_out, sep='\t', index=False)

    return df



def res_info(fpath_res, fpath_match, modi_one, dct_chrc):
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

    df['modi_site'] = df['modi_dct'].apply(get_modi_site)  # 计数从1开始


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
            ls.append(0) # 不存在，则强度为0
            ls.append(0) # 不存在，则强度为0

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
        if ls_by_inten:
            match_median = np.median(ls_by_inten)
        else: # 如果没有匹配上任何谱峰
            match_median = 1
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
    """统计特征离子出现的比例
    
    根据能量划分
    """

    df = pd.read_csv(fpath_res, sep='\t')
    ls_nce = df['NCE'].unique()
    ls_nce = sorted(ls_nce, reverse=False)

    cols_q = ['cov_b', 'cov_y', 'cov_by'
              , 'cov_inten_b', 'cov_inten_y', 'cov_inten_by']

    def stat_num(df):
        num = df.shape[0]
        print('---'*3, f'All num {num}')

        for c in dct_chrc.keys():
            print('---'*3, c)

            num = df[f'is_{c}'].sum()
            ratio = num/df.shape[0]
            print(f'---num {num}, ratio {ratio:.4f}')

            s_p = '------NCE n r, \t inten: q1 q2 q3 \t q1m q2m q3m \t\t'
            s_p += '\t\t  '.join(cols_q)
            print(s_p)

            for nce in ls_nce:
                _df1 = df[df['NCE'] == nce]
                n = sum(_df1[f'is_{c}'])
                r = n/_df1.shape[0]

                _df2 = _df1.copy()
                q1 = _df2[f'inten_{c}'].quantile(0.25)
                q2 = _df2[f'inten_{c}'].quantile(0.5)
                q3 = _df2[f'inten_{c}'].quantile(0.75)
            
                q1m = _df2[f'inten_{c}_match_median'].quantile(0.25)
                q2m = _df2[f'inten_{c}_match_median'].quantile(0.5)
                q3m = _df2[f'inten_{c}_match_median'].quantile(0.75)

                print(f'---{nce} {n} {r:.3f} | {q1:.3f} {q2:.3f} {q3:.3f} | {q1m:.3f} {q2m:.3f} {q3m:.3f}', end='|')

                q_ls = [0.25, 0.5, 0.75]
                p_ls = ''
                for col in cols_q:
                    _ls = []

                    # q1 = _df2[col].quantile(0.25)
                    # q2 = _df2[col].quantile(0.5)
                    # q3 = _df2[col].quantile(0.75)

                    # # 1.5倍四分位距法, 最小值和最大值
                    # k = 1.5
                    # iqr = q3 - q1
                    # q0 = q1 - k*iqr
                    # q4 = q3 + k*iqr

                    # mi = _df2[col].min()
                    # ma = _df2[col].max()
                    
                    # _ls.append(mi)
                    # _ls.append(q2)
                    # _ls.append(ma)

                    for q in q_ls:
                        _ls.append(_df2[col].quantile(q))
                    p_ls += ' '.join([f'{x:.3f}' for x in _ls]) + ' | '
                print(f'{p_ls}')

                    

    stat_num(df)


def _main():
    ## 路径

    # mgf文件夹
    # dpath_mgf = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw'
    # dpath_mgf = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20240409_NCEshuffle'
    # dpath_mgf = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20250409_NCEsequential'
    # dpath_mgf = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20250418_stepHCD_35trigger'
    # dpath_mgf = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20250418_stepHCD_40trigger'
    dpath_mgf = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20250623_stepHCD_35trigger'

    # pfind结果路径
    # fpath_res = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\result\search_constrain\result\pFind-Filtered.spectra'
    # fpath_res = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\result\20240409_NCEshuffle_search_constrain\result\pFind-Filtered.spectra'
    # fpath_res = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\result\20250409_NCEsequential_search_constrain\result\pFind-Filtered.spectra'
    # fpath_res = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\result\20250418_stepHCD_35trigger_constrain\result\pFind-Filtered.spectra'
    # fpath_res = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\result\20250418_stepHCD_40trigger_constrain\result\pFind-Filtered.spectra'
    fpath_res = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\result\20250623_stepHCD_35trigger_constrain\result\pFind-Filtered.spectra'

    # 输出路径（非必须修改）
    fpath_out_flt = fpath_res.replace('.spectra', '_flt.spectra') # 过滤仅修饰结果
    fpath_out_NCE = fpath_res.replace('.spectra', '_NCE.spectra') # NCE对齐结果

    # 输出碎片离子匹配路径（非必须修改）
    fpath_out_match = str(Path(fpath_out_NCE).parent/(Path(fpath_out_NCE).stem+'_fragMatch.csv'))

    # 目标修饰，以条肽仅一个
    modi_one = set(['Xlink_BS2G[114][K](Gluratylation[K])'])

    # 除目标修饰外，仅考虑的修饰
    # 若考虑全部开放式修饰, 则modi_ls = None
    modi_ls = set(['Carbamidomethyl[C]', 'Oxidation[M]', 'Acetyl[ProteinN-term]'])
    # modi_ls = None

    # 酶切类型
    enz_type = [3] # 仅特异
    # enz_type = [0, 1, 2, 3] # 0非特异 1仅C端特异，2仅N端特异，3特异


    # 触发的能量
    nce_trigger = '35'
    # nce_trigger = '40'
    

    # 连续采集的能量
    # nce_ls = ['21', '24', '27', '30', '33', '36', '39']
    # nce_ls = ['24,39', '27', '30,39', '33,39', '24,42', '27,42', '30', '39']
    nce_ls = ['27', '30', '27,39', '30,39', '33,39', '27,42', '39', '42']
    # nce_ls = ['24,39', '27,39', '30,39', '33,39', '24,42', '27,42', '30', '39']

    # # 序列覆盖度过滤
    ls_filter_cov = [] # 不过滤
    # ls_filter_cov = [0, 0, 0.7] # b离子单独、y离子单独、by离子合并

    # 碎片离子匹配误差
    frag_ppm = 20 # ppm

    # 匹配是使用进程数, 速度慢可以调大一点, 但不要超过cpu总线程数
    n_process = 10

    # 特征离子信息Characteristic ion
    # 注意：质量不带电荷，即不加质子质量
    dct_chrc = {'Linlm':214.1317315, 'Cyclm':197.1051845} # Xlink_BS2G[114][K](Gluratylation[K])

    # 碎片离子类型模板路径（非必须修改）
    fpath_iontype_temp = str(Path.cwd()/'utils2'/'iontype_ini_used_linear_template.csv')

    # 碎片离子类型生成路径（非必须修改）
    fpath_iontype = fpath_iontype_temp.replace('_template.csv', '_use.csv')

    # ----------------------------------------------------
    df_iontype = generate_iontype(fpath_iontype_temp, fpath_iontype, dct_chrc)

    print(f"{'='*18} filter result {'='*18}")
    df_res = filter_csv(fpath_res, fpath_out_flt, modi_one, modi_ls, enz_type)

    print(f"{'='*18} generate nce align {'='*18}")
    df_nce = generate_nce_align(dpath_mgf, fpath_out_flt, fpath_out_NCE, nce_trigger=nce_trigger, nce_ls=nce_ls)

    print(f"{'='*18} match fragment {'='*18}")
    df_match = match_frag(dpath_mgf, fpath_out_NCE, fpath_iontype, fpath_out=fpath_out_match, progress_visual=True, n_process=n_process, match_ppm=frag_ppm)

    print(f"{'='*18} calc coverage {'='*18}")
    df_label = calc_coverage(fpath_out_NCE, fpath_out_match, dpath_mgf, ls_filter_cov)

    print(f"{'='*18} add info {'='*18}")
    df_label = res_info(fpath_out_NCE, fpath_out_match, modi_one, dct_chrc)

    print(f"{'='*18} stats res {'='*18}")
    stats_res(fpath_out_NCE, dct_chrc)



    

if __name__ == '__main__':
    _main()
