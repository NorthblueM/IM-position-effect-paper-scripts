#-*- coding: UTF-8 -*-

# @Date     : 30th Nov, 2023
# @Author   : Northblue

"""蛋白质谱数据处理, pFind相关, 工具函数
"""

import os
import shutil
import pandas as pd
from pathlib import Path

MASS_PROTON = 1.007276466621

def get_4rawname_fpath_pparse(rawname):
    return rawname+'_HCDFT.mgf'


def get_plink_param_fpath(dpath):
    """获取pLink参数文件路径

    Args:
        dpath: pLink结果文件夹路径
    """

    fnames_plink = []
    for name in os.listdir(dpath):
        if name.endswith('.plink'):
            fnames_plink.append(name)
    assert len(fnames_plink) == 1, '.plink超过一个'
    
    return str(Path(dpath)/fnames_plink[0])

def get_plink_report_xl_csm_fpath(dpath, dmid='reports', endswith_str='.filtered_cross-linked_spectra.csv'):
    """获取pLink的CSM结果文件csv路径

    Args:
        dpath: pLink结果文件夹路径, 不包括reports
    """

    if dmid:
        dpath = str(Path(dpath)/dmid)

    fnames_plink = []
    for fname in [p.name for p in Path(dpath).iterdir() if not p.is_dir()]:
        if fname.endswith(endswith_str):
            fnames_plink.append(fname)
    assert len(fnames_plink) == 1, f'Not unique {endswith_str} file'
    
    return str(Path(dpath)/fnames_plink[0])


def get_plink_report_total_fpath(dpath, dmid='reports', endswith_str='.filtered_cross-linked_spectra.csv'):
    """获取pLink的总表结果文件csv路径

    Args:
        dpath: pLink结果文件夹路径, 不包括reports
    """

    fpath_xl_csm =  get_plink_report_xl_csm_fpath(dpath, dmid=dmid, endswith_str=endswith_str)
    fpath_total = fpath_xl_csm[:-1*len(endswith_str)]+'.csv'

    if not Path(fpath_total).exists():
        assert False, 'Not exist '+fpath_total

    return fpath_total

def get_plink_report_csv_fpath(dpath, mid=r'search', suf=r'filtered_cross-linked_spectra.csv'):
    """获取pLink的csv结果文件路径

    Args:
        dpath: pLink结果文件夹路径, 不包括reports
        mid: 中间文件夹名
    """

    if mid:
        dpath = os.path.join(dpath, mid)
    dpath = os.path.join(dpath, 'reports')
    fnames = []
    for name in os.listdir(dpath):
        if name.endswith(suf):
            fnames.append(name)
    assert len(fnames) == 1, '超过一个'
    return str(Path(dpath)/fnames[0])


def get_plink_pfd_rank1_fpaths(dpath):
    """获取pLink结果tmps中的rank1的pfd路径集合

    Args:
        dpath: pLink结果文件夹路径
    """

    dpath_tmp = Path(dpath)/'tmps'
    
    fpaths_pfd_rank1 = []
    for name in os.listdir(dpath_tmp):
        segs = name.split('.')
        if len(segs)>=2 and segs[-1] == 'pfd' and segs[-2].startswith('File'):
            fpaths_pfd_rank1.append(str(Path(dpath_tmp)/name))
    
    return fpaths_pfd_rank1

def get_title_rawname(title, mode='pParse'):
    if mode == 'pParse':
        return '.'.join(title.split('.')[:-5])
    if mode == 'pXtract':
        return '.'.join(title.split('.')[:-4])
    else:
        assert False, 'missing mode'

def get_title_scannum(title, mode='pParse'):
    if mode == 'pParse':
        return int(title.split('.')[-4])
    if mode == 'pXtract':
        return int(title.split('.')[-3])
    else:
        assert False, 'missing mode'

def get_title_charge(title, mode='pParse'):
    if mode == 'pParse':
        return int(title.split('.')[-3])
    if mode == 'pXtract':
        return int(title.split('.')[-2])
    else:
        assert False, 'missing mode'

def get_title_pparseid(title):
    return int(title.split('.')[-2])


def get_title_split_pparse(title):
    segs = title.split('.')
    return ['.'.join(segs[:-5])
            , int(segs[-4])
            , int(segs[-3])
            , int(segs[-2])
            ]

def get_title_split_pxtract(title):
    segs = title.split('.')
    return ['.'.join(segs[:-4])
            , int(segs[-3])
            , int(segs[-2])
            ]

def split_title(df, col_title='Title', mode='pParse'):

    if mode == 'pParse':
        cols = ['rawname', 'scannum', 'charge_title', 'pparseid']
        df[cols] = pd.DataFrame(df[col_title].apply(get_title_split_pparse).to_list(), index=df.index)

    else:
        df['rawname'] = df[col_title].apply(lambda x:get_title_rawname(x, mode=mode))
        df['scannum'] = df[col_title].apply(lambda x:get_title_scannum(x, mode=mode))
        df['charge_title'] = df[col_title].apply(lambda x:get_title_charge(x, mode=mode))
        # df['pparseid'] = df[col_title].apply(lambda x:get_title_pparseid(x))

    return df


def set_rawscan(df):
    df['rawscan'] = df['rawname'].str.cat(df['scannum'].astype(str),sep='.')
    return df


def set_rawscan_charge(df):
    df['rawscan_c'] = df['rawname'].str.cat(df['scannum'].astype(str), sep='.').str.cat(df['charge_title'].astype(str), sep='.')
    return df


def set_split_tilte(df, is_print=False, **kwargs):

    df = split_title(df, **kwargs)
    df = set_rawscan(df)
    df = set_rawscan_charge(df)
    if is_print:
        cols = ['rawname', 'rawscan', 'rawscan_c', kwargs['col_title']]
        for col in cols:
            print(col, df[col].nunique())
    return df


def get_mh2mz(mh, charge):
    return (mh+(charge-1)*MASS_PROTON)/charge

def get_mz2mass(mz, charge, mass_proton=1.007276466621):
    """中性质量"""
    return (mz-mass_proton)*charge

def set_mz2mass(df, col_mz, col_charge, mass_proton=1.007276466621):
    df['mass'] = (df[col_mz]-mass_proton)*df[col_charge]
    return df

def cal_match_error(theo, exp):
    return abs(exp-theo)/theo*1e6
