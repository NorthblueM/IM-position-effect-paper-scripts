# -*- coding: UTF-8 -*-
# @Date     : 15th Jul, 2022
# @Author   : Northblue

"""plink xl results ion match
search not 15N labeled
"""

from tqdm import tqdm
import sys
import os
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import spec_exp_file
import spec_theo_file
import ion_match

from multiprocessing import Process, Queue

sys.path.append(os.path.dirname(os.getcwd()))
from tools import df_apply_mutli as dfm

from pathlib import Path

def get_rawname_pstudio(fpath_mgf):
    """raw file name"""
    fname_mgf = os.path.split(fpath_mgf)[1]
    if '_HCDFT.' in fname_mgf:
        return os.path.split(fpath_mgf)[1][:-10]
    if '.HCD.FTMS.' in fname_mgf: # pXtract
        return os.path.split(fpath_mgf)[1][:-13]


def get_title_rawname_pparse(title):
    """raw file name from pParse format title
    """
    return '.'.join(title.split('.')[:-5])


def build_mgf_idx_title2line(path):

    if path.endswith('.mgf'):
        return build_mgf_idx_title2line_fpath(path)
    else:
        return build_mgf_idx_title2line_dpath(path)


def build_mgf_idx_title2line_dpath(dpath):
    """
    Args:
        dpath_mgf

    Return:
        {raw_name:{title:[start_id, end_id],},}
    """

    raw_idx_title2line = {}
    fnames_mgf = [x for x in list(os.walk(dpath))[0][2] if os.path.splitext(x)[1] == '.mgf']
    for fname in fnames_mgf:
        fpath = os.path.join(dpath, fname)
        raw_idx_title2line = {**raw_idx_title2line, **build_mgf_idx_title2line_fpath(fpath)}

    return raw_idx_title2line


def build_mgf_idx_title2line_fpath(fpath):
    """
    Args:
        fpath_mgf

    Returns:
        {raw_name:{title:[start_id, end_id],},}
    """

    raw_idx_title2line = {}
    raw_name = get_rawname_pstudio(fpath)
    idx_title2line = spec_exp_file.get_mgf_title2line_idx(fpath)
    raw_idx_title2line[raw_name] = idx_title2line

    return raw_idx_title2line


def add_mgf_idx_title2line(fpath_mgf, raw_idx_title2line):
    """
    """

    raw_name = get_rawname_pstudio(fpath_mgf)
    idx_title2line = spec_exp_file.get_mgf_title2line_idx(fpath_mgf)
    raw_idx_title2line[raw_name] = idx_title2line

    return raw_idx_title2line


def build_res_mgf_idx_title2line(df_res, mgf2fpath):
    """
    The title is all in the result
    """
    raw_idx_title2line = {}
    rawname_ls = set(list(df_res['Title'].apply(lambda x:get_title_rawname_pparse(x))))
    for rawname in rawname_ls:
        if rawname not in raw_idx_title2line:
            print('===index:', rawname)
            raw_idx_title2line = add_mgf_idx_title2line(mgf2fpath[rawname], raw_idx_title2line)

    return raw_idx_title2line


def add_mgf_idx_title2offset(fpath_mgf, raw_idx_title2line):
    """
    """

    raw_name = get_rawname_pstudio(fpath_mgf)
    idx_title2line = spec_exp_file.get_mgf_title2offset_idx(fpath_mgf)
    raw_idx_title2line[raw_name] = idx_title2line

    return raw_idx_title2line


def build_res_mgf_idx_title2offset(df_res, mgf2fpath, n_process=0):
    """
    The title is all in the result
    """
    raw_idx_title2line = {}

    col_title = 'Title'
    if 'Title' not in df_res.columns:
        col_title = 'File_Name'
    rawname_ls = set(list(df_res[col_title].apply(lambda x:get_title_rawname_pparse(x))))
        

    # 单进程
    if n_process <= 1:
        for rawname in rawname_ls:
            if rawname not in raw_idx_title2line:
                print('===index:', rawname)
                raw_idx_title2line = add_mgf_idx_title2offset(mgf2fpath[rawname], raw_idx_title2line)
    else:
        rawname_ls = [x for x in rawname_ls if x not in raw_idx_title2line]
        param_ls = []
        for rawname in rawname_ls:
            print('===index:', rawname)
            param_ls.append([mgf2fpath[rawname]])
        ls = dfm.multi_pool(spec_exp_file.get_mgf_title2offset_idx, param_ls, n_process=n_process)
        raw_idx_title2line = dict(zip(rawname_ls, ls))

    return raw_idx_title2line


def build_mgf2fpath(dpath):
    """
    Return:
        {raw_name:fpath_mgf}
    """

    mgf2fpath = {}
    fnames_mgf = [x for x in list(os.walk(dpath))[0][2] if os.path.splitext(x)[1] == '.mgf']
    for fname in fnames_mgf:
        fpath = os.path.join(dpath, fname)
        raw_name = get_rawname_pstudio(fpath)
        mgf2fpath[raw_name] = fpath

    return mgf2fpath


def match_spec_xl_15n_one(ser, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath):
    """ser: plink one line results
    Parameter for object, mainly save df_iontype generation time
    """

    title = ser['Title']
    raw_name = get_title_rawname_pparse(title)
    fpath_mgf = mgf2fpath[raw_name]

    if raw_name not in raw_idx_title2line:
        raw_idx_title2line = add_mgf_idx_title2line(fpath_mgf, raw_idx_title2line)
    idx_title2line = raw_idx_title2line[raw_name]
    
    spec_exp.set_spec_mgf(fpath_mgf, title, idx_title2line)
    spec_exp.set_inten_relat()
    df_peak = spec_exp.df_peak

    pep = ser['Peptide']
    mods = ser['Modifications']
    linker = ser['Linker']
    prec_charge = ser['Charge']
    
    label_id = ser['LabelID']
    spec_theo.update_mass_labeled(label_id)

    spec_theo.set_pep_plink(pep, mods, linker)
    spec_theo.cal_theo_mz()
    spec_theo.set_df_iontype_str()
    spec_theo.set_df_iontype_charge()
    df_ion = spec_theo.df_ion

    df_peak_match = ion_match.match_spec(df_peak, df_ion, keep_unmatched=False, iontype_str=True, prec_charge=prec_charge)

    return df_peak_match


def match_plink_xl_15n_1p(df_res, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath, q, progress_visual=False):
    """1 process
    """


    if progress_visual:
        iterm = tqdm(list(df_res.index))
    else:
        iterm = list(df_res.index)
    for i in iterm:
        ser = df_res.loc[i]
        df_peak_match = match_spec_xl_15n_one(ser, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath)
        df_peak_match['res_id'] = i
        q.put(df_peak_match)


def match_plink_xl_15n(dpath_mgf, fpath_res, fpath_out='', progress_visual=False, n_process=1, head=0):
    """
    """

    df_res = pd.read_csv(fpath_res)
    if head:
        df_res = df_res.iloc[:head] # for testing

    mgf2fpath = build_mgf2fpath(dpath_mgf)

    # raw_idx_title2line = build_mgf_idx_title2line(dpath_mgf)
    raw_idx_title2line = build_res_mgf_idx_title2line(df_res, mgf2fpath)
    
    spec_exp = spec_exp_file.SpecExp()
    spec_theo = spec_theo_file.SpecTheoXl()

    if progress_visual:
        print('===log: spec ion peak matching...')
    df_res = df_res.reset_index(drop=True)
    ls = []
    # ======one process
    if n_process == 1:
        if progress_visual:
            iterm = tqdm(range(df_res.shape[0]))
        else:
            iterm = range(df_res.shape[0])
        for i in  iterm:
            ser = df_res.iloc[i]
            df_peak_match = match_spec_xl_15n_one(ser, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath)

            df_peak_match['res_id'] = i
            ls.append(df_peak_match)
    
    # ======multi process
    else:
        q = Queue()
        interval = int(df_res.shape[0]/n_process)
        process_ls = []
        for i in range(n_process):
            start = i*interval
            if i != n_process-1: # not last process
                end = (i+1)*interval
            else:
                end = df_res.shape[0]
            p = Process(target=match_plink_xl_15n_1p, args=(df_res.iloc[start:end], spec_exp, spec_theo, raw_idx_title2line, mgf2fpath, q), kwargs={'progress_visual':progress_visual})
            process_ls.append(p)
            p.start()
        
        # must be consumed in advance
        for p in process_ls:
            while(p.is_alive()):
                while(not q.empty()):
                    ls.append(q.get())
    
        for p in process_ls:
            p.join()
    
        ls.extend([q.get() for i in range(q.qsize())])
    
    if len(ls) != df_res.shape[0]:
        print('match_plink_xl_15n() error!!!')
    
    df_match = pd.concat(ls)
    df_match = df_match.sort_values(['res_id', 'mz'])

    if fpath_out:
        df_match.to_csv(fpath_out, index=False)
    return df_match



def match_spec_xl_iontype_1mgf_one(ser, spec_exp, spec_theo, idx_title2line, fpath_mgf):
    """ser: plink one line results
    Parameter for object, mainly save df_iontype generation time

    iontype from read .csv
    one mgf include multi raw
    """

    # exp spec
    title = ser['Title']
    
    spec_exp.set_spec_mgf(fpath_mgf, title, idx_title2line)
    spec_exp.set_inten_relat()
    df_peak = spec_exp.df_peak

    # theo spec
    pep = ser['Peptide']
    mods = ser['Modifications']
    prec_charge = ser['Charge']

    spec_theo.set_pep_plink(pep, mods)
    spec_theo.cal_theo_mz()
    spec_theo.set_df_iontype_str()
    spec_theo.set_df_iontype_charge()
    df_ion = spec_theo.df_ion

    # match
    df_peak_match = ion_match.match_spec(df_peak, df_ion, keep_unmatched=False, iontype_str=True, prec_charge=prec_charge)

    return df_peak_match


def match_plink_xl_iontype_1mgf(fpath_mgf, fpath_res, fpath_iontype, fpath_out='', progress_visual=False, n_process=1, head=0):
    """

    iontype from read .csv
    one mgf include multi raw
    """

    df_res = pd.read_csv(fpath_res)
    if head:
        df_res = df_res.iloc[:head] # for testing

    # custom mgf
    raw_idx_title2line = build_mgf_idx_title2line(fpath_mgf)
    idx_title2line = list(raw_idx_title2line.values())[0]
    
    spec_exp = spec_exp_file.SpecExp()
    spec_theo = spec_theo_file.SpecTheoXl()

    # custom iontype
    spec_theo.read_df_iontype(fpath_iontype)

    if progress_visual:
        print('===log: spec ion peak matching...')
    df_res = df_res.reset_index(drop=True)
    ls = []
    # ======one process
    if n_process == 1:
        if progress_visual:
            iterm = tqdm(range(df_res.shape[0]))
        else:
            iterm = range(df_res.shape[0])
        for i in  iterm:
            ser = df_res.iloc[i]
            df_peak_match = match_spec_xl_iontype_1mgf_one(ser, spec_exp, spec_theo, idx_title2line, fpath_mgf)

            df_peak_match['res_id'] = i
            ls.append(df_peak_match)
    
    # ======multi process
    else:
        assert False, 'not support multi process'
    
    if len(ls) != df_res.shape[0]:
        print('match_plink_xl_iontype() error!!!')
    
    df_match = pd.concat(ls)
    df_match = df_match.sort_values(['res_id', 'mz'])

    if fpath_out:
        df_match.to_csv(fpath_out, index=False)
    return df_match


def match_spec_pfind_linear_iontype_one(ser, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath):
    """ser: plink one line results
    Parameter for object, mainly save df_iontype generation time

    iontype from read .csv
    """

    # exp spec
    title = ser['Title']

    if isinstance(mgf2fpath, dict): # 多个mgf，文件名符合pFind直接导出的格式
        raw_name = get_title_rawname_pparse(title)
        fpath_mgf = mgf2fpath[raw_name]

        # to do (done√): first build all mgf_idx_title2line
        if raw_name not in raw_idx_title2line:
            raw_idx_title2line = add_mgf_idx_title2offset(fpath_mgf, raw_idx_title2line)
        idx_title2line = raw_idx_title2line[raw_name]
    else: # 只有一个mgf文件
        fpath_mgf = mgf2fpath
        idx_title2line = raw_idx_title2line

    spec_exp.set_spec_mgf_offset(fpath_mgf, title, idx_title2line)
    # spec_exp.set_inten_relat(mode='sum') # 相对强度
    spec_exp.set_inten_relat(mode='max') # 相对强度
    df_peak = spec_exp.df_peak

    # theo spec, pfind linear
    pep = ser['Sequence']
    mods = ser['Modification']
    prec_charge = ser['Charge']

    spec_theo.set_pep_pfind(pep, mods)
    spec_theo.cal_theo_mz()
    spec_theo.set_df_iontype_str()
    spec_theo.set_df_iontype_charge()
    df_ion = spec_theo.df_ion

    # match
    df_peak_match = ion_match.match_spec(df_peak, df_ion, keep_unmatched=False, iontype_str=True, prec_charge=prec_charge)

    return df_peak_match


def match_pfind_linear_iontype_1p(df_res, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath, q, progress_visual=False):
    """1 process
    """

    if progress_visual:
        iterm = tqdm(list(df_res.index))
    else:
        iterm = list(df_res.index) 
    for i in iterm: # 索引，而非行顺序
        ser = df_res.loc[i]
        df_peak_match = match_spec_pfind_linear_iontype_one(ser, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath)
        df_peak_match['res_id'] = i
        q.put(df_peak_match)


def match_pfind_linear_iontype(dpath_mgf, fpath_res, fpath_iontype, fpath_out='', progress_visual=False, n_process=1, head=0):
    """

    iontype from read .csv
    """

    if isinstance(fpath_res, str): # 路径
        df_res = pd.read_csv(fpath_res, sep='\t')
    else: # 已经是dataframe
        df_res = fpath_res
    
    col_title = 'Title'
    if 'Title' not in df_res.columns:
        col_title = 'File_Name'
        df_res['Title'] = df_res['File_Name']


    if head:
        df_res = df_res.iloc[:head] # for testing

    if Path(dpath_mgf).is_dir(): # 一个文件夹，mgf文件名符合pFind直接导出的格式
        mgf2fpath = build_mgf2fpath(dpath_mgf)

        # raw_idx_title2line = build_mgf_idx_title2line(dpath_mgf)
        raw_idx_title2line = build_res_mgf_idx_title2offset(df_res, mgf2fpath, n_process=n_process)
    else: # 只有一个mgf文件
        mgf2fpath = dpath_mgf
        print('===index:', Path(mgf2fpath).name)
        raw_idx_title2line = spec_exp_file.get_mgf_title2offset_idx(dpath_mgf)
    

    spec_exp = spec_exp_file.SpecExp()
    spec_theo = spec_theo_file.SpecTheoLinear()

    # custom iontype
    if isinstance(fpath_iontype, str): # 路径
        spec_theo.read_df_iontype(fpath_iontype)
    else: # 已经是dataframe
        spec_theo.df_iontype = fpath_iontype

    if progress_visual:
        print('===log: spec ion peak matching...')
    df_res = df_res.reset_index(drop=True)
    ls = []
    # ======one process
    if n_process == 1:
        if progress_visual:
            iterm = tqdm(range(df_res.shape[0]))
        else:
            iterm = range(df_res.shape[0])
        for i in  iterm:
            ser = df_res.iloc[i]
            df_peak_match = match_spec_pfind_linear_iontype_one(ser, spec_exp, spec_theo, raw_idx_title2line, mgf2fpath)

            df_peak_match['res_id'] = i
            ls.append(df_peak_match)
    
    # ======multi process
    else:
        q = Queue()
        interval = int(df_res.shape[0]/n_process)
        process_ls = []
        for i in range(n_process):
            start = i*interval
            if i != n_process-1: # not last process
                end = (i+1)*interval
            else:
                end = df_res.shape[0]
            p = Process(target=match_pfind_linear_iontype_1p, args=(df_res.iloc[start:end], spec_exp, spec_theo, raw_idx_title2line, mgf2fpath, q), kwargs={'progress_visual':progress_visual})
            process_ls.append(p)
            p.start()
        
        # must be consumed in advance
        for p in process_ls:
            while(p.is_alive()):
                while(not q.empty()):
                    ls.append(q.get())
    
        for p in process_ls:
            p.join()
    
        ls.extend([q.get() for i in range(q.qsize())])
    
    if len(ls) != df_res.shape[0]:
        s = 'match_pfind_linear_iontype() error!!!'
        print(s)
        assert False, s

    df_match = pd.concat(ls)
    df_match = df_match.sort_values(['res_id', 'mz'])

    if fpath_out:
        df_match.to_csv(fpath_out, index=False)
    return df_match


def _test():
    dpath_mgf = r'F:\DataSet\pLink2_test\Leiker_15N_labeling\data'
    fpath_res = r'F:\DataSet\pLink2_test\Leiker_15N_labeling\mpz\pLink_task_2025.01.09.12.12.44\reports\uniprot-ecoli-20171023_2025.01.09.filtered_cross-linked_spectra.csv'
    match_plink_xl_15n(dpath_mgf, fpath_res, head=200, progress_visual=True)


def _main():
    dpath_mgf = r'F:\DataSet\pLink2_test\Leiker_15N_labeling\data'
    fpath_res = r'F:\DataSet\pLink2_test\Leiker_15N_labeling\mpz\pLink_task_2025.01.09.12.12.44\reports\uniprot-ecoli-20171023_2025.01.09.filtered_cross-linked_spectra.csv'
    fpath_out = fpath_res[:-4]+'_ionmatch.csv'

    # match_plink_xl_15n(dpath_mgf, fpath_res, progress_visual=True)
    match_plink_xl_15n(dpath_mgf, fpath_res, fpath_out, n_process=2, progress_visual=True)


if __name__ == '__main__':
    _main()
