# -*- coding: UTF-8 -*-
# Date      : Mar 26, 2022
# Author    : Northblue

"""pandas多进程"""

import multiprocessing as mp
from multiprocessing import Pool, Process
from functools import partial
from re import X
from tabnanny import verbose
import numpy as np
import pandas as pd
import joblib
from joblib import Parallel, delayed
from joblib import wrap_non_picklable_objects
import psutil


def is_joblib_parallel_subprocess(pid):
    """判断是joblib.parallel产生的子进程"""

    # # joblib.Parallel, 实际子进程
    # # ['c:\\ProgramStudy\\Anaconda3\\envs\\py3\\python.exe', '-c', 'from joblib.externals.loky.backend.popen_loky_win32 import main; main()', '--multiprocessing-fork', '2112']
    s ='from joblib.externals.loky.backend.popen_loky_win32 import main; main()'
    
    # # joblib.Parallel进程
    # # ['c:\\ProgramStudy\\Anaconda3\\envs\\py3\\python.exe', '-c', 'from joblib.externals.loky.backend.resource_tracker import main; main(500, False)']
    # s = 'from joblib.externals.loky.backend.resource_tracker import main; main('

    try:
        cmdline = psutil.Process(pid).cmdline()
        for param in cmdline:
            if s in param and psutil.Process(pid).name() == 'python.exe':
                return 1
    except:
        return 0
    return 0


def terminate_joblib_parallel_subprocess(subproc):

    for pid in subproc:
        if is_joblib_parallel_subprocess(pid):
            psutil.Process(pid).terminate()



def split_id(n_proess, num):
    """总共需要计算num次, 平均分配给n_proess线程
    [(start, end)]
        start从0开始
        end索引时不被选择
    """
    ls = []
    interval = int(num/n_proess)
    rem = num%n_proess #余数
    # 数据划分样例：012|345|67|89
    for i in range(n_proess):
        start = i*interval
        end = (i+1)*interval
        if i < rem:
            start += i
            end += (i+1)
        else:
            start += rem
            end += rem
        ls.append((start,end))
    return ls

def apply_1col(df, col, fun, args=[], kwargs={}, n_process=mp.cpu_count()-4, terminate_subproc=True):
    """df的apply变换, 多进程, 针对一列
    Parameter:
        col: 变换的列
        fun: 变换函数fun(x, *args)
        args: 函数参数
        n_process: 进程数
    Return:
        df[col_new]
    """
    current_process = psutil.Process()
    subproc_before = set([p.pid for p in current_process.children(recursive=True)])


    if n_process < 0: n_process = 1
    col_new = '__!#%@$^__'
    # col_i = '__@$^!#%__'
    # df[col_i] = list(range(df.shape[0]))
    def key_fun(df):
        df[col_new] = df[col].apply(lambda x:fun(x, *args, **kwargs))
        return df
    
    df_ls = []
    id_ls = split_id(n_process, df.shape[0])
    for idx in id_ls:
        df_ls.append(df[idx[0]:idx[1]])

    df_deals = Parallel(n_jobs=n_process)(delayed(key_fun)(x) for x in df_ls)
    df = pd.concat(df_deals)

    subproc_after = set([p.pid for p in current_process.children(recursive=True)])
    if terminate_subproc:
        terminate_joblib_parallel_subprocess(subproc_after-subproc_before)
    return df[col_new]


def apply_raw(df, fun, args=[], kwargs={}, n_process=mp.cpu_count()-4, terminate_subproc=True):
    """df的apply变换, 多进程, 针对行
    Parameter:
        fun: 变换函数fun(x, *args)
        args: 函数参数
        n_process: 进程数
    Return:
        df[col_new]
    """
    current_process = psutil.Process()
    subproc_before = set([p.pid for p in current_process.children(recursive=True)])

    if n_process < 0: n_process = 1
    col_new = '__!#%@$^__'
    # col_i = '__@$^!#%__'
    # df[col_i] = list(range(df.shape[0]))
    def key_fun(df):
        df[col_new] = df.apply(lambda x:fun(x, *args, **kwargs), axis=1)
        return df

    df_ls = []
    id_ls = split_id(n_process, df.shape[0])
    for idx in id_ls:
        df_ls.append(df[idx[0]:idx[1]])



    df_deals = Parallel(n_jobs=n_process)(delayed(key_fun)(x) for x in df_ls)
    df = pd.concat(df_deals)

    subproc_after = set([p.pid for p in current_process.children(recursive=True)])
    if terminate_subproc:
        terminate_joblib_parallel_subprocess(subproc_after-subproc_before)
    return df[col_new]

# def apply_raw_id(df, fun, args=[], kwargs={}, n_process=mp.cpu_count()-4):
    # """
    # 分割df，变换为分割df index，当df特别大的时候费时间
    # """


def apply_group(df, col, fun, t='series', args=[], kwargs={}, n_process=mp.cpu_count()-4, terminate_subproc=True):
    """df的apply变换, 多进程, 针对groupby后
    Parameter:
        col: groupby列
        fun: 变换函数fun(x, *args)
        args: 函数参数
        n_process: 进程数
    Return:
        df
    注意：
        df可以执行series操作，但选series更快
        只有当fun特别复杂时，才有加速效果，否则慢2~5倍
    """
    current_process = psutil.Process()
    subproc_before = set([p.pid for p in current_process.children(recursive=True)])


    if n_process < 0: n_process = 1
    df_gs = df.groupby(col)

    def get_idx(name, s, t='1col'):
        if isinstance(s, pd.Series) or isinstance(s, pd.DataFrame):
            idx_len = s.shape[0]
        else:
            idx_len = 1
        
        # 多级索引
        if isinstance(name, tuple):
            if t == 'df' and len(s.shape) == 1:
                idx_len = 1
            return pd.MultiIndex.from_arrays([[x]*idx_len for x in name])
        else:
            return [name]*idx_len
    
    def key_fun_series(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        s.name = name
        return s

    def key_fun_df(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        if len(s.shape) == 1:
            cols = list(s.index)
            v = [s.values]
        else:
            cols = list(s.columns)
            v = s.values
        s = pd.DataFrame(v, index=get_idx(name, s, t), columns=cols)
        return s

    def key_fun_1col(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        cols = [s.name]
        s = pd.DataFrame(s.values, index=get_idx(name, s), columns=cols)
        return s

    def key_fun_value(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        s = pd.DataFrame([s], index=get_idx(name, s))
        return s

    if t=='series':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_series)((name, group)) for name, group in df_gs)
        df_deal = pd.DataFrame(df_deals)
    if t=='1col':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_1col)((name, group)) for name, group in df_gs)
        df_deal = pd.concat(df_deals)
    if t=='df':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_df)((name, group)) for name, group in df_gs)
        df_deal = pd.concat(df_deals)
    if t=='value':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_value)((name, group)) for name, group in df_gs)
        df_deal = pd.concat(df_deals)


    subproc_after = set([p.pid for p in current_process.children(recursive=True)])
    if terminate_subproc:
        terminate_joblib_parallel_subprocess(subproc_after-subproc_before)
    return df_deal

def groupby_fast(df, col, fun, t='series', args=[], kwargs={}, n_process=mp.cpu_count()-4, terminate_subproc=True):
    """快速groupby
    实现：
        只是更改了apply_group的groupby环节
    注意：
        比apply_group慢很多，约5~10倍
    """
    current_process = psutil.Process()
    subproc_before = set([p.pid for p in current_process.children(recursive=True)])


    if n_process < 0: n_process = 1

    col_id = '__!#%@$^__'
    col_str = '__@$^!#%__'
    if len(col) == 1:
        group_set = set(df[col[0]])
        group_dict = dict(zip(group_set, range(len(group_set))))
        id_dict = dict(zip(range(len(group_set)), group_set))
        df[col_id] = df[col[0]].apply(lambda x:group_dict[x])
    else:
        # df[col_str] = df.apply(lambda x:'_'.join(x[col])) # 拼接成字符串，又要考虑int等转换，速度慢
        df[col_str] = df.apply(lambda x:tuple(x[col]), axis=1) # tuple hashable
        group_set = set(df[col_str])
        group_dict = dict(zip(group_set, range(len(group_set))))
        id_dict = dict(zip(range(len(group_set)), group_set))
        df[col_id] = df[col_str].apply(lambda x:group_dict[x])
    gs = df.groupby(col_id)
    df_gs = []
    for name, group in gs:
        del group[col_id]
        if col_str in group.columns: del group[col_str]
        df_gs.append((id_dict[name], group))

    # df_gs = df.groupby(col)

    def get_idx(name, s, t='1col'):
        if isinstance(s, pd.Series) or isinstance(s, pd.DataFrame):
            idx_len = s.shape[0]
        else:
            idx_len = 1
        
        # 多级索引
        if isinstance(name, tuple):
            if t == 'df' and len(s.shape) == 1:
                idx_len = 1
            return pd.MultiIndex.from_arrays([[x]*idx_len for x in name])
        else:
            return [name]*idx_len
    
    def key_fun_series(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        s.name = name
        return s

    def key_fun_df(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        if len(s.shape) == 1:
            cols = list(s.index)
            v = [s.values]
        else:
            cols = list(s.columns)
            v = s.values
        s = pd.DataFrame(v, index=get_idx(name, s, t), columns=cols)
        return s

    def key_fun_1col(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        cols = [s.name]
        s = pd.DataFrame(s.values, index=get_idx(name, s), columns=cols)
        return s

    def key_fun_value(z):
        name, group = z
        s = fun(group, *args, **kwargs)
        s = pd.DataFrame([s], index=get_idx(name, s))
        return s

    # def deal_fun(t, )

    if t=='series':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_series)((name, group)) for name, group in df_gs)
        df_deal = pd.DataFrame(df_deals)
    if t=='1col':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_1col)((name, group)) for name, group in df_gs)
        df_deal = pd.concat(df_deals)
    if t=='df':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_df)((name, group)) for name, group in df_gs)
        df_deal = pd.concat(df_deals)
    if t=='value':
        df_deals = Parallel(n_jobs=n_process)(delayed(key_fun_value)((name, group)) for name, group in df_gs)
        df_deal = pd.concat(df_deals)


    subproc_after = set([p.pid for p in current_process.children(recursive=True)])
    if terminate_subproc:
        terminate_joblib_parallel_subprocess(subproc_after-subproc_before)
    return df_deal


def joblib_parallel(fun, params, n_process=joblib.cpu_count()-4, terminate_subproc=True):
    """多进程，joblab库
    Arg:
        fun: 一个函数名
        params：参数列表[[,],...]
    Return:
        res_ls: 结果列表[,...]
    """
    current_process = psutil.Process()
    subproc_before = set([p.pid for p in current_process.children(recursive=True)])

    if n_process < 0: n_process = 1

    def key_fun(param):
        return fun(*param)

    # res_ls = Parallel(n_jobs=n_process, verbose=1)(delayed(key_fun)(param) for param in params)
    # res_ls = Parallel(n_jobs=n_process)(delayed(key_fun)(param) for param in params)
    with Parallel(n_jobs=n_process) as parallel:
        res_ls = parallel(delayed(key_fun)(param) for param in params)

    subproc_after = set([p.pid for p in current_process.children(recursive=True)])
    if terminate_subproc:
        terminate_joblib_parallel_subprocess(subproc_after-subproc_before)
    return res_ls

# @wrap_non_picklable_objects
def fun_1(param, fun):
    """不能放在multi_pool里面
    AttributeError: Can't pickle local object 'multi_pool.<locals>.fun_1'
    @wrap_non_picklable_objects装饰后可放在函数中定义
    """
    return fun(*param)

def multi_pool(fun, params, n_process=mp.cpu_count()-4):
    """多进程，multiprocessing.Pool库
    Arg:
        fun: 一个函数名
        params：参数列表[[,],...]
    Return:
        res_ls: 结果列表[,...]
    """
    if n_process < 0: n_process = 1
    if n_process > len(params): n_process = len(params)

    pool = Pool(processes=n_process)
    # res_ls = pool.map(fun, params)

    partial_fun = partial(fun_1, fun=fun)
    res_ls = pool.map(partial_fun, params)
    pool.close()
    pool.join()
    # with Pool(10) as pool:
        # res = pool.map(fun1, params)

    return res_ls

# def multi_process(funs, params, n_process=mp.cpu_count()-4):
    # """multiprocessing.Process库
    # Arg:
    #     funs: 一个函数名列表[fun1,...]
    #     params：参数列表[[,],...]
    # Return:
    #     # res_ls: 结果列表[,...]
    # """

    


