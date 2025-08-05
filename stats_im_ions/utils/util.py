# -*- coding: UTF-8 -*-

# @Date     : Sep 12, 2020
# @Author   : Northblue
"""工具函数"""

import pandas as pd

def print_sep_line():
    """输出分割行"""

    print('=================================')

def print_map_column(dct):
    """字典值两列输出

        例如：{'a':1,'b':2}
        输出：  a   1
                b   2
    """

    for k, value in dct.items():
        print("%s\t%s" %(k, value))

def is_num(string):
    """判断一个字符串是否为一个数

    Args:
        string: 字符串
    Returns:
        True or False
    """

    try:
        float(string)
        return True
    except ValueError:
        pass
    return False

def is_int(string):
    """判断一个字符串是否为整数

    Args:
        string: 字符串
    Returns:
        True or False
    """

    try:
        int(string)
        return True
    except ValueError:
        pass
    return False


def convert_ls2ser(ls):
    """ls数据结构变为ser

    主要是为了使用ser的value_counts功能
    """

    ser = pd.Series(ls)
    return ser


def print_dct_head(dct, n=5, is_print=True):
    """输出字典的前n个元素"""

    _dct = {}
    for i, (k, v) in enumerate(dct.items()):
        if i >= n:
            break
        _dct[k] = v
        if is_print:
            print(f'{k}: {v}')
    return _dct


def get_dct_v_idx(dct, idx=0):
    """获取字典的第idx个元素"""

    k = list(dct.keys())[idx]
    v = list(dct.values())[idx]
    return k, v
