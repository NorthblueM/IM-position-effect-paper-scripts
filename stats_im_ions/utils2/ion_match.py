# -*- coding: UTF-8 -*-
# @Date     : 15th Jul, 2022
# @Author   : Northblue

import bisect
from ion_match_param import Param
import spec_theo_file

# 中间Param()改变, MS2_MATCH_TOL不会变，因为只import了一次
MS2_MATCH_WIND = Param().ms2_match_wind
MS2_MATCH_TOL = Param().ms2_match_tolerance
FRAG_CHARGE_LOWER_PREC = Param().matching_fragment_charge_lower_precursor

def cal_match_error(theo, exp):
    return abs(exp-theo)/theo*1e6


def match_ion(mz, mz_ls):
    """match one peak
    bisection method
    mz_ls: sort from smallest to biggest
    """

    ls = [] # [[idx, match_error],]

    mz_low, mz_high = mz-MS2_MATCH_WIND , mz+MS2_MATCH_WIND # 0.5Da must more than 20ppm
    idx_low, idx_high = bisect.bisect_left(mz_ls, mz_low), bisect.bisect_right(mz_ls, mz_high)
    n = len(mz_ls)

    for i in range(idx_low, idx_high+1):
        if i >= 0 and i < n:
            match_error = cal_match_error(mz_ls[i], mz)
            if match_error <= MS2_MATCH_TOL:
                ls.append([i, match_error])
    return ls


def match_spec(df_peak, df_ion, keep_unmatched=True, iontype_str=False, prec_charge=100, match_type=''):
    """
    keep_unmatched: whether to keep unmatched peaks
    """

    if FRAG_CHARGE_LOWER_PREC:
        df_ion = df_ion[df_ion['charge'] <= prec_charge] # 允许precursor peak
        # df_ion = df_ion[df_ion['charge'] < prec_charge]
    df_ion = df_ion.sort_values('mz').reset_index(drop=True)
    
    ls = []
    mz_theo_ls = list(df_ion['mz'])
    for i in range(df_peak.shape[0]): # each experimental peak
        mz_exp = df_peak.iloc[i]['mz']

        # [['ion_theo_idx', 'match_error'], ]
        match_ls = match_ion(mz_exp, mz_theo_ls)

        # [['tp_ion_theo', 'mz_theo', 'match_error'], ]
        match_info = []

        for match in match_ls:
            # 'tp_ion_theo': see spec_theo_file.py get_iontheo_str()
            tp_ion_theo = spec_theo_file.get_iontheo_str(df_ion.iloc[match[0]], iontype_str=iontype_str, match_type=match_type)

            mz_theo = df_ion.iloc[match[0]]['mz']
            match_info.append([tp_ion_theo, mz_theo, match[1]])
        ls.append(match_info)

    df_peak['match_info'] = ls

    if not keep_unmatched:
        df_peak = df_peak[[True if x else False for x in df_peak['match_info']]]

    return df_peak