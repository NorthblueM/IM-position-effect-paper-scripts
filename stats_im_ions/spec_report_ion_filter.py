#-*- coding: UTF-8 -*-

# @Date     : 18th Jun, 2025
# @Author   : Northblue

"""pf2, 根据report ion过滤谱图
"""

import os
import sys
from pathlib import Path
sys.path.append(str(Path.cwd().parent))

from tools import df_apply_mutli as dfm

import pandas as pd
import bisect
import struct


def write_pf2_from_mp(df_head, scan2peak, rawname, fpath_out, fpath_out_idx=''):

    if not os.path.exists(os.path.dirname(fpath_out)):
        os.makedirs(os.path.dirname(fpath_out))

    # https://blog.51cto.com/u_16099279/6358555
    def x_mz_c_tolist(x):
        """pepmass and charge to list
        Args:
            x: col; pepmass, charge;
        Returns:
            ls: [[pepmass,charge],]
        """
        ls = list(zip(*[x['pepmass'].tolist(), x['charge'].tolist()]))
        # ls = list(map(list, ls))
        return ls

    # {scan1:[[pepmass, charge],], scan2:}
    col_scan = 'scannum'
    if col_scan not in df_head: col_scan = 'scanid'
    scan2mz_c = df_head.groupby([col_scan]).apply(x_mz_c_tolist, include_groups=False).to_dict()

    scan2offset = {}

    offset = 0

    # write pf2
    with open(fpath_out, 'wb') as f:
        f.write(struct.pack("2i", len(scan2mz_c), len(rawname)))
        f.write(struct.pack("%ds"%len(rawname), bytes(rawname, 'utf-8')))
        offset += 2*4 + len(rawname)*1

        for scannum, precursors in scan2mz_c.items():
            scan2offset[scannum] = offset

            num_peak = len(scan2peak[scannum][0])

            # to [mz,c,mz,c]
            ls_peak = [y for x in list(zip(*scan2peak[scannum])) for y in x]

            f.write(struct.pack("i", scannum)) # scan_no
            f.write(struct.pack("i", num_peak)) # nPeak
            f.write(struct.pack("%dd"%(num_peak*2), *ls_peak)) # peaks
            f.write(struct.pack("i", len(precursors)))  # nMix
            offset += 3*4 + num_peak*2*8

            for precursor in precursors:
                f.write(struct.pack("d", precursor[0]))  # pepmass
                f.write(struct.pack("i", precursor[1]))  # charge
            offset += len(precursors) * (4+8)


    if fpath_out_idx:
        # write pf2idx
        with open(fpath_out_idx, 'wb') as f:
            for scannum in scan2mz_c:
                f.write(struct.pack("i", scannum))  # scan_no
                f.write(struct.pack("q", scan2offset[scannum]))  # offset, long long


def write_mgf_from_ms2_mp(fpath_out, df_head, ms2_scan2peak):

    with open(fpath_out, 'w', encoding='utf-8') as fo:
        for i, row in enumerate(df_head.itertuples()):

            lines = []
            lines.append('BEGIN IONS')
            lines.append('TITLE=%s'%(getattr(row, 'title')))
            lines.append('CHARGE=%d+'%(getattr(row, 'charge')))
            lines.append('RTINSECONDS=%lf'%(getattr(row, 'rt')))
            lines.append('PEPMASS=%lf'%(getattr(row, 'pepmass')))

            ls_mz, ls_inten = ms2_scan2peak[getattr(row, 'scannum')]
            for i, mz in enumerate(ls_mz):
                lines.append('%.6f %.1f'%(mz, ls_inten[i]))
            
            lines.append('END IONS\n')

            fo.write('\n'.join(lines))


def cal_match_error(theo, exp):
    return abs(exp-theo)/theo*1e6


def match_ion_isin(mz, mz_ls, match_tol=20, match_wind=0.5):
    """match one peak
    bisection method
    mz_ls: sort from smallest to biggest
    """

    mz_low, mz_high = mz-match_wind , mz+match_wind # 0.5Da must more than 20ppm
    idx_low, idx_high = bisect.bisect_left(mz_ls, mz_low), bisect.bisect_right(mz_ls, mz_high)
    n = len(mz_ls)

    for i in range(idx_low, idx_high+1):
        if i >= 0 and i < n:
            match_error = cal_match_error(mz_ls[i], mz)
            if match_error <= match_tol:
                return 1
    return 0


def is_spec_report_ion(ls_mz, ls_report_ions, match_tol=20, match_wind=0.5):
    """
    判断谱图是否满足report ion的要求
    默认: match_tol=20ppm, match_wind=0.5Da
    """

    for ls_ion in ls_report_ions:
        flag = 1
        for mz in ls_ion:
            flag &= match_ion_isin(mz, ls_mz, match_tol=match_tol, match_wind=match_wind)
        if flag:
            return 1
    return 0


def load_pf2_and_filter(fpath_pf2, ls_report_ions, fpath_out):

    print(f'------------{fpath_pf2}')

    scan2peak = {}
    cols_head = ['title', 'charge', 'rt', 'pepmass', 'peak_num', 'inten_sum', 'inten_max']
    ls_head = []

    # title = '_'.join(Path(fpath_pf2).stem.split('_')[:-1]) # 去掉_HCDFT.pf2
    # total_num = 0

    with open(fpath_pf2, 'rb') as pf2_file:
        total_num = struct.unpack('i', pf2_file.read(4))[0]  # 读取完整谱图数量total_num(int)
        title_char_len = struct.unpack('i', pf2_file.read(4))[0]  # 读取title字符串长度(int) 
        title = pf2_file.read(title_char_len).decode('utf-8').strip('\x00')  # 读取并处理title, ChenRF
        # title = struct.unpack("%ds"%title_char_len, pf2_file.read(title_char_len)) # .decode('utf-8')

        # 遍历所有谱图
        for idx in range(total_num):
            scan_no = struct.unpack('i', pf2_file.read(4))[0]  # 读取scan_no
            peak_num = struct.unpack('i', pf2_file.read(4))[0]  # 读取谱峰数量

            # [mz,c,mz,c]
            ls_mz_inten = struct.unpack(str(peak_num*2)+"d", pf2_file.read(peak_num*2*8))
            # pf2_file.read(2*8*peak_num)  # 跳过谱峰内容（谱峰需要从ms2中读取并写入）

            ls_mz = ls_mz_inten[0::2]
            ls_inten = ls_mz_inten[1::2]
            assert len(ls_mz) == len(ls_inten) == peak_num, 'peak_num error'
            inten_sum = sum(ls_inten)

            if peak_num == 0: inten_max = 0.0
            else: inten_max = max(ls_inten)

            # 过滤谱图
            flag = 0
            flag = is_spec_report_ion(ls_mz, ls_report_ions)

            if flag:
                scan2peak[scan_no] = [ls_mz, ls_inten]

            precursor_num = struct.unpack('i', pf2_file.read(4))[0]  # 读取母离子数量
            for i in range(precursor_num):
                pepmass = struct.unpack('d', pf2_file.read(8))[0]
                charge = struct.unpack('i', pf2_file.read(4))[0]

                if flag:
                    ls_head.append([f'{title}.{scan_no}.{scan_no}.{charge}.{i}.dta', charge, 0.0, pepmass, peak_num, inten_sum, inten_max])

    print(f'------{title}: #scan {total_num} -> {len(scan2peak)}')

    df_head = pd.DataFrame(ls_head, columns=cols_head)
    df_head['scanid'] = df_head['title'].apply(lambda x:int(x.split('.')[-4]))

    rawname = title
    write_pf2_from_mp(df_head, scan2peak, rawname, fpath_out, fpath_out_idx='')


def load_mgf_and_filter(fpath_mgf, ls_report_ions, fpath_out):

    print(f'------------{fpath_mgf}')

    scan2peak = {}
    cols = ['title', 'pepmass', 'charge', 'rt', 'scannum']
    ls_ls = []

    spec_total = 0
    spec_filtered = 0
    with open(fpath_mgf, 'r', encoding='utf-8') as f:
        title = ''
        scannum = 0
        pepmass = 0.0
        charge = 0
        rt = 0.0
        ls_mz = []
        ls_inten = []
        for line in f:
            if line.isspace():
                continue
            if line[0].isalpha():
                if line.startswith('TITLE='):
                    title = line.strip().split('=')[1]
                    scannum = int(title.split('.')[-4]) # pParse
                    ls_mz = []
                    ls_inten = []
                    continue
                if line.startswith('CHARGE='):
                    charge = int(line.strip().split('=')[1][:-1])
                    continue
                if line.startswith('RTINSECONDS='):
                    rt = float(line.strip().split('=')[1])
                    continue
                if line.startswith('PEPMASS='):
                    pepmass = float(line.strip().split('=')[1])
                    continue
                if line.startswith('END IONS'):
                    assert charge > 0, 'charge must be > 0'
                    assert pepmass > 0, 'pepmass must be > 0'
                    assert title, 'title must not be empty'
                    spec_total += 1
                    # 过滤谱图
                    # print(f'  {title}: scannum={scannum}, pepmass={pepmass}, charge={charge}, rt={rt}, peak_num={len(ls_mz)}')
                    if is_spec_report_ion(ls_mz, ls_report_ions):
                        spec_filtered += 1
                        ls_ls.append([title, pepmass, charge, rt, scannum])
                        if scannum not in scan2peak: # 同scan num的谱图只保留一个
                            scan2peak[scannum] = [ls_mz, ls_inten]
            if title and line.split()[0][0].isdigit():
                segs = line.strip().split()
                assert len(segs) >= 2, 'line: mz intensity'
                mz = float(segs[0].strip())
                inten = float(segs[1].strip())
                ls_mz.append(mz)
                ls_inten.append(inten)
                continue

    print(f'------{Path(fpath_mgf).stem}: #spec {spec_total} -> {spec_filtered}')
    
    df_head = pd.DataFrame(ls_ls, columns=cols)

    write_mgf_from_ms2_mp(fpath_out, df_head, scan2peak)


def bat_folder(dpath_spec, ls_report_ions, dpath_out, n_process=8, spec_format='pf2'):
    """一个文件夹的pf2文件批量处理"""

    ls_param = []
    if spec_format == 'pf2':
        for fpath_spec in Path(dpath_spec).glob('*.pf2'):
            fpath_out = str(Path(dpath_out)/fpath_spec.name)
            ls_param.append([str(fpath_spec), ls_report_ions, fpath_out])
        print(f'============批量处理{len(ls_param)}个文件')
        dfm.multi_pool(load_pf2_and_filter, ls_param, n_process=n_process)
    elif spec_format == 'mgf':
        for fpath_spec in Path(dpath_spec).glob('*.mgf'):
            fpath_out = str(Path(dpath_out)/fpath_spec.name)
            ls_param.append([str(fpath_spec), ls_report_ions, fpath_out])
        print(f'============批量处理{len(ls_param)}个文件')
        dfm.multi_pool(load_mgf_and_filter, ls_param, n_process=n_process)
    else:
        raise ValueError(f'Unsupported spec_format: {spec_format}')


def _main():
    
    # [[]], 内层list是且的关系, 外层list是或的关系
    # # Propionyl[K](Delta_H(4)C(3)O(1)[K]), C(8)H(16)N(2)O(1)H+, C(8)H(13)N(1)O(1)H+
    # ls_report_ions = [[157.133531], [140.106984]] #+质子
    # ls_report_ions = [[140.106984]]
    # ls_report_ions = [[157.133531]]

    # Xlink_BS2G[114][K](Gluratylation[K]) CycIm ion_C10H16O3N_m/z= 198.11302, LinIm ion_ C10H19O3N2_m/z= 215.13957 # 一价+氢原子
    # ls_report_ions = [[215.139008], [198.112461]] #+质子
    ls_report_ions = [[198.112461]]
    # ls_report_ions = [[215.139008]]
    # ls_report_ions = [[215.139008, 198.112461]]


    print(f'============ls_report_ions: {ls_report_ions}')

    # 单个文件
    # path_spec = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20240409_NCEshuffle\ZZY_BS2G_d0_Trypsin_multi_nce_disordering_0408_HCDFT.pf2'

    # 批量文件夹
    path_spec = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw'
    # path_spec = r'\\10.68.0.17\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20250418_stepHCD_40trigger'

    spec_format = 'mgf'  # pf2 or mgf
    print(f'============spec format: {spec_format}')

    n_process = 8
    dname_out = 'spec_report_ion_filtered'
    
    if os.path.isdir(path_spec): # 如果是文件夹, 批量
        dpath_out = str(Path(path_spec)/dname_out)
        bat_folder(path_spec, ls_report_ions, dpath_out, n_process=n_process, spec_format=spec_format)
    else:
        fpath_out = str(Path(path_spec).parent/dname_out/Path(path_spec).name)
        if spec_format == 'pf2':
            load_pf2_and_filter(path_spec, ls_report_ions, fpath_out)
        elif spec_format == 'mgf':
            load_mgf_and_filter(path_spec, ls_report_ions, fpath_out)
        else:
            raise ValueError(f'Unsupported spec_format: {spec_format}')
    

if __name__ == "__main__":
    _main()
