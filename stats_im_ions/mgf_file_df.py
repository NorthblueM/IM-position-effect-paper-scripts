#-*- coding: UTF-8 -*-

# @Date     : 14th Sep, 2021
# @Author   : Northblue

"""mgf文件的读写、修改
    pandas存一张谱图信息
"""

import os
import pandas as pd
import numpy as np
from multiprocessing import Process, Queue
import struct
from pathlib import Path

class MgfFile(object):
    """"MGF文件"""

    def __init__(self, fpath):
        self.fpath = fpath
        self.df_head = pd.DataFrame()
        self.df_peak = pd.DataFrame()

        self.fname_suf_head = '_head.csv'
        self.fname_suf_peak = '_peak.csv'

        self.title2scanid = {}
        self.scanid2title = {}
    
    def read_mgf(self):
        """"读入一个mgf
            谱图头信息（title索引）、谱峰信息（title索引）两张表
        """

        self.read_mgf_head()
        self.set_scanid_pparse()
        self.read_mgf_peak()
        return self.df_head, self.df_peak
    

    def read_pf2(self):
        """"读入一个mgf
            谱图头信息（title索引）、谱峰信息（title索引）两张表
        """

        self.read_pf2_head()
        self.set_scanid_pparse()
        self.read_pf2_peak()
        return self.df_head, self.df_peak


    def save_head_csv(self):
        self.df_head.to_csv(self.get_fpath_head(), index=False)
    
    def save_peak_csv(self):
        self.df_peak.to_csv(self.get_fpath_peak(), index=False)
    
    def read_mgf_head(self):
        """读取谱图头信息
        """

        cols = ['title', 'charge', 'rt', 'pepmass', 'peak_num', 'inten_sum']
        ls = []

        with open(self.fpath, 'r') as fin:
            for line in fin:
                if 'BEGIN IONS' in line:
                    peak_num = 0
                    inten_sum = 0.0
                    continue
                if 'END IONS' in line:
                    one = [title, charge, rt, pepmass, peak_num, inten_sum]
                    ls.append(one)
                    continue
                if 'TITLE=' in line:
                    title = line.split('=')[1].strip()
                    continue
                if 'CHARGE=' in line:
                    charge = int(line.split('=')[1].strip()[:-1])
                    continue
                if 'RTINSECONDS=' in line:
                    rt = float(line.split('=')[1].strip())
                    continue
                if 'PEPMASS=' in line: # 此处质量为M/Z，pLink CSV中为中性质量+H
                    pepmass = float(line.split('=')[1].strip())
                    continue
                if not line.isspace() and line.split()[0][0].isdigit():
                    peak_num += 1
                    inten_sum += float(line.split()[1].strip())
                    continue
            self.df_head = pd.DataFrame(ls, columns=cols)

        return self.df_head

    def read_mgf_peak(self):
        """读取谱峰，以Title作为索引
        """

        cols = ['scanid', 'mz', 'inten']
        set_scanid = set()
        b_scanid = 0
        ls = []
        with open(self.fpath, 'r') as fin:
            for line in fin:
                if 'TITLE=' in line:
                    title = line.split('=')[1].strip()
                    scanid = self.title2scanid[title]
                    if scanid in set_scanid:
                        b_scanid = 0
                    else:
                        b_scanid = 1
                        set_scanid.add(scanid)
                    continue
                if b_scanid and not line.isspace() and line.split()[0][0].isdigit():
                    segs = line.strip().split()
                    mz = float(segs[0].strip())
                    inten = float(segs[1].strip())
                    one = [scanid, mz, inten]
                    ls.append(one)
                    continue
            self.df_peak = pd.DataFrame(ls, columns=cols)

        return self.df_peak


    def read_pf2_head(self):
        """读取谱图头信息
        """

        cols = ['title', 'charge', 'rt', 'pepmass', 'peak_num', 'inten_sum']
        ls = []

        with open(self.fpath, 'rb') as pf2_file:
            total_num = struct.unpack('i', pf2_file.read(4))[0]  # 读取完整谱图数量total_num(int)
            title_char_len = struct.unpack('i', pf2_file.read(4))[0]  # 读取title字符串长度(int) 
            title = pf2_file.read(title_char_len).decode('utf-8').strip('\x00')  # 读取并处理title, ChenRF
            # title = struct.unpack("%ds"%title_char_len, pf2_file.read(title_char_len)) # .decode('utf-8')

            # 遍历所有谱图
            for i in range(total_num):
                scan_no = struct.unpack('i', pf2_file.read(4))[0]  # 读取scan_no
                peak_num = struct.unpack('i', pf2_file.read(4))[0]  # 读取谱峰数量

                # [mz,c,mz,c]
                ls_mz_inten = struct.unpack(str(peak_num*2)+"d", pf2_file.read(peak_num*2*8))
                # pf2_file.read(2*8*peak_num)  # 跳过谱峰内容（谱峰需要从ms2中读取并写入）

                # ls_mz = ls_mz_inten[0::2]
                # ls_inten = ls_mz_inten[1::2]
                # assert len(ls_mz) == len(ls_inten) == peak_num, 'peak_num error'
                # inten_sum = sum(ls_inten)

                precursor_num = struct.unpack('i', pf2_file.read(4))[0]  # 读取母离子数量
                for i in range(precursor_num):
                    pepmass = struct.unpack('d', pf2_file.read(8))[0]
                    charge = struct.unpack('i', pf2_file.read(4))[0]
                    # ls.append([f'{title}.{scan_no}.{scan_no}.{charge}.{i}.dta', charge, 0.0, pepmass, peak_num, inten_sum])
                    ls.append([f'{title}.{scan_no}.{scan_no}.{charge}.{i}.dta', charge, 0.0, pepmass, peak_num, 0.0])

        self.df_head = pd.DataFrame(ls, columns=cols)

        return self.df_head

    def read_pf2_peak(self):
        """读取谱峰，以Title作为索引
        """

        cols = ['scanid', 'mz', 'inten']
        ls = []
        with open(self.fpath, 'rb') as pf2_file:
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

                ls_tmp = [scan_no]
                for i in range(peak_num):
                    _ls_tmp = ls_tmp.copy()
                    _ls_tmp.append(ls_mz[i])
                    _ls_tmp.append(ls_inten[i])
                    ls.append(_ls_tmp)

                precursor_num = struct.unpack('i', pf2_file.read(4))[0]  # 读取母离子数量
                for i in range(precursor_num):
                    pepmass = struct.unpack('d', pf2_file.read(8))[0]
                    charge = struct.unpack('i', pf2_file.read(4))[0]

        self.df_peak = pd.DataFrame(ls, columns=cols)

        return self.df_peak


    def set_title2scanid(self):
        """根据df_head, 设置对象中的title2scanid成员
        """

        self.title2scanid = dict(zip(list(self.df_head['title']), list(self.df_head['scanid'])))

        return self.title2scanid
    
    
    def set_scanid2title(self):
        """根据df_head, 设置对象中的scanid2title成员
        注意: scanid可能重复
        """

        self.scanid2title = dict(zip(list(self.df_head['scanid']), list(self.df_head['title'])))

        return self.scanid2title


    def set_scanid_pparse(self):
        """按照pParse设置Title的规则，头文件加scan id
        """

        self.df_head['scanid'] = self.df_head['title'].apply(lambda x:int(x.split('.')[-4]))

        self.set_title2scanid()
        return self.title2scanid

    def set_scanid_pxtract(self):
        """按照pXtract设置Title的规则，头文件加scan id
        """

        self.df_head['scanid'] = self.df_head['title'].apply(lambda x:int(x.split('.')[-3]))

        self.set_title2scanid()
        return self.title2scanid


    def get_fpath_head(self):
        """"获取谱图头信息文件路径"""
        # return self.fpath[:-4]+self.fname_suf_head
        fparent = Path(self.fpath).parent
        fstem = Path(self.fpath).stem
        return str(fparent/fstem)+self.fname_suf_head

    def get_fpath_peak(self):
        """"获取谱峰信息文件路径"""
        # return self.fpath[:-4]+self.fname_suf_peak
        fparent = Path(self.fpath).parent
        fstem = Path(self.fpath).stem
        return str(fparent/fstem)+self.fname_suf_peak
    

    def set_inten_max(self):
        """头文件设置谱峰最大强度"""
        self.df_head = self.df_head.set_index(['scanid'])
        self.df_head[['inten_max']] = self.df_peak.groupby(['scanid']).agg({'inten':['max']})

        self.df_head = self.df_head.reset_index()

        # peak_num = 0, 强度为0
        self.df_head.loc[self.df_head['peak_num'] == 0, 'inten_max'] = 0.0

        return self.df_head
    
    def set_inten_sum(self):
        """头文件设置谱峰强度和"""
        self.df_head = self.df_head.set_index(['scanid'])
        self.df_head[['inten_sum']] = self.df_peak.groupby(['scanid']).agg({'inten':['sum']})

        self.df_head = self.df_head.reset_index()

        # peak_num = 0, 强度为0
        self.df_head.loc[self.df_head['peak_num'] == 0, 'inten_sum'] = 0.0

        return self.df_head

    def set_inten_relat(self):
        """谱峰文件设置相对谱峰强度"""
        self.df_peak['inten_relat'] = self.df_peak.groupby(['scanid'])[['inten']].transform(lambda x:x/x.max())

        return self.df_peak
    
    def save_inten_relat_flow(self):
        
        self.read_mgf()

        self.set_inten_max()
        self.set_inten_relat()

        self.save_head_csv()
        self.save_peak_csv()
    
    def save_inten_relat_flow_pf2(self):
        
        self.read_pf2()

        self.set_inten_max()
        self.set_inten_sum()
        self.set_inten_relat()

        self.save_head_csv()
        self.save_peak_csv()


def bat_folder_mgf_save_df_1(processor_datas, start, end, file_type='mgf'):
    """单线程"""
    for i in range(start, end):
        data = processor_datas[i]

        fpath = data
        print('===', fpath)
        mgf = MgfFile(fpath)
        if file_type in ['mgf']:
            mgf.save_inten_relat_flow()
        elif file_type in ['pf1', 'pf2']:
            mgf.save_inten_relat_flow_pf2()


def bat_folder_mgf_save_df_multiprocess(dpath, file_type='mgf', n_process=8):
    """处理一个文件夹的mgf转换为df"""
    
    fpaths = []
    for fname in list(os.walk(dpath))[0][2]:
        if fname.endswith('.'+file_type):
            fpath = os.path.join(dpath, fname)
            fpaths.append(fpath)

    processor_datas = fpaths
    processor_num = n_process
    interval = int(len(processor_datas)/processor_num)
    rem = len(processor_datas)%processor_num #余数
    process_ls = []
    # 数据划分样例：012|345|67|89
    for i in range(processor_num):
        start = i*interval
        end = (i+1)*interval
        if i < rem:
            start += i
            end += (i+1)
        else:
            start += rem
            end += rem
        print('===data split:\t', start, end)

        p = Process(target=bat_folder_mgf_save_df_1, args=(processor_datas, start, end, file_type))
        process_ls.append(p)
        p.start()

    for p in process_ls:
        p.join()

def _main():
    
    # # === mgf读入，写出df测试
    # fpath = r'.mgf'
    # mgf = MgfFile(fpath)
    # mgf.read_mgf()
    # # ===

    # # ===一个文件
    # fpath = r'.pf1'
    
    
    # file_type = fpath.split('.')[-1]
    # mgf = MgfFile(fpath)
    # if file_type in ['mgf']:
    #    mgf.save_inten_relat_flow()
    #elif file_type in ['pf1', 'pf2']:
    #    mgf.save_inten_relat_flow_pf2()
    # # ===

    # # # ===处理一个文件夹的mgf转换为df
    dpath = r'D:\users\pzmao\pFindCooperation\202503_imide_position\20250405_Multi-HCD\data\raw_20250418_stepHCD_40trigger'
    
    
    # bat_folder_mgf_save_df_multiprocess(dpath, file_type='mgf', n_process=8)
    # # bat_folder_mgf_save_df_multiprocess(dpath, file_type='pf1', n_process=8)
    bat_folder_mgf_save_df_multiprocess(dpath, file_type='pf2', n_process=4)
    # # # ===


if __name__ == '__main__':
    _main()





