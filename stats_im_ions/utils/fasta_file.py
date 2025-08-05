# -*- coding: UTF-8 -*-

# @Date    : Sep 12, 2020
# @Author  : Nrothblue
"""fasta文件的一些操作"""

import os
import sys
sys.path.append(os.path.dirname(os.getcwd()))

from utils import util
import random
import pandas as pd



def read_fasta_file(path):
    """读取fasta文件

    Args:
        path: fasta文件路径
    Returns:
        字典{蛋白质名称：蛋白质序列}
    """

    map_protein_seq = {}

    with open(path, 'r') as fin:
        i_total = 0
        prtn_name = ''
        prtn_seq = ''
        prtn_seq_ls = []
        for line in fin:
            if line[0] == '>':
                i_total += 1
                if i_total != 1:
                    prtn_seq = ''.join(prtn_seq_ls)
                    map_protein_seq[prtn_name] = prtn_seq
                    prtn_name = ''
                    prtn_seq = ''
                    prtn_seq_ls = []
                prtn_name = line.strip()[1:]
            else:
                prtn_seq_ls.append(line.strip())

        prtn_seq = ''.join(prtn_seq_ls)
        map_protein_seq[prtn_name] = prtn_seq
    return map_protein_seq




def covert_fasta_shortac(prtac2seq):
    """蛋白质名称，用短名，以空格分割"""

    prtac2seq_new = {}
    for prtac, seq in prtac2seq.items():
        prtac2seq_new[prtac.split()[0]] = seq
    return prtac2seq_new

def covert_fasta_i2l(prtac2seq):
    """蛋白质序列i2l"""

    prtac2seq_i2l = {}
    for prtac, seq in prtac2seq.items():
        prtac2seq_i2l[prtac] = seq.replace('I', 'L')
    return prtac2seq_i2l

def read_fasta_i2l(fpath, is_i2l=False, is_shortac=False):

    prtac2seq = read_fasta_file(fpath)
    if is_i2l:
        prtac2seq = covert_fasta_i2l(prtac2seq)
        if is_shortac:
            prtac2seq = covert_fasta_shortac(prtac2seq)
    elif is_shortac:
        prtac2seq = covert_fasta_shortac(prtac2seq)
    return prtac2seq


def stats_aa_freq(fpath, i2l=False):
    """统计氨基酸在fasta中的频率
    [表示蛋白质N端
    """

    prtac2seq = read_fasta_i2l(fpath, is_i2l=i2l, is_shortac=False)
    aa2freq = {}
    for ac, seq in prtac2seq.items():
        if '[' not in aa2freq: aa2freq.setdefault('[', 0)
        aa2freq['['] += 1
        for aa in seq:
            if aa not in  aa2freq: aa2freq.setdefault(aa, 0)
            aa2freq[aa] += 1
    return aa2freq


def write_fasta_file(path, map_prtn_seq):
    """写出fasta文件

    Args:
        path: fasta文件路径
        map_prtn_seq: 字典{蛋白质名称：蛋白质序列}
    """

    prtn_line_len = 60 # 每行蛋白质序列最长

    with open(path, 'w') as f:
        for k, v in map_prtn_seq.items():
            f.write('>' + k + '\n')
            for i in range(int(len(v)/prtn_line_len)):
                f.write(v[i*prtn_line_len:(i+1)*prtn_line_len] + '\n')
            if len(v)%prtn_line_len != 0:
                i = int(len(v)/prtn_line_len)
                f.write(v[i*prtn_line_len:] + '\n')

def merge_fasta_file(map_prtn1, map_prtn2):
    """合并两个fasta文件

    注意：两个fasta没有重名ac
    """

    map_prtn = {}

    for ac, seq in map_prtn1.items():
        map_prtn[ac] = seq

    for ac, seq in map_prtn2.items():
        map_prtn[ac] = seq

    return map_prtn


def extract_top_fasta(map_prtn, num):
    """提取前top个蛋白"""

    map_prtn_top = {}
    i = 0
    for ac, seq in map_prtn.items():
        if i < num:
            map_prtn_top[ac] = seq
            i += 1

    return map_prtn_top

def extract_shuffle_fasta(map_ac2seq, num, seed):
    """提取固定数目的蛋白，先打乱"""

    map_ac2seq_num = {}
    acs = list(map_ac2seq.keys())
    random.seed(seed)
    random.shuffle(acs)
    i = 0
    while(i<num):
        map_ac2seq_num[acs[i]] = map_ac2seq[acs[i]]
        i += 1
    return map_ac2seq_num

def enzyme_specific(map_ac2seq, max_miss_site, min_pep_len, cleavs, term):
    """进行特异酶切
    max_miss_site=3, min_pep_len=5, cleavs=['K','R'], term='N'
    """

    map_ac2pep_seqs = {}

    # N端酶切
    for ac, seq in map_ac2seq.items():
        pep_seqs = []
        cleav_sites = [0] # 酶切位点坐标，从0开始
        for i, aa in enumerate(seq):
            if aa in cleavs:
                cleav_sites.append(i)

        if cleav_sites[-1] != len(seq)-1:
            cleav_sites.append(len(seq)-1)

        for i, site_start in enumerate(cleav_sites):
            for j in range(max_miss_site+1):
                k = i+1+j
                if k < len(cleav_sites):
                    site_end = cleav_sites[k]
                else:
                    continue

                if i == 0:
                    pep_seq = seq[0: site_end+1]
                else:
                    pep_seq = seq[site_start+1: site_end+1]

                if len(pep_seq) >= min_pep_len:
                    pep_seqs.append(pep_seq)
        map_ac2pep_seqs[ac] = pep_seqs

    return map_ac2pep_seqs


def enzyme_specific2(ac2seq, digest, miss_site_max=3, pep_len_min=6, pep_len_max=60):
    """蛋白特异酶切

    Args:
        ac2seq: key,蛋白名称; seq,蛋白质序列
        digest: 酶切位点，[[['K','R']['P']],],前面为C端切[[切割位点],[忽略位点]],后面为N端切
        miss_site_max: 最大遗漏酶切位点
        pep_len_min: 肽段最小长度
        pep_len_max: 肽段最大长度
    """

    ac2pep_seqs = {}
    for ac, seq in ac2seq.items():
        
        digest_sites = [] # 酶切位点, 当前氨基酸“前”切
        for i, aa in enumerate(seq):
            if i == 0: # 第一个氨基酸
                digest_sites.append(i)
                continue
            if i == len(seq)-1: # 最后一个氨基酸
                digest_sites.append(i+1)
                continue
            if aa in digest[0][0]: # C端切
                if seq[i+1] not in digest[0][1]: # 忽略位点
                    digest_sites.append(i+1)
            if aa in digest[1][0]: # N端切
                if seq[i-1] not in digest[1][1]: # 忽略位点
                    # TODO: 判断上一个酶切氨基酸不是当前氨基酸，防止相邻氨基酸C端和N端酶切同时合法，实则是同一个酶切位点
                    digest_sites.append(i)
        print(len(digest_sites))
        
        pep_seqs = []
        for i in range(len(digest_sites)-1):
            for j in range(1, miss_site_max+1+1):
                if i+j >= len(digest_sites):
                    continue
                pep_seq = seq[digest_sites[i]:digest_sites[i+j]]
                if len(pep_seq) >= pep_len_min and len(pep_seq) <= pep_len_max:
                    pep_seqs.append(pep_seq)
        # pep_seqs = list(set(pep_seqs))
        ac2pep_seqs[ac] = pep_seqs

    return ac2pep_seqs



def get_reverse_rename(ac2seq):
    """获得REVERSE序列并重命名，去掉REVERSE_"""

    ac2seq_reverse = {}
    for ac, seq in ac2seq.items():
        if 'REVERSE_' in ac:
            ac2seq_reverse[ac[8:]] = seq
    return ac2seq_reverse

def rename_ac(ac2seq, k_str):
    """重命名蛋白的ac字段，以'k_str_'+从1开始的编号"""

    ac2seq_rename = {}
    for i, (ac, seq) in enumerate(ac2seq.items()):
        ac_rename = '_'.join([k_str, str(i+1)])
        ac2seq_rename[ac_rename] = seq
    return ac2seq_rename


def batch_merge_fasta(org_paths, trap_paths, outpath):
    """原始库与陷阱库合并

    Args:
        org_paths: 原始库完整路径列表
        trap_paths: 陷阱库完整路径列表
        output: 合并fasta输出文件夹路径

    """

    def get_file_name(inpath):
        return inpath.split('\\')[-1]

    def get_merge_name(name1, name2):
        [name1_str1, name_ends1] = os.path.splitext(name1)
        [name1_str2, name_ends2] = os.path.splitext(name2)
        return '_'.join([name1_str1, name1_str2]) + name_ends1

    print('===@()', sys._getframe().f_code.co_name)

    for org_path in org_paths:
        for trap_path in trap_paths:

            map_prtn_org = read_fasta_file(org_path)
            map_prtn_trap = read_fasta_file(trap_path)
            map_prtn_merge = merge_fasta_file(map_prtn_org, map_prtn_trap)

            # print(map_prtn_org)
            # print(map_prtn_merge)

            merge_name = get_merge_name(get_file_name(org_path), get_file_name(trap_path))
            merge_path = os.path.join(outpath, merge_name)

            write_fasta_file(merge_path, map_prtn_merge)
            print('======org_name, trap_num, merge_name:\t', get_file_name(org_path), get_file_name(trap_path), get_file_name(merge_path))
            print('======org_num, trap_num, merge_num:\t', len(map_prtn_org), len(map_prtn_trap), len(map_prtn_merge))


def extract_ac(ac2seq, ac_ls, cmp_shortac=False):
    """根据蛋白名称提取蛋白
        cmp_shortac: 通过短文件名比较
    """
    if cmp_shortac:
        ac_ls = [x.split()[0] for x in ac_ls]

    ac_set = set(ac_ls)
    ac2seq_new = {}
    for ac, seq in ac2seq.items():
        ac_str = ac
        if cmp_shortac:
            ac_str = ac_str.split()[0]
        if ac_str in ac_set:
            ac2seq_new[ac] = seq
    print('===original protein: %d\textract protein: %d'%(len(ac2seq),len(ac2seq_new)))
    return ac2seq_new


def extract_ac_pfind_spectra(fpath):
    """从pFind的spectra鉴定文件提取蛋白质名称
        去掉反库
    """

    df = pd.read_csv(fpath, sep='\t')
    ls_ac_str = list(df['Proteins'])
    ls_ac = []
    for ac_str in ls_ac_str:
        ls_ac.extend(ac_str.strip('/').split('/'))
    
    ls_ac_target = []
    for ac in ls_ac:
        if not ac.startswith('REV_'):
            ls_ac_target.append(ac)
    
    print('===total protein: %d\t target: %d'%(len(set(ls_ac)), len(set(ls_ac_target))))

    return set(ls_ac_target)

def build_fasta_pfind_spectra(fpath_fasta, fpath_pfind, fpath_out):
    """从pFind的spectra鉴定文件提取蛋白质名称
        建新库
    """

    ac2seq = read_fasta_file(fpath_fasta)
    set_ac = extract_ac_pfind_spectra(fpath_pfind)
    ac2seq = extract_ac(ac2seq, list(set_ac), cmp_shortac=True)

    write_fasta_file(fpath_out, ac2seq)



def _main():
    # # # # # # # 测试fasta文件的读入
    # fasta_path = r'\\10.68.0.17\dataset\mpz_xl\Mechtler_NC_2020_PXD014337\fasta\db_org'
    # fasta_name = r'cas9_human20220419.fasta'
    # fasta_file = os.path.join(fasta_path, fasta_name)

    # map_protein_seq = read_fasta_file(fasta_file)
    # # util.print_map_column(map_protein_seq)
    # print(len(map_protein_seq))
    # # # # # # # ===


    # # # 测试fasta文件的写出
    # fasta_name = r'cas9_pep_bkp.fasta'
    # fasta_file = os.path.join(fasta_path, fasta_name)

    # write_fasta_file(fasta_file, map_protein_seq)
    # # ===


    # # # === 测试提取前top个蛋白
    # fasta_fpath = r'F:\pFindWork\2103_Rappsilber-PPI-Review\E.coli-Dong-Leiker\fasta\uniprot-human-20171023.fasta'
    # num = 4489
    # fasta_top_fpath = fasta_fpath[:-6] + '_' + str(num) + '.fasta'
    # map_prt_seq = read_fasta_file(fasta_fpath)
    # map_prt_seq = extract_top_fasta(map_prt_seq, num)
    # write_fasta_file(fasta_top_fpath, map_prt_seq)
    # # # ===

    # # === 测试提取固定数目的蛋白，先打乱
    # fasta_fpath = r'F:\pFindResearch\pDeepXL\Data_pLink_Leiker\fasta\OrgDB\uniprot-human-20171023.fasta'
    # num = 20000
    # seed = 6
    # fasta_shuffle_fpath = fasta_fpath[:-6] + '_' + str(num) + '_s' + str(seed) + '.fasta'
    # map_prt_seq = read_fasta_file(fasta_fpath)
    # map_prt_seq = extract_shuffle_fasta(map_prt_seq, num, seed)
    # write_fasta_file(fasta_shuffle_fpath, map_prt_seq)
    # # ===


    # # # === 只获得Decoy序列并重命名
    # fpath = r'F:\pFindWork\2103_Rappsilber-PPI-Review\E.coli-Dong-Leiker\fasta\ecoli_humanNonHomo_4489x2KRp1_R.fasta'
    # fpath = r'F:\pFindWork\2103_Rappsilber-PPI-Review\PXD017620_Yeast-Urlaub-DSS\fasta\yeast_human_6729x2KR_R.fasta'
    # ac2seq = read_fasta_file(fpath)
    # print(len(ac2seq))
    # ac2seq_reverse = get_reverse_rename(ac2seq)
    # print(len(ac2seq_reverse))
    # write_fasta_file(fpath[:-6] + '_R.fasta', ac2seq_reverse)

    # # # ===

    # # # # # # # 重命名fasta的ac
    # fasta_path = r'\\10.29.0.89\dataset\mpz_xl\Mechtler_NC_2020_PXD014337\fasta\db_org'
    # fasta_name = r'uniprot-human-20171023.fasta'
    # fasta_file = os.path.join(fasta_path, fasta_name)

    # k_str = 'human'
    # map_protein_seq = read_fasta_file(fasta_file)
    # map_protein_seq = rename_ac(map_protein_seq, k_str)

    # dpath_out = r'\\10.29.0.89\dataset\mpz_xl\Mechtler_NC_2020_PXD014337\fasta\db_org_rename'
    # fname_out = k_str + '_rename.fasta'
    # fpath_out = os.path.join(dpath_out, fname_out)
    # write_fasta_file(fpath_out, map_protein_seq)

    # util.print_map_column(dict(zip(list(map_protein_seq.keys())[:3], list(map_protein_seq.values())[:3])))
    # print(len(map_protein_seq))
    # # # # # # # ===

    # # # # # # # # 数据库合并
    # dpath_fasta = r'\\10.29.0.89\dataset\mpz_xl\Mechtler_NC_2020_PXD014337\fasta\db_org_ac_rename'
    # fname_fasta1 = r'ecoli.fasta'
    # fname_fasta2 = r'worm.fasta'
    # fname_out = '-'.join([os.path.splitext(fname_fasta1)[0], fname_fasta2])
    # fpath_fasta1 = os.path.join(dpath_fasta, fname_fasta1)
    # fpath_fasta2 = os.path.join(dpath_fasta, fname_fasta2)
    # fpath_out = os.path.join(dpath_fasta, fname_out)
  
    # ac2seq1 = read_fasta_file(fpath_fasta1)
    # ac2seq2 = read_fasta_file(fpath_fasta2)
    # ac2seq_merge = merge_fasta_file(ac2seq1, ac2seq2)

    # write_fasta_file(fpath_out, ac2seq_merge)

    # util.print_map_column(dict(zip(list(ac2seq_merge.keys())[:3], list(ac2seq_merge.values())[:3])))
    # print(len(ac2seq_merge))
    # # # # # # # # ===

    # # # # # # # # 数据库I2L
    # fpath = r'F:\pFindWork\2106_pLink-DSSO_PPI\rerank\blast\test\cas9pep_human.fasta'
    # ac2seq = read_fasta_i2l(fpath, is_i2l=True, is_shortac=False)
    # write_fasta_file(fpath[:-6]+'_i2l.fasta', ac2seq)
    # # # # # # # # ===

    # # # # # # # 数据库合并，批量
    # org_paths = [r'\\10.29.0.89\dataset\mpz_xl\Mechtler_NC_2020_PXD014337\fasta\db_org\cas9pep.fasta']
    # trap_path = r'\\10.29.0.89\dataset\mpz_xl\Mechtler_NC_2020_PXD014337\fasta\db_org_rename'
    # trap_paths = [os.path.join(trap_path, 'plus10.fasta'), os.path.join(trap_path, 'crapome.fasta'), os.path.join(trap_path, 'ecoli1k.fasta'), os.path.join(trap_path, 'ecoli.fasta'), os.path.join(trap_path, 'worm.fasta'), os.path.join(trap_path, 'human.fasta'), os.path.join(trap_path, 'crapome-ecoli.fasta'), os.path.join(trap_path, 'crapome-human.fasta'), os.path.join(trap_path, 'ecoli-human.fasta'), os.path.join(trap_path, 'ecoli-worm-human.fasta')]
    # outpath = r'\\10.29.0.89\dataset\mpz_xl\Mechtler_NC_2020_PXD014337\fasta\db_run_rename'

    # batch_merge_fasta(org_paths, trap_paths, outpath)
    # # # # # # # ===


    # # # # # # # 根据蛋白名称提取蛋白
    fpath_fasta = r'\\10.68.0.17\users\pzmao\pFindWork\2206_CaoYongSTY\data_Ecoli-Leiker-15N\fasta\runDB\Ecoli.fasta'
    # fpath_ac_ls = r'\\10.68.0.17\users\pzmao\pFindWork\2206_CaoYongSTY\SpacePfind\Ecoli-Leiker-5raw_leiker-mono\result\pFind_protein_list_target_ProteinGroupsFDR1.csv'
    fpath_ac_ls = r'\\10.68.0.17\users\pzmao\pFindWork\2206_CaoYongSTY\SpacePfind\Ecoli-Leiker-F2_leiker-mono\result\pFind_protein_list_target_ProteinGroupsFDR1.csv'

    # fpath_fasta_out = r'\\10.68.0.17\users\pzmao\pFindWork\2206_CaoYongSTY\data_Ecoli-Leiker-15N\fasta\runDB\Ecoli-Reduce5Raw.fasta'
    fpath_fasta_out = r'\\10.68.0.17\users\pzmao\pFindWork\2206_CaoYongSTY\data_Ecoli-Leiker-15N\fasta\runDB\Ecoli-ReduceF2.fasta'

    ac_ls = list(pd.read_csv(fpath_ac_ls)['protein'])
    print('===要提取的蛋白数：%d'%(len(ac_ls)))

    ac2seq = read_fasta_file(fpath_fasta)
    ac2seq_new = extract_ac(ac2seq, ac_ls, cmp_shortac=True)

    write_fasta_file(fpath_fasta_out, ac2seq_new)

    # # # # # # # ===


if __name__ == "__main__":
    _main()
