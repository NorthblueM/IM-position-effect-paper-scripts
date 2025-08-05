#-*- coding: UTF-8 -*-

# @Date     : 4th Aug, 2024
# @Author   : Northblue

"""get modification mass
"""

import os
from pathlib import Path
MODIFICATION_FILE_DIR = os.path.dirname(os.path.abspath(__file__))

def read_mod2mass(fpath_mod=str(Path(MODIFICATION_FILE_DIR)/'modification.ini')):
    """read modification.ini

    just read mass
    """

    mod2mass = {}
    mod_name = ''
    with open(fpath_mod, 'r') as fin:
        for line in fin:
            if line.startswith('name'):
                mod_name = line.split('=')[1].strip().split()[0]
                continue
            if mod_name and line.startswith(mod_name):
                segs = line.split('=')[1].strip().split()
                mod_mass_mono = float(segs[2])
                mod2mass[mod_name] = mod_mass_mono
                continue
    return mod2mass