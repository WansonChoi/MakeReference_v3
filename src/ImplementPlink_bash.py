#-*- coding: utf-8 -*-

import os, sys, re
from shutil import which


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


PLINK = "dependency/plink" if os.path.exists("dependency/plink") else which("plink")


"""

Python interface for bash execution of Plink.

"""

def make_bed(*args, **kwargs):

    _arguments = ["_file", "_bfile", "_out"]
    _arguments2 = ["_missing_genotype", "_exclude"]


    keys = kwargs.keys()

    # "--file" or "--bfile"
    if "_file" in keys:
        _file = kwargs["_file"]
        _file_c = "--file"
    elif "_bfile" in keys:
        _file = kwargs["_bfile"]
        _file_c = "--bfile"
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Error in file or bfile arguments.")
        return -1

    # "--out"
    if "_out" in keys:
        _out = kwargs["_out"]
    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Error in out arguments.")
        return -1


    # Optional arguments.

    command_optionals = []

    for item in _arguments2:
        if item in kwargs.keys():

            _arguments2_c = re.sub('_', '-', re.sub(r'(^_)', '--', item))

            command_optionals.append(_arguments2_c)
            command_optionals.append(str(kwargs[item]))



    command = ' '.join([PLINK, "--noweb", "--make-bed",
                        _file_c, _file,
                        "--out", _out,
                        ' '.join(command_optionals)])

    print(command)

    if not os.system(command):
        return _out
    else:
        return -1



