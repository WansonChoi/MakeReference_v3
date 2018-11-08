#-*- coding: utf-8 -*-

import os, sys, logging,re
import argparse, textwrap
from src.HLAtoSequences_forOld import HLAtoSequences
from src.encodeHLA_forOld import encodeHLA
from src.encodeVariants_forOld import encodeVariants


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]


def MakeReference_v3(**kwargs):


    ########## < Argument Checking > ##########

    # indispensable arguments
    _arguments = ["_hped", "_input", "_out", "_dict_AA", "_dict_SNPS"]

    # optional arguments
    _arguments2 = ["_p_plink", "_p_beagle", "_p_linkage2beagle"]


    # Checking indispensable arguments.
    if not (kwargs and (set(_arguments).issubset(set(kwargs.keys())))):
        print(std_ERROR_MAIN_PROCESS_NAME + "Not all of required arguments were given. Please check them again.\n")
        sys.exit()


    # Checking optional arguments.
    for item in _arguments2:
        if item not in kwargs.keys():
            kwargs[item] = None

    print(kwargs)






    ########## < Core Variables > ##########


    ### Path variables.
    _p_plink = kwargs["_p_plink"] if kwargs["_p_plink"] else "dependency/plink"
    _p_beagle = kwargs["_p_beagle"] if kwargs["_p_beagle"] else "dependency/beagle.jar"
    _p_linkage2beagle = kwargs["_p_linkage2beagle"] if kwargs["_p_linkage2beagle"] else "dependency/linkage2beagle.jar"



    ### Dictionary files
    _dictionary_AA_seq = kwargs["_dict_AA"] + ".txt" # From now on, official extension of HLA sequence information dictionary is ".txt". (2018. 9. 25.)
    _dictionary_AA_map = kwargs["_dict_AA"] + ".map"

    _dictionary_SNPS_seq = kwargs["_dict_SNPS"] + ".txt"
    _dictionary_SNPS_map = kwargs["_dict_SNPS"] + ".map"



    ### Intermediate path.
    _out = kwargs["_out"] if not kwargs["_out"].endswith('/') else kwargs["_out"].rstrip('/')
    if bool(os.path.dirname(_out)): os.makedirs(os.path.dirname(_out), exist_ok=True)




    ########## < Dependency Checking > ##########

    # Plink
    if not os.path.exists(_p_plink):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"Plink(1.07)\" not found. Please check it again(\"{0}\").\n".format(_p_plink))
        sys.exit()
    # Beagle
    if not os.path.exists(_p_beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"beagle.jar\" not found. Please check it again(\"{0}\").\n".format(_p_beagle))
        sys.exit()
    # linkage2beagle.jar
    if not os.path.exists(_p_linkage2beagle):
        print(std_ERROR_MAIN_PROCESS_NAME + "\"linkage2beagle.jar\" not found. Please check it again(\"{0}\").\n".format(_p_linkage2beagle))
        sys.exit()


    # dictionaries
    if not os.path.exists(_dictionary_AA_seq):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA AA sequence file not found. Please check it again(\"{0}\").\n".format(_dictionary_AA_seq))
        sys.exit()
    if not os.path.exists(_dictionary_AA_map):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA AA map file not found. Please check it again(\"{0}\").\n".format(_dictionary_AA_map))
        sys.exit()
    if not os.path.exists(_dictionary_SNPS_seq):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA SNPS sequence file not found. Please check it again(\"{0}\").\n".format(_dictionary_SNPS_seq))
        sys.exit()
    if not os.path.exists(_dictionary_SNPS_map):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA SNPS map file not found. Please check it again(\"{0}\").\n".format(_dictionary_SNPS_map))
        sys.exit()









    ########## < Main > ##########

    ### [1] ENCODE_AA


    ### [2] ENCODE_HLA


    ### [3] ENCODE_SNPS


    ### [4] EXTRACT_FOUNDERS


    ### [5] MERGE


    ### [6] QC


    ### [7] PREPARE


    ### [8] PHASE


    ### [9] CLEANUP




    return 0






def b_MARKER_HLA(**kwargs):

    if "_type" not in kwargs.keys():
        print(std_ERROR_MAIN_PROCESS_NAME + "\"_type\" argument must be passed.\n")
        return -1
    else:
        _type = kwargs["_type"] if kwargs["_type"] == "AA" or kwargs["_type"] == "SNPS" or kwargs["_type"] == "HLA" else None





    return 0



if __name__ == "__main__":


    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #################################################################################################

        MakeReference_v3.py

        Generating reference panel for SNP2HLA based on HLA sequence information dictionaries.

        (ex.)
        : python3 MakeReference_v3.py  
            -hped ./data/b_MarkerGenerator.py/HAPMAP_CEU_HLA.4field.ped 
            -i ./Trial_HAPMAP_CEU
            -dict-AA ./data/b_MarkerGenerator.py/HLA_DICTIONARY_AA.hg18.imgt370
            -dict-SNPS ./data/b_MarkerGenerator.py/HLA_DICTIONARY_SNPS.hg18.imgt370 

        HLA PED file should contain HLA alleles in the following (alphabetical) order:
        HLA-A, B, C, DPA1, DPB1, DQA1, DQB1, DRB1

    #################################################################################################
                                     '''),
                                     # epilog="-*- Recoded to Python script by Wansun Choi in Han lab. at Asan Medical Center -*-",
                                     add_help=False)

    parser._optionals.title = "OPTIONS"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("--input", "-i", help="\nInput Data file(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("-hped", help="\nHLA Type Data(.ped)\n\n", required=True)

    parser.add_argument("-dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", required=True)
    parser.add_argument("-dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", required=True)
    # parser.add_argument("-b", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n")





    ##### <for Test> #####
    args = parser.parse_args(["-hped", "data/HAPMAP_CEU_HLA.ped",
                              "-i", "data/HAPMAP_CEU",
                              "-dict-AA", "data/HLA_DICTIONARY_AA",
                              "-dict-SNPS", "data/HLA_DICTIONARY_SNPS",
                              "-o", "test/20181108/Trial_1"])


    ##### <for Publications> #####
    # args = parser.parse_args()
    print(args)




    MakeReference_v3(_hped=args.hped, _input=args.input, _out=args.out,
                     _dict_AA=args.dict_AA, _dict_SNPS=args.dict_SNPS)

    # MakeReference_v3(_hped=args.hped, _input=args.input, _out=args.out,
    #                  _dict_AA=args.dict_AA, _dict_SNPS=args.dict_SNPS,
    #                  _p_plink = "./dependency/plink", _p_beagle="./dependency/beagle.jar", _p_linkage2beagle="lalal")