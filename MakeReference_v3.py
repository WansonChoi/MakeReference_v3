#-*- coding: utf-8 -*-

import os, sys, re
import argparse, textwrap
import pandas as pd
from numpy import concatenate
from shutil import which

from src.HLAtoSequences_forOld import HLAtoSequences
from src.encodeHLA_forOld import encodeHLA
from src.encodeVariants_forOld import encodeVariants
import src.ImplementPlink_bash as Plink


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))


HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

# Files to be removed in the end.
__REMOVE__ = []


def MakeReference_v3(**kwargs):


    ########## < Argument Checking > ##########

    # indispensable arguments
    _arguments = ["_hped", "_input", "_out", "_dict_AA", "_dict_SNPS"]

    # optional arguments
    _arguments2 = ["_p_plink", "_p_JAVA", "_p_beagle", "_p_linkage2beagle"]


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
    _p_JAVA =  kwargs["_p_JAVA"] if kwargs["_p_JAVA"] else which("java")
    _p_beagle = kwargs["_p_beagle"] if kwargs["_p_beagle"] else "dependency/beagle.jar"
    _p_linkage2beagle = kwargs["_p_linkage2beagle"] if kwargs["_p_linkage2beagle"] else "dependency/linkage2beagle.jar"


    ### HLA type data.
    _hped = kwargs["_hped"]

    ### Traditional SNP markers
    _SNPs = kwargs["_input"]

    _SNPs_dirname = os.path.dirname(_SNPs)
    _SNPs_basename = os.path.basename(_SNPs)


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


    # HLA type data
    if not os.path.exists(_hped):
        print(std_ERROR_MAIN_PROCESS_NAME + "HLA type data can't be found. Please check it again(\"{0}\").\n".format(_hped))

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

    ##### [1] ENCODE_AA
    print(std_MAIN_PROCESS_NAME + "[1] Generating HLA AA binary markers.\n")

    __HLA_AA__ = b_MARKER_HLA(_type="AA", _hped=_hped, _dict=kwargs["_dict_AA"], _out=_out)



    ##### [2] ENCODE_HLA
    print(std_MAIN_PROCESS_NAME + "[2] Generating HLA binary markers.\n")

    __HLA__ = b_MARKER_HLA(_type="HLA", _hped=_hped, _out=_out)



    ##### [3] ENCODE_SNPS
    print(std_MAIN_PROCESS_NAME + "[3] Generating HLA SNPS binary markers.\n")

    __HLA_SNPS__ = b_MARKER_HLA(_type="SNPS", _hped=_hped, _dict=kwargs["_dict_SNPS"], _out=_out)



    ##### [4] EXTRACT_FOUNDERS
    print(std_MAIN_PROCESS_NAME + "[4] Extracting Founders in traditional SNP markers.\n")

    # Filtering for traditional SNP markers.("--input")

    # making bed file.
    __SNP_FOUNDERS__ = Plink.make_bed(_bfile=_SNPs, _out=(_out+"."+_SNPs_basename+".FOUNDERS"), _filter_founders=True, _mind=0.3, _alleleACGT=True)

    # QC("--hardy", "--freq", "--missing")
    __hardy__ = Plink.Quality_Control("--hardy", _bfile=__SNP_FOUNDERS__, _out=__SNP_FOUNDERS__)
    __freq__ = Plink.Quality_Control("--freq", _bfile=__SNP_FOUNDERS__, _out=__SNP_FOUNDERS__)
    __missing__ = Plink.Quality_Control("--missing", _bfile=__SNP_FOUNDERS__, _out=__SNP_FOUNDERS__)

    print(__hardy__)
    print(__freq__)
    print(__missing__)


    __df_hardy__ = pd.read_table(__hardy__, sep='\t|\s+', engine='python', header=0)

    flag_p_value = __df_hardy__.iloc[:, 8] < 0.000001
    p_value = __df_hardy__.loc[flag_p_value, :].iloc[:, 1].unique()
    p_value.sort()
    # print(p_value)


    __df_freq__ = pd.read_table(__freq__, sep='\t|\s+', engine='python', header=0)

    flag_MAF = __df_freq__.iloc[:, 4] < 0.01
    MAF = __df_freq__.loc[flag_MAF, :].iloc[:, 1].unique()
    MAF.sort()
    # print(MAF)


    __df_missing__ = pd.read_table(__missing__, sep='\t|\s+', engine='python', header=0)

    flag_F_MISS = __df_missing__.iloc[:, 4] > 0.05
    F_MISS = __df_missing__.loc[flag_F_MISS, :].iloc[:, 1].unique()
    F_MISS.sort()
    # print(F_MISS)

    target_remove =  list(set(concatenate((p_value, MAF, F_MISS), axis=0)))
    target_remove.sort()

    pd.Series(target_remove).to_csv(_out+".unqualified.txt", sep='\t', header=False, index=False)


    # excluding
    __SNP_FOUNDERS_2__ = Plink.make_bed(_bfile=__SNP_FOUNDERS__, _out=__SNP_FOUNDERS__+".QC",
                         _exclude=_out+".unqualified.txt", _allow_no_sex=True)




    # Filtering generated HLA binary markers(AA, SNPS, HLA)

    __HLA_AA_FOUNDERS__ = Plink.make_bed(_bfile=__HLA_AA__, _out=__HLA_AA__+".FOUNDERS",
                                         _filter_founders=True, _maf=0.0001)
    __HLA_SNPS_FOUNDERS__ = Plink.make_bed(_bfile=__HLA_SNPS__, _out=__HLA_SNPS__+".FOUNDERS",
                                           _filter_founders=True, _maf=0.0001)
    __HLA_HLA_FOUNDERS__ = Plink.make_bed(_bfile=__HLA__, _out=__HLA__+".FOUNDERS",
                                         _filter_founders=True, _maf=0.0001)


    # Garbage Collecting
    __REMOVE__.append(_out + ".unqualified.txt")




    ##### [5] MERGE
    print(std_MAIN_PROCESS_NAME + "[5] Merging whole binary markers.\n")

    bed_suffix = (".bed", ".bim", ".fam")
    l_temp = (__HLA_AA_FOUNDERS__, __HLA_SNPS_FOUNDERS__, __HLA_HLA_FOUNDERS__)

    target_merge = [[item+suffix for suffix in bed_suffix] for item in l_temp]
    print(target_merge)
    pd.DataFrame(target_merge).to_csv(_out+".mergelist.txt", sep='\t', header=False, index=False)


    __MERGED_HLA__ = Plink.make_bed(_bfile=__SNP_FOUNDERS_2__, _out=(_out+"."+_SNPs_basename+".MERGED.FOUNDERS"),
                                    _merge_list = _out+".mergelist.txt")
    print(__MERGED_HLA__)


    # Garbage Collecting
    __REMOVE__.extend([__HLA_AA__+".*", __HLA__+".*", __HLA_SNPS__+".*", _out+".mergelist.txt"])




    ##### [6] QC
    print(std_MAIN_PROCESS_NAME + "[6] Quality Control on MERGED binary markers.\n")

    __HLA_MERGED_freq__ = Plink.Quality_Control("--freq", _bfile=__MERGED_HLA__, _out=__MERGED_HLA__+".FRQ")

    __df_HLA_MERGED_freq__ = pd.read_table(__HLA_MERGED_freq__, sep='\t|\s+', engine='python', header=0)
    # print(__df_HLA_MERGED_freq__.head())


    # Filtering MAF.
    flag_MAF = (__df_HLA_MERGED_freq__.iloc[:, 4] < 0.0001) | (__df_HLA_MERGED_freq__.iloc[:, 4] > 0.9999)
    MAF = __df_HLA_MERGED_freq__.loc[flag_MAF, :].iloc[:, 1]

    MAF.to_csv(_out+".unqualified2.txt", sep='\t', header=False, index=False)


    # Filtering "P", "A" allele order.
    sr_A1 = __df_HLA_MERGED_freq__.iloc[:, 2]
    sr_A2 = __df_HLA_MERGED_freq__.iloc[:, 3]

    # flag_allele_order = ((sr_A1 == "A") & (sr_A2 == "P")) | ((sr_A1 == "P") & (sr_A2 == "A"))
    flag_allele_order = ((sr_A1 == "A") & (sr_A2 == "P"))
    sr_Labels = __df_HLA_MERGED_freq__.loc[flag_allele_order, "SNP"].reset_index(drop=True)
    sr_temp = pd.Series(["P" for i in range(0, sr_Labels.shape[0])])

    pd.concat([sr_Labels, sr_temp], axis=1).to_csv(_out+".allele.order.txt", sep='\t', header=False, index=False)


    # Filtered MERGED.FRQ(*Final Output*)
    __HLA_MERGED_output__ = Plink.make_bed(_bfile=__MERGED_HLA__, _out=_out, _geno=0.5, _exclude=_out+".unqualified2.txt",
                                           _reference_allele=_out+".allele.order.txt")

    __HLA_MERGED_freq_2__ = Plink.Quality_Control("--freq", _bfile=__HLA_MERGED_output__, _out=__HLA_MERGED_output__+".FRQ",
                                                  _keep_allele_order=True)

    # Garbage Collecting
    __REMOVE__.extend([_out+".unqualified2.txt", __MERGED_HLA__+".*", _out+".allele.order.txt", __SNP_FOUNDERS__+".*"])



    ##### [7] PREPARE
    print(std_MAIN_PROCESS_NAME + "[7] Preparing Beagle files.\n")

    df_markers = pd.read_table(__HLA_MERGED_output__+".bim", sep='\t|\s+', engine='python', header=0, dtype=str).iloc[:, [1, 3, 4, 5]]
    df_markers.to_csv(_out+".markers", sep='\t', header=False, index=False)

    # final output in .ped/.map format.
    __HLA_MERGED_output_2__ = Plink.recode(_bfile=__HLA_MERGED_output__, _out=_out,
                                           _keep_allele_order=True, _alleleACGT=True)

    df_dat = pd.read_table(__HLA_MERGED_output_2__+".map", sep='\t|\s+', engine='python', header=None, dtype=str, usecols=[1])
    df_dat.index = pd.Index(["M" for i in range(0, df_dat.shape[0])])
    df_dat.to_csv(_out+".dat", sep='\t', header=False, index=True)

    df_nopheno = pd.read_table(__HLA_MERGED_output_2__+".ped", sep='\t|\s+', engine='python', header=None, dtype=str).drop(5, axis=1)
    df_nopheno.to_csv(_out+".nopheno.ped", sep='\t', header=False, index=False)

    # linkage2beagle.jar
    command = ' '.join([_p_JAVA, "-Xmx500m", "-jar", _p_linkage2beagle, "pedigree={0}".format(_out+".nopheno.ped"), "data={0}".format(_out+".dat"),
                        "beagle={0}".format(_out+".bgl"), "standard=true", ">", "{0}".format(_out+".bgl.log")])
    print(command)

    if bool(os.system(command)):
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to implement linkage2beagle.jar\n")
        sys.exit()


    # Garbage Collecting
    __REMOVE__.extend([(_out + item) for item in [".nopheno.ped", ".ped", ".map", ".dat"]])



    ##### [8] PHASE
    print(std_MAIN_PROCESS_NAME + "[8] Beagle Phasing.\n")

    # beagle.jar
    command = ' '.join([_p_JAVA, "-Xmx500m", "-jar", _p_beagle, "maxwindow=1000", "verbose=true", "missing=0", "nsamples=4", "niterations=10",
                        "unphased={0}".format(_out+".bgl"), "log={0}".format(_out+".phasing"), ">>", "{0}".format(_out+".bgl.log")])
    print(command)

    if bool(os.system(command)):
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to implement beagle.jar\n")
        sys.exit()


    # Garbage Collecting
    __REMOVE__.extend([(_out + item) for item in [".bgl.gprobs", ".bgl.r2", ".bgl", ".phasing.log"]])




    ##### [9] CLEANUP
    print(std_MAIN_PROCESS_NAME + "[9] Clean-up.\n")

    command = ' '.join(["rm", ' '.join(__REMOVE__)])
    print(command)

    if bool(os.system(command)):
        print(std_ERROR_MAIN_PROCESS_NAME + "Failed to remove target files successfully.\n")
        return -1



    return 0






def b_MARKER_HLA(**kwargs):

    """
    - _type (AA, SNPS, HLA)
    - _dictionary_seq (AA, SNPS)
    - _dictionary_map (AA, SNPS)
    - _out

    """

    if "_type" not in kwargs.keys():
        print(std_ERROR_MAIN_PROCESS_NAME + "\"_type\" argument must be passed.\n")
        return -1
    else:
        TYPE = kwargs["_type"] if kwargs["_type"] == "AA" or kwargs["_type"] == "SNPS" or kwargs["_type"] == "HLA" else None



    if TYPE == "AA" or TYPE == "SNPS":

        # Passed arguments
        HPED = kwargs["_hped"]
        OUT = kwargs["_out"]
        _dictionary_seq = kwargs["_dict"] + ".txt"
        _dictionary_map = kwargs["_dict"] + ".map"


        df_ped = HLAtoSequences(_hped=HPED, _dictionary_seq=_dictionary_seq, _type=TYPE, _out=OUT, _return_as_dataframe=True)
        # df_ped = HLAtoSequences(_hped=HPED, _dictionary_seq=_dictionary_seq, _type=TYPE, _out=OUT)

        __ENCODED_1__ = encodeVariants(_p_ped=df_ped, _p_map=_dictionary_map, _out=OUT+".{0}.ENCODED".format(TYPE))


        # Plink : --missing-genotype
        __ENCODED_2__ = Plink.make_bed(_file=__ENCODED_1__, _out=OUT+".{0}.RAW".format(TYPE), _missing_genotype=0)
        print(__ENCODED_2__)


        # Plink : --exclude
        df_temp = pd.read_table(__ENCODED_2__+".bim", sep='\t|\s+', engine='python', header=None, dtype=str)
        print(df_temp.head())

        flag_5_0 = (df_temp.iloc[:, 4] == '0')
        flag_5_x = (df_temp.iloc[:, 4] == 'x')
        flag_6_x = (df_temp.iloc[:, 5] == 'x')
        flag_INS = (df_temp.iloc[:, 1].str.match(r'^INS'))

        flag_union = (flag_5_0 | flag_5_x | flag_6_x) & ~flag_INS

        df_temp2 = df_temp.loc[flag_union, :].iloc[:, 1]
        df_temp2.to_csv(OUT+".{0}.noPolymorphic.txt".format(TYPE), sep='\t', header=False, index=False)


        __ENCODED_3__ = Plink.make_bed(_bfile=__ENCODED_2__, _out=OUT+".{0}.CODED".format(TYPE), _exclude=OUT+".{0}.noPolymorphic.txt".format(TYPE))
        print(__ENCODED_3__)


        # Clean-up
        rm_pattern = ["*.ENCODED.*", "*.RAW.*", "*.noPolymorphic.*"]
        rm_list = [os.path.join(os.path.dirname(OUT), rm_pattern[i]) for i in range(0, len(rm_pattern))]

        command = ' '.join(["rm", ' '.join(rm_list)])
        print(command)
        os.system(command)


        return __ENCODED_3__



    elif TYPE == "HLA":

        # Passed arguments
        HPED = kwargs["_hped"]
        OUT = kwargs["_out"]

        __ENCODED_1__ = encodeHLA(HPED, OUT+".HLA")
        __ENCODED_2__ = Plink.make_bed(_file=__ENCODED_1__, _out=OUT+".HLA")

        return __ENCODED_2__

    else:
        print(std_ERROR_MAIN_PROCESS_NAME + "Something wrong with \"_type\" argument of b_MARKER_HLA.")
        return -1



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

    parser.add_argument("--input", "-i", help="\nInput SNP data file(.bed/.bim/.fam)\n\n", required=True)
    parser.add_argument("--out", "-o", help="\nOutput file prefix\n\n", required=True)

    parser.add_argument("-hped", help="\nHLA Type Data(.ped)\n\n", required=True)

    parser.add_argument("-dict-AA", help="\nPrefix of AA HLA Dictionary file(*.txt, *.map).\n\n", required=True)
    parser.add_argument("-dict-SNPS", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n", required=True)
    # parser.add_argument("-b", help="\nPrefix of SNP HLA Dictionary file(*.txt, *.map).\n\n")





    ##### <for Test> #####
    # args = parser.parse_args(["-hped", "data/HAPMAP_CEU_HLA.ped",
    #                           "-i", "data/HAPMAP_CEU",
    #                           "-dict-AA", "data/HLA_DICTIONARY_AA",
    #                           "-dict-SNPS", "data/HLA_DICTIONARY_SNPS",
    #                           "-o", "test/20181121_1/20181121_SNPfixed2"])


    ##### <for Publications> #####
    args = parser.parse_args()
    print(args)




    MakeReference_v3(_hped=args.hped, _input=args.input, _out=args.out,
                     _dict_AA=args.dict_AA, _dict_SNPS=args.dict_SNPS)

    # MakeReference_v3(_hped=args.hped, _input=args.input, _out=args.out,
    #                  _dict_AA=args.dict_AA, _dict_SNPS=args.dict_SNPS,
    #                  _p_plink = "./dependency/plink", _p_beagle="./dependency/beagle.jar", _p_linkage2beagle="lalal")