# -*- coding: utf-8 -*-
import os, sys, re
import pandas as pd
import argparse, textwrap


########## < Core Global Variables > ##########

std_MAIN_PROCESS_NAME = "\n[%s]: " % (os.path.basename(__file__))
std_ERROR_MAIN_PROCESS_NAME = "\n[%s::ERROR]: " % (os.path.basename(__file__))
std_WARNING_MAIN_PROCESS_NAME = "\n[%s::WARNING]: " % (os.path.basename(__file__))

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]
HLA_names2 = ["A", "C", "B", "DRB1", "DQA1", "DQB1", "DPA1", "DPB1"]
isREVERSE = {'A': False, 'C': True, 'B': True, 'DRB1': True, 'DQA1': False, 'DQB1': True, 'DPA1': True, 'DPB1': False}


def HLAtoSequences(_hped, _dictionary_seq, _type, _out,
                   _return_as_dataframe=False):



    ########## < Core Variables > ##########


    HLA_DICTIONARY = pd.DataFrame()
    HLA_DICTIONARY2 = {}
    ALLELES_SEQ_LENGTH = {}




    ########## < Argument checking > ##########

    # (1) ped file existence
    if not os.path.exists(_hped):
        print(_hped)
        print(std_MAIN_PROCESS_NAME + "Given ped file dosen't exist. Please check it again.\n")
        sys.exit()

    # (2) HLA DICTIONARY file
    if not os.path.exists(_dictionary_seq):
        print(std_MAIN_PROCESS_NAME + "Given dictionary file dosen't exist. Please check it again.\n")
        sys.exit()

    # (3) Chekcing `_type`
    if not (_type == "AA" or _type == "SNPS"):
        print(std_MAIN_PROCESS_NAME + "Given value for argument `_type` has wrong value. Please check it again.\n")
        sys.exit()




    ########## < Control Flags > ##########

    LOADING_DICTIONARY = 1
    LOADING_HPED = 1
    CLEANING_HPED = 1
    BRINGING_SEQUENCES = 1
    EXPORTING_OUTPUT_PED = 1



    if LOADING_DICTIONARY:

        ########## <1. Dictionary Preparation> ##########

        # print("\n[1] Loading \"Old\" Dictionary Data(created by Sherman Jia.")

        HLA_DICTIONARY = pd.read_table(_dictionary_seq, sep='\t', header=None, names=["Alleles", "Seqs", "INS"], index_col=0).fillna("")

        # Dividing `HLA_DICTIONARY` in advance.
        HLA_DICTIONARY2 = {HLA_names[i]: HLA_DICTIONARY.filter(regex= ''.join(["^", HLA_names[i], "\:"]), axis=0) for i in range(0, len(HLA_names))}

        # for i in range(0, len(HLA_names)):
        #     print("\nSequence Information of %s" % HLA_names[i])
        #     print(HLA_DICTIONARY2[HLA_names[i]].head())

        ALLELES_SEQ_LENGTH = {HLA_names[i] : (len(HLA_DICTIONARY2[HLA_names[i]].iat[0, 0]) + len(HLA_DICTIONARY2[HLA_names[i]].iat[0, 1])) for i in range(0, len(HLA_names))}



    if LOADING_HPED:

        ########## <2. Loading Coated PED(Input PED) file> ##########

        # print("\n[2] Loading Input PED file.")

        HPED = pd.read_table(_hped, sep='\t', header=None, index_col=[0, 1, 2, 3, 4, 5], dtype=str)
        HPED.columns = pd.Index([name + '_' + str(j + 1) for name in HLA_names for j in range(0, 2)])

        # print(HPED.head())



    if CLEANING_HPED:

        ########## <2. Loading Coated PED(Input PED) file> ##########

        # print("\n[3] Cleaning HLA allele names.")

        CHPED = pd.concat([HPED.filter(regex='_'.join([HLA_names[i], '\d{1}']), axis=1).applymap(lambda x : NomenCleaner_forOld(x, HLA_names[i], HLA_DICTIONARY2[HLA_names[i]].index.to_series()) if x != "0" else x) for i in range(0, len(HLA_names))], axis=1)
        # print(CHPED.head())

        # Exporting
        # CHPED.to_csv(_out+".chped", sep='\t', header=False, index=True)



    if BRINGING_SEQUENCES:

        ########## <3. Bringing Sequences> ##########

        # print("\n[4]: Transforming Allele names to Sequences.")

        l_FOR_OUTPUT = []

        for i in range(0, len(HLA_names2)):

            ### Bringing Sequences

            curr_dict = HLA_DICTIONARY2[HLA_names2[i]]

            # HLA columns(each 2 columns)
            df_temp = CHPED.filter(regex='_'.join([HLA_names2[i], '\d{1}']), axis=1).applymap(lambda x : BringSequence(x, HLA_names2[i], curr_dict) if x != "0" else x)
            # print(df_temp.head())



            ### Chop up

            COLUMNS_BY_HLA = []

            for j in range(0, len(df_temp)):

                if df_temp.iat[j, 0] != '0' and df_temp.iat[j, 1] != '0':

                    # (Case1) Most normal case - when allele_name is found as pair.

                    # seq1 = df_temp.iat[j, 0] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 0][::-1]
                    # seq2 = df_temp.iat[j, 1] if not isREVERSE[HLA_name[i]] else df_temp.iat[j, 1][::-1]

                    seq1 = df_temp.iat[j, 0]
                    seq2 = df_temp.iat[j, 1]

                    PAIRED = [value for item in zip(seq1, seq2) for value in item]

                else:

                    # (Case2) when not found as a pair of alleles, but as a only single allele, => 0 is given
                    # (0, 0) will be compensated as length of `HLA_seq`, Due to here, I need to prepared `len(HLA_seq)` beforehand.

                    seq1 = ['0' for z in range(0, ALLELES_SEQ_LENGTH[HLA_names2[i]])]

                    PAIRED = [value for item in zip(seq1, seq1) for value in item]

                COLUMNS_BY_HLA.append(PAIRED)

            l_FOR_OUTPUT.append(pd.DataFrame(COLUMNS_BY_HLA))



    if EXPORTING_OUTPUT_PED:

        ########## <4. Exporting OUTPUT PED file> ##########

        # print("\n[4]: Exporting OUTPUT PED file.")

        df_OUTPUT = pd.concat(l_FOR_OUTPUT, axis=1)
        df_OUTPUT.index = HPED.index
        df_OUTPUT.columns = pd.Index(range(0, df_OUTPUT.shape[1]))
        # print(df_OUTPUT.head())


        ### Final Output ped file.

        if not _return_as_dataframe:

            if _type == 'AA':
                df_OUTPUT.to_csv(_out + '.AA.ped', sep='\t', header=False, index=True)
                return _out+".AA.ped"

            elif _type == 'SNPS':
                df_OUTPUT.to_csv(_out + '.SNPS.ped', sep='\t', header=False, index=True)
                return _out+".SNPS.ped"

        else:
            return df_OUTPUT




def BringSequence(_single_allele, _hla, _dict):

    # try:
    #     Seq = _dict.loc[_single_allele, "Seqs"]
    # except KeyError:
    #     Seq = "0"
    #
    # return Seq

    df_temp = _dict.filter(regex=re.escape(_single_allele), axis=0)

    if df_temp.shape[0] > 0:
        Seq = df_temp.iat[0, 0] if not isREVERSE[_hla] else df_temp.iat[0, 0][::-1]
        Ins = df_temp.iat[0, 1]

        return (Seq + Ins)
    else:
        return "0"


def NomenCleaner_forOld(_single_allele, _hla, _sr_alleles):


    if len(_single_allele) == 2 or len(_single_allele) == 3:
        return ':'.join([_hla, _single_allele])


    if len(_single_allele) == 4:

        p = re.compile('\:'.join([_hla, _single_allele[0:2], _single_allele[2:4]]))
        flag = _sr_alleles.str.match(p)

        if flag.any():
            return ':'.join([_hla, _single_allele[0:2], _single_allele[2:4]])
        else:
            return "0"


    elif len(_single_allele) == 5:

        p_tag = re.compile(r'.+[A-Z]$')
        hasTag = p_tag.match(_single_allele)

        if hasTag:
            t_name = _single_allele[:-1]
            t_name_tag = _single_allele[-1]

            return ':'.join([_hla, t_name[0:2], t_name[2:4]]) + t_name_tag

        else:
            p_2_3 = re.compile('\:'.join([_hla, _single_allele[0:2], _single_allele[2:5]]))
            p_3_2 = re.compile('\:'.join([_hla, _single_allele[0:3], _single_allele[3:5]]))

            flag_2_3 = _sr_alleles.str.match(p_2_3)
            flag_3_2 = _sr_alleles.str.match(p_3_2)

            if flag_2_3.any():
                return ':'.join([_hla, _single_allele[0:2], _single_allele[2:5]])
            elif flag_3_2.any():
                return ':'.join([_hla, _single_allele[0:3], _single_allele[3:5]])
            else:
                return "0"

    # elif len(_single_allele) == 6:
    #     p_3_3 = re.compile('\:'.join([_hla, _single_allele[0:3], _single_allele[3:6]]) + "[A-Z]?$")
    #     flag_3_3 = _sr_alleles.str.match(p_3_3)
    #
    #     if flag_3_3.any():
    #         return p_3_3
    #     else:
    #         return "-1"

    else:
        return "0"






if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description=textwrap.dedent('''\
    #########################################################################################

        HLAtoSequences.py (2018. 11. 8.)
        
        - This version of "HLAtoSequences.py" is the module for compatitability with old 
        version of MakeReference(check the project "MakeReference_v3") 
        - This script Converts HLA alleles (in .ped file format) to amino acid or DNA sequences

        Input file should contain: FID, IID, pID, mID, sex, pheno, HLA-A (2), B (2), C (2), 
        DPA1 (2), DPB1 (2), DQA1 (2), DQB1 (2), DRB1 (2) ... Broad Order
        
        
        - return type : (1) written file name or (2) pandas.DataFrame object(*.ped file)
        

    #########################################################################################
                                     '''),
                                     add_help=False
                                     )

    parser._optionals.title = "OPTIONS"
    # parser._optionals.description = "- Necessary main options.\n"

    parser.add_argument("-h", "--help", help="\nShow this help message and exit\n\n", action='help')

    parser.add_argument("-hped", help="\nHLA Type Data(Standard 4-field allele \"*.ped\" file).\n\n", required=True)
    parser.add_argument("-dict", help="\nHLA dictonary file name(ex. 'HLA_DICTIONARY_AA.txt')\n\n", required=True)
    parser.add_argument("-type", help="\nAA(for Amino Acid) or SNP(for SNPs)\n\n", choices=["AA", "SNPS"], required=True)
    parser.add_argument("-o", help="\nOutput file prefix.\n\n", required=True)




    ##### <for Test> #####

    # # (2018. 8. 27.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Data/HATK/data/HLA_Analysis/AssociationTest/CancerResearch_Example/data_Rev_merged.4field.ped",
    #                           "-dict", "/Users/wansun/Data/HATK/data/MakeReference/HLA_DICTIONARY_AA.hg19.imgt3320.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/CHECK_HLAtoSeuqences.txt"])

    # # (2018. 8. 28.)
    # args = parser.parse_args(["-ped", "/Users/wansun/Data/HATK/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/IMGT370_fixed/HLA_DICTIONARY_AA.hg18.imgt370.txt",
    #                           "-type", "AA",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK"])

    # args = parser.parse_args(["-ped", "/Users/wansun/Data/HATK/data/MakeReference/HAPMAP_CEU_HLA.4field.ped",
    #                           "-dict", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/IMGT370_fixed/HLA_DICTIONARY_SNPS.hg18.imgt370.txt",
    #                           "-type", "SNPS",
    #                           "-o", "/Users/wansun/Git_Projects/HATK_2nd/hatk_2nd/TEST_HATK"])


    ##### <for Publication> #####

    args = parser.parse_args()


    print(args)


    # Implementing Main Function
    HLAtoSequences(args.hped, args.dict, args.type, args.o)


