# MakeReference_v3
Third generation of MakeReference.


to Seungho.

(Requirements.)

1. make "dependency/" folder in the project foler and prepare 

    (1) plink(1.07)
    (2) beagle.jar
    (3) linkage2beagle.jar

2. Implement

ex) 
python3 MakeReference_v3.py \
        -hped data/HAPMAP_CEU_HLA.ped \
        -i data/HAPMAP_CEU \
        -dict-AA data/HLA_DICTIONARY_AA \
        -dict-SNPS data/HLA_DICTIONARY_SNPS \
        -o Ref_Panel_v3


