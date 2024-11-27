# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 11:01:54 2022
@author: ZR48SA
"""

#%% change directory to script directory
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd()))
print(os.getcwd())
from Download import *
from setup_NCBI_taxonomy import *
from prep_Database import *
from make_DIAMOND_database import *

#%% Set execution


download_DIAMOND=True
download_ncbi_taxonomy=True
download_database=True
construct_DIAMOND_database=True

#Database choices: "RefSeq", "NCBI_NR", "Swiss-Prot", "TrEMBL", "UniProt", "UniRef100", "UniRef90", "UniRef50"
db="Swiss-Prot" 



#%% Step 1: Download DIAMOND

if download_DIAMOND:
    print("downloading DIAMOND")
    diamond_path=download_diamond(path=False)


#%% Step 2: Setup NCBI taxonomy

if download_ncbi_taxonomy:
    print("Setting up NCBI taxonomy")
    names,nodes=download_ncbi_taxdump(path=False)
    ncbi_taxonomy_path=parse_NCBI_taxonomy(names,nodes)

    
#%% Step 3: Setup Diamond database

if download_database:
    print("Downloading database : "+db)
    fasta_files=download_db(db)

if construct_DIAMOND_database:
    print("constructing DIAMOND database : " +db)
    merged_database=merge_files(fasta_files,delete_old=False)
    prepped_database=prep_database(merged_database,delete_old=False)
    diamond_database_path=make_diamond_database(diamond_path,database_path=prepped_database,output_path=False,delete_old=False)

#%% Step 4: Write filepaths to file

path_df=pd.DataFrame([["diamond_path",                  diamond_path],
                      ["ncbi_taxonomy_path",      ncbi_taxonomy_path],
                      ["diamond_database_path",diamond_database_path]],columns=["name","path"]).set_index("name")

path_df.to_csv(Path(basedir,"setup_filepaths.tsv"),sep="\t")

print("done")