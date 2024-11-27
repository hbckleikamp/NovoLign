# -*- coding: utf-8 -*-
"""
Created on Fri Oct 29 15:55:59 2021
@author: hbckleikamp | updated mpabst August 2023
"""

# =============================================================================
# IMPORT MODULES, DEFINE PATHS AND PARAMETERS
# =============================================================================

import glob, os
from pathlib import Path
from inspect import getsourcefile
import pandas as pd
import time
from bar_graphs import graphs
start=time.time()
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
basedir=os.getcwd()

# =============================================================================
# BEGIN DEFINE PARAMETERS
# =============================================================================

##### Input filepaths ##### 

#input_files=glob.glob("".join((basedir,"\*Input_*"))) # location of input folders
# fasta_database_path  =str(Path(basedir,"Setup","Swiss-Prot\\uniprot_sprot_NoAmb_IJeqL.fasta")) 
# diamond_database_path=str(Path(basedir,"Setup","Swiss-Prot\\uniprot_sprot_NoAmb_IJeqL.dmnd"))



input_folders=["F:/Phages in anaerobic digesters/Denovo/Sample 1.denovo.csv"]
fasta_database_path  ="F:/nr/nr/nr_NoAmb_IJeqL.fasta" 
diamond_database_path="F:/nr/nr/nr_NoAmb_IJeqL.dmnd"

diamond_path=str(Path(basedir,"Setup","diamond","diamond.exe"))
diamond_folder=str(Path(basedir,"Setup","diamond"))
ncbi_taxonomy_path=   str(Path(basedir,"Setup","ncbi_taxonomy","parsed_ncbi_taxonomy.tsv"))      # placeholder path to parsed ncbi taxonomy
Temporary_directory=basedir     # DMD temmporary folder

#### Performance parameters ####
min_ALC_score=70                # minimum ALC(%)
bit=25                          # minimum bitscore
freq_cut=5                      # lineage frequency filter for composition and DB creation
Write_to_database="Proteins"    # Options (False, "Proteins","Taxids"): do not make a database(False), use aligned proteins ("Proteins") use aligned taxids("Taxids").
DB_rank="genus"                 # "OX" "species" "genus" or "family"

# =============================================================================
# END DEFINE PARAMETERS
# =============================================================================    
    
#load NCBI parsed taxonomy
ranks=["superkingdom","phylum","class","order","family","genus","species"]
ncbi_taxdf=pd.read_csv(ncbi_taxonomy_path,sep="\t").fillna("")

#default diamond output columns
output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]
if Write_to_database: output_columns+=["full_sseq"]

# =============================================================================
# START NOVOLIGN PIPELINE
# =============================================================================

from write_to_fasta import write_to_fasta
from diamond_alignment import align
from process_alignment import *
from experiment_qc import Plot_high_scoring
from database_qc import *

for infolder in input_folders:
    print("analysing: "+infolder)

    # define input and NovoLign modules
    de_novo_file=''.join(glob.glob("".join((infolder,"\*de novo* *peptides.csv"))))
    database_searching_file=''.join(glob.glob("".join((infolder,"\*psm.csv"))))
    fasta_database=''.join(glob.glob("".join((infolder,"\*.fasta"))))
    taxa=''.join(glob.glob("".join((infolder,"\*TaxIDs.txt"))))

    """
    Prepare folders and paths
    """
    Output_directory=str(Path(Path(infolder).parents[0],Path(infolder).name.replace("Input_","Output_")))
    if not os.path.exists(Output_directory): os.mkdir(Output_directory)
    if not os.path.exists(Temporary_directory): os.mkdir(Temporary_directory)

    """
    1/5 DMD alignment
    """
    print("Step 1 of 5: Start diamond alignment")
    fasta_file,denovo_peptides=write_to_fasta(de_novo_file=de_novo_file,
                                              database_searching_file=database_searching_file,
                                              Output_directory=Output_directory,base_ALC_Score_cut=40,add_decoy=True)
    dmd_file=str(Path(Output_directory,"diamond_alignments"+"\\de novo peptides.tsv"))
    if os.path.exists(dmd_file) == True: 
        print("Diamond alignment exists")
        alignments=dmd_file
    else:
        alignments=align(fasta_file,Output_directory,diamond_path=diamond_path,database_path=diamond_database_path,output_columns=output_columns)
    target_decoys=Process_diamond_alignment(alignments,Output_directory,minimum_bitscore=bit,min_ALC_score=min_ALC_score,output_columns=output_columns)

    """
    2/5 LCAs for alignments
    """
    print("Step 2 of 5: Construct LCAs")
    for target_decoy in target_decoys: 
        # 1. Conventional lca
        denovo_peptides_lca=lca(target_decoy,Output_directory,
                                denovo_peptides=denovo_peptides,
                                method="conventional",
                                filter_cutoff=freq_cut,
                                minimum_rank=DB_rank,
                                write_database=Write_to_database)
        # 2. Bitsore lca 
        denovo_peptides_blca=lca(target_decoy,Output_directory,
                                 denovo_peptides=denovo_peptides,
                                 method="bitscore",weight_column="bitscore",
                                 filter_cutoff=freq_cut,
                                 minimum_rank=DB_rank,
                                 write_database=Write_to_database)
        # # 3. Weighted lca
        denovo_peptides_wlca=lca(target_decoy,Output_directory,
                                 denovo_peptides=denovo_peptides,
                                 method="weighted",weight_column="weights",
                                 filter_cutoff=freq_cut,
                                 minimum_rank=DB_rank,
                                 write_database=Write_to_database)

    """
    3/5 Grouped taxonomy report
    """
    print("Step 3 of 5: Get DN compositions")
    files=glob.glob("".join((Output_directory,"\lca\*.tsv")))
    for file in files: graphs(file,Output_directory,ranks,freq_cut)

    """
    4/5 Check quality of fragmentation spectra
    """  
    print("Step 4 of 5: Check spectral quality")
    Plot_high_scoring(denovo_peptides,target_decoys[0],database_searching_file,Output_directory,de_novo_file)
    

    """
    5/5 Check coverage obtained by database searching
    """
    print("Step 5 of 5: Check DB searching output")
    Compare_Bar(denovo_peptides_wlca,database_searching_file,fasta_database,taxa,Output_directory,ncbi_taxdf)

# =============================================================================
# DONE
# =============================================================================

end=time.time()
print("NovoLign completed: Run time",round((end-start)/60,2),"minutes")
globals().clear()
