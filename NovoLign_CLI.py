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
# Default arguments
# =============================================================================
input_file=""
diamond_datbase=""

# =============================================================================
# Parse arguments
# =============================================================================



import argparse

parser = argparse.ArgumentParser(
                    prog='NovoLign',
                    description='DIAMOND alignment of De Novo Sequenced peptides')

#required arguments
parser.add_argument("-i", "--input_file",       required = True, 
                    help="Required: input peptide file, should be tabular and contain the column 'Peptide'.\
                        if a folder is supplied, it should have the filstructure as described on https://github.com/hbckleikamp/NovoLign")
parser.add_argument("-d", "--diamond_database_path", required = True, help="Required: path to diamond database")                   

#Output_folder
parser.add_argument("-o", "--Output_directory", default="", required = False, help="Output folder")      

#Performance parameters
parser.add_argument("-ALC", "--min_ALC_score", required = False, default=70,  help="Minimum ALC score (PEAKS specific score)") 
parser.add_argument("--bit", "-min_bit_score", required = False, default=25,  help="Minimum bitscore (DIAMOND alignment score)") 
parser.add_argument("-freq", "--freq_cut",     required = False, default=5,   help="Minium taxa frequency for denoising") 
parser.add_argument("-lcas", required = False, default=['lca','bitlca','wcla'],   help="Wich LCA algorithm to use") 


#Optional: Database construction
parser.add_argument("-f", "--fasta_database_path", required = False, default="", help="Required for database construction: Path to fasta database")  
parser.add_argument("-DBwrite", "--Write_to_database", required = False, default="Proteins", 
                    help="Required for database construction. Options: (False, 'Proteins','Taxids'): do not make a database(False), use aligned proteins ('Proteins') use aligned taxids('Taxids').")  

parser.add_argument("-DB_rank", required = False, default="genus",
                    help="Used in database construction with -DBwrite 'Taxids', selects the taxonomic rank for database construction, Options: 'OX' 'species' 'genus' or 'family'")


#Optional: Comparison with database searching file
parser.add_argument("-di", "--database_searching_file", required = False, default="",
                    help="Optional: database searching file, for comparison with de novo sequenced peptides, should be tabular and contain the column 'Peptide'")
parser.add_argument("-taxa",   required = False, default="", help="Specified taxa for the comparison of database searching and de novo sequencing outputs")


#default filepaths
parser.add_argument("-diamond_path",   required = False, help="path to diamond executable", default=str(Path(basedir,"Setup","diamond","diamond.exe"))) 
parser.add_argument("-diamond_folder", required = False, help="path to diamond exeutable folder", default=str(Path(basedir,"Setup","diamond")))
parser.add_argument("-ncbi_taxonomy_path", required = False, help="path to ncbi taxonomy file", default=str(Path(basedir,"Setup","ncbi_taxonomy","parsed_ncbi_taxonomy.tsv")))  
parser.add_argument("-Temporary_directory", required = False, help="Temporary directory for writing DIAMOND index files", default=basedir) 

args = parser.parse_args()
args = {k:v for k,v in vars(parser.parse_args()).items() if v is not None}

print("")
print(args) #debug
print("")
locals().update(args)

# =============================================================================
# Load arguments
# =============================================================================    


ncbi_taxdf=pd.read_csv(ncbi_taxonomy_path,sep="\t").fillna("")
ranks=["superkingdom","phylum","class","order","family","genus","species"]
output_columns=["qseqid","sseqid","stitle","pident","bitscore","qseq","sseq"]  
if Write_to_database: output_columns+=["full_sseq"]

def read_list(x): 
    if type(x) is list: return x
    return x.strip("[]").split(",")
    
de_novo_files,database_searching_files,lcas=[read_list(i) for i in [input_file,database_searching_file,lcas]]



if os.path.isdir(input_file):
    
    if not len(Output_directory): Output_directory=str(Path(Path(input_file).parents[0],"Output_"+Path(input_file).name))
    

    de_novo_files=''.join(glob.glob("".join((input_file,"\*de novo* *peptides.csv"))))
    database_searching_files=''.join(glob.glob("".join((input_file,"\*psm.csv"))))
    fasta_database=''.join(glob.glob("".join((input_file,"\*.fasta"))))
    taxa=''.join(glob.glob("".join((input_file,"\*TaxIDs.txt"))))

    de_novo_files.sort()
    database_searching_files.sort()


print("De Novo Files found:")
print(de_novo_files)
print("")



# =============================================================================
# START NOVOLIGN PIPELINE
# =============================================================================


from write_to_fasta import write_to_fasta
from diamond_alignment import align
from process_alignment import *
from experiment_qc import Plot_high_scoring
from database_qc import *

toggle_usefileame=False
for ix,de_novo_file in enumerate(de_novo_files):
    print("")
    print("analysing: "+de_novo_file)

    if len(database_searching_files)>ix:
        database_searching_file=database_searching_files[ix]
    else:
        database_searching_file=""
        
    if database_searching_file=="":
        print("no database searching file found, proceeding without comparison!")
        print("")
        

    """
    Prepare folders and paths
    """
    if not len(Output_directory): toggle_usefilename=True 
    if toggle_usefilename: Output_directory=str(Path(Path(de_novo_file).parents[0],"Output_"+Path(de_novo_file).name))
    if not os.path.exists(Output_directory): os.mkdir(Output_directory)
    if not os.path.exists(Temporary_directory): os.mkdir(Temporary_directory)

    """
    1/5 DMD alignment
    """
    print("Step 1 of 5: Start diamond alignment")
    fasta_file,denovo_peptides=write_to_fasta(de_novo_file=de_novo_file,Output_directory=Output_directory,base_ALC_Score_cut=40,add_decoy=True)
    dmd_file=str(Path(Output_directory,"diamond_alignments"+"\\de novo peptides.tsv"))
    if os.path.exists(dmd_file) == True: 
        print("Diamond alignment exists")
        alignments=dmd_file
    else:
        alignments=align(fasta_file,Output_directory,diamond_path=diamond_path,database_path=diamond_database_path,output_columns=output_columns)
    target_decoys=Process_diamond_alignment(alignments,Output_directory=Output_directory,minimum_bitscore=bit,min_ALC_score=min_ALC_score,output_columns=output_columns)
    
    if not len(target_decoys):
        print("Quitting, no alignments found that pass the threshold")
        continue


    """
    2/5 LCAs for alignments
    """
    print("Step 2 of 5: Construct LCAs")
    for target_decoy in target_decoys: 

        if "lca" in lcas:   # 1. Conventional lca
            denovo_peptides_lca=lca(target_decoy,Output_directory=Output_directory,
                                    denovo_peptides=denovo_peptides,
                                    method="conventional",
                                    filter_cutoff=freq_cut,
                                    minimum_rank=DB_rank,
                                    write_database=Write_to_database)
            
            
            
        if "bitlca" in lcas: # 2. Bitsore lca 
        
            denovo_peptides_blca=lca(target_decoy,Output_directory=Output_directory,
                                     denovo_peptides=denovo_peptides,
                                     method="bitscore",weight_column="bitscore",
                                     filter_cutoff=freq_cut,
                                     minimum_rank=DB_rank,
                                     write_database=Write_to_database)
        
        if "wlca" in lcas: # # 3. Weighted lca
        
            denovo_peptides_wlca=lca(target_decoy,Output_directory=Output_directory,
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

    if len(database_searching_file):

        """
        4/5 Check quality of fragmentation spectra
        """  
        print("Step 4 of 5: Check spectral quality")
        Plot_high_scoring(denovo_peptides,target_decoys[0],database_searching_file,Output_directory,de_novo_file)
        
        if len(fasta_database_path):
    
            """
            5/5 Check coverage obtained by database searching
            """
            print("Step 5 of 5: Check DB searching output")
            
            if   "wlca" in lcas: Compare_Bar(denovo_peptides_wlca,database_searching_file,fasta_database_path,taxa,Output_directory,ncbi_taxdf)
            elif "blca" in lcas: Compare_Bar(denovo_peptides_blca,database_searching_file,fasta_database_path,taxa,Output_directory,ncbi_taxdf)
            else:                Compare_Bar(denovo_peptides_lca,database_searching_file,fasta_database_path,taxa,Output_directory,ncbi_taxdf)
            
        else:
            print("No fasta database found, skipping step 5!")

    else:
        print("No database searching file found, skipping steps 4 & 5 !")

# =============================================================================
# DONE
# =============================================================================

end=time.time()
print("NovoLign completed: Run time",round((end-start)/60,2),"minutes")
globals().clear()
