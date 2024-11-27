# -*- coding: utf-8 -*-
"""
Created on Tue Jan 10 12:20:18 2023

@author: hugokleikamp
"""

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
basedir=str(Path(os.getcwd()))
print(os.getcwd())

#%% import 

import shutil
import Bio
from Bio import SeqIO
import pandas as pd
import itertools
import random


#%%
def merge_files(files,path=False,delete_old=False):
    
    if type(files)==str:
        files=[files]
        
    if not path:
        ppath=Path(files[0]).parents[0]
        outpath=str(Path(ppath,Path(ppath.name+".fa")))
        
    if len(files)>1:
    
        with open(outpath,'wb') as wfd:
            for f in files:
                print(f)
                with open(f,'rb') as fd:
                    shutil.copyfileobj(fd, wfd)
                    
                if delete_old:
                    os.remove(f)
        
        return outpath
    
    else:
        return files[0]


def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g

def prep_database(
                  Path_to_db,
                  Equate_IL=True,                                 # Equates I and J to L
                  Remove_ambiguous=["J","B","X","Z","[","(",">"], # removes uncertain aas: (could also include "O","U"), does not remove J in case of Equate_IL
                  No_Fragments=False,                             # Remove fragmented sequences (contain "Fragment" in header)
                  
                  Prokaryote_only=False,                          # Only retain prokaryotes in database: requires taxonomy file
                  No_Dump=False,                                  # Remove dump taxa from database: requires taxonomy file
                  Path_to_taxonomy="",                            # path to parsed NCBI taxonomy (as created by https://github.com/hbckleikamp/NCBI2Lineage)
                  
                  Add_taxid=False,                                # add taxonomy id to accession, useful since database searching tools often only annotate accession
                  Taxid_delimiter="OX=",                          # left delimiter of taxid in header (right delimiter =space), used to (Uniref would have Taxid=)  
                                                                  # is used to parse taxonomy in Prokaryotes_only, No_Dump and Add_taxid
                  Add_decoy=False,                                # options: False (no decoy), "reverse" (reversed) or "scramble" (fully randomized)
                  decoy_delimiter="decoy_",                       # gets appended to front of headers of decoy sequences            
                  
                  Output_path=False,                              # full path to database output
                  delete_old=False,
                  
                  ):

    if Prokaryote_only or No_Dump:
        ranks=["superkingdom","phylum","class","order","family","genus","species"] 
        taxdf=pd.read_csv(Path_to_taxonomy,sep="\t")
        taxdf.columns=['OX']+ranks+["OS"]
        
        if Prokaryote_only:                           taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)
        if No_Dump and "Dump_taxid" in taxdf.columns: taxdf=taxdf[taxdf["Dump_taxid"].astype(str)=="False"] 
    
    #parse output_path
    if not Output_path:
        Output_path=Path_to_db
        if Prokaryote_only:   Output_path=Output_path.replace(".fasta","_BacArch.fasta")
        if Remove_ambiguous:  Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
        if No_Dump:           Output_path=Output_path.replace(".fasta","_NoDump.fasta")
        if No_Fragments:      Output_path=Output_path.replace(".fasta","_NoFrag.fasta")
        if Equate_IL:         Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
        if Add_decoy:         Output_path=Output_path.replace(".fasta","_Decoy.fasta")
        if Add_taxid:         Output_path=Output_path.replace(".fasta","_taxid.fasta")

    #read    
    recs=SeqIO.parse(Path_to_db,format="fasta")
    chunks=chunk_gen(recs)
    
    #write 
    once=True
    print("writing "+Path(Output_path).stem)
    with open(Output_path,"w") as f: #clears file if exists
        pass
    with open(Output_path,"a") as f:
    
        for ic,c in enumerate(chunks):
            print("chunk "+ str(ic))
    
            chunk_df=pd.DataFrame([[str(r.seq),r.description] for r in c],columns=["seq","description"])

            #attempt to identify Taxid delimiter
            if once:
                test=chunk_df.head(10)["description"]
                if not test.str.contains(Taxid_delimiter).any():
                    if   test.str.contains("OX="   ).any(): Taxid_delimiter="OX=" 
                    elif test.str.contains("Taxid=").any(): Taxid_delimiter="Taxid=" 
                once=False


            if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
            if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Remove_ambiguous],axis=1).any(axis=1)]
            if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)]

            if Prokaryote_only or No_Dump: chunk_df=chunk_df[chunk_df.description.str.split(Taxid_delimiter).apply(lambda x:x[-1]).str.split(" ").apply(lambda x: x[0]).isin(taxdf["OX"])]
           
            if Add_taxid:        chunk_df["id"]=chunk_df["id"]+"|"+chunk_df.description.str.split(Taxid_delimiter).apply(lambda x: x[-1]).str.split(" ").apply(lambda x:x[0])
            if Add_decoy:
                decoy=chunk_df.copy()
                if Add_decoy=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
                if Add_decoy=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
                decoy["description"]=decoy_delimiter+decoy["description"]
                chunk_df=pd.concat([chunk_df,decoy])
                
            f.write("\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"]))

    if delete_old:
        os.remove(Path_to_db)
    return Output_path
