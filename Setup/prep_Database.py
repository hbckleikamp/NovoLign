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
import numpy as np
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
                  
                  Path_to_taxonomy="",                            # path to parsed NCBI taxonomy
                  Path_to_synonyms="",                            # path to ncbi taxdump names.dmp file (used for renameing OS-style headers to OX-style)
                  
                  Add_taxid=False,                                # add taxonomy id to accession, useful since database searching tools often only annotate accession
                  Taxid_delimiter_left="OX=",                     # left delimiter of taxid in header (right delimiter =space), used to (Uniref would have Taxid=)  
                  Taxid_delimiter_right=" ",
                  db=False,                                       #predefined delimiters
                  #Database choices: "RefSeq", "NCBI_NR", "Swiss-Prot", "TrEMBL", "UniProt", "UniRef100", "UniRef90", "UniRef50"
                  
                  Add_decoy=False,                                # options: False (no decoy), "reverse" (reversed) or "scramble" (fully randomized)
                  decoy_delimiter="decoy_",                       # gets appended to front of headers of decoy sequences            
                  
                  Output_path=False,                              # full path to database output
                  delete_old=False,
                  ranks=["superkingdom","phylum","class","order","family","genus","species"]
   
                  ):


    #pre-defined taxa delimiters
    if db in ["NCBI_NR","RefSeq"]:                  Taxid_delimiter_left,Taxid_delimiter_right=" [","] " #NCBI-style header
    if db in ["Swiss-Prot", "TrEMBL", "UniProt"]:   Taxid_delimiter_left,Taxid_delimiter_right="OX="," " #UniprotKB-style header
    if db in ["UniRef100", "UniRef90", "UniRef50"]: Taxid_delimiter_left,Taxid_delimiter_right="Taxid="," " #UniRef-style header
    suf=Path(Path_to_db).suffix
    
    #read taxonomy 
    if Prokaryote_only or No_Dump:
        taxdf=pd.read_csv(Path_to_taxonomy,sep="\t",index_col=[0])
    
        
        if Prokaryote_only:                           taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)
        if No_Dump and "Dump_taxid" in taxdf.columns: taxdf=taxdf[taxdf["Dump_taxid"].astype(str)=="False"] 
    
    #parse output_path
    if not Output_path:
        
        if len(suf): Output_path=Path_to_db.replace(suf,".fasta")
        else:        Output_path=Path_to_db+".fasta"
        
        if Prokaryote_only:   Output_path=Output_path.replace(".fasta","_BacArch.fasta")
        if Remove_ambiguous:  Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
        if No_Dump:           Output_path=Output_path.replace(".fasta","_NoDump.fasta")
        if No_Fragments:      Output_path=Output_path.replace(".fasta","_NoFrag.fasta")
        if Equate_IL:         Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
        if Add_decoy:         Output_path=Output_path.replace(".fasta","_Decoy.fasta")
    
    #read    
    recs=SeqIO.parse(Path_to_db,format="fasta")
    chunks=chunk_gen(recs)
    
    #write 
    once=True
    rename=True
    print("writing "+Path(Output_path).stem)
    with open(Output_path,"w") as f: #clears file if exists
        pass
    with open(Output_path,"a") as f:
       
        for ic,c in enumerate(chunks):
            print("chunk "+ str(ic))
       
            chunk_df=pd.DataFrame([[str(r.seq),r.description] for r in c],columns=["seq","description"])
            chunk_df.description=chunk_df.description.str.replace("\x01"," ")
            
            ### dynamic taxid delmiter detection
            if once:
                test=chunk_df.head(10)["description"]
                if not test.str.contains(Taxid_delimiter_left).any(): 
                    if   test.str.contains("OX="   ).any(): Taxid_delimiter_left,Taxid_delimiter_right ="OX="," " 
                    elif test.str.contains("Taxid=").any(): Taxid_delimiter_left,Taxid_delimiter_right ="Taxid="," "
                    else:                                   Taxid_delimiter_left,Taxid_delimiter_right =" [","] "
                once=False 
       
                #test if taxonomies are OS or OX? 
                taxids=test.str.split(Taxid_delimiter_left,regex=False).apply(lambda x: x[-1]).str.split(Taxid_delimiter_right,regex=False).apply(lambda x: x[0])
                if not taxids.str.isnumeric().all():
                    taxdf=pd.read_csv(Path_to_taxonomy,sep="\t",index_col=[0])
                    
                    with open(Path_to_synonyms,"r") as syns: lines=pd.DataFrame([l.split("\t")[0:7] for l in syns.readlines()])
                    namesdf=lines.iloc[:,[0,2,4,6]]
                    namesdf.columns=["taxid","name","x","type"]
                    names2tax=namesdf[["name","taxid"]].drop_duplicates().set_index("name").astype(np.int64)
                    rename=True #set rename flag
       
            #sequence related database editing
            if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
            if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Remove_ambiguous],axis=1).any(axis=1)].reset_index(drop=True)
            if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)].reset_index(drop=True)
       
            #taxonomy related database editing
            
            if Prokaryote_only or No_Dump or rename:
                taxids=(chunk_df.description+" ").str.split(Taxid_delimiter_left,regex=False)
                taxids=taxids.apply(lambda x: x[1::2]).explode().dropna().str.split(Taxid_delimiter_right,regex=False).apply(lambda x: x[0])
                
                
                if rename: #turn OS into OX
                    q=taxids.isin(names2tax.index)
                    taxids[~q]="root" #placeholder (or taxids=taxids[q])
                    
                    t=taxids.reset_index().merge(names2tax,left_on="description",right_on="name",how="left").drop_duplicates()
                    lins=taxdf.set_index("OX").loc[t.taxid].reset_index()
                    lins.index=t["index"]
                    lins["OX"]=lins["OX"].astype(str)
                    lca=lins.groupby(lins.index)[ranks+["OX"]].nth(0)[(lins.groupby(lins.index)[ranks+["OX"]].nunique()==1)]
            
                    taxa=[]
                    counter=0
                    while lca.OX.isnull().sum():
                        taxa.append(lca.OX.dropna())
                        
                        t=lca.loc[lca["OX"].isnull(),ranks].ffill(axis=1).species.fillna("root")
                        t.name="OS"
                        lins=taxdf.merge(t.reset_index(),how="right").set_index("index")
                        lins["OX"]=lins["OX"].astype(str)
                        lca=lins.groupby(lins.index)[ranks+["OX"]].nth(0)[(lins.groupby(lins.index)[ranks+["OX"]].nunique()==1)]
                        
                        counter+=1
                        if counter>5:
                            break
                
                
                    taxa=pd.concat(taxa).reset_index().astype(int)
                    taxids=np.ones(len(chunk_df),dtype=int)
                    taxids[taxa["index"]]=taxa["OX"]
                    
                    schunk=chunk_df.description.str.split(" ",n=1) 
                    chunk_df["description"]=schunk.apply(lambda x: x[0])+" OX="+taxids.astype(str)+" "+schunk.apply(lambda x: x[1]) #add to description
                    
                 
       
                if Prokaryote_only or No_Dump: chunk_df=chunk_df[taxids.isin(taxdf["OX"])].reset_index(drop=True)
       
       
       
            if Add_decoy:
                decoy=chunk_df.copy()
                if Add_decoy=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
                if Add_decoy=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
                decoy["description"]=decoy_delimiter+decoy["description"]
                chunk_df=pd.concat([chunk_df,decoy])
            
      
            
            f.write("\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"]+"\n"))
       



    if delete_old:
        os.remove(Path_to_db)
    return Output_path
