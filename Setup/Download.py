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
script_dir=os.getcwd()
print(os.getcwd())
basedir=os.getcwd()

#%% import 
import urllib, gzip, zipfile, shutil, tarfile, time, requests
from sys import platform

def download(url,filename):

    #download
    if not os.path.exists(filename): #check if file is already downloaded
        while True:
            try:
                urllib.request.urlretrieve(url,filename)
                break
            except:
                time.sleep(2)
                print("retry")
        return filename
        

def extract(path): #path to folder

    #recursive extraction            
    while any([f.endswith((".zip",".gz",".tar")) for f in os.listdir(path)]):
        
        for f in os.listdir(path):
            
            i=str(Path(path,Path(f)))
            o=str(Path(path,Path(f).stem))
            
            if f[0].isalnum() and f.endswith(".zip"):
                print("extracting "+f)
                with zipfile.ZipFile(i, 'r') as zip_ref:
                    zip_ref.extractall(path)
                if os.path.exists(i): os.remove(i)

            if f[0].isalnum() and f.endswith(".gz"):
                print("extracting "+f)
                with gzip.open(i,'rb') as f_in:
                    with open(o,'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                if os.path.exists(i): os.remove(i)
                
            if f[0].isalnum() and f.endswith(".tar"):
                print("extracting "+f)
                tar = tarfile.open(i, "r:")
                tar.extractall(path)
                tar.close()
                if os.path.exists(i): os.remove(i)
                    
        
def extract_subfolders(folder):

    for subdir, dirs, files in os.walk(folder):
        for s in subdir():

            extract(s)
   
    
def download_extract(urls,path,
                     fasta_exts=[".faa",".fasta",".fa"] #accepted file extenstions for fasta files that are returned as output
                     ):

    if type(urls)==str:
        urls=[urls]
    
    for url in urls:
        print(url)
        
        filename = str(Path(path,url.split("/")[-1]))
        if not os.path.exists(path): os.mkdir(path)
        
        #Add prompt if exists!
        if os.path.exists(filename):
            w = input("File "+filename+" already exists, delete? (y/n): ")
            if w == "y":
                if os.path.isfile(filename):
                    os.remove(filename)
                elif os.path.isdir(filename):
                    shutil.rmtree(filename)
            if w == "n":
                continue
        
        while True:
            try:
                file=download(url,filename) #download
                extract(path)               #extract folder contents
                break
                
            except:
                #in case of incomplete download, decompression can fail, so retry
                os.remove(str(Path(path,filename)))

        for root, dirs, files in os.walk(path):
            for d in dirs:
                extract(str(Path(root,d)))

        #return all fasta files
        fastas=[]
        for root, dirs, files in os.walk(path):
            for file in files:
                for ext in fasta_exts:
                    if file.endswith(ext):
                        fastas.append(str(Path(root,file)))
 
    return list(set(fastas))
        
#%%

#Wrapper for downloading db
def download_db(DB,path=False):

    #RefSeq
    if DB=="RefSeq": #Refseq select
        urls=requests.get("https://ftp.ncbi.nlm.nih.gov/blast/db/refseq_protein-prot-metadata.json").json().get("files")
        if not path: path=str(Path(basedir,"RefSeq"))
    
    #NCBI NR
    if DB=="NCBI_NR":
        urls=requests.get("https://ftp.ncbi.nlm.nih.gov/blast/db/nr-prot-metadata.json").json().get("files")
        if not path: path=str(Path(basedir,"NCBI_NR"))
    
    #Swiss-Prot
    if DB=="Swiss-Prot":
        urls="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
        if not path: path=str(Path(basedir,"Swiss-Prot"))
    
    #TrEMBL
    if DB=="TrEMBL":
        urls="https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
        if not path: path=str(Path(basedir,"TrEMBL"))
    
    #UniProt (Swiss-Prot+TrEBML)
    if DB=="UniProt":
        urls=["https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz",
              "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"]
        if not path: path=str(Path(basedir,"UniProt"))
    
    #Uniref 100
    if DB=="UniRef100":
        urls="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz"
        if not path: path=str(Path(basedir,"UniRef100"))
    
    #Uniref 90
    if DB=="UniRef90":
        urls="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
        if not path: path=str(Path(basedir,"UniRef90"))
    
    #Uniref 50
    if DB=="UniRef50":
        urls="https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz"
        if not path: path=str(Path(basedir,"UniRef50"))
    
    fasta_files=download_extract(urls,path)
    
    
    return fasta_files

def download_ncbi_taxdump(path=False):

    urls="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"    
    if not path: path=str(Path(basedir,"taxdump"))
    download_extract(urls,path)
    
    #cleanup
    for i in os.listdir(path):
        if i=="names.dmp":
            names=str(Path(path,i))
        elif i=="nodes.dmp":
            nodes=str(Path(path,i))

    return names,nodes


def download_diamond(path=False):


    if platform == "linux" or platform == "linux2":
        # linux
        urls="https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-linux64.tar.gz"
    elif platform == "win32":
        # Windows
        urls="https://github.com/bbuchfink/diamond/releases/download/v2.1.8/diamond-windows.zip"    
    else:
        print("this pipeline only supports widows or linux systems, since DIAMOND installation on OSX requires Miniconda")
        return
    
    if not path: path=str(Path(basedir,"diamond"))
    download_extract(urls,path)
    
    for i in os.listdir(path):
        if i.startswith("diamond"):
            return str(Path(path,i))
    
