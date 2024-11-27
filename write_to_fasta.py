"""
module for converting PEAKS and DeepNovo outputs into target and decoy peak lists
Author: HK / mpabst
Date: 26-06-2022
"""
from pathlib import Path
import pandas as pd
import numpy as np
import re, random
from collections import Counter
import os
import warnings

# add decoy peptides        
def scramble(x):
    return ''.join(random.sample(x, len(x)))
def scramble_C(x):
    return ''.join(random.sample(x[:len(x)-1], len(x)-1)+[x[-1]])
def reverse(x):
    return x[::-1]
def reverse_C(x):
    return x[-2::-1]+x[-1]

# simple conversion
def con_csv(peaksfiles):

    # convert files into target/decoy format
    for p in peaksfiles:
        df=pd.read_csv(p)
        df['Peptide']=df['Peptide'].apply(lambda x: re.sub("[\(\[].*?[\)\]]", "", x).replace(",","")) #remove ptms in peptides
        heads=">"+df["Peptide"]+"_"+df["ALC (%)"].astype(str)+"_"+df["ppm"].astype(str)
        with open("target_"+Path(p).stem+".fa","w") as f: # create target
            f.write("\n".join([heads[ix]+"\n"+peptide for ix,peptide in enumerate(df["Peptide"])])+"\n")

        with open("decoy_"+Path(p).stem+".fa","w") as f: # create decoy
            f.write("\n".join([heads[ix]+"\n"+peptide[::-1]  for ix,peptide in enumerate(df["Peptide"])])+"\n")
    return(df)


#read table with dynamic delmiter detection
def read_table(tabfile, *,
               Keyword="Peptide",
               ):

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        try:
            tab = pd.read_excel(tabfile, engine='openpyxl')
        except:
            with open(tabfile, "r") as f:
                tab = pd.DataFrame(f.read().splitlines())

        # dynamic delimiter detection: if file delimiter is different, split using different delimiters until the desired column name is found
        if Keyword:
            if Keyword not in tab.columns:
                delims = [i[0] for i in Counter(
                    [i for i in str(tab.iloc[0]) if not i.isalnum()]).most_common()]
                for delim in delims:
                    if delim == " ":
                        delim = "\s"
                    try:
                        tab = pd.read_csv(tabfile, sep=delim)
                        if Keyword in tab.columns:
                            return tab
                    except:
                        pass

    return tab


# universl peaks and deepnovo conversion
def write_to_fasta(de_novo_file,
                   Output_directory,          # str: one  novo sequencing output file (full filepath), in .csv or .tsv format
                   database_searching_file="",
                   base_ALC_Score_cut=40,             # numeric: mininum required score cutoff (PEAKS score, ALC(%))
                   Deepnovo_Score_cutoff=-0.1,         # numeric: mininum required score cutoff (DeepNovo score)
                   ppm_cutoff=None,                    # numeric: maximum ppm tolerance (absolute)
                   length_cutoff=None,                 # numeric: mininum required peptide length
                   Area_cutoff=None,                   # numeric: mininum required peak area
                   Intensity_cutoff=None,              # bool: mininum required intensity
                   equate_IL=True,                     # bool: change I and J to L
                   matrix_path=None,                   # str: location of scoring matrix, used to remove petpides with amino acids that are not in matrix
                   add_decoy=True,                     # wether decoy sequences should be included. (Necessary for some filtering steps and correcting weights, in corrected weighted LCA)
                   output_folder="diamond_fasta",
                   decoy_mode="scramble",              # scramble, reverse or scramble_C, reverse_C     
                   header_info=["Scan","alignment_Target_Decoy",
                                "ALC (%)","predicted_score"] # which data should be retained in the fasta output [Scan and Target_Decoy are minimum requirements for the script to run]
                   ): 

 
    df=read_table(de_novo_file)
    
    # convert Deepnovo output to PEAKS output format    
    if "predicted_sequence" in df.columns: df=format_DeepNovo(df)
    
    # clean sequences
    df['Peptide']=df['Peptide'].apply(lambda x: re.sub("[\(\[].*?[\)\]]", "", x).replace(",","")) #remove ptms in peptides, everything in between ()
    df['Peptide']=df['Peptide'].apply(lambda x: re.sub("[\[\[].*?[\]\]]", "", x))                 #remove ptms in peptides, everything in between []
    
    # center ppm around 0, then make absolute
    if 'ppm' in df.columns:
        hist, bin_edges=np.histogram(df['ppm'].tolist(),bins=100)
        df['ppm']=abs(df['ppm']-(bin_edges[np.where(hist==np.max(hist))[0]]
                                       +bin_edges[np.where(hist==np.max(hist))[0]+1])/2).round(1)
    
    if 'Tag Length' in df.columns: df["Tag Length"]=df["Peptide"].apply(len)   # add tag length
    for i in ['Tag Length','ALC (%)','predicted_score','ppm','Area','Intensity']:
                if i in df.columns:
                    df[i]=df[i].astype(float) 
    
    if (type(base_ALC_Score_cut   )!=type(None)) & ('ALC (%)'         in df.columns): df=df[df['ALC (%)'         ]>=base_ALC_Score_cut]   # filter on scoring Peaks 
    if (type(Deepnovo_Score_cutoff)!=type(None)) & ('predicted_score' in df.columns): df=df[df['predicted_score' ]<=Deepnovo_Score_cutoff] # filter on scoring Deepnovo 
    if (type(ppm_cutoff           )!=type(None)) & ('ppm'             in df.columns): df=df[df['ppm'             ]<=ppm_cutoff]            # filter on ppm
    if (type(length_cutoff        )!=type(None)) & ('Area'            in df.columns): df=df[df['Area'            ]>=Area_cutoff]           # filter on area
    if (type(Area_cutoff          )!=type(None)) & ('Intensity'       in df.columns): df=df[df['Intensity'       ]>=Intensity_cutoff]      # filter on intensity
    if (type(Intensity_cutoff     )!=type(None)) & ('Tag Length'      in df.columns): df=df[df['Tag Length'      ]>=length_cutoff]         # filter on tag length

    if "Scan" not in df.columns: df["Scan"]=np.arange(len(df)) # add scan 
    if equate_IL: df["Peptide"]=df["Peptide"].str.replace("I","L").str.replace("J","L") # equate I and J to L

    # remove any peptide with characters not present in the scoring matrix
    if type(matrix_path)!=type(None):
        m_aas=set((pd.read_csv(matrix_path).iloc[0].values[0]).split())
        u_aas=set("".join(df["Peptide"].values.tolist())) 
        d=[df["Peptide"].str.contains(i) for i in list(u_aas-m_aas)]
        if d:
            df=df[~df[pd.concat(d,axis=1).any(axis=1)]]

    #df=df.sample(500).reset_index() # optional subsampling for method development purpose
    peptide_df=df

    if add_decoy: # add decoy alignment peptides
        df.loc[:,"alignment_Target_Decoy"]="Target"
        decoy=df.copy()
        decoy.loc[:,"alignment_Target_Decoy"]="Decoy"

        if decoy_mode=="scramble_C": decoy["Peptide"]=decoy["Peptide"].apply(scramble_C)
        if decoy_mode=="scramble":   decoy["Peptide"]=decoy["Peptide"].apply(scramble)
        if decoy_mode=="reverse_C":  decoy["Peptide"]=decoy["Peptide"].apply(reverse_C)
        if decoy_mode=="reverse":    decoy["Peptide"]=decoy["Peptide"].apply(reverse)

        df=pd.concat([df,decoy]).astype(str).reset_index()

    # prepare fasta output [Scan and Target_Decoy are minimum requirements for the script to run]
    header_info=[i for i in header_info if i in df.columns]
    df=df[["Peptide"]+header_info].drop_duplicates().reset_index(drop=True)
    hdict=df[header_info].T.to_dict()
    shdict=[str(hdict[ix]).replace(" ","") for ix in range(len(df))] #spaces are removed from header info, otherwise eval wont work
    heads=">"+df["Peptide"]+";"+shdict
    
    # write complete DN output
    if not os.path.exists(str(Path(Output_directory,output_folder))): os.mkdir(str(Path(Output_directory,output_folder)))
    out_path_all=str(Path(Output_directory,output_folder,Path(de_novo_file).stem+".fa"))
    with open(out_path_all,"w") as f:
        f.write("\n".join([heads[ix]+"\n"+peptide for ix,peptide in enumerate(df["Peptide"])])+"\n")

    if len(database_searching_file):
        db_df=pd.read_csv(database_searching_file)
        
        #de novo only scans, de novo only precursors
        peptide_df["DeNovo_only_Scan"]=     ~peptide_df["Scan"].isin(db_df["Scan"])
        peptide_df["DeNovo_only_Precursor"]=~peptide_df[["Scan","Mass"]].astype(str).sum(axis=1).isin(db_df[["Scan","Mass"]].astype(str).sum(axis=1))
        peptide_df["DeNovo_only_PSM"]=      ~peptide_df[["Scan","Mass","Peptide"]].astype(str).sum(axis=1).isin(db_df[["Scan","Mass","Peptide"]].astype(str).sum(axis=1))
        
    return out_path_all,peptide_df


def format_DeepNovo(df):

    df["Peptide"]=df["predicted_sequence"]         
    df=df[df["Peptide"].notnull()]

    # calculate ppm shift from m/z
    mass=np.zeros((1,len(df)))
    # currently only carbamidomethylation and oxidation are implemented
    mass+=df["Peptide"].str.count("(Carbamidomethylation)").fillna(0)*57.021463 #add possible modifications
    mass+=df["Peptide"].str.count("(Oxidation)").fillna(0)*15.994915

    df['Peptide']=df['Peptide'].apply(lambda x: re.sub("[\(\[].*?[\)\]]", "", x).replace(",","")) #remove ptms in peptides

    # calculate peptide mass foor DeepNovo  
    df['calculated_mass']=mass[0]+df['Peptide'].apply(lambda x: mass_calc(x)).values
    df['precursor_mass']=df['precursor_mz']*df['precursor_charge']-df['precursor_charge']*1.007277                
    df["ppm"]=(1000000/df['calculated_mass'])*(df['calculated_mass']-df['precursor_mass'])

    # rename columns
    if "feature_id" in df.columns: df=df.rename(columns={"feature_id":"Scan"})  
    if "feature_area" in df.columns: df=df.rename(columns={"feature_area":"Area"})  
    if "feature_intensity" in df.columns: df=df.rename(columns={"feature_intensity":"Intensity"})  

    return df

# mass calculation for DeepNovo            
std_aa_mass = {'G': 57.02146, 'A': 71.03711, 'S': 87.03203, 'P': 97.05276, 'V': 99.06841,
               'T': 101.04768,'C': 103.00919,'L': 113.08406,'I': 113.08406,'J': 113.08406,
               'N': 114.04293,'D': 115.02694,'Q': 128.05858,'K': 128.09496,'E': 129.04259,
               'M': 131.04049,'H': 137.05891,'F': 147.06841,'U': 150.95364,'R': 156.10111,
               'Y': 163.06333,'W': 186.07931,'O': 237.14773}

def mass_calc(x,std_aa_mass=std_aa_mass):
    return sum(std_aa_mass.get(aa) for aa in x)+18.01056
