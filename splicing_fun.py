#/usr/bin/python

# splicing_fun.py
# Caleb Matthew Radens
# 2016_4_7

### A bunch of majiq-related functions

import os
import sys
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt

DPSI_HEADER = ["Gene Name",
               "Gene ID",
               "LSV ID",
               "E(dPSI) per LSV junction",
               "P(|E(dPSI)|>=0.20) per LSV junction",
               "E(PSI)",                                # This header also has the 1st condition name
               "E(PSI)",                                # This header also has the 2nd condition name
               "LSV Type",
               "A5SS",
               "A3SS",
               "ES",
               "Num. Junctions",
               "Num. Exons",
               "De Novo Junctions?",
               "chr",
               "strand",
               "Junctions coords",
               "Exons coords",
               "Exons Alternative Start",
               "Exons Alternative End",
               "IR coords",
               "Voila link"]

def import_dpsi(File_path):
    """
        Given a file path pointing at voila .txt output,
        return a dictionary with LSV_ID->data about that LSV.
    """
    
    lsv_dictionary = dict()
    
    with open(File_path, "rb") as handle:
        line_i = 0
        for line in handle:
            line_split = line.rstrip("\r\n").split("\t")
            
            if line_i == 0:
                condition_1_name = line_split[5].split(" ")[0]
                condition_2_name = line_split[6].split(" ")[0]
                line_i += 1
                continue

            gene_name = str(line_split[0])
            
            gene_id = str(line_split[1])
            
            dPSIs_floated = map(float,str(line_split[3]).split(";"))
            
            prob_dPSIs_floated = map(float,str(line_split[4]).split(";"))
            
            ePSI_1_floated = map(float,str(line_split[5]).split(";"))
            
            ePSI_2_floated = map(float,str(line_split[6]).split(";"))
            
            lsv_type = str(line_split[7])
                      
            A5SS = bool(line_split[8])
            
            A3SS = bool(line_split[9])
            
            ES = bool(line_split[10])
            
            n_junctions = int(line_split[11])
            
            n_exons = int(line_split[12])
            
            de_novo_junct = int(line_split[13])
            
            chr = str(line_split[14])
            
            strand = str(line_split[15])
            
            junct_coord = str(line_split[16])
            
            exon_coord = str(line_split[17]).split(";")
            
            exon_alt_start = str(line_split[18])
            
            exon_alt_end = str(line_split[19])
            
            ir_coords = str(line_split[20])
            
            voila_link = str(line_split[21])
            
            LSV_ID = str(line_split[2])
            lsv_dictionary[LSV_ID] = dict({
                                                DPSI_HEADER[0]:gene_name,
                                                DPSI_HEADER[1]:gene_id,
                                                DPSI_HEADER[3]:dPSIs_floated,
                                                DPSI_HEADER[4]:prob_dPSIs_floated,
                                                DPSI_HEADER[5]:ePSI_1_floated,
                                                DPSI_HEADER[6]:ePSI_2_floated,
                                                DPSI_HEADER[7]:lsv_type,
                                                DPSI_HEADER[8]:A5SS,
                                                DPSI_HEADER[9]:A3SS,
                                                DPSI_HEADER[10]:ES,
                                                DPSI_HEADER[11]:n_junctions,
                                                DPSI_HEADER[12]:n_exons,
                                                DPSI_HEADER[13]:de_novo_junct,
                                                DPSI_HEADER[14]:chr,
                                                DPSI_HEADER[15]:strand,
                                                DPSI_HEADER[16]:junct_coord,
                                                DPSI_HEADER[17]:exon_coord,
                                                DPSI_HEADER[18]:exon_alt_start,
                                                DPSI_HEADER[19]:exon_alt_end,
                                                DPSI_HEADER[20]:ir_coords,
                                                DPSI_HEADER[21]:voila_link})
            
            line_i += 1
        
        lsv_dictionary["condition_1_name"] = condition_1_name
        lsv_dictionary["condition_2_name"] = condition_2_name
        
        print str(line_i)+ " LSVs extracted from "+os.path.basename(File_path)
        return lsv_dictionary
 
def get_name_set(LSV_dict, Cutoff = 0):
    """
        Given LSV dictionary, return set of unique LSV IDs over cutoff
        
        Arguments:
            LSV_dict: output of import_dpsi
            Cutoff:   Only return LSV IDs with at least 1 junction dPSI >= Cutoff
            
        Return:
            Set
    """
    
    # names AKA LSV IDs
    names = LSV_dict.keys()
    names_over_cutoff = set()
    
    for name in names:
        if name == "condition_1_name" or name == "condition_2_name":
            continue
        dPSIs = LSV_dict[name]["E(dPSI) per LSV junction"]
        for dPSI in dPSIs:
            if dPSI > Cutoff:
                names_over_cutoff.add(name)
    
    return names_over_cutoff

def no_intron_retention(LSV_dict):
    """
        Given LSV dictionary, return LSV dictionary without intron-retention LSVs.
    """
    # names AKA LSV IDs
    names = LSV_dict.keys()
    non_intron_names = list()
    intron_names = list()
    
    for name in names:
        if name == "condition_1_name" or name == "condition_2_name":
            non_intron_names.append(name)
            intron_names.append(name)
            continue
        # If LSV type is intron:
        if LSV_dict[name]["LSV Type"][-1:] == "i":
            # Save these names, too, just in case I want em later
            intron_names.append(name)
        else:
            non_intron_names.append(name)
    
    # Copy subset of dictionary using found names.
    #     The None is not nec., because I know all keys will be in this
    #     Dict, but for future Caleb/aliens modifying this code, I'm keeping it.
    new_dict = {k: LSV_dict.get(k, None) for k in non_intron_names}
    
    print str(len(non_intron_names))+ " out of " + str(len(non_intron_names)+len(intron_names)-3) + " non-intronic LSVs found."
    
    return new_dict

def plot_venn(List_of_sets,
              Set_labels,
              Main = "I forgot to give this plot a name.",
              Out_File = ""):
    """
        Given a list of sets, generate a venn diagram in Out_Dir.
        
        Arguments:
            List_of_sets (two or three only!)
            Set_labels: Label for each circle
            Main: Title of plot
            Out_File: Where should plot be saved? And what should the file be named?
                Parent directory expected to already exist...
                This will overwrite plots if they already exist
    """
    if not os.path.isdir(os.path.dirname(Out_File)):
        raise ValueError(os.path.dirname(Out_File)+" <--- PATH DOESN'T EXIST")
    if len(List_of_sets) == 2:
        if len(Set_labels) != 2:
            raise ValueError("Set_labels needs to be the same length as the number of sets...")
        # Default figure dimensions...
        plt.figure()
        venn2(List_of_sets,Set_labels)
        plt.title(Main)
        plt.savefig(Out_File)
        
    elif len(List_of_sets) == 3:
        if len(Set_labels) != 3:
            raise ValueError("Set_labels needs to be the same length as the number of sets...")
        # Default figure dimensions...
        plt.figure()
        venn3(List_of_sets,Set_labels)
        plt.title(Main)
        plt.savefig(Out_File)
    else:
        raise ValueError("List_of_sets needs to be of length 2 or 3.")
        
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            