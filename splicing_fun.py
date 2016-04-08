#/usr/bin/python

# splicing_fun.py
# Caleb Matthew Radens
# 2016_4_7

### A bunch of majiq-related functions

import os
import sys

DPSI_HEADER = ["Gene Name",
               "Gene ID",
               "LSV ID",
               "(dPSI) per LSV junction",
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
        
        print str(line_i)+ " LSVs extracted from "+os.path.basename(File_path)
        return lsv_dictionary

def top_lsvs(Cutoff = 0.95):
    """
        Given an imported dpsi file, return names of LSVs that contain
        at least one junction with P(E(dPSI)) >= Cutoff.
    """
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            