#/usr/bin/python

# splicing_fun.py
# Caleb Matthew Radens
# 2016_4_7

### A bunch of majiq-related functions

import os
import sys
import itertools
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt

import pdb

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

CUTOFF_DPSI = 0.2
CUTOFF_PROB = 0.95

def quick_import(Directory, Cutoff_dPSI = CUTOFF_DPSI, Cutoff_prob = CUTOFF_PROB):
    """
        Given a directory with '*_quantify_deltapsi' files, import all dPSI
            text files, throwing away intron events, return dictionary as follows.
            Key : Value
            NAME_OF_COMPARISON : Non-intron containing imported LSV dictionary
        
        Writes a result file to the directory.
        
        Assumptions:
            Directory contains files that end in ".deltapsi_quantify_deltapsi.txt"
            
            Expected format of input files:
                ~/path/NAME_OF_COMPARISON.deltapsi_quantify_deltapsi.txt
            
    """
    if not (os.path.isdir(Directory)):
        raise ValueError(Directory+" not found.")
#     pdb.set_trace()
    
    dpsi_comparison_name = list()
    dpsi_files = list()
    
    for root, subdirs, files in os.walk(Directory):
        # ID files that contain the Pattern
        for f in files:
            if ".deltapsi_quantify_deltapsi.txt" in f:
                split_file_name = f.split(".")
                dpsi_comparison_name.append(split_file_name[0])
                dpsi_files.append(os.path.join(root,f))
    
    print "Found "+str(len(dpsi_comparison_name))+ " dPSI text files in directory..."
    
    imported_dpsi_dicts = list()
    
    for f in dpsi_files:
        imported_dpsi_dicts.append(import_dpsi(f,Cutoff_dPSI,Cutoff_prob))
        
    print "Imported "+str(len(imported_dpsi_dicts))+ " dpSI text files..."
        
    imported_dpsi_dicts_no_introns = list()
    
    for dpsi_dict in imported_dpsi_dicts:
        imported_dpsi_dicts_no_introns.append(no_intron_retention(dpsi_dict))
        
    print "Threw away intronic LSVs from "+str(len(imported_dpsi_dicts_no_introns))+ " imported dPSI files..."
    
    import_dictionary = dict()
    
    if len(dpsi_comparison_name) != len(imported_dpsi_dicts_no_introns):
        raise ValueError("Something is probably screwy with the names of the dPSI text files...")
    
    for name, imported_file in itertools.izip(dpsi_comparison_name, imported_dpsi_dicts_no_introns):
        # Minus two because two of the keys are not LSVs
        print "Successfully imported "+name+" with "+ str(len(imported_file)-2)+" LSVs!"
        import_dictionary[name]=imported_file
    
    return import_dictionary

def import_dpsi(File_path, Cutoff_dPSI = CUTOFF_DPSI, Cutoff_prob = CUTOFF_PROB):
    """
        Given a file path pointing at voila .txt output,
        return a dictionary with LSV_ID->data about that LSV.
        
        Cutoff_dPSI: Which junctions in each LSV make this dPSI cutoff (+/-)?
        Cutoff_prob: Which junctions in each LSV make the dPSI *and* the prob(dPSI) cutoffs?
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
            
            # Saves the junction # and directionality of change
            over_cutoff = list()
            
            # Just saves which junction index for those that are significant
            sig_junctions = list()
            
            # Junction index
            i = 0
            for dPSI, prob_dPSI in zip(dPSIs_floated, prob_dPSIs_floated):
                # If the dPSI is greater than 0 and over the cutoff dPSI and probability
                if ((dPSI > Cutoff_dPSI) and (prob_dPSI > Cutoff_prob)):
                    # junctionIndex_over
                    over_cutoff.append(str(i)+"_over")
                    sig_junctions.append(i)
                # Else if the dPSI is less than 0 and over the cutoff dPSI and probability
                elif ((dPSI < -Cutoff_dPSI) and (prob_dPSI > Cutoff_prob)):
                    # junctionIndex_under
                    over_cutoff.append(str(i)+"_under")
                    sig_junctions.append(i)
                i+=1
            
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
                                                "dPSI_over_cutoff":over_cutoff,
                                                "sig_junctions":sig_junctions,
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
 
def get_name_set(LSV_dict, Cutoff = 0, probability_dPSI = 0):
    """
        Given LSV dictionary, return set of unique LSV IDs over cutoff
        
        Arguments:
            LSV_dict: output of import_dpsi
            Cutoff:   Only return LSV IDs with at least 1 junction abs(dPSI) >= Cutoff
            
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
        prob_dPSIs = LSV_dict[name]["P(|E(dPSI)|>=0.20) per LSV junction"]
        
        for dPSI, prob_dPSI in zip(dPSIs, prob_dPSIs):
            if (abs(dPSI) >= Cutoff) and (prob_dPSI >= probability_dPSI):
                names_over_cutoff.add(name)
    
    return names_over_cutoff

def no_intron_retention(LSV_dict):
    """
        Given LSV dictionary, return LSV dictionary without LSVs containing 
            significantly changing intron-retention junctions.
            
        If dPSI of the intron retention junction is < 0.05, the LSV IS 
            included in the return dictionary. Using 0.05 instead of 0.2
            because, by definition, at this point the LSV is complex. 
            Thus, a dPSI of 0.2 is potentially *very* significant. 
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
        
    for intron in intron_names:
        if intron == "condition_1_name" or intron == "condition_2_name":
            continue
        dPSI_intron = LSV_dict[intron]["E(dPSI) per LSV junction"][-1]
        # IF the dPSI of the intron junction is not changing significantly, keep it
        #     0.05 was suggested by Yoseph Barash
        if abs(dPSI_intron) < 0.05:
            non_intron_names.append(intron)
    
    # Copy subset of dictionary using found names.
    #     The None is not nec., because I know all keys will be in this
    #     Dict, but for future Caleb/aliens modifying this code, I'm keeping it.
    new_dict = {k: LSV_dict.get(k, None) for k in non_intron_names}
    
    print str(len(non_intron_names))+ " non-intronic LSVs out of " + str(len(non_intron_names)+len(intron_names)-3) + " total LSVs found."
    
    return new_dict

def intersect_LSVs(LSV_dict_1, LSV_dict_2, Cutoff_dPSI = CUTOFF_DPSI, Cutoff_prob= CUTOFF_PROB):
    """
        Given two LSV dictionaries from the same majiq build, return LSV IDs that are
            In both
            Unique to the first
            Unique to the second
            
        Whereby LSVs are compared with respect to specific junctions. If both LSVs have
            at least one junction (same junction) with dPSIs significantly changing,
            then the LSV is "In both." 
            Otherwise, the LSV is unique to one or the other as long as it has a 
            significantly changing dPSI. Significantly changing means it was over both of
            the cutoffs specified by the function import_dpsi().
            
        If a LSV from both dictionaries shares even 1 junction, any other info in that LSV 
            is lost. For example:
                LSV_yflsv_1 = [0.4,0.0,-0.8,0.5] (probabilities are >0.95)
                LSV_yflsv_2 = [0.3,-0.4,0.0,-0.7] (probabilities are >0.95)
            Would be a shared LSV, b/c the first and last junctions are significantly shared.
            The 2nd+3rd junctions, despite being significant changes in their own right,
            will not count towards "Unique to the first" or "Unique to the second"...
                See function: find_disagreements() if you want to find these bad boys.
                
        If both LSVs have junctions with significant dPSI, but not associated with the same 
            junctions, then the LSV will be added to both Unique to First and Unique to Second.
            A print message will state that this happened, because I expect it to be 
            rare and thus interesting.
            
        Cutoff_dPSI: Which junctions in each LSV make this dPSI cutoff (+/-)?
        Cutoff_prob: Which junctions in each LSV make the dPSI *and* the prob(dPSI) cutoffs?
    """
    lsv_1_ids = get_name_set(LSV_dict_1, Cutoff_dPSI, Cutoff_prob)
    lsv_2_ids = get_name_set(LSV_dict_2, Cutoff_dPSI, Cutoff_prob)
    # Matched IDs in both (not necessarily significantly matched, yet: just by ID)
    shared_lsvs = list(lsv_1_ids & lsv_2_ids)
    # Significantly in both LSVs
    in_both_sig_ids = list()
    in_both_sig_junctions = list()
    
    # Which LSV IDs are unique to first dictionary
    only_1_ids = list()
    # Set function: I know these will be unique to 1
    only_1_ids.extend(lsv_1_ids.difference(lsv_2_ids))
    
    # Which junctions are unique to first dictionary
    only_1_junctions = list()
    
    # Add all the significant junctions I know are Unique to First
    for lsv_id in only_1_ids:
        significant_dPSI_1 = set(LSV_dict_1[lsv_id]["sig_junctions"])
        for dPSI in significant_dPSI_1:
            # Re-naming the junction index to something a bit more meaningful
            lsv_junction_name = lsv_id+"_"+str(dPSI)
            only_1_junctions.append(lsv_junction_name)
    
    # Which LSV IDs are unique to second dictionary
    only_2_ids = list()
    # Set function: I know these will be unique to 2
    only_2_ids.extend(lsv_2_ids.difference(lsv_1_ids))
    
    # Which junctions are unique to second dictionary
    only_2_junctions = list()
    
    # Add all the significant junctions I know are Unique to Second
    for lsv_id in only_2_ids:
        significant_dPSI_2 = set(LSV_dict_2[lsv_id]["sig_junctions"])
        for dPSI in significant_dPSI_2:
            # Re-naming the junction index to something a bit more meaningful
            lsv_junction_name = lsv_id+"_"+str(dPSI)
            only_2_junctions.append(lsv_junction_name)
    
    # Loop over all lsvs that are potentially shared
    for lsv_id in shared_lsvs:
        # See import_dPSI() for what 'sig_functions' means
        # In short, each junction index in each LSV is added to this list if
        #  and only if it is significantly changing
        significant_dPSI_1 = set(LSV_dict_1[lsv_id]["sig_junctions"])
        significant_dPSI_2 = set(LSV_dict_2[lsv_id]["sig_junctions"])
        
        in_both = significant_dPSI_1 & significant_dPSI_2
        
        # If the LSVs share significant junctions
        if len(in_both) > 0:
            # Add the LSV ID to the "In both" list
            in_both_sig_ids.append(lsv_id)
            # Also add the junctions to the in both junctions list
            for junction in list(in_both):
                # Re-naming the junction index to something a bit more meaningful
                lsv_junction_name = lsv_id+"_"+str(junction)
                in_both_sig_junctions.append(lsv_junction_name)
                
        # Else both LSVs contain significant junctions, but none are shared in the same direction
        else:
            only_1_ids.append(lsv_id)
            for dPSI in significant_dPSI_1:
                # Re-naming the junction index to something a bit more meaningful
                lsv_junction_name = lsv_id+"_"+str(dPSI)
                only_1_junctions.append(lsv_junction_name)
            only_2_ids.append(lsv_id)
            for dPSI in significant_dPSI_2:
                # Re-naming the junction index to something a bit more meaningful
                lsv_junction_name = lsv_id+"_"+str(dPSI)
                only_2_junctions.append(lsv_junction_name) 
    
    all_possible_lsvs = lsv_1_ids | lsv_2_ids     
    print "Out of "+ str(len(all_possible_lsvs))+" unique LSV IDs in both dictionaries:"
    print str(len(in_both_sig_junctions)) +" junctions out of "+ str(len(in_both_sig_ids))+" LSVs in both"               
    print str(len(only_1_junctions)) +" junctions out of "+ str(len(only_1_ids))+" LSVs in first only"               
    print str(len(only_2_junctions)) +" junctions out of "+ str(len(only_2_ids))+" LSVs in second only"   
                
    return set(in_both_sig_ids), only_1_ids, only_2_ids
   
def find_extreme_diffs(LSV_dict_1, LSV_dict_2, Cutoff_dPSI = CUTOFF_DPSI, Cutoff_prob= CUTOFF_PROB):
    """
        Given two LSV dictionaries, return list of LSV IDs containing at least 
            one identical junction that is significantly changing in both, 
            but in opposite directions.
    """
    lsv_ids_1 = get_name_set(LSV_dict_1, Cutoff = Cutoff_dPSI, probability_dPSI = Cutoff_prob)
    lsv_ids_2 = get_name_set(LSV_dict_2, Cutoff = Cutoff_dPSI, probability_dPSI = Cutoff_prob)
    
    in_both = lsv_ids_1.intersection(lsv_ids_2)
    
    extreme_diffs = list()
    
    # For each LSV containing a highly significant dPSI in both dictionaries
    for lsv_id in in_both:
#         LSV_dict_1[lsv_id]["dPSI_over_cutoff"]
        sig_indeces_from_1 =  LSV_dict_1[lsv_id]["sig_junctions"]
        sig_indeces_from_2 =  LSV_dict_2[lsv_id]["sig_junctions"]
        
        sig_juntion_indeces_in_both = set(sig_indeces_from_1).intersection(set(sig_indeces_from_2))
        
        for index in sig_juntion_indeces_in_both:   
            dPSI_1 = LSV_dict_1[lsv_id]["E(dPSI) per LSV junction"][index]
            dPSI_2 = LSV_dict_2[lsv_id]["E(dPSI) per LSV junction"][index]
            # Both positive?
            if ((dPSI_1 > 0) and (dPSI_2 > 0)):
                continue
            
            # Both negative?
            if ((dPSI_1 < 0) and (dPSI_2 < 0)):
                continue  
                      
            # 1st is positive, second is negative?
            if ((dPSI_1 > 0) and (dPSI_2 < 0)):
                extreme_diffs.append(lsv_id)
                
            # 1st is negative, second is 0 or positive?
            if ((dPSI_1 < 0) and (dPSI_2 > 0)):
                extreme_diffs.append(lsv_id)
    
    print "Found " + str(len(extreme_diffs))+ " junctions showing opposite dPSI corresponding to "+str(len(set(extreme_diffs)))+" LSVs."
    
    return list(set(extreme_diffs))
                
def plot_venn(List_of_sets,
              Set_labels,
              Main = "I forgot to give this plot a name.",
              Out_File = "",
              Custom_overlap_numbers = []):
    """
        Given a list of sets, generate a venn diagram in Out_Dir.
        
        Arguments:
            List_of_sets (two or three only!)
            Set_labels: Label for each circle
            Main: Title of plot
            Out_File: Where should plot be saved? And what should the file be named?
                Parent directory expected to already exist...
                This will overwrite plots if they already exist
            Custom_overlap_numbers: optional. If you want to supply your own 3 overlap sets:
                [# in first, # in second, # in both]
    """
    if not os.path.isdir(os.path.dirname(Out_File)):
        raise ValueError(os.path.dirname(Out_File)+" <--- PATH DOESN'T EXIST")
    
    if len(Custom_overlap_numbers) != 0 and len(Custom_overlap_numbers) != 3:
        raise ValueError("Custom overlap only works for 2 circle venn diagrams at the moment...")
    
    if len(Custom_overlap_numbers) == 3:
        plt.figure()
        venn2(subsets={'10': Custom_overlap_numbers[0], '01': Custom_overlap_numbers[1], '11': Custom_overlap_numbers[2]}, set_labels = Set_labels)
        plt.title(Main)
        plt.savefig(Out_File)
        return
                
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
    
def all_pairwise(List_1, List_2, Order_Matters = False):
    """
        Given two lists, return list of lists of all
            possible pairwise comparisons.
            
        If List_1 == List_2, this function is effectively performing:
            List CHOOSE/COMBINATION 2, if Order_Matters = False
            List PICK/PERMUTATION 2, if Order_Matters = True
            
        Returns [[a,b],[a,c],[etc,etc]] of combinations/permutations
    """
    if type(List_1) is not list:
        raise ValueError("List_1 needs to be a list...")
    
    if type(List_2) is not list:
        raise ValueError("List_2 needs to be a list...")
    
    if type(Order_Matters) is not bool:
        raise ValueError("Order_Matters needs to be True or False...")                        
                         
    tuples = [(x,y) for x in List_1 for y in List_2 if x != y]
    for entry in tuples:
        if not Order_Matters:
            if (entry[1], entry[0]) in tuples:
                tuples.remove((entry[1],entry[0]))
    
    matrix = list()
    
    for t in tuples:
        matrix.append(list(t))
    return matrix         
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            