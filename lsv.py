#/usr/bin/python

# lsv.py
# Caleb Matthew Radens
# 2016_5_14

##### API for analyzing voila text files.

import os
import pdb

class LSV(object):
    """
    Data container for one or more voila dpsi text files that comprises methods
        to import/export, analyze, and compare local splice variants.
    """
    
    # Default E(dPSI) and Prob(E(dPSI)) values.
    _CUTOFF_DPSI = 0.2
    _CUTOFF_PROB = 0.95
    
    # Expected order of columns from a voila text file:
    _DPSI_HEADER = ["Gene Name",
                   "Gene ID",
                   "LSV ID",
                   "E(dPSI) per LSV junction",
                   "P(|E(dPSI)|>=0.20) per LSV junction",
                   "E(PSI)",     # This header also has the 1st condition name
                   "E(PSI)",     # This header also has the 2nd condition name
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
    
    # TODO: should these be properties?
    n_lsvs = 0
    condition_1_name = ""
    condition_2_name = ""
    lsv_ids = []
    lsv_dictionary = dict()
    
    
    def __init__(self,match_ids="all", file_path = "", user_friendly = False):
        if not isinstance(file_path, str):
            raise ValueError("File path needs to be a string.")
        if user_friendly:
            if len(file_path) > 0:
                raise IOError("User friendly mode is set to True, but a "
                    "non-empty string argument was provided for file_path.")
            self.file_path = self.user_friendly_mode()
        else:
            self.is_this_a_file(file_path)
            self.is_this_a_voila_file(file_path)
            self.__file_path = file_path
        
        self.import_voila(match_lsv_ids=match_ids)
        
    def import_voila(self, match_lsv_ids):
        
        try:
            assert isinstance(match_lsv_ids, list)
            
            assert len(match_lsv_ids)>0
            
            try_to_match = True
        except:
            if not match_lsv_ids=="all":
                raise ValueError("If not a list, match_ids needs to be 'all'")
            
            try_to_match = False
            
        with open(self.__file_path, "rb") as handle:
            line_i = 0
            for line in handle:
                line_split = line.rstrip("\r\n").split("\t")
                
                if line_i == 0:
                    condition_one_name = line_split[5].split(" ")[0]
                    condition_two_name = line_split[6].split(" ")[0]
                    line_i += 1
                    continue
                
                LSV_ID = str(line_split[2])
                
                
                if try_to_match:
                    if LSV_ID not in match_lsv_ids:
                        
                        line_i+=1
                        continue
                    print "Found",LSV_ID
    
                gene_name = str(line_split[0])
                
                gene_id = str(line_split[1])
                
                dPSIs_floated = map(float,str(line_split[3]).split(";"))
                
                prob_dPSIs_floated = map(float,str(line_split[4]).split(";"))
                
                # Saves the junction # and directionality of change
                over_cutoff = list()
                
                # Just saves which junction index for those 
                #     that are significant
                sig_junctions = list()
                
                # Junction index
                i = 0
                for dPSI, prob_dPSI in zip(dPSIs_floated, prob_dPSIs_floated):
                    # If the dPSI is greater than 0 
                    #     and over the cutoff dPSI and probability
                    if ((dPSI > self._CUTOFF_DPSI) and 
                        (prob_dPSI > self._CUTOFF_PROB)):
                        # junctionIndex_over
                        over_cutoff.append(str(i)+"_over")
                        sig_junctions.append(i)
                    # Else if the dPSI is less than 0
                    #     and over the cutoff dPSI and probability
                    elif ((dPSI < -self._CUTOFF_DPSI) 
                          and (prob_dPSI > self._CUTOFF_PROB)):
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
                
                chrom = str(line_split[14])
                
                strand = str(line_split[15])
                
                junct_coord = str(line_split[16])
                
                exon_coord = str(line_split[17]).split(";")
                
                exon_alt_start = str(line_split[18])
                
                exon_alt_end = str(line_split[19])
                
                ir_coords = str(line_split[20])
                
                voila_link = str(line_split[21])
                
                self.lsv_dictionary[LSV_ID] = dict({
                                    self._DPSI_HEADER[0]:gene_name,
                                    self._DPSI_HEADER[1]:gene_id,
                                    self._DPSI_HEADER[3]:dPSIs_floated,
                                    self._DPSI_HEADER[4]:prob_dPSIs_floated,
                                    "dPSI_over_cutoff":over_cutoff,
                                    "sig_junctions":sig_junctions,
                                    self._DPSI_HEADER[5]:ePSI_1_floated,
                                    self._DPSI_HEADER[6]:ePSI_2_floated,
                                    self._DPSI_HEADER[7]:lsv_type,
                                    self._DPSI_HEADER[8]:A5SS,
                                    self._DPSI_HEADER[9]:A3SS,
                                    self._DPSI_HEADER[10]:ES,
                                    self._DPSI_HEADER[11]:n_junctions,
                                    self._DPSI_HEADER[12]:n_exons,
                                    self._DPSI_HEADER[13]:de_novo_junct,
                                    self._DPSI_HEADER[14]:chrom,
                                    self._DPSI_HEADER[15]:strand,
                                    self._DPSI_HEADER[16]:junct_coord,
                                    self._DPSI_HEADER[17]:exon_coord,
                                    self._DPSI_HEADER[18]:exon_alt_start,
                                    self._DPSI_HEADER[19]:exon_alt_end,
                                    self._DPSI_HEADER[20]:ir_coords,
                                    self._DPSI_HEADER[21]:voila_link})
                
                line_i += 1
            
            self.condition_1_name = condition_one_name
            self.condition_2_name = condition_two_name
            
            if try_to_match:
#                 pdb.set_trace()
                if len(self.lsv_dictionary.keys()) < len(match_lsv_ids):
                    raise ValueError("Not all IDs found in file...")
                #print str(len(match_lsv_ids))+ " LSVs extracted from "\
                #                +os.path.basename(self.__file_path)
                pass
            else:
                #print str(line_i-1)+ " LSVs extracted from "\
                #                +os.path.basename(self.__file_path)
                pass
                             
        self.lsv_ids = self.lsv_dictionary.keys()
        self.n_lsvs = len(self.lsv_ids)
                             
    def lookup_gene(self, Gene_name, Ensembl_id=False):
        """
            Look up Gene name OR Ensembl ID
            
            If Ensembl_id = True, then use ENSG00... instead.
            
            Returns new LSV object with only LSVs in the gene.
        """
        if not isinstance(Gene_name, str):
            raise ValueError("Gene_name needs to be a string.")
        if not isinstance(Ensembl_id, bool):
            raise ValueError("Ensembl_id needs to be True/False.")
        
        matched_ids = list()
        for lsvid in self.lsv_ids:
            if not Ensembl_id:
                if self.lsv_dictionary[lsvid]["Gene Name"] == Gene_name:
                    matched_ids.append(lsvid)
            elif Ensembl_id:
                if self.lsv_dictionary[lsvid]["Gene ID"] == Gene_name:
                    matched_ids.append(lsvid)   
                    
        if len(matched_ids) == 0:
            raise ValueError(Gene_name+" wasn't found.")
        
        #print "Found "+str(len(matched_ids))+" LSV in gene: "+Gene_name
        return self.subset(matched_ids)
    
    def lsv_id_check(self, LSV_ID):
        """
            Check if an LSV ID exists in this LSV object.
            
            Return True or False.
        """
        if not isinstance(LSV_ID, str):
            raise ValueError("LSV_ID needs to be a string.")
        if not LSV_ID in self.lsv_ids:
            return False
        else:
            return True
    
    def subset(self, subset_lsv_ids = "all"):
        """
        Return a subset of this LSV, based on provided lsvids.
        """
        if subset_lsv_ids != "all":
            for lsv_id in subset_lsv_ids:
                found_id = self.lsv_id_check(lsv_id)
                if not found_id:
                    raise ValueError("LSV ID "+lsv_id+" not found in file.")
            
        return LSV(file_path=self.__file_path,
                   user_friendly=False,
                   match_ids=subset_lsv_ids)
    
    def user_friendly_mode(self):
        """
            Interact with user in terminal. Allows user to drag n' drop voila
                text file. 
        """
        user_arg = input("Please enter a voila file path:")
        self.is_this_a_file(user_arg)
        self.is_this_a_voila_file(user_arg)
        self.__file_path = user_arg
    
    def is_this_a_file(self,path):
        if not isinstance(path, str):
            raise ValueError("Expected a string.")
        if len(path) == 0:
            raise IOError("You forgot to provide a file path.")
        if not (os.path.isdir(path) or os.path.isfile(path)):
            raise IOError("Please only provide an extant file path.")
                
    def is_this_a_voila_file(self, path):
        if not ".deltapsi_quantify_deltapsi.txt" in path:
            raise ValueError("That's not a voila text file...")
        with open(path, "rb") as handle:
            for line in handle:
                n_cols = len(line.rstrip("\n\r").split("\t"))
                if not n_cols == len(self._DPSI_HEADER):
                    raise IOError("Provided file didn't have the expected"
                                  " number of columns for a voila txt file...")
                break
    
    
    # TODO: Should these be set/get types of things or properties or something?
    def __get_path__(self):
        print "This LSV was constructed from: \n   "+self.__file_path
        
    def __len__(self):
        return self.n_lsvs
    
    