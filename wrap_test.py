import os
import sys

sys.path.append("/home/cradens/home_base/splicing/")
from splicing_fun import *

afp="/home/cradens/home_base/splicing/spurlock_dpsi/Th1_e_A_Th2_e_A.deltapsi_quantify_deltapsi.txt"

bfp="/home/cradens/home_base/splicing/spurlock_dpsi/Th1_e_A_Th17_e_A.deltapsi_quantify_deltapsi.txt"

a=import_dpsi(afp)
a_label = a["condition_1_name"]+"_vs_"+a["condition_2_name"]

b=import_dpsi(bfp)
b_label = b["condition_1_name"]+"_vs_"+b["condition_2_name"]

cutoff = 0.0
an = get_name_set(a, cutoff)
bn = get_name_set(b, cutoff)

list_of_sets = list()
list_of_sets.append(an)
list_of_sets.append(bn)

print a_label
print b_label

plot_venn(list_of_sets, ["TH1vTH2", "TH1vTH17"],
		 "test_venn",
		 "/home/cradens/home_base/splicing/test_venn")
