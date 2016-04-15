import os
import sys

sys.path.append("/home/cradens/home_base/splicing/")
from splicing_fun import *

# One junction in LSV must be great than or equal to:
cutoff = 0.2

TH1e_TH17e="/home/cradens/home_base/splicing/spurlock_dpsi/Th1_e_A_Th17_e_A.deltapsi_quantify_deltapsi.txt"

TH1e_TH17e_imp=import_dpsi(TH1e_TH17e)
TH1e_TH17e_no_int = no_intron_retention(TH1e_TH17e_imp)
TH1e_TH17e_names = get_name_set(TH1e_TH17e_no_int, cutoff)
TH1e_TH17e_label = TH1e_TH17e_no_int["condition_1_name"]+"_vs_"+TH1e_TH17e_no_int["condition_2_name"]
print "TH1e_TH17e_label label is: "+TH1e_TH17e_label

N_1 = "/home/cradens/home_base/splicing/ranzani_dpsi/T4_N_T4_1.deltapsi_quantify_deltapsi.txt"
N_2 = "/home/cradens/home_base/splicing/ranzani_dpsi/T4_N_T4_2.deltapsi_quantify_deltapsi.txt"
N_17 = "/home/cradens/home_base/splicing/ranzani_dpsi/T4_N_T4_17.deltapsi_quantify_deltapsi.txt"
one_2 = "/home/cradens/home_base/splicing/ranzani_dpsi/T4_1_T4_2.deltapsi_quantify_deltapsi.txt"
one_17 = "/home/cradens/home_base/splicing/ranzani_dpsi/T4_1_T4_17.deltapsi_quantify_deltapsi.txt"
two_17 = "/home/cradens/home_base/splicing/ranzani_dpsi/T4_2_T4_17.deltapsi_quantify_deltapsi.txt"

N_1_imp = import_dpsi(N_1)
N_1_no_int = no_intron_retention(N_1_imp)
N_1_names = get_name_set(N_1_no_int, cutoff)

N_2_imp = import_dpsi(N_2)
N_2_no_int = no_intron_retention(N_2_imp)

N_17_imp = import_dpsi(N_17)
N_17_no_int = no_intron_retention(N_17_imp)
N_17_names = get_name_set(N_17_no_int, cutoff)

one_2_imp = import_dpsi(one_2)
one_2_no_int = no_intron_retention(one_2_imp)

one_17_imp = import_dpsi(one_17)
one_17_no_int = no_intron_retention(one_17_imp)
one_17_names = get_name_set(one_17_no_int, cutoff)
one_17_label = one_17_no_int["condition_1_name"]+"_vs_"+one_17_no_int["condition_2_name"]

two_17_imp = import_dpsi(two_17)
two_17_no_int = no_intron_retention(two_17_imp)


# print "Overlap of sets: "+str(len(N_1_names & N_17_names))
# for in_both in (N_1_names & N_17_names):
# 	print N_1_no_int[in_both]["Gene Name"]

print "Only N_1: " + str(len(N_1_names.difference(N_17_names)))
for n1 in (N_1_names.difference(N_17_names)):
	pass
	#print N_1_no_int[n1]["Gene Name"]

print "Only N_17: " + str(len(N_17_names.difference(N_1_names)))
for n17 in (N_17_names.difference(N_1_names)):
	pass
	#print N_17_no_int[n17]["Gene Name"]

list_of_sets = list()
list_of_sets.append(N_1_names)
list_of_sets.append(N_17_names)

plot_venn(list_of_sets, ["Naive<->TH1", "Naive<->TH17"],
		 "Naive to TH1 or TH17 (Ranzani)",
		 "/home/cradens/home_base/splicing/ranzani_N_1or17")

print "Overlap of TH1<->TH17: "+str(len(one_17_names & TH1e_TH17e_names))
for in_both in (one_17_names & TH1e_TH17e_names):
	print one_17_no_int[in_both]["Gene Name"]

print "Only TH1_TH17 Ranzani: " + str(len(one_17_names.difference(TH1e_TH17e_names)))
for one17_ranz in (one_17_names.difference(TH1e_TH17e_names)):
	pass
	#print one_17_no_int[one17_ranz]["Gene Name"]

print "Only TH1_TH17 Spurlock: " + str(len(TH1e_TH17e_names.difference(one_17_names)))
for one17_spur in (TH1e_TH17e_names.difference(one_17_names)):
	pass
	#print TH1e_TH17e_no_int[one17_spur]["Gene Name"]

list_of_sets = list()
list_of_sets.append(one_17_names)
list_of_sets.append(TH1e_TH17e_names)

plot_venn(list_of_sets, ["Ranzani", "Spurlock"],
		 "Ranzani vs. Spurlock : [TH1 vs. TH17]",
		 "/home/cradens/home_base/splicing/ranzani_vs_spurlock_1v17")

consistently_different = (one_17_names & TH1e_TH17e_names)
surely_TH1 = (consistently_different & N_1_names)
print "surely_TH1: "+str(surely_TH1)
surely_TH17= (consistently_different & N_17_names)
print "surely_TH17: "+str(surely_TH17)

print one_17_no_int['ENSG00000060237:990637-990955:target']["Gene Name"]
print TH1e_TH17e_no_int['ENSG00000060237:990637-990955:target']["Gene Name"]
print N_17_no_int['ENSG00000060237:990637-990955:target']["Gene Name"]

print one_17_no_int['ENSG00000060237:990637-990955:target']["E(dPSI) per LSV junction"]
print TH1e_TH17e_no_int['ENSG00000060237:990637-990955:target']["E(dPSI) per LSV junction"]
print N_17_no_int['ENSG00000060237:990637-990955:target']["E(dPSI) per LSV junction"]

print one_17_no_int['ENSG00000060237:990637-990955:target']["P(|E(dPSI)|>=0.20) per LSV junction"]
print TH1e_TH17e_no_int['ENSG00000060237:990637-990955:target']["P(|E(dPSI)|>=0.20) per LSV junction"]
print N_17_no_int['ENSG00000060237:990637-990955:target']["P(|E(dPSI)|>=0.20) per LSV junction"]
