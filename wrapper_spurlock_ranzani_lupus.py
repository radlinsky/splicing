import sys

folder = "/home/cradens/script/"

sys.path.append(folder)

from get_exec_lines import *

# bash file that will run majiq is being written to (and will be run from)  this location:
fw = open('/data/THelper/spurlock_ranzani_sullivan/run_MAJIQ_all_pairwise.bash','w')

fw.write('majiq build /data/DB/hg19/ensembl.hg19.gff3 -conf /data/THelper/spurlock_ranzani_sullivan/settings.txt --output ./build --nthreads 16\n')

print "RANZANI"

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_TN','Ra_T1','./build','./majiq/dpsi/TNvsT1','./voila/dpsi/TNvsT1',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_TN','Ra_T2','./build','./majiq/dpsi/TNvsT2','./voila/dpsi/TNvsT2',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_TN','Ra_T17','./build','./majiq/dpsi/TNvsT17','./voila/dpsi/TNvsT17',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_TN','Ra_TR','./build','./majiq/dpsi/TNvsTR','./voila/dpsi/TNvsTR',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_TN','Ra_TCM','./build','./majiq/dpsi/TNvsTCM','./voila/dpsi/TNvsTCM',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_TN','Ra_TEM','./build','./majiq/dpsi/TNvsTEM','./voila/dpsi/TNvsTEM',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_T1','Ra_T2','./build','./majiq/dpsi/T1vsT2','./voila/dpsi/T1vsT2',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_T1','Ra_T17','./build','./majiq/dpsi/T1vsT17','./voila/dpsi/T1vsT17',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt','Ra_T2','Ra_T17','./build','./majiq/dpsi/T2vsT17','./voila/dpsi/T2vsT17',show_all=False)

fw.write(m1+v1)

print "SPURLOCK"

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T1_e_A','Sp_T2_e_A','./build','./majiq/dpsi/T1e_T2e','./voila/dpsi/T1e_T2e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T1_e_A','Sp_T17_e_A','./build','./majiq/dpsi/T1e_T17e','./voila/dpsi/T1e_T17e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T2_e_A','Sp_T17_e_A','./build','./majiq/dpsi/T2e_T17e','./voila/dpsi/T2e_T17e',show_all=False)

fw.write(m1+v1)

print "Primary: T1<->T2, T1<->T17, T2<->T17"

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T1_p_A','Sp_T2_p_A','./build','./majiq/dpsi/T1p_T2p','./voila/dpsi/T1p_T2p',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T1_p_A','Sp_T17_p_A','./build','./majiq/dpsi/T1p_T17p','./voila/dpsi/T1p_T17p',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T2_p_A','Sp_T17_p_A','./build','./majiq/dpsi/T2p_T17p','./voila/dpsi/T2p_T17p',show_all=False)

fw.write(m1+v1)

print "Primary <-> Effector: T1, T2, T17"

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T1_p_A','Sp_T1_e_A','./build','./majiq/dpsi/T1p_T1e','./voila/dpsi/T1p_T1e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T2_p_A','Sp_T2_e_A','./build','./majiq/dpsi/T2p_T2e','./voila/dpsi/T2p_T2e',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T17_p_A','Sp_T17_e_A','./build','./majiq/dpsi/T17p_T17e','./voila/dpsi/T17p_T17e',show_all=False)

fw.write(m1+v1)

print "PolyA <-> Total: primary <-> primary"

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T1_p_A','Sp_T1_p_r','./build','./majiq/dpsi/T1p_PolyAvTotal','./voila/dpsi/T1p_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T2_p_A','Sp_T2_p_r','./build','./majiq/dpsi/T2p_PolyAvTotal','./voila/dpsi/T2p_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T17_p_A','Sp_T17_p_r','./build','./majiq/dpsi/T17p_PolyAvTotal','./voila/dpsi/T17p_PolyAvTotal',show_all=False)

fw.write(m1+v1)

print "PolyA <-> Total: effector <-> effector"

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T1_e_A','Sp_T1_e_r','./build','./majiq/dpsi/T1e_PolyAvTotal','./voila/dpsi/T1e_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T2_e_A','Sp_T2_e_r','./build','./majiq/dpsi/T2e_PolyAvTotal','./voila/dpsi/T2e_PolyAvTotal',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Sp_T17_e_A','Sp_T17_e_r','./build','./majiq/dpsi/T17e_PolyAvTotal','./voila/dpsi/T17e_PolyAvTotal',show_all=False)

fw.write(m1+v1)

print "SULLIVAN"

m1, v1 = get_EXECdpsi_line('./settings.txt','Su_T_ND','Su_T_SLE','./build','./majiq/dpsi/T_NDvSLE','./voila/dpsi/T_NDvSLE',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Su_B_ND','Su_B_SLE','./build','./majiq/dpsi/B_NDvSLE','./voila/dpsi/B_NDvSLE',show_all=False)

fw.write(m1+v1)

m1, v1 = get_EXECdpsi_line('./settings.txt','Su_M_ND','Su_M_SLE','./build','./majiq/dpsi/M_NDvSLE','./voila/dpsi/M_NDvSLE',show_all=False)

fw.write(m1+v1)

print "SULLIVAN T vs RANZANI T NAIVE"

m1, v1 = get_EXECdpsi_line('./settings.txt','Su_T_ND','Ra_TN','./build','./majiq/dpsi/NDvsNAIVE','./voila/dpsi/NDvsNAIVE',show_all=False)

fw.write(m1+v1)


samps = ["Ra_TN","Ra_T1","Ra_T2","Ra_T17","Ra_TR","Ra_TCM","Ra_TEM","Sp_T1_e_r","Sp_T1_e_A","Sp_T2_e_r","Sp_T2_e_A", "Sp_T17_e_r","Sp_T17_e_A", "Su_T_ND","Su_T_SLE","Su_B_ND","Su_B_SLE","Su_M_ND","Su_M_SLE"]

# Psi calculations
for s in samps:
    m1,v1 = get_EXECpsi_line('/data/THelper/spurlock_ranzani_sullivan/settings.txt',s,'./build/','./majiq/psi/%s'%s,'./voila/psi/%s'%s)
    fw.write(m1+v1)
   
fw.close()
